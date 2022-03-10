# The functions
# > "_kelley",
# > "_solve_Lagrangian_relaxation",
# are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and Lea Kapelevich released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson, Lea Kapelevich

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

################################################################################

# The function
# > "_bundle_level"
# is derived from a similar named function in the 'SDDiP.jl' package by
# Lea Kapelevich released under the MIT Expat License.
# This specific function is also relased under MIT Expat License.

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2017: LEAXPS-15\lkape.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

################################################################################

"""
Solving the Lagrangian relaxation problem, i.e. the inner problem of the
Lagrangian dual
"""

function _solve_Lagrangian_relaxation!(
    node::SDDP.Node,
    π_k::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Set the Lagrangian relaxation of the objective in the primal model
    old_obj = JuMP.objective_function(model)
    JuMP.set_objective_function(model, JuMP.@expression(model, old_obj - π_k' * h_expr))

    # Optimization
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == MOI.OPTIMAL

    # Update the correct values
    L_k = JuMP.objective_value(model)
    if update_subgradients
        h_k .= -JuMP.value.(h_expr)
    end

    # Reset old objective
    JuMP.set_objective_function(model, old_obj)

    return L_k
end


"""
Kelley's method to solve Lagrangian dual
"""
function solve_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    primal_obj::Float64,
    π_k::Vector{Float64},
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.Kelley
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)
    # The original value for x (former old_rhs)
    x_in_value = zeros(number_of_states)
    # The current estimate for π (in our case determined in initialization)
    # π_k
    # The best estimate for π (former best_mult)
    π_star = zeros(number_of_states)
    # The best estimate for the dual objective value (ignoring optimization sense)
    L_star = -Inf
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states)
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = algo_params.duality_regime.atol
    rtol = algo_params.duality_regime.rtol
    iteration_limit = algo_params.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, algo_params.state_approximation_regime)
    end
    node.ext[:backward_data][:old_rhs] = x_in_value
    @infiltrate

    ############################################################################
    # LOGGING OF LAGRANGIAN DUAL
    ############################################################################
    #lag_log_file_handle = open("C:/Users/cg4102/Documents/julia_logs/Lagrange.log", "a")
    #print_helper(print_lagrange_header, lag_log_file_handle)

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_approx_model" begin
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        approx_model = JuMP.Model(GAMS.Optimizer)
        set_solver!(approx_model, algo_params, applied_solvers, :kelley)

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        set_objective_bound!(approx_model, s, bound_results.obj_bound)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Magnanti Wong part later
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        set_multiplier_bounds!(node, approx_model, number_of_states, bound_results.dual_bound, algo_params)
    end

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = 0 # why zero?
    #-inf is not possible, since then the while loop would not start at all

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        end

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
        end

        ########################################################################
        # ADD CUTTING PLANE
        ########################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "add_cut" begin
            JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))
        end

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        TimerOutputs.@timeit DynamicSDDiP_TIMER "outer_sol" begin
            JuMP.optimize!(approx_model)
        end
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        π_k .= JuMP.value.(π)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        #print("UB: ", f_approx, ", LB: ", f_actual)
        #println()

        ########################################################################
        if L_star > t_k + atol/10.0
            error("Could not solve for Lagrangian duals. LB > UB.")
        end
    end

    ############################################################################
    # CONVERGENCE ANALYSIS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "convergence_check" begin
        if isapprox(L_star, t_k, atol = atol, rtol = rtol)
            # CONVERGENCE ACHIEVED
            if isapprox(L_star, s * primal_obj, atol = atol, rtol = rtol)
                # CONVERGENCE TO TRUE OPTIMUM (APPROXIMATELY)
                lag_status = :opt
            else
                # CONVERGENCE TO A SMALLER VALUE THAN THE PRIMAL OBJECTIVE
                # sometimes this occurs due to numerical issues
                # still leads to a valid cut
                lag_status = :conv
            end
        elseif all(h_k .== 0)
            # NO OPTIMALITY ACHIEVED, BUT STILL ALL SUBGRADIENTS ARE ZERO
            # may occur due to numerical issues
            lag_status = :sub
        elseif iter == iteration_limit
            # TERMINATION DUE TO ITERATION LIMIT
            # stil leads to a valid cut
            lag_status = :iter
        end
    end

    ############################################################################
    # APPLY MAGNANTI AND WONG APPROACH IF INTENDED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "magnanti_wong" begin
        magnanti_wong!(node, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_star,
            iteration_limit, atol, rtol, algo_params.duality_regime.dual_choice_regime, iter)
    end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, algo_params.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

    ############################################################################
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k = π_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(state)
        # Store expression for slack
        h_expr[i] = JuMP.@expression(node.subproblem, state - x_in_value[i])
        # Relax copy constraint (i.e. z does not have to take the value of ̄x anymore)
        JuMP.unfix(state)

        # Set bounds to ensure that inner problems are feasible
        # As we use binary approximation, 0 and 1 can be used
        JuMP.set_lower_bound(state, 0)
        JuMP.set_upper_bound(state, 1)

    end

    return
end

function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    for (i, (_, state)) in enumerate(node.states)
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(state.in)
        # Store expression for slack
        h_expr[i] = JuMP.@expression(node.subproblem, state.in - x_in_value[i])
        # Relax copy constraint (i.e. z does not have to take the value of ̄x anymore)
        JuMP.unfix(state.in)

        # Set bounds to ensure that inner problems are feasible
        # Bound shouldn't be too tight, so use 1e9 as a default if nothing
        # else is specified
        lb = has_lower_bound(state.out) ? lower_bound(state.out) : -1e9
        ub = has_upper_bound(state.out) ? upper_bound(state.out) : 1e9
        JuMP.set_lower_bound(state.in, lb)
        JuMP.set_upper_bound(state.in, ub)
    end

    return
end

function restore_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    for (i, (_, bin_state)) in enumerate(node.ext[:backward_data][:bin_states])
        # prepare_state_fixing!(node, state_comp)
        JuMP.fix(bin_state, x_in_value[i], force = true)
    end

    return
end

function restore_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    for (i, (_, state)) in enumerate(node.states)
        # prepare_state_fixing!(node, state_comp)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end

    return
end

function set_objective_bound!(approx_model::JuMP.Model, s::Int, obj_bound::Float64)

    JuMP.set_upper_bound(approx_model[:t], s * obj_bound)

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model, number_of_states::Int, dual_bound::Float64, algo_params::DynamicSDDiP.AlgoParams)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]

    for (i, (key, value)) in enumerate(node.ext[:backward_data][:bin_states])
    	associated_original_state = node.ext[:backward_data][:bin_x_names][key]
    	beta = algo_params.state_approximation_regime.binary_precision[associated_original_state]
    	associated_k = node.ext[:backward_data][:bin_k][key]
    	bound = dual_bound * 2^(associated_k-1) * beta
    	JuMP.set_upper_bound(π⁺[i], bound)
        JuMP.set_upper_bound(π⁻[i], bound)
	end
end

"""
Given the optimal dual objective value from the Kelley's method, try to find
the optimal dual multipliers with the smallest L1-norm.

Note that this is only done until the maximum number of iterations is
achieved in total.
"""

function magnanti_wong!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    t_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    s::Int,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.MagnantiWongChoice,
    iter::Int,
    )

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    t = approx_model[:t]
    π = approx_model[:π]

    # Reset objective
    JuMP.@objective(approx_model, Min, sum(π⁺) + sum(π⁻))
    JuMP.set_lower_bound(t, t_k)

    # The worst-case scenario in this for-loop is that we run through the
    # iterations without finding a new dual solution. However if that happens
    # we can just keep our current λ_star.
    for _ in (iter+1):iteration_limit
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        π_k .= value.(π)
        L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        if isapprox(L_star, L_k, atol = atol, rtol = rtol)
            # At this point we tried the smallest ‖π‖ from the cutting plane
            # problem, and it returned the optimal dual objective value. No
            # other optimal dual vector can have a smaller norm.
            π_star = π_k
            return
        end
        JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))

    end

    return
end


function magnanti_wong!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    t_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    s::Int,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.StandardChoice,
    iter::Int,
    )

    return
end



"""
Level Bundle method to solve Lagrangian dual
"""
function solve_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    primal_obj::Float64,
    π_k::Vector{Float64},
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.LevelBundle
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)
    # The original value for x (former old_rhs)
    x_in_value = zeros(number_of_states)
    # The current estimate for π (in our case determined in initialization)
    # π_k
    # The best estimate for π (former best_mult)
    π_star = zeros(number_of_states)
    # The best estimate for the dual objective value (ignoring optimization sense)
    L_star = -Inf
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states)
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = algo_params.duality_regime.atol
    rtol = algo_params.duality_regime.rtol
    iteration_limit = algo_params.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax)

    # Set bundle_parameters
    #---------------------------------------------------------------------------
    level_factor = dual_solution_regime.level_factor

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, algo_params.state_approximation_regime)
    end
    node.ext[:backward_data][:old_rhs] = x_in_value

    ############################################################################
    # LOGGING OF LAGRANGIAN DUAL
    ############################################################################
    #lag_log_file_handle = open("C:/Users/cg4102/Documents/julia_logs/Lagrange.log", "a")
    #print_helper(print_lagrange_header, lag_log_file_handle)

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_approx_model" begin
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        approx_model = JuMP.Model(GAMS.Optimizer)

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        set_objective_bound!(approx_model, s, bound_results.obj_bound)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Magnanti Wong part later
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        set_multiplier_bounds!(approx_model, number_of_states, bound_results.dual_bound, algo_params)
    end

    ############################################################################
    # BUNDLE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = 0

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        end
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]
        #@infiltrate

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
        end

        ########################################################################
        # ADD CUTTING PLANE
        ########################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "add_cut" begin
            JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))
        end

        ########################################################################
        # RESET OBJECTIVE FOR APPROX_MODEL AFTER NONLINEAR MODEL
        ########################################################################
        JuMP.@objective(approx_model, Max, t)
        set_solver!(approx_model, algo_params, applied_solvers, :kelley)

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        TimerOutputs.@timeit DynamicSDDiP_TIMER "outer_sol" begin
            JuMP.optimize!(approx_model)
        end
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        #print("UB: ", f_approx, ", LB: ", f_actual)
        #println()

        ########################################################################
        # COMPUTE GAP AND FORM A NEW LEVEL
        ########################################################################
        f_up = t_k # objective value of the approx_model
        f_down = L_star # best lower bound so far
        gap = f_up - f_down

        """
        We use a convex combination of f_up and f_down for the new level.
        level = (1-a) f_up + a f_down = f_up - a (f_up - f_down) = f_up - a gap

        For level_factor = 0 => level = f_up, i.e. t_k, that means, this
        reduces to the basic cutting-plane method.

        For level_factor = 1 => level = f_down, i.e. L_star. That means that
        we search for the closest multiplier, for which the level is larger than
        L_star. At least for L_star = L_k then the multiplier does not change
        anymore, because the condition is trivially satisfied for the current
        multiplier. Otherwise, this should yield one of the previous multipliers,
        so again no improvement.
        """

        # TODO: POSSIBLE IMPROVEMENTS/CHANGES
        """
        1.) STRUCTURE
        In the literature often a different order of the steps is proposed,
        so maybe we should change this.

        In particular, the gap is often determined before the approx_model
        is solved with the new multipliers. That means that the gap is determined
        using f_down and f_up = t_{k-1}. In this case, t_0 has to be determined
        appropriately, e.g. we can just choose primal_obj.

        2.) STOPPING
        In the literature it is often proposed to stop if the gap is sufficiently
        small. In our case, this doesn't make a difference, but since the gap
        can be determined differently as well, then our stopping criterion
        would change.

        3.) STABILITY CENTER
        We always choose the new multiplier obtained by solving the proximal
        problem as new stability center. It is also possible to change the
        stability center only if a sufficiently large improvement is achieved.

        4.) DETERMINING f_up
        Instead of solving approx_model, we could also determine f_up in the
        following way:
        > If the proximal problem is infeasible, set f_up = level.
        > If the proximal problem is feasible, do not change f_up.

        5.) DETERMINE gap
        Right now, the gap is determined as f_up - f_down with
        f_up = t_k
        f_down = L_star

        We could also choose f_up = primal_obj, since we already know that
        this is the best possible upper bound. Is this beneficial?
        """

        level = f_up - gap * level_factor
        # - atol/10.0 for numerical issues?
        JuMP.set_lower_bound(t, level)

        ########################################################################
        # DETERMINE NEXT ITERATION USING PROXIMAL PROBLEM
        ########################################################################
        # Objective function of approx model has to be adapted to new center
        # TODO: Does this work with π[i]?
        JuMP.@objective(approx_model, Min, sum((π_k[i] - π[i])^2 for i in 1:number_of_states))
        set_solver!(approx_model, algo_params, applied_solvers, :level_bundle)
        TimerOutputs.@timeit DynamicSDDiP_TIMER "bundle_sol" begin
            JuMP.optimize!(approx_model)
        end

        """ This is removed, as sometimes the QCP is solved to optimality, but
        Gurobi is not able to get an optimality certificate
        (QCP status(13): Unable to satisfy optimality tolerances;
        a sub-optimal solution is available), even though other solvers prove
        that the solution is indeed optimal.
        Even if the solution is suboptimal, this should not be a problem, as long
        as it is feasible. Therefore, I think the algorithm should not
        stop here.
        """
        #@assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL

        π_k .= JuMP.value.(π)
        #@infiltrate
        #print(L_k, ", ", t_k, ", ", level)

        ########################################################################
        if L_star > t_k + atol/10.0
            error("Could not solve for Lagrangian duals. LB > UB.")
        end
    end

    ############################################################################
    # CONVERGENCE ANALYSIS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "convergence_check" begin
        if isapprox(L_star, t_k, atol = atol, rtol = rtol)
            # CONVERGENCE ACHIEVED
            if isapprox(L_star, s * primal_obj, atol = atol, rtol = rtol)
                # CONVERGENCE TO TRUE OPTIMUM (APPROXIMATELY)
                lag_status = :opt
            else
                # CONVERGENCE TO A SMALLER VALUE THAN THE PRIMAL OBJECTIVE
                # sometimes this occurs due to numerical issues
                # still leads to a valid cut
                lag_status = :conv
            end
        elseif all(h_k .== 0)
            # NO OPTIMALITY ACHIEVED, BUT STILL ALL SUBGRADIENTS ARE ZERO
            # may occur due to numerical issues
            lag_status = :sub
        elseif iter == iteration_limit
            # TERMINATION DUE TO ITERATION LIMIT
            # stil leads to a valid cut
            lag_status = :iter
        end
    end

    ############################################################################
    # APPLY MAGNANTI AND WONG APPROACH IF INTENDED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "magnanti_wong" begin
        magnanti_wong!(node, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_star,
            iteration_limit, atol, rtol, algo_params.duality_regime.dual_choice_regime, iter)
    end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, algo_params.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

    ############################################################################
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k = π_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


"""
Solve lagrangian relaxation to obtain intercept for strengthened Benders cuts
"""
function _getStrengtheningInformation(
    node::SDDP.Node,
    π_k::Vector{Float64},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)
    # The original value for x (former old_rhs)
    x_in_value = zeros(number_of_states)
    # The current estimate for π (in our case determined in initialization)
    # π_k
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states)

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    relax_copy_constraints(node, x_in_value, h_expr, algo_params.state_approximation_regime)

    ########################################################################
    # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
    ########################################################################
    # Evaluate the inner problem and determine a subgradient
    Lag_obj = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, false)
    @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    restore_copy_constraints!(node, x_in_value, algo_params.state_approximation_regime)

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

    return lag_obj
end
