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

    #Infiltrator.@infiltrate
    # Optimization
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_Lag_inner" begin
        JuMP.optimize!(model)
    end
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
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
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
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
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
    atol = cut_generation_regime.duality_regime.atol
    rtol = cut_generation_regime.duality_regime.rtol
    iteration_limit = cut_generation_regime.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    # Augmented Lagrangian dual using 1-norm?
    #---------------------------------------------------------------------------
    augmented = cut_generation_regime.duality_regime.augmented

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.copy_regime)
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
        approx_model = JuMP.Model()
        JuMP.set_optimizer(approx_model, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4,"TimeLimit"=>300,"NumericFocus"=>algo_params.numerical_focus
        ))
        JuMP.set_silent(approx_model)

        approx_model.ext[:sddp_policy_graph] = node.subproblem.ext[:sddp_policy_graph]
        set_solver!(approx_model, algo_params, applied_solvers, :kelley, algo_params.solver_approach)

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        set_objective_bound!(approx_model, s, bound_results.obj_bound)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Norm Minimization part later
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        set_multiplier_bounds!(node, approx_model, number_of_states, bound_results.dual_bound,
            algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
            cut_generation_regime.duality_regime)
    end

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none
    feas_flag = false

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
            if !augmented
                L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
            else
                L_k = _augmented_Lagrangian_relaxation!(node, node_index, π_k, h_expr, h_k, algo_params.regularization_regime, true)
            end
        end
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

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

            # Try recovering from numerical issues
            if (JuMP.termination_status(approx_model) != MOI.OPTIMAL)
                #elude_numerical_issues!(approx_model, algo_params)
                feas_flag = true
                break
            end
        end
        # @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        π_k .= JuMP.value.(π)
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        if L_star > t_k + atol/10.0
            #error("Could not solve for Lagrangian duals. LB > UB.")
            break
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
        elseif L_star > t_k + atol/10.0
            # NUMERICAL ISSUES, LOWER BOUND EXCEEDS UPPER BOUND
            lag_status = :bound_issues
        elseif feas_flag
            # NUMERICAL ISSUES, SOME SUBPROBLEM BECAME INFEASIBLE (OR UNBOUNDED)
            lag_status = :feas_issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    if lag_status == :opt || lag_status == :unbounded
    # In other cases we do not have an optimal solution from Kelley's method,
    # so finding the minimal norm optimal solution does not make sense.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "minimal_norm" begin
            mn_results = minimal_norm_choice!(node, node_index, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_star,
                iteration_limit, atol, rtol, cut_generation_regime.duality_regime.dual_choice_regime, iter, lag_status, augmented, algo_params)

            iter = mn_results.iter
            lag_status = mn_results.lag_status
        end
    # elseif isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
        # println("Proceeding without minimal norm choice.")
    end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, cut_generation_regime.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass, algo_params.solver_approach)

    ############################################################################
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k .= π_star

    #Infiltrator.@infiltrate

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


"""
Given the optimal dual objective value from the Kelley's method, try to find
the optimal dual multipliers with the smallest L1-norm.

Note that this is only done until the maximum number of iterations is
achieved in total.
"""

function minimal_norm_choice!(
    node::SDDP.Node,
    node_index::Int64,
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
    dual_choice_regime::DynamicSDDiP.MinimalNormChoice,
    iter::Int,
    lag_status::Symbol,
    augmented::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
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
    for it in (iter+1):iteration_limit
        JuMP.optimize!(approx_model)

        try
            @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        catch err
            return (iter=it, lag_status=:mn_issue)
        end

        π_k .= JuMP.value.(π)

        if !augmented
            L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        else
            L_k = _augmented_Lagrangian_relaxation!(node, node_index, π_k, h_expr, h_k, algo_params.regularization_regime, true)
        end
        if isapprox(L_star, L_k, atol = atol, rtol = rtol)
            # At this point we tried the smallest ‖π‖ from the cutting plane
            # problem, and it returned the optimal dual objective value. No
            # other optimal dual vector can have a smaller norm.
            π_star .= π_k
            return (iter=it, lag_status=:mn_opt)
        end
        JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))

    end

    return (iter=iteration_limit, lag_status=:mn_iter)
end


function minimal_norm_choice!(
    node::SDDP.Node,
    node_index::Int64,
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
    lag_status::Symbol,
    augmented::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    )

    return (iter=iter, lag_status=lag_status)
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
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
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
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
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
    atol = cut_generation_regime.duality_regime.atol
    rtol = cut_generation_regime.duality_regime.rtol
    iteration_limit = cut_generation_regime.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    # Set bundle_parameters
    #---------------------------------------------------------------------------
    level_factor = dual_solution_regime.level_factor

    # Augmented Lagrangian dual using 1-norm?
    #---------------------------------------------------------------------------
    augmented = cut_generation_regime.duality_regime.augmented

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.copy_regime)
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
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        approx_model = JuMP.Model()
        JuMP.set_optimizer(approx_model, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4,"TimeLimit"=>300,"NumericFocus"=>algo_params.numerical_focus
        ))
        JuMP.set_silent(approx_model)

        approx_model.ext[:sddp_policy_graph] = node.subproblem.ext[:sddp_policy_graph]
        set_solver!(approx_model, algo_params, applied_solvers, :kelley, algo_params.solver_approach)

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        set_objective_bound!(approx_model, s, bound_results.obj_bound)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Norm Minimization part later
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        set_multiplier_bounds!(node, approx_model, number_of_states, bound_results.dual_bound,
            algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
            cut_generation_regime.duality_regime)

    end

    ############################################################################
    # BUNDLE METHOD
    ############################################################################
    iter = 0
    lag_status = :none
    feas_flag = false
    π_k_dummy = zeros(length(π_k))

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
            if !augmented
                L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
            else
                L_k = _augmented_Lagrangian_relaxation!(node, node_index, π_k, h_expr, h_k, algo_params.regularization_regime, true)
            end
        end
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]
        #Infiltrator.@infiltrate

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
        #set_solver!(approx_model, algo_params, applied_solvers, :kelley, algo_params.solver_approach)

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        TimerOutputs.@timeit DynamicSDDiP_TIMER "outer_sol" begin
            JuMP.optimize!(approx_model)

            # Try recovering from numerical issues
            if (JuMP.termination_status(approx_model) != MOI.OPTIMAL)
                #Infiltrator.@infiltrate
                #elude_numerical_issues!(approx_model, algo_params)
                feas_flag = true
                break
            end
        end
        #@assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        π_k_dummy .= JuMP.value.(π)
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

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
        #set_solver!(approx_model, algo_params, applied_solvers, :level_bundle, algo_params.solver_approach)
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
        if (JuMP.termination_status(approx_model) != JuMP.MOI.OPTIMAL && JuMP.termination_status(approx_model) != JuMP.MOI.LOCALLY_SOLVED)
            if dual_solution_regime.switch_to_kelley
                # in case of an error, we proceed with the Kelley step, by using the multipliers obtained from the inner problem
                π_k .= π_k_dummy
            else
                # in case of an error, we leave the bundle method and use the current multipliers to construct a cut
                feas_flag = true
                break
            end
        else
            π_k .= JuMP.value.(π)
        end

        π_k .= JuMP.value.(π)

        # Delete the level lower bound for the original approx_model again
        JuMP.delete_lower_bound(t)

        #Infiltrator.@infiltrate
        #print(L_k, ", ", t_k, ", ", level)

        ########################################################################
        if L_star > t_k + atol/10.0
            #error("Could not solve for Lagrangian duals. LB > UB.")
            break
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
        elseif L_star > t_k + atol/10.0
            # NUMERICAL ISSUES, LOWER BOUND EXCEEDS UPPER BOUND
            lag_status = :bound_issues
        elseif feas_flag
            # NUMERICAL ISSUES, SOME SUBPROBLEM BECAME INFEASIBLE (OR UNBOUNDED)
            lag_status = :feas_issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    if lag_status == :opt || lag_status == :unbounded
    # In other cases we do not have an optimal solution from Kelley's method,
    # so finding the minimal norm optimal solution does not make sense.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "minimal_norm" begin
            mn_results = minimal_norm_choice!(node, node_index, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_star,
                iteration_limit, atol, rtol, cut_generation_regime.duality_regime.dual_choice_regime, iter, lag_status, augmented, algo_params)

            iter = mn_results.iter
            lag_status = mn_results.lag_status
        end
    # elseif isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
        # println("Proceeding without minimal norm choice.")
    end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, cut_generation_regime.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass, algo_params.solver_approach)

    ############################################################################
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k .= π_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


"""
Solve lagrangian relaxation to obtain intercept for strengthened Benders cuts
"""
function _getStrengtheningInformation(
    node::SDDP.Node,
    π_k::Vector{Float64},
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
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
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    relax_copy_constraints!(node, x_in_value, h_expr, cut_generation_regime.state_approximation_regime, DynamicSDDiP.StateSpaceCopy())

    ########################################################################
    # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
    ########################################################################
    # Evaluate the inner problem and determine a subgradient
    lag_obj = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, false)
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    restore_copy_constraints!(node, x_in_value, cut_generation_regime.state_approximation_regime)

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass, algo_params.solver_approach)

    return lag_obj
end


"""
Subgradient method to solve Lagrangian dual
"""
function solve_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    primal_obj::Float64,
    π_k::Vector{Float64},
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.Subgradient
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
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
    atol = cut_generation_regime.duality_regime.atol
    rtol = cut_generation_regime.duality_regime.rtol
    iteration_limit = cut_generation_regime.duality_regime.iteration_limit

    # Set subgradient parameters
    #---------------------------------------------------------------------------
    # As in SDDiP.jl we check whether the lower bound has improved in an iteration.
    # If it hasn't, we reduce the step-size parameter γ. If it hasn't for a
    # predefined number of times, the algorithm stops due to bounds stalling.
    times_unchanged = 0
    wait = cut_generation_regime.duality_regime.dual_solution_regime.wait
    max_times_unchanged = cut_generation_regime.duality_regime.dual_solution_regime.max_times_unchanged
    gamma_step = cut_generation_regime.duality_regime.dual_solution_regime.gamma

    # Cache for lower bound for these checks
    cached_bound = -Inf

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    # Augmented Lagrangian dual using 1-norm?
    #---------------------------------------------------------------------------
    # NOTE: We only allow the augmented Lagrangian for Kelley's method and the Bundle method so far.
    @assert cut_generation_regime.duality_regime.augmented != true

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.copy_regime)
    end
    node.ext[:backward_data][:old_rhs] = x_in_value

    ############################################################################
    # SET-UP THE PROJECTION MODEL
    ############################################################################
    """ Note that if we use bounds on the dual multipliers (e.g. when a regularization
    is applied or if we have sign conditions on the multipliers), we have to make
    sure that the new incumbent reached by the subgradient step is feasible.
    If it isn't, we have to project it to the feasible set Π.
    To this end, we solve a projection problem, which finds the point in Π
    that is closest to the computed point π_k.

    Even if this problem is quadratic, it should be really easy to solve in
    most cases, and even trivially if we do not consider any bounds on π.

    Note that the objective is created in each iteration.
    """

    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_proj_model" begin
        # Optimizer is re-set anyway
        proj_model = JuMP.Model()
        proj_model.ext[:sddp_policy_graph] = node.subproblem.ext[:sddp_policy_graph]

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Norm Minimization part later
        JuMP.@variable(proj_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(proj_model, π⁻[1:number_of_states] >= 0)
        JuMP.@expression(proj_model, π, π⁺ .- π⁻) # not required to be a constraint
        set_multiplier_bounds!(node, proj_model, number_of_states, bound_results.dual_bound,
            algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
            cut_generation_regime.duality_regime)

        # Set solver
        JuMP.set_optimizer(proj_model, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4,"TimeLimit"=>300,"NumericFocus"=>algo_params.numerical_focus
        ))
        JuMP.set_silent(proj_model)
        #set_solver!(proj_model, algo_params, applied_solvers, :level_bundle, algo_params.solver_approach)
    end

    ############################################################################
    # SUBGRADIENT METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # In this case, we initialize t_k directly as the primal objective
    t_k = s * primal_obj
    # TODO: Maybe change this later.

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol) && times_unchanged <= max_times_unchanged
        iter += 1

        ########################################################################
        # CHECK FOR TIMES UNCHANGED
        ########################################################################
        # We check if the lower bound has improved every wait-th iteration
        if mod(iter, wait) == 0
            if cached_bound < L_star
                cached_bound = L_star
                times_unchanged = 0
            else
                gamma_step = gamma_step / 2
                times_unchanged += 1
            end
        end

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        end
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
        end

        ########################################################################
        # GET A NEW INCUMBENT
        ########################################################################
        # t_k = JuMP.objective_value(approx_model)

        # Calculate the new step-size
        if sum(h_k.^2) == 0
            # If the subgradients are zero already, then we can set the step to
            # 1 as we will not move anyway. Otherwise, in the below formula
            # we would divide by zero.
            step = 1
        else
            step = gamma_step * (t_k - L_star) / sum(h_k.^2)
        end

        # Update the multipliers by doing a subgradient step
        π_k .+= step * h_k

        # Create the objective
        JuMP.@objective(proj_model, Min, sum((π_k[i] - π[i])^2 for i in 1:number_of_states))

        # Solve the projection problem
        TimerOutputs.@timeit DynamicSDDiP_TIMER "proj_sol" begin
            JuMP.optimize!(proj_model)
        end

        @assert JuMP.termination_status(proj_model) == JuMP.MOI.OPTIMAL

        # Get the new incumbent as the projection of the original point
        π_k .= JuMP.value.(π)

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        if L_star > t_k + atol/10.0
            #error("Could not solve for Lagrangian duals. LB > UB.")
            break
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
                # NOTE: In the current setting this can't happen, as we set
                # t_k = primal_obj initially and never change it
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
        elseif L_star > t_k + atol/10.0
            # NUMERICAL ISSUES, LOWER BOUND EXCEEDS UPPER BOUND
            lag_status = :issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    if isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
        lag_status = :mn_issue
        #println("Proceeding without minimal norm choice.")
    end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, cut_generation_regime.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass, algo_params.solver_approach)

    # Set dual_vars (here π_k) to the optimal solution
    π_k .= π_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end
