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
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    h_k::Vector{Float64},
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Set the Lagrangian relaxation of the objective in the primal model
    old_obj = JuMP.objective_function(model)
    JuMP.set_objective_function(model, @expression(old_obj - π_k' * h_expr))

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
    bound_results::Tuple{Float64,Float64}
    algo_params::NCNBD.AlgoParams,
    applied_solvers::NCNBD.AppliedSolvers,
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
    h_expr = Vector{AffExpr}(undef, number_of_states)
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = algo_params.duality_regime.atol
    rtol = algo_params.duality_regime.rtol
    iteration_limit = algo_params.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver(node.subproblem, algo_params, applied_solvers, :lagrange_relax)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    relax_copy_constraints(node, x_in_value, h_expr, algo_params.state_approximation_regime)

    ############################################################################
    # LOGGING OF LAGRANGIAN DUAL
    ############################################################################
    #lag_log_file_handle = open("C:/Users/cg4102/Documents/julia_logs/Lagrange.log", "a")
    #print_helper(print_lagrange_header, lag_log_file_handle)

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    # Approximation of Lagrangian dual by cutting planes
    # Optimizer is re-set anyway
    approx_model = JuMP.Model(GLPK.Optimizer)
    set_solver(approx_model, algo_params, applied_solvers, :kelley)

    # Create the objective
    # Note that it is always formulated as a maximization problem, but that
    # s modifies the sense appropriately
    @variable(approx_model, t)
    set_objective_bound(approx_model, s, bound_results.obj_bound)
    @objective(approx_model, Max, t)

    # Create the dual variables
    # Note that the real dual multipliers are split up into two non-negative
    # variables here, which is required for the Magnanti Wong part later
    @variable(approx_model, π⁺[1:number_of_states] >= 0)
    @variable(approx_model, π⁻[1:number_of_states] >= 0)
    @expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
    set_multiplier_bounds(approx_model, number_of_states bound_results.dual_bound)

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = -Inf

    while iter <= iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # ADD CUTTING PLANE
        ########################################################################
        JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
        end

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        #print("UB: ", f_approx, ", LB: ", f_actual)
        #println()

        ########################################################################
        # PREPARE NEXT ITERATION
        ########################################################################
        # dual_vars .= value.(x)
        # can be deleted with the next update of GAMS.jl
        # replace!(dual_vars, NaN => 0)

        if L_star > t_k + atol/10.0
            error("Could not solve for Lagrangian duals. LB > UB.")
        end
    end

    ############################################################################
    # CONVERGENCE ANALYSIS
    ############################################################################
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
    elseif all(h_k .== 0)
        # NO OPTIMALITY ACHIEVED, BUT STILL ALL SUBGRADIENTS ARE ZERO
        # may occur due to numerical issues
        lag_status = :sub
    elseif iter == iteration_limit
        # TERMINATION DUE TO ITERATION LIMIT
        # stil leads to a valid cut
        lag_status = :iter
    end

    ############################################################################
    # APPLY MAGNANTI AND WONG APPROACH IF INTENDED
    ############################################################################
    magnanti_wong(node, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_k, L_star,
        iteration_limit, atol, rtol, algo_params.dual_choice_regime)

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    restore_copy_constraints(node, x_in_value, algo_params.state_approximation_regime)

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver(node.subproblem, algo_params, applied_solvers, :forward_pass)

    ############################################################################
    # LOGGING
    ############################################################################
    print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k = π_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(bin_state)
        # Store expression for slack
        h_expr[i] = @expression(node.subproblem, bin_state - x_in_value[i])
        # Relax copy constraint (i.e. z does not have to take the value of ̄x anymore)
        JuMP.unfix(bin_state)

        # Set bounds to ensure that inner problems are feasible
        # As we use binary approximation, 0 and 1 can be used
        JuMP.set_lower_bound(bin_state, 0)
        JuMP.set_upper_bound(bin_state, 1)

    end

    return
end

function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    for (i, (_, state)) in enumerate(node.states)
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(state.in)
        # Store expression for slack
        h_expr[i] = @expression(node.subproblem, state.in - x_in_value[i])
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

function restore_copy_constraints(
    node::SDDP.Node,
    x_in_value::Vector{Float64}
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    for (i, (_, bin_state)) in enumerate(node.ext[:backward_data][:bin_states])
        #prepare_state_fixing!(node, state_comp)
        JuMP.fix(bin_state, x_in_value[i], force = true)
    end
end

function restore_copy_constraints(
    node::SDDP.Node,
    x_in_value::Vector{Float64}
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    for (i, (_, state)) in enumerate(node.states)
        #prepare_state_fixing!(node, state_comp)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end
end

function set_objective_bound(approx_model::JuMP.Model, s::Int, obj_bound::Float64)

    JuMP.set_upper_bound(approx_model[:t], s * obj_bound)

end

function set_multiplier_bounds(approx_model::JuMP.Model, number_of_states::Int, dual_bound::Float64)

    for i in 1:number_of_states
        JuMP.set_upper_bound(π⁺[i], dual_bound)
        JuMP.set_upper_bound(π⁻[i], dual_bound)
    end
end

"""
Given the optimal dual objective value from the Kelley's method, try to find
the optimal dual multipliers with the smallest L1-norm.

Note that this is only done until the maximum number of iterations is
achieved in total.
"""

function magnanti_wong(
    node::SDDP.Node,
    approx_model::JuMP.model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    t_k::Float64,
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    h_k::Vector{Float64},
    s::Int,
    L_k::Float64,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.MagnantiWongChoice,
    )

    # Reset objective
    @objective(approx_model, Min, sum(π⁺) + sum(π⁻))
    JuMP.set_lower_bound(t, t_k)

    # The worst-case scenario in this for-loop is that we run through the
    # iterations without finding a new dual solution. However if that happens
    # we can just keep our current λ_star.
    for _ in (iter+1):iteration_limit
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        π_k .= value.(π)
        L_k = _solve_Lagrangian_relaxation(node.subproblem, π_k, h_expr, h_k)
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


function magnanti_wong(
    node::SDDP.Node,
    approx_model::JuMP.model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    t_k::Float64,
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    h_k::Vector{Float64},
    s::Int,
    L_k::Float64,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.StandardChoice,
    )

    return
end



"""
Level bundle method to solve the Lagrangian duals.
"""
function _bundle_level(
    node::SDDP.Node,
    node_index::Int64,
    obj::Float64,
    dual_vars::Vector{Float64},
    integrality_handler::SDDP.SDDiP,
    algo_params::NCNBD.AlgoParams,
    applied_solvers::NCNBD.AppliedSolvers,
    dual_bound::Union{Float64,Nothing}
    )

    # INITIALIZATION
    ############################################################################
    atol = integrality_handler.atol # corresponds to deltabar
    rtol = integrality_handler.rtol # corresponds to deltabar
    model = node.ext[:linSubproblem]
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal duals

    # initialize bundle parameters
    level_factor = algo_params.level_factor

    for (i, (name, bin_state)) in enumerate(node.ext[:backward_data][:bin_states])
        integrality_handler.old_rhs[i] = JuMP.fix_value(bin_state)
        integrality_handler.slacks[i] = bin_state - integrality_handler.old_rhs[i]
        JuMP.unfix(bin_state)
        #JuMP.unset_binary(state_comp.in) # TODO: maybe not required
        JuMP.set_lower_bound(bin_state, 0)
        JuMP.set_upper_bound(bin_state, 1)
    end

    # LOGGING OF LAGRANGIAN DUAL
    ############################################################################
    lag_log_file_handle = open("C:/Users/cg4102/Documents/julia_logs/Lagrange.log", "a")
    print_helper(print_lagrange_header, lag_log_file_handle)

    # SET-UP APPROXIMATION MODEL
    ############################################################################
    # Subgradient at current solution
    subgradients = integrality_handler.subgradients
    # Best multipliers found so far
    best_mult = integrality_handler.best_mult
    # Dual problem has the opposite sense to the primal
    dualsense = (
        JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? JuMP.MOI.MAX_SENSE :
            JuMP.MOI.MIN_SENSE
    )

    # Approximation of Lagrangian dual as a function of the multipliers
    approx_model = JuMP.Model(Gurobi.Optimizer)
    # even if objective is quadratic, it should be possible to use Gurobi
    if applied_solvers.Lagrange == "CPLEX"
        set_optimizer(approx_model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0, "numericalemphasis"=>0))
        set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0, "numericalemphasis"=>0))
    elseif applied_solvers.Lagrange == "Gurobi"
        set_optimizer(approx_model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0, "NumericFocus"=>1))
        set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0, "numericalemphasis"=>0))
    else
        set_optimizer(approx_model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0))
        set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.Lagrange, "optcr"=>0.0, "numericalemphasis"=>0))
    end

    # Define Lagrangian dual multipliers
    @variables approx_model begin
        θ
        x[1:length(dual_vars)]
    end

    if dualsense == MOI.MIN_SENSE
        JuMP.set_lower_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (Inf, Inf, -Inf)
    else
        #JuMP.set_upper_bound(θ, 10000.0)

        JuMP.set_upper_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (-Inf, -Inf, Inf)
    end

    # BOUND DUAL VARIABLES IF INTENDED
    ############################################################################
    if !isnothing(dual_bound)
        for i in 1:length(dual_vars)
            JuMP.set_lower_bound(x[i], -dual_bound)
            JuMP.set_upper_bound(x[i], dual_bound)
        end
    end

    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none
    while iter < integrality_handler.iteration_limit
        iter += 1

        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the real function and determine a subgradient
        f_actual = _solve_Lagrangian_relaxation!(subgradients, node, dual_vars, integrality_handler.slacks, :yes)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange] #|| model.ext[:sddp_policy_graph].ext[:iteration] == 12

        # ADD CUTTING PLANE TO APPROX_MODEL
        ########################################################################
        # Update the model and update best function value so far
        if dualsense == MOI.MIN_SENSE
            JuMP.@constraint(
                approx_model,
                θ >= f_actual + LinearAlgebra.dot(subgradients, x - dual_vars)
                # TODO: Reset upper bound to inf?
            )
            if f_actual <= best_actual
                best_actual = f_actual
                best_mult .= dual_vars
            end
        else
            JuMP.@constraint(
                approx_model,
                θ <= f_actual + LinearAlgebra.dot(subgradients, x - dual_vars)
                # TODO: Reset lower boumd to -inf?
            )
            if f_actual >= best_actual
                # bestmult is not simply getvalue.(x), since approx_model may just haven gotten lucky
                # same for best_actual
                best_actual = f_actual
                best_mult .= dual_vars
            end
        end

        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Define objective for approx_model
        JuMP.@objective(approx_model, dualsense, θ)

        # Get an upper bound from the approximate model
        # (we could actually also use obj here)
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        f_approx = JuMP.objective_value(approx_model)

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange] #|| model.ext[:sddp_policy_graph].ext[:iteration] == 12

        # Construct the gap (not directly used for termination, though)
        #gap = abs(best_actual - f_approx)
        gap = abs(best_actual - obj)

        print("UB: ", f_approx, ", LB: ", f_actual, best_actual)
        println()

        # CONVERGENCE CHECKS AND UPDATE
        ########################################################################
        # convergence achieved
        if isapprox(best_actual, f_approx, atol = atol, rtol = rtol)
            # convergence to obj -> tight cut
            if isapprox(best_actual, obj, atol = atol, rtol = rtol)
                lag_status = :aopt
            # convergence to a smaller value than obj
            # maybe possible due to numerical issues
            # -> valid cut
            else
                lag_status = :conv
            end

        # zero subgradients (and no further improvement), despite no convergence
        # maybe possible due to numerical issues
        # -> valid cut
        elseif all(subgradients.== 0)
            lag_status = :sub

        # lb exceeds ub: no convergence
        elseif best_actual > f_approx + atol/10.0
            error("Could not solve for Lagrangian duals. LB > UB.")
        end

        # return
        if lag_status == :sub || lag_status == :aopt || lag_status == :conv
            dual_vars .= best_mult
            if dualsense == JuMP.MOI.MIN_SENSE
                dual_vars .*= -1
            end

            for (i, (name, bin_state)) in enumerate(node.ext[:backward_data][:bin_states])
                #prepare_state_fixing!(node, state_comp)
                JuMP.fix(bin_state, integrality_handler.old_rhs[i], force = true)
            end

            if applied_solvers.MILP == "CPLEX"
                set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.MILP, "optcr"=>0.0, "numericalemphasis"=>0))
            elseif applied_solvers.MILP == "Gurobi"
                set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.MILP, "optcr"=>0.0, "NumericFocus"=>1))
            else
                set_optimizer(model, optimizer_with_attributes(GAMS.Optimizer, "Solver"=>applied_solvers.MILP, "optcr"=>0.0))
            end

            return (lag_obj = best_actual, iterations = iter, lag_status = lag_status)
        end

        # FORM A NEW LEVEL
        ########################################################################
        if dualsense == :Min
            level = f_approx + gap * level_factor
            #TODO: + atol/10.0 for numerical issues?
            JuMP.setupperbound(θ, level)
        else
            level = f_approx - gap * level_factor
            #TODO: - atol/10.0 for numerical issues?
            JuMP.setlowerbound(θ, level)
        end

        # DETERMINE NEXT ITERATE USING PROXIMAL PROBLEM
        ########################################################################
        # Objective function of approx model has to be adapted to new center
        JuMP.@objective(approx_model, Min, sum((dual_vars[i] - x[i])^2 for i=1:length(dual_vars)))
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL

        # Next iterate
        dual_vars .= value.(x)
        # can be deleted with the next update of GAMS.jl
        replace!(dual_vars, NaN => 0)

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange] #|| model.ext[:sddp_policy_graph].ext[:iteration] == 12

        # Logging
        print_helper(print_lag_iteration, lag_log_file_handle, iter, f_approx, best_actual, f_actual)

    end

    lag_status = :iter
    #error("Could not solve for Lagrangian duals. Iteration limit exceeded.")
    return (lag_obj = best_actual, iterations = iter, lag_status = lag_status)

end
