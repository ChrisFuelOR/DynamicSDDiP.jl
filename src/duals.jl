# The function
# > "get_dual_solution",
# is derived from the same function in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################


"""
Solving the dual problem to obtain cut information - using LP relaxation
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    solver_obj::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.LinearDuality,
    )

end

"""
Solving the dual problem to obtain cut information - using LP relaxation
and strengthening by Lagrangian relaxation
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    solver_obj::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.StrengthenedDuality,
    )

end

"""
Solving the dual problem to obtain cut information - using Lagrangian dual
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    solver_obj::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.LagrangianDuality,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    subproblem = node.subproblem

    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)

    # storages for information on Lagrangian dual
    lag_obj = 0
    lag_iterations = 0
    lag_status = :none

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "dual_initialization" begin
        dual_vars = initialize_duals(node, subproblem, algo_params, applied_solvers, duality_regime.dual_initialization_regime)
    end

    ############################################################################
    # GET PRIMAL SOLUTION TO BOUND LAGRANGIAN DUAL, TRACK CONVERGENCE AND DEBUGGING
    ############################################################################
    # REGULARIZE PROBLEM IF REGULARIZATION IS USED
    node.ext[:regularization_data] = Dict{Symbol,Any}()
    regularize_binary!(node, node_index, subproblem, algo_params.regularization_regime)

    # SOLVE PRIMAL PROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    primal_obj = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    # ADAPT TOTAL SOLVES APPROPRIATELY
    if haskey(model.ext, :total_solves)
        model.ext[:total_solves] += 1
    else
        model.ext[:total_solves] = 1
    end

    # DEREGULARIZE PROBLEM IF REQUIRED
    deregularize_binary!(node, subproblem, algo_params.regularization_regime)

    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # GET BOUNDS FOR LAGRANGIAN DUAL
    ############################################################################
    bound_results = get_dual_bounds(node, algo_params, primal_obj, duality_regime.dual_bound_regime)

    @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    try
        ########################################################################
        # CALL SOLUTION METHOD
        ########################################################################
        # Solve dual and return a dict with the multiplier of the copy constraints.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
            results = solve_lagrangian_dual(
                    node,
                    node_index,
                    primal_obj,
                    dual_vars,
                    bound_results,
                    integrality_handler,
                    algo_params,
                    applied_solvers,
                    duality_regime.dual_solution_regime
                    )
        end

        lag_obj = results.lag_obj
        lag_iterations = results.iterations
        lag_status = results.lag_status

        ########################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ########################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(lag_status, duality_regime.dual_status_regime)

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    catch e
        SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    store_dual_values!(node, dual_vars, binary_state, integrality_handler, algo_params.state_approximation_regime)

    return (
        dual_values=dual_values,
        bin_state=bin_state,
        intercept=lag_obj,
        iterations=lag_iterations,
        lag_status=lag_status,
    )


"""
Determining objective and/or variable bounds for the Lagrangian dual
if ValueBound is used.

Note that we always solve the primal problem, even if we do not use its
objective value as the objective bound, as this is more convenient for
debugging purposes.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.ValueBound,
    )

    return (
        obj_bound = primal_obj,
        dual_bound = Inf
    )

end

"""
Determining objective and/or variable bounds for the Lagrangian dual
if NormBound is used.

Note that we always solve the primal problem, even if we do not use its
objective value as the objective bound, as this is more convenient for
debugging purposes.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.NormBound,
    )

    return (
        obj_bound = Inf,
        dual_bound = get_norm_bound(node, node_index, algo_params)
    )

end

"""
Determining objective and/or variable bounds for the Lagrangian dual
if BothBound is used.

Note that we always solve the primal problem, even if we do not use its
objective value as the objective bound, as this is more convenient for
debugging purposes.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.BothBounds,
    )

    return (
        obj_bound = primal_obj,
        dual_bound = get_norm_bound(node, node_index, algo_params)
    )

end

"""
Determining the norm bound to be used for the Lagrangian dual.

Actually, we attempt to calculate a norm of B (coefficient matrix of the
binary expansion), e.g. the column sum norm. However, we can also use an
overestimator by taking the maximum upper bound of all state variables.
"""
function get_norm_bound(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    )

    B_norm_bound = 0
    for (name, state_comp) in node.states
        if state_comp.info.in.upper_bound > B_norm_bound
            B_norm_bound = state_comp.info.in.upper_bound
        end
    end
    dual_bound = algo_params.regularization_regime.sigma[node_index] * B_norm_bound

    return dual_bound
end


"""
Checking the status of the Lagrangian dual solution and throw an error if required
under rigorous regime.
"""
function lagrangian_status_check(
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Rigorous,
    )

    if lag_status == :conv
        error("Lagrangian dual converged to value < solver_obj.")
    elseif lag_status == :sub
        error("Lagrangian dual had subgradients zero without LB=UB.")
    elseif lag_status == :iter
        error("Solving Lagrangian dual exceeded iteration limit.")
    end

    return
end


"""
Trivial check of the status of the Lagrangian dual solution under lax regime.
"""
function lagrangian_status_check(
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Lax,
    )

    return
end


"""
Initializing duals with zero vector.
"""
function initialize_duals(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_initalization_regime::DynamicSDDiP.ZeroDuals,
)

    # Get number of states and create zero vector for duals
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)
    dual_vars_initial = zeros(number_of_states)

    return

end


"""
Initializing duals by solving LP relaxation.
"""
function initialize_duals(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_initalization_regime::DynamicSDDiP.LPDuals,
)

    # Get number of states and create zero vector for duals
    number_of_states = get_number_of_states(node, algo_params.state_approximation_regime)
    dual_vars_initial = zeros(number_of_states)

    # Create LP Relaxation
    undo_relax = JuMP.relax_integrality(subproblem);

    # Define appropriate solver
    set_solver(subproblem, algo_params, applied_solvers, :LP_relax)

    # Solve LP Relaxation
    JuMP.optimize!(subproblem)

    # Get dual values (reduced costs) for binary states as initial solution #TODO
    get_and_set_dual_values!(node, dual_vars_initial, algo_params.state_approximation_regime)

    # Undo relaxation
    undo_relax()

    return

end

function get_and_set_dual_values!(node::SDDP.Node, dual_vars_initial::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation)

    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
       variable_name = node.ext[:backward_data][:bin_states][name]
       reference_to_constr = FixRef(variable_name)
       dual_vars_initial[i] = JuMP.getdual(reference_to_constr)
    end

    return
end

function get_and_set_dual_values!(node::SDDP.Node, dual_vars_initial::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.NoState)

    for (i, name) in enumerate(keys(node,states))
        reference_to_constr = FixRef(name)
        dual_vars_initial[i] = JuMP.getdual(reference_to_constr)
    end

    return
end

function store_dual_values!(node::SDDP.Node, dual_vars::Vector{Float64}, bin_state::Dict{Symbol, BinaryState},
    integrality_handler::SDDP.SDDiP, state_approximation_regime::DynamicSDDiP.BinaryApproximation)

    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
        dual_values[name] = dual_vars[i]

        value = integrality_handler.old_rhs[i]
        x_name = node.ext[:backward_data][:bin_x_names][name]
        k = node.ext[:backward_data][:bin_k][name]
        bin_state[name] = BinaryState(value, x_name, k)
    end

    return
end

function store_dual_values!(node::SDDP.Node, dual_vars::Vector{Float64}, bin_state::Dict{Symbol, BinaryState},
    integrality_handler::SDDP.SDDiP, state_approximation_regime::DynamicSDDiP.NoStateApproximation)

    for (i, name) in enumerate(keys(node.states))
        dual_values[name] = dual_vars[i]
    end

    return
end
