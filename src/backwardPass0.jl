# The functions
# > "solve_all_children",
# > "solve_subproblem_backward",
# > "get_dual_variables_backward",
# > "calculate_bound",
# > "solve_first_stage_problem"
# are derived from similar named functions (backward_pass,
# solve_all_children, solve_subproblem, get_dual_variables, calculate_bound,
# solve_first_stage_problem) in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Executing the backward pass of a loop of DynamicSDDiP.
"""
function backward_pass(
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    scenario_path::Vector{Tuple{T,NoiseType}},
    sampled_states::Vector{Dict{Symbol,Float64}},
    objective_states::Vector{NTuple{N,Float64}},
    belief_states::Vector{Tuple{Int,Dict{T,Float64}}}) where {T,NoiseType,N}

    ####################################################################
    # INITIALIZATION
    ####################################################################

    # storage for cuts
    cuts = Dict{T,Vector{Any}}(index => Any[] for index in keys(model.nodes))

    # storage for data on solving Lagrangian dual
    model.ext[:lag_iterations] = Int[]
    model.ext[:lag_status] = Symbol[]

    ############################################################################
    # Traverse backwards through the stages
    ############################################################################
    for index in length(scenario_path):-1:1
        outgoing_state = sampled_states[index]
        items = BackwardPassItems(T, SDDP.Noise)

        node_index, _ = scenario_path[index]
        node = model[node_index]
        if length(node.children) == 0
            continue
        end

        # Dict to store values of binary approximation of the state
        # Note that we could also retrieve this from the actual trial point
        # (outgoing_state) or from its approximation via binexpand. However,
        # this collection is not only important to obtain the correct values,
        # but also to store them together with the symbol/name of the variable.
        node.ext[:binary_state_values] = Dict{Symbol, Float64}()

        ####################################################################
        # SOLVE ALL CHILDREN PROBLEMS
        ####################################################################
        solve_all_children(
            model,
            node,
            node_index,
            items,
            1.0,
            belief_state,
            objective_state,
            outgoing_state,
            algo_params.backward_sampling_scheme,
            scenario_path[1:index],
            algo_params,
            applied_solvers
        )

        # RECONSTRUCT ANCHOR POINTS IN BACKWARD PASS
        ####################################################################
        anchor_points = Dict{Symbol,Float64}()
        for (name, value) in outgoing_state
            state_comp = node.states[name]
            epsilon = algo_params.binary_precision[name]
            (approx_state_value, )  = determine_anchor_states(state_comp, value, epsilon)
            anchor_points[name] = approx_state_value
        end

        @infiltrate algo_params.infiltrate_state in [:all]

        # REFINE BELLMAN FUNCTION BY ADDING CUTS
        ####################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "update_bellman" begin
            new_cuts = refine_bellman_function(
                model,
                node,
                node_index,
                node.bellman_function,
                options.risk_measures[node_index],
                outgoing_state,
                anchor_points,
                items.bin_state,
                items.duals,
                items.supports,
                items.probability,
                items.objectives,
                algo_params,
                applied_solvers
            )
        end
        push!(cuts[node_index], new_cuts)
        #NOTE: Has to be adapted for stochastic case
        push!(model.ext[:lag_iterations], sum(items.lag_iterations))
        push!(model.ext[:lag_status], items.lag_status[1])

        #TODO: Implement cut-sharing as in SDDP

    end
    return cuts
end


"""
Solving all children within the backward pass.
"""
function solve_all_children(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    items::BackwardPassItems,
    belief::Float64,
    belief_state,
    objective_state,
    outgoing_state::Dict{Symbol,Float64},
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    scenario_path,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers
) where {T}
    length_scenario_path = length(scenario_path)
    for child in node.children
        if isapprox(child.probability, 0.0, atol = 1e-6)
            continue
        end
        child_node = model[child.term]
        for noise in SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node)
            if length(scenario_path) == length_scenario_path
                push!(scenario_path, (child.term, noise.term))
            else
                scenario_path[end] = (child.term, noise.term)
            end
            ####################################################################
            # IF SOLUTIONS FOR THIS NODE ARE CACHED ALREADY, USE THEM
            ####################################################################
            if haskey(items.cached_solutions, (child.term, noise.term))
                sol_index = items.cached_solutions[(child.term, noise.term)]
                push!(items.duals, items.duals[sol_index])
                push!(items.supports, items.supports[sol_index])
                push!(items.nodes, child_node.index)
                push!(items.probability, items.probability[sol_index])
                push!(items.objectives, items.objectives[sol_index])
                push!(items.belief, belief)
                push!(items.bin_state, items.bin_state[sol_index])
                push!(items.lag_iterations, items.lag_iterations[sol_index])
                push!(items.lag_status, items.lag_status[sol_index])
            else
                ################################################################
                # SOLVE THE BACKWARD PASS PROBLEM
                ################################################################
                TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_BP" begin
                    subproblem_results = solve_subproblem_backward(
                        model,
                        child_node,
                        node_index+1,
                        outgoing_state,
                        noise.term,
                        scenario_path,
                        algo_params,
                        applied_solvers
                    )
                end
                push!(items.duals, subproblem_results.duals)
                push!(items.supports, noise)
                push!(items.nodes, child_node.index)
                push!(items.probability, child.probability * noise.probability)
                push!(items.objectives, subproblem_results.objective)
                push!(items.belief, belief)
                push!(items.bin_state, subproblem_results.bin_state)
                push!(items.lag_iterations, subproblem_results.iterations)
                push!(items.lag_status, subproblem_results.lag_status)
                items.cached_solutions[(child.term, noise.term)] = length(items.duals)
            end
        end
    end
    if length(scenario_path) == length_scenario_path
        # No-op. There weren't any children to solve.
    else
        # Drop the last element (i.e., the one we added).
        pop!(scenario_path)
    end
end


"""
Solving the backward pass problem for one specific child
"""
# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise. If require_duals=true, also return the dual variables
# associated with the fixed constraint of the incoming state variables.
function solve_subproblem_backward(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    state::Dict{Symbol,Float64},
    noise,
    scenario_path::Vector{Tuple{T,S}},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers;
) where {T,S}

    ############################################################################
    # MODEL PARAMETRIZATION
    ############################################################################
    subproblem = node.subproblem

    # Storage for backward pass data
    node.ext[:backward_data] = Dict{Symbol,Any}()

    # Parameterize the model. Fix the value of the incoming state variables.
    # Then parameterize the model depending on `noise` and set the objective.
    set_incoming_state(node, state)
    parameterize(node, noise)

    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # CHANGE STATE SPACE IF BINARY APPROXIMATION IS USED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        changeStateSpace!(node, subproblem, state, algo_params.state_approximation_regime)
    end

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "dual_initialization" begin
        dual_vars_initial = initialize_duals(node, subproblem, algo_params.dual_initialization_method)
    end

    # RESET SOLVER (as it may have been changed in between for some reason)
    ########################################################################
    DynamicSDDiP.set_solver(node.subproblem, algo_params, applied_solvers, :backward_pass)

    # GET PRIMAL SOLUTION TO BOUND LAGRANGIAN DUAL
    ############################################################################
    @infiltrate algo_params.infiltrate_state in [:all]

    # REGULARIZATION
    if algo_params.regularization
        node.ext[:regularization_data] = Dict{Symbol,Any}()
        regularize_backward!(node, subproblem, algo_params.sigma[node_index])
    end

    # SOLVE PRIMAL PROBLEM (REGULARIZED OR NOT)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    if haskey(model.ext, :total_solves)
        model.ext[:total_solves] += 1
    else
        model.ext[:total_solves] = 1
    end

    # NOTE: TO-DO: Attempt numerical recovery as in SDDP

    solver_obj = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    @infiltrate algo_params.infiltrate_state in [:all]

    # PREPARE ACTUAL BACKWARD PASS METHOD BY DEREGULARIZATION
    ############################################################################
    if algo_params.regularization
        deregularize_backward!(node, subproblem)
    end

    # DUAL SOLUTION
    ############################################################################
    # Solve dual and return a dict with the multiplier of the copy constraints.
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
        lagrangian_results = get_dual_variables_backward(node, node_index, solver_obj, algo_params, applied_solvers, dual_vars_initial)
    end

    dual_values = lagrangian_results.dual_values
    bin_state = lagrangian_results.bin_state
    objective = lagrangian_results.intercept
    iterations = lagrangian_results.iterations
    lag_status = lagrangian_results.lag_status

    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # REGAIN ORIGINAL MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        rechangeStateSpace!(node, subproblem, state, algo_params.state_approximation_regime)
    end

    return (
        duals = dual_values,
        bin_state = bin_state,
        objective = objective,
        iterations = iterations,
        lag_status = lag_status,
    )
end


"""
Calling the Lagrangian dual solution method and determining dual variables
required to construct a cut
"""
function get_dual_variables_backward(
    node::SDDP.Node,
    node_index::Int64,
    solver_obj::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_vars_initial::Vector{Float64}
    )

    # storages for return of dual values and binary state values (trial point)
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()

    # TODO: implement smart choice for initial duals
    number_of_states = length(node.ext[:backward_data][:bin_states])
    # dual_vars = zeros(number_of_states)
    #solver_obj = JuMP.objective_value(node.ext[:linSubproblem])
    dual_vars = dual_vars_initial

    lag_obj = 0
    lag_iterations = 0
    lag_status = :none

    # Create an SDDiP integrality_handler here to store the Lagrangian dual information
    #TODO: Store tolerances in algo_params
    integrality_handler = SDDP.SDDiP(iteration_limit = algo_params.lagrangian_iteration_limit, atol = algo_params.lagrangian_atol, rtol = algo_params.lagrangian_rtol)
    integrality_handler = SDDP.update_integrality_handler!(integrality_handler, applied_solvers.MILP, number_of_states)
    node.ext[:lagrange] = integrality_handler

    # DETERMINE AND ADD BOUNDS FOR DUAL VARIABLES
    ############################################################################
    #TODO: Determine a norm of B (coefficient matrix of binary expansion)
    # We use the column sum norm here
    # But instead of calculating it exactly, we can also use the maximum
    # upper bound of all state variables as a bound

    B_norm_bound = 0
    for (name, state_comp) in node.states
        if state_comp.info.in.upper_bound > B_norm_bound
            B_norm_bound = state_comp.info.in.upper_bound
        end
    end
    dual_bound = algo_params.sigma[node_index] * B_norm_bound

    @infiltrate algo_params.infiltrate_state in [:all, :lagrange]
    #|| node.ext[:linSubproblem].ext[:sddp_policy_graph].ext[:iteration] == 12

    try
        # SOLUTION WITHOUT BOUNDED DUAL VARIABLES (BETTER TO OBTAIN BASIC SOLUTIONS)
        ########################################################################
        if algo_params.lagrangian_method == :kelley
            results = _kelley(node, node_index, solver_obj, dual_vars, integrality_handler, algo_params, applied_solvers, nothing)
            lag_obj = results.lag_obj
            lag_iterations = results.iterations
            lag_status = results.lag_status
        elseif algo_params.lagrangian_method == :bundle_level
            results = _bundle_level(node, node_index, solver_obj, dual_vars, integrality_handler, algo_params, applied_solvers, nothing)
            lag_obj = results.lag_obj
            lag_iterations = results.iterations
            lag_status = results.lag_status
        end

        # OPTIMAL VALUE CHECKS
        ########################################################################
        if algo_params.lag_status_regime == :rigorous
            if lag_status == :conv
                error("Lagrangian dual converged to value < solver_obj.")
            elseif lag_status == :sub
                error("Lagrangian dual had subgradients zero without LB=UB.")
            elseif lag_status == :iter
                error("Solving Lagrangian dual exceeded iteration limit.")
            end

        elseif algo_params.lag_status_regime == :lax
            # all cuts will be used as they are valid even though not necessarily tight
        end

        # DUAL VARIABLE BOUND CHECK
        ########################################################################
        # if one of the dual variables exceeds the bounds (e.g. in case of an
        # discontinuous value function), use bounded version of Kelley's method
        bound_check = true
        for dual_var in dual_vars
            if dual_var > dual_bound
                bound_check = false
            end
        end

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    catch e
        SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end

    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
        # TODO (maybe) change dual signs inside kelley to match LP duals
        dual_values[name] = -dual_vars[i]

        value = integrality_handler.old_rhs[i]
        x_name = node.ext[:backward_data][:bin_x_names][name]
        k = node.ext[:backward_data][:bin_k][name]
        bin_state[name] = BinaryState(value, x_name, k)
    end

    return (
        dual_values=dual_values,
        bin_state=bin_state,
        intercept=lag_obj,
        iterations=lag_iterations,
        lag_status=lag_status,
    )
end


"""
Calculate the lower bound (if minimizing, otherwise upper bound) of the problem
model at the point state.
"""
function calculate_bound(
    model::SDDP.PolicyGraph{T},
    root_state::Dict{Symbol,Float64} = model.initial_root_state;
    risk_measure = SDDP.Expectation(),
) where {T}

    # Note that here all children of the root node are solved, since the root
    # node is not node 1, but node 0.
    # In our case, this means that only stage 1 problem is solved again,
    # using the updated Bellman function from the backward pass.

    # Initialization.
    noise_supports = Any[]
    probabilities = Float64[]
    objectives = Float64[]
    problem_size = Dict{Symbol,Int64}[]
    current_belief = SDDP.initialize_belief(model)

    # Solve all problems that are children of the root node.
    for child in model.root_children
        if isapprox(child.probability, 0.0, atol = 1e-6)
            continue
        end
        node = model[child.term]
        for noise in node.noise_terms
            subproblem_results = solve_first_stage_problem(
                model,
                node,
                root_state,
                noise.term,
                Tuple{T,Any}[(child.term, noise.term)],
            )
            push!(objectives, subproblem_results.objective)
            push!(probabilities, child.probability * noise.probability)
            push!(noise_supports, noise.term)
            push!(problem_size, subproblem_results.problem_size)
        end
    end
    # Now compute the risk-adjusted probability measure:
    risk_adjusted_probability = similar(probabilities)
    offset = SDDP.adjust_probability(
        risk_measure,
        risk_adjusted_probability,
        probabilities,
        noise_supports,
        objectives,
        model.objective_sense == MOI.MIN_SENSE,
    )
    # Finally, calculate the risk-adjusted value.
    return (bound = sum(obj * prob for (obj, prob) in zip(objectives, risk_adjusted_probability)) +
           offset, problem_size = problem_size)
end


"""
Solving the first-stage problem to determine a lower bound
"""
function solve_first_stage_problem(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    state::Dict{Symbol,Float64},
    noise,
    scenario_path::Vector{Tuple{T,S}};
) where {T,S}

    ############################################################################
    # MODEL PARAMETRIZATION
    ############################################################################
    subproblem = node.subproblem

    # Parameterize the model. First, fix the value of the incoming state
    # variables. Then parameterize the model depending on `noise`. Finally,
    # set the objective.
    set_incoming_state(node, state)
    parameterize(node, noise)

    ############################################################################
    # SOLUTION
    ############################################################################
    JuMP.optimize!(subproblem)

    if haskey(model.ext, :total_solves)
        model.ext[:total_solves] += 1
    else
        model.ext[:total_solves] = 1
    end

    # Maybe attempt numerical recovery as in SDDP

    state = get_outgoing_state(node)
    stage_objective = JuMP.value(node.stage_objective)
    objective = JuMP.objective_value(subproblem)
    dual_values = get_dual_variables(node, node.integrality_handler)

    ############################################################################
    # DETERMINE THE PROBLEM SIZE
    ############################################################################
    problem_size = Dict{Symbol,Int64}()
    problem_size[:total_var] = size(JuMP.all_variables(subproblem),1)
    problem_size[:bin_var] = JuMP.num_constraints(subproblem, VariableRef, MOI.ZeroOne)
    problem_size[:int_var] = JuMP.num_constraints(subproblem, VariableRef, MOI.Integer)
    problem_size[:total_con] = JuMP.num_constraints(subproblem, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
                                + JuMP.num_constraints(subproblem, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64})
                                + JuMP.num_constraints(subproblem, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64})

    return (
        state = state,
        duals = dual_values,
        objective = objective,
        stage_objective = stage_objective,
        problem_size = problem_size
    )
end
