# The functions
# > "solve_all_children",
# > "solve_subproblem_backward",
# > "calculate_bound",
# > "solve_first_stage_problem"
# are derived from similar named functions (backward_pass,
# solve_all_children, solve_subproblem, calculate_bound,
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
    # objective_states::Vector{NTuple{N,Float64}},
    # belief_states::Vector{Tuple{Int,Dict{T,Float64}}}) where {T,NoiseType,N}
    ) where {T,NoiseType}

    ############################################################################
    # INITIALIZATION
    ############################################################################

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

        ########################################################################
        # SOLVE ALL CHILDREN PROBLEMS
        ########################################################################
        solve_all_children(
            model,
            node,
            node_index,
            items,
            1.0,
            # belief_state,
            # objective_state,
            outgoing_state,
            algo_params.backward_sampling_scheme,
            scenario_path[1:index],
            algo_params,
            applied_solvers
        )

        ########################################################################
        # RECONSTRUCT ANCHOR POINTS IN BACKWARD PASS
        ########################################################################
        variable_info = node.ext[:state_info_storage][name].out
        anchor_states = determine_anchor_states(node, outgoing_state, algo_params.state_approximation_regime, variable_info)
        @infiltrate algo_params.infiltrate_state in [:all]

        ########################################################################
        # REFINE BELLMAN FUNCTION BY ADDING CUTS
        ########################################################################
        """
        Note that for all backward openings, bin_state is the same,
        so we can just use bin_state[1] in the following.
        Maybe this should be changed later.
        """

        TimerOutputs.@timeit DynamicSDDiP_TIMER "update_bellman" begin
            new_cuts = refine_bellman_function(
                model,
                node,
                node_index,
                node.bellman_function,
                options.risk_measures[node_index],
                outgoing_state,
                anchor_states,
                items.bin_state[1],
                items.duals,
                items.supports,
                items.probability,
                items.objectives,
                algo_params,
                applied_solvers
            )
        end

        push!(cuts[node_index], new_cuts)
        # NOTE: This has to be adapted for stochastic case
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
    # belief_state,
    # objective_state,
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

    return
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
    set_incoming_state!(node, state)
    parameterize(node, noise)

    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # CHANGE STATE SPACE IF BINARY APPROXIMATION IS USED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        changeStateSpace!(node, subproblem, state, algo_params.state_approximation_regime)
    end

    ############################################################################
    # SOLVE DUAL PROBLEM TO OBTAIN CUT INFORMATION
    ############################################################################
    # Solve dual and return a dict with the multipliers of the copy constraints.
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_dual" begin
        dual_results = get_dual_solution(node, node_index, algo_params, applied_solvers, algo_params.duality_regime)
    end

    ############################################################################
    # REGAIN ORIGINAL MODEL IF BINARY APPROXIMATION IS USED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        rechangeStateSpace!(node, subproblem, state, algo_params.state_approximation_regime)
    end

    @infiltrate algo_params.infiltrate_state in [:all]

    return (
        duals = dual_results.dual_values,
        bin_state = dual_results.bin_state,
        objective = dual_results.intercept,
        iterations = dual_results.iterations,
        lag_status = dual_results.lag_status,
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
    set_incoming_state!(node, state)
    parameterize(node, noise)

    ############################################################################
    # SOLUTION
    ############################################################################
    JuMP.optimize!(subproblem)

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
