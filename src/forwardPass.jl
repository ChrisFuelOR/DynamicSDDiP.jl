# The functions
# > "forward_pass",
# > "solve_subproblem_forward",
# are derived from similar named functions (forward_pass, solve_subproblem)
# in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Executing the forward pass of an inner loop iteration for DynamicSDDiP.
"""
function forward_pass(model::SDDP.PolicyGraph{T}, options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams, applied_solvers::DynamicSDDiP.AppliedSolvers,
    ::SDDP.DefaultForwardPass) where {T}

    ############################################################################
    # SAMPLING AND INITIALIZATION (SIMLAR TO SDDP)
    ############################################################################
    # First up, sample a scenario. Note that if a cycle is detected, this will
    # return the cycle node as well.
    scenario_path, terminated_due_to_cycle = SDDP.sample_scenario(model, algo_params.sampling_scheme)

    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol,Float64}[]
    # Our initial incoming state.
    incoming_state_value = copy(model.initial_root_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0

    ############################################################################
    # ACTUAL ITERATION
    ########################################################################
    # Iterate down the scenario tree.
    for (depth, (node_index, noise)) in enumerate(scenario_path)
        node = model[node_index]

        ########################################################################
        # SET SOLVER
        ########################################################################
        DynamicSDDiP.set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

        ########################################################################
        # SUBPROBLEM SOLUTION
        ########################################################################
        # Solve the subproblem, note that `require_duals = false`.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_FP" begin
            subproblem_results = solve_subproblem_forward(
                model,
                node,
                node_index,
                incoming_state_value, # only values, no State struct!
                noise,
                scenario_path[1:depth],
                algo_params,
                algo_params.regularization_regime,
            )
        end
        # Cumulate the stage_objective.
        cumulative_value += subproblem_results.stage_objective
        # Set the outgoing state value as the incoming state value for the next
        # node.
        incoming_state_value = copy(subproblem_results.state)
        # Add the outgoing state variable to the list of states we have sampled
        # on this forward pass.
        push!(sampled_states, incoming_state_value)

    end

    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        # objective_states = objective_states,
        # belief_states = belief_states,
        cumulative_value = cumulative_value,
    )
end


"""
Solving the subproblem within the forward pass of a loop of DynamicSDDiP.
"""
# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise. If require_duals=true, also return the dual variables
# associated with the fixed constraint of the incoming state variables.
function solve_subproblem_forward(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    state::Dict{Symbol,Float64},
    noise,
    scenario_path::Vector{Tuple{T,S}},
    algo_params::DynamicSDDiP.AlgoParams,
    regularization_regime::DynamicSDDiP.AbstractRegularizationRegime;
) where {T,S}

    subproblem = node.subproblem

    ############################################################################
    # MODEL PARAMETRIZATION
    ############################################################################
    # Parameterize the model. First, fix the value of the incoming state
    # variables. Then parameterize the model depending on `noise`. Finally,
    # set the objective.
    set_incoming_state!(node, state)
    parameterize(node, noise)

    ############################################################################
    # REGULARIZE SUBPROBLEM IF REQUIRED
    ############################################################################
    if node_index > 1
        node.ext[:regularization_data] = Dict{Symbol,Any}()
        regularize_subproblem!(node, node_index, subproblem, regularization_regime)
    end

    ############################################################################
    # SOLUTION
    ############################################################################
    @infiltrate algo_params.infiltrate_state in [:all]
    JuMP.optimize!(subproblem)

    # Maybe attempt numerical recovery as in SDDP

    state = get_outgoing_state(node)
    objective = JuMP.objective_value(subproblem)
    stage_objective = objective - JuMP.value(bellman_term(node.bellman_function))
    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # DE-REGULARIZE SUBPROBLEM IF REQUIRED
    ############################################################################
    if node_index > 1
        deregularize_subproblem!(node, subproblem, regularization_regime)
    end

    return (
        state = state,
        objective = objective,
        stage_objective = stage_objective,
    )
end
