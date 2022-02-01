# Some functions are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Executing the backward pass of a loop of DynamicSDDiP if the unified framework
is used, that is, an aggregated Lagrangian dual is solved for all backward
openings together.
"""
function backward_pass(
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    scenario_path::Vector{Tuple{T,NoiseType}},
    sampled_states::Vector{Dict{Symbol,Float64}},
    epi_states::Vector{Float64},
    framework_regime::DynamicSDDiP.UnifiedFramework,
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
        epi_state = epi_states[index]
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
        # SOLVE ALL CHILDREN PROBLEMS TO GET A BOUND FOR THE LAGRANGIAN DUAL
        ########################################################################
        primal_obj = solve_all_children_primal(
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
        # SOLVE THE LAGRANGIAN DUAL
        ########################################################################
        solve_aggregated_dual(
            model,
            node,
            node_index,
            items,
            1.0,
            # belief_state,
            # objective_state,
            outgoing_state,
            epi_state,
            primal_obj,
            algo_params.backward_sampling_scheme,
            scenario_path[1:index],
            algo_params,
            applied_solvers
        )

        ########################################################################
        # RECONSTRUCT ANCHOR POINTS IN BACKWARD PASS
        ########################################################################
        anchor_states = determine_anchor_states(node, outgoing_state, algo_params.state_approximation_regime)
        @infiltrate algo_params.infiltrate_state in [:all]

        ########################################################################
        # REFINE BELLMAN FUNCTION BY ADDING CUTS
        ########################################################################
        """
        Note that for all backward openings, bin_state is the same,
        so we can just use bin_state[1] in the following.
        Maybe this should be changed later.
        """

        @infiltrate
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

        # NOTE: Implement using cuts for similar nodes as in SDDP

    end

    return cuts
end


"""
Solving all children within the backward pass to determine a primal solution.
"""
function solve_all_children_primal(
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

    primal_obj = 0

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
            ################################################################
            # SOLVE THE BACKWARD PASS PRIMAL PROBLEM
            ################################################################
            TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_BP" begin
                primal_obj_scenario = solve_subproblem_backward_primal(
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

            primal_obj += noise.probability * primal_obj_scenario
        end
    end

    return primal_obj
end


"""
Solving the backward pass problem for one specific child to obtain
    primal_obj bound for Lagrangian dual
"""
# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise. If require_duals=true, also return the dual variables
# associated with the fixed constraint of the incoming state variables.
function solve_subproblem_backward_primal(
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
    # CHANGE STATE SPACE AND REGULARIZE IF REQUIRED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        changeStateSpace!(node, subproblem, state, algo_params.state_approximation_regime)
    end

    # REGULARIZE PROBLEM IF REGULARIZATION IS USED
    node.ext[:regularization_data] = Dict{Symbol,Any}()
    regularize_bw!(node, node_index, subproblem, algo_params.regularization_regime, algo_params.state_approximation_regime)

    # RESET SOLVER (as it may have been changed in between for some reason)
    DynamicSDDiP.set_solver!(subproblem, algo_params, applied_solvers, :backward_pass)

    ############################################################################
    # SOLVE PRIMAL PROBLEM
    ############################################################################
    # SOLVE PRIMAL PROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    primal_obj_scenario = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    ############################################################################
    # REGAIN UNREGULARIZED MODEL IF REQUIRED
    ############################################################################
    # DEREGULARIZE PROBLEM IF REQUIRED
    deregularize_bw!(node, subproblem, algo_params.regularization_regime, algo_params.state_approximation_regime)

    @infiltrate algo_params.infiltrate_state in [:all]

    return primal_obj_scenario

end
