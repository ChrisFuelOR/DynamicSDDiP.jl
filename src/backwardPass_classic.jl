"""
Executing the backward pass for a given stage in DynamicSDDiP
if the classical cut generation method is used (for multi-cut or single-cut).
"""
function backward_pass_node(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    index::Int64,
    items::BackwardPassItems,
    # belief_state,
    # objective_state,
    outgoing_state::Dict{Symbol,Float64},
    epi_states::Vector{Float64},
    scenario_path,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::Union{DynamicSDDiP.LagrangianDuality, DynamicSDDiP.LinearDuality, DynamicSDDiP.StrengthenedDuality},
    cut_aggregation_regime::Union{DynamicSDDiP.SingleCutRegime, DynamicSDDiP.MultiCutRegime}
    ) where {T,NoiseType}

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
        epi_states,
        algo_params.backward_sampling_scheme,
        scenario_path[1:index],
        add_cut_flag,
        algo_params,
        cut_generation_regime,
        applied_solvers
    )

    return
end


"""
Executing the backward pass for a given stage in DynamicSDDiP
if the unified framework cut generation method is used (for multi-cut).
"""
function backward_pass_node(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    index::Int64,
    items::BackwardPassItems,
    # belief_state,
    # objective_state,
    outgoing_state::Dict{Symbol,Float64},
    epi_states::Vector{Float64},
    scenario_path,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    cut_aggregation_regime::DynamicSDDiP.MultiCutRegime,
    ) where {T,NoiseType}

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
        epi_states,
        algo_params.backward_sampling_scheme,
        scenario_path[1:index],
        add_cut_flag,
        algo_params,
        cut_generation_regime,
        applied_solvers
    )

    return
end


"""
Solving all children on this stage in the backward pass.
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
    epi_states::Vector{Float64},
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    scenario_path,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
) where {T}
    length_scenario_path = length(scenario_path)
    for child in node.children
        if isapprox(child.probability, 0.0, atol = 1e-6)
            continue
        end
        child_node = model[child.term]
        for i in 1:length(SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node))
            noise = SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node)[i]

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
                push!(items.dual_0_var, items.dual_0_var[sol_index])
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
                # DETERMINE ASSOCIATED EPI_STATE
                ################################################################
                if algo_params.cut_type == SDDP.SINGLE_CUT
                    #epi_state = epi_states[1]
                    epi_state = Inf
                elseif algo_params.cut_type == SDDP.MULTI_CUT
                    epi_state = epi_states[i]
                end

                ################################################################
                # SOLVE THE BACKWARD PASS PROBLEM
                ################################################################
                TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_BP" begin
                    subproblem_results = solve_subproblem_backward(
                        model,
                        child_node,
                        node_index+1,
                        outgoing_state,
                        epi_state,
                        noise.term,
                        i,
                        scenario_path,
                        add_cut_flag,
                        algo_params,
                        cut_generation_regime,
                        applied_solvers
                    )
                end
                push!(items.duals, subproblem_results.duals)
                push!(items.dual_0_var, subproblem_results.dual_0_var)
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
    epi_state::Float64,
    noise,
    i::Int64,
    scenario_path::Vector{Tuple{T,S}},
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
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

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # CHANGE STATE SPACE IF BINARY APPROXIMATION IS USED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        changeStateSpace!(node, subproblem, state, cut_generation_regime.state_approximation_regime)
    end

    ############################################################################
    # SOLVE DUAL PROBLEM TO OBTAIN CUT INFORMATION
    ############################################################################
    # Solve dual and return a dict with the multipliers of the copy constraints.
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_dual" begin
        dual_results = get_dual_solution(node, node_index, i, epi_state, add_cut_flag, algo_params, cut_generation_regime, applied_solvers, cut_generation_regime.duality_regime)
    end

    if haskey(model.ext, :total_solves)
        model.ext[:total_solves] += 1
    else
        model.ext[:total_solves] = 1
    end

    ############################################################################
    # REGAIN ORIGINAL MODEL IF BINARY APPROXIMATION IS USED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        rechangeStateSpace!(node, subproblem, state, cut_generation_regime.state_approximation_regime)
    end

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    return (
        duals = dual_results.dual_values,
        dual_0_var = dual_results.dual_0_var,
        bin_state = dual_results.bin_state,
        objective = dual_results.intercept,
        iterations = dual_results.iterations,
        lag_status = dual_results.lag_status,
    )

end
