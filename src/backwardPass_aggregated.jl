"""
Executing the backward pass of a loop of DynamicSDDiP if the unified framework
is used, that is, an aggregated Lagrangian dual is solved for all backward
openings together.
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
    cut_aggregation_regime::DynamicSDDiP.SingleCutRegime,
    ) where {T,NoiseType}

    ############################################################################
    # SOLVE PRIMAL PROBLEMS FOR ALL CHILDREN
    ############################################################################
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
        cut_generation_regime,
        applied_solvers
    )
    #TODO: Is this required? We do not really use the primal_obj value in the
    #dual of the unified framework anyway.

    ############################################################################
    # SOLVE THE AGGREGATED DUAL PROBLEM
    ############################################################################
    solve_aggregated_dual(
            model,
            node,
            node_index,
            index,
            items,
            1.0,
            # belief_state,
            # objective_state,
            outgoing_state,
            epi_states[1], # there is only one epi_state for Single_cut
            primal_obj,
            add_cut_flag,
            algo_params.backward_sampling_scheme,
            scenario_path[1:index],
            algo_params,
            cut_generation_regime,
            applied_solvers,
            duality_regime.dual_solution_regime,
        )

    return
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
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
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
                    cut_generation_regime,
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
    # CHANGE STATE SPACE AND REGULARIZE IF REQUIRED
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        changeStateSpace!(node, subproblem, state, cut_generation_regime.state_approximation_regime)
    end

    # REGULARIZE PROBLEM IF REGULARIZATION IS USED
    node.ext[:regularization_data] = Dict{Symbol,Any}()
    regularize_bw!(node, node_index, subproblem, cut_generation_regime, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    # RESET SOLVER (as it may have been changed in between for some reason)
    set_solver!(subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)

    ############################################################################
    # SOLVE PRIMAL PROBLEM
    ############################################################################
    # SOLVE PRIMAL PROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    primal_obj_scenario = JuMP.objective_value(subproblem)
    # Try recovering from numerical issues
    if (JuMP.termination_status(subproblem) != MOI.OPTIMAL)
        elude_numerical_issues!(subproblem, algo_params)
    end

    ############################################################################
    # REGAIN UNREGULARIZED MODEL IF REQUIRED
    ############################################################################
    # DEREGULARIZE PROBLEM IF REQUIRED
    deregularize_bw!(node, subproblem, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    return primal_obj_scenario

end
