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
    epi_states::Dict{Symbol,Vector{Float64}},
    # objective_states::Vector{NTuple{N,Float64}},
    # belief_states::Vector{Tuple{Int,Dict{T,Float64}}}) where {T,NoiseType,N}
    ) where {T,NoiseType}

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # storage for data on solving Lagrangian dual
    model.ext[:agg_lag_iterations] = 0
    model.ext[:corr_lag_iterations] = 0
    model.ext[:corr_realizations] = 0
    model.ext[:lag_iterations] = Float64[]

    ############################################################################
    # TRAVERSE THE STAGES BACKWARDS
    ############################################################################
    for index in length(scenario_path):-1:1
        # Determine trial state
        outgoing_state = sampled_states[index]

        # Determine current node (stage)
        node_index, _ = scenario_path[index]
        node = model[node_index]
        if length(node.children) == 0
            continue
        end

        # Reset cut counter
        node.ext[:total_cuts] = 0
        node.ext[:active_cuts] = 0

        ########################################################################
        # ITERATE OVER CUT GENERATION REGIMES
        ########################################################################
        for cut_generation_regime in algo_params.cut_generation_regimes

            # New cuts for cut_generation regime are only generated and added
            # if we are in the right iterations
            if (model.ext[:iteration] >= cut_generation_regime.iteration_to_start
                && model.ext[:iteration] <= cut_generation_regime.iteration_to_stop)

                items = BackwardPassItems(T, SDDP.Noise)

                # Determine required epi_states
                epi_states_node = epi_states[Symbol(node_index)]

                # Add flag to decide whether a cut is added
                add_cut_flag = true

                # Dict to store values of binary approximation of the state
                # Note that we could also retrieve this from the actual trial point
                # (outgoing_state) or from its approximation via binexpand. However,
                # this collection is not only important to obtain the correct values,
                # but also to store them together with the symbol/name of the variable.
                node.ext[:binary_state_values] = Dict{Symbol, Float64}()

                ################################################################
                # Solve backward pass problems for all realizations depending on framework
                ################################################################
                backward_pass_node(
                    model,
                    node,
                    node_index,
                    index,
                    items,
                    outgoing_state,
                    epi_states_node,
                    scenario_path[1:index],
                    add_cut_flag,
                    algo_params,
                    cut_generation_regime,
                    applied_solvers,
                    cut_generation_regime.duality_regime,
                    algo_params.cut_aggregation_regime,
                )

                ################################################################
                # RECONSTRUCT ANCHOR POINTS IN BACKWARD PASS
                ################################################################
                anchor_state = determine_anchor_state(node, outgoing_state, cut_generation_regime.state_approximation_regime)
                Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

                ################################################################
                # REFINE BELLMAN FUNCTION BY ADDING CUTS
                ################################################################
                """
                Note that for all backward openings, bin_state is the same,
                so we can just use bin_state[1] in the following.
                Maybe this should be changed later.
                """
                TimerOutputs.@timeit DynamicSDDiP_TIMER "update_bellman" begin
                    refine_bellman_function(
                        model,
                        node,
                        node_index,
                        node.bellman_function,
                        options.risk_measures[node_index],
                        outgoing_state,
                        epi_states_node,
                        anchor_state,
                        items.bin_state[1],
                        items.duals,
                        items.dual_0_var,
                        items.supports,
                        items.probability,
                        items.objectives,
                        items.add_cut_flags,
                        algo_params,
                        cut_generation_regime,
                        applied_solvers
                    )
                end

                ################################################################
                #NOTE: I did not include the similar node thing from SDDP.jl.
                #Not really sure what it means anyway.

                ################################################################
                # UPDATE LIST OF BENDERS CUTS FOR CHEN & LUEDTKE APPROACH
                ################################################################
                """ Note that we store the indices for all Benders cuts instead
                of only the last K, since it is not known here which value K has."""

                if isa(cut_generation_regime.duality_regime, DynamicSDDiP.LinearDuality) || isa(cut_generation_regime.duality_regime, DynamicSDDiP.StrengthenedDuality)
                    update_Benders_cut_list!(node, items.add_cut_flags, algo_params.cut_aggregation_regime, cut_generation_regime.state_approximation_regime)
                end
                # TODO: A similar approach can be used for Lagrangian cuts

                ################################################################
                # LOGGING
                ################################################################
                # Logging of lag_iterations
                if isa(cut_generation_regime.duality_regime, Union{DynamicSDDiP.LagrangianDuality,DynamicSDDiP.UnifiedLagrangianDuality})
                    push!(model.ext[:lag_iterations], Statistics.mean(items.lag_iterations))
                end

            end
        end

        # Update cut count
        count_cuts(node, node.bellman_function.global_theta, 1)
        for V in node.bellman_function.local_thetas
            count_cuts(node, V)
        end

    end

    return
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
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_first" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    state = get_outgoing_state(node)
    stage_objective = JuMP.value(node.stage_objective)
    objective = JuMP.objective_value(subproblem)

    ############################################################################
    # DETERMINE THE PROBLEM SIZE
    ############################################################################
    problem_size = Dict{Symbol,Int64}()
    problem_size[:total_var] = size(JuMP.all_variables(subproblem),1)
    problem_size[:bin_var] = JuMP.num_constraints(subproblem, JuMP.VariableRef, MOI.ZeroOne)
    problem_size[:int_var] = JuMP.num_constraints(subproblem, JuMP.VariableRef, MOI.Integer)
    problem_size[:total_con] = JuMP.num_constraints(subproblem, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, MOI.LessThan{Float64})
                                + JuMP.num_constraints(subproblem, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, MOI.GreaterThan{Float64})
                                + JuMP.num_constraints(subproblem, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, MOI.EqualTo{Float64})

    return (
        state = state,
        objective = objective,
        stage_objective = stage_objective,
        problem_size = problem_size
    )
end


"""
Updating the list of Benders cuts for the approach by Chen & Luedtke
for current cut_generation_regime without state approximation
and for SingleCutRegime.
"""
function update_Benders_cut_list!(
    node::SDDP.Node,
    add_cut_flags::Vector{Bool},
    cut_aggregation_regime::DynamicSDDiP.SingleCutRegime,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    if any(add_cut_flags)
        cut_index = lastindex(node.bellman_function.global_theta.cuts)
        push!(node.ext[:Benders_cuts_original], (cut_index, :all))
    end
end

"""
Updating the list of Benders cuts for the approach by Chen & Luedtke
for current cut_generation_regime without state approximation
and for MultiCutRegime.
"""
function update_Benders_cut_list!(
    node::SDDP.Node,
    add_cut_flags::Vector{Bool},
    cut_aggregation_regime::DynamicSDDiP.MultiCutRegime,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    for i = 1:length(add_cut_flags)
        if add_cut_flags[i]
            cut_index = lastindex(node.bellman_function.local_thetas[i].cuts)
            push!(node.ext[:Benders_cuts_original], (cut_index, Symbol(i)))
        end
    end
end

"""
Updating the list of Benders cuts for the approach by Chen & Luedtke
for current cut_generation_regime with binary state approximation
and for SingleCutRegime.
"""
function update_Benders_cut_list!(
    node::SDDP.Node,
    add_cut_flags::Vector{Bool},
    cut_aggregation_regime::DynamicSDDiP.SingleCutRegime,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    if any(add_cut_flags)
        cut_index = lastindex(node.bellman_function.global_theta.cuts)
        push!(node.ext[:Benders_cuts_binary], (cut_index, :all))
    end
end

"""
Updating the list of Benders cuts for the approach by Chen & Luedtke
for current cut_generation_regime with binary state approximation
and for MultiCutRegime.
"""
function update_Benders_cut_list!(
    node::SDDP.Node,
    add_cut_flags::Vector{Bool},
    cut_aggregation_regime::DynamicSDDiP.MultiCutRegime,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    for i = 1:length(add_cut_flags)
        if add_cut_flags[i]
            cut_index = lastindex(node.bellman_function.local_thetas[i].cuts)
            push!(node.ext[:Benders_cuts_binary], (cut_index, Symbol(i)))
        end
    end
end
