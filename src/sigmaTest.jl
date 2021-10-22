# The functions
# > "forward_sigma_test"
# > "solve_subproblem_sigma_test"
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
Executing a forward pass to check if parameter sigma is chosen sufficiently
large in DynamicSDDiP.
"""
function forward_sigma_test(
    model::SDDP.PolicyGraph{T}, options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams, applied_solvers::DynamicSDDiP.AppliedSolvers,
    scenario_path::Vector{Tuple{T,S}}, sigma_increased::Bool) where {T,S}

    # INITIALIZATION (NO SAMPLING HERE!)
    ############################################################################
    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol,Float64}[]
    # Our initial incoming state.
    incoming_state_value = copy(options.initial_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0

    # ACTUAL ITERATION
    ########################################################################
    # Iterate down the scenario tree.
    for (depth, (node_index, noise)) in enumerate(scenario_path)
        node = model[node_index]

        # SET SOLVER
        ########################################################################
        DynamicSDDiP.set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

        # SOLVE REGULARIZED PROBLEM
        ############################################################################
        # This has to be done in this sigma test again, since the current upper bound
        # may not be the best upper bound and thus, just comparing the bounds is not
        # sufficient. Moreover, we do it in this loop to not increase sigma for
        # all stages, but only where it is needed

        # Solve the subproblem, note that `require_duals = false`.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_sigma_test_reg" begin
            reg_results = solve_subproblem_forward(
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

        # SOLVE NON-REGULARIZED PROBLEM
        ########################################################################
        # Solve the subproblem, note that `require_duals = false`.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_sigma_test" begin
            non_reg_results = solve_subproblem_forward(
                model,
                node,
                node_index,
                incoming_state_value, # only values, no State struct!
                noise,
                scenario_path[1:depth],
                algo_params,
                DynamicSDDiP.NoRegularization(),
            )
        end

        # COMPARE SOLUTIONS
        ########################################################################
        if !isapprox(reg_results.objective, non_reg_results.objective)
            sigma = algo_params.regularization_regime.sigma[node_index]
            sigma_factor = algo_params.regularization_regime.sigma_factor
            # if stage objectives are not approximately equal, then the
            # regularization is not exact and sigma should be increased
            sigma = sigma * sigma_factor
            # marker if new inner loop iteration should be started
            # instead of heading to outer loop
            sigma_increased = true
        end

        # Cumulate the stage_objective. (NOTE: not really required anymore)
        cumulative_value += non_reg_results.stage_objective
        # Set the outgoing state value as the incoming state value for the next node.
        # NOTE: We use the states determined by the regularized problem here
        #incoming_state_value = copy(subproblem_results.state)
        incoming_state_value = copy(reg_results.state)
        # Add the outgoing state variable to the list of states we have sampled
        # on this forward pass.
        push!(sampled_states, incoming_state_value)

    end

    # ===== End: drop off starting state if terminated due to cycle =====
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        # objective_states = objective_states,
        # belief_states = belief_states,
        cumulative_value = cumulative_value,
        sigma_increased = sigma_increased,
    )
end
