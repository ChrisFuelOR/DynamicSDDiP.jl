# The functions
# > "solve",
# > "master_loop"
# > "inner_loop"
# > "iteration"
# are derived from similar named functions (solve, master_loop, iteration) in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Solves the `model`. In contrast to SDDP.jl, all parameters configuring
the algorithm are given (and possibly pre-defined) in algo_params.
"""
function solve(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers
)

    ############################################################################
    # INITIALIZATION (SIMILAR TO SDDP.jl)
    ############################################################################
    # Reset the TimerOutput
    TimerOutputs.reset_timer!(DynamicSDDiP_TIMER)

    # Prepare logging
    log_file_handle = open(algo_params.log_file, "a")
    log = Log[]

    if print_level > 0
        print_helper(print_banner, log_file_handle)
    end

    if print_level > 1
        print_helper(print_parameters, log_file_handle, algo_params, applied_solvers)
    end

    # Maybe add run_numerical_stability_report as in SDDP.jl later

    # Prepare stopping rules
    #---------------------------------------------------------------------------
    # Convert the vector to an AbstractStoppingRule. Otherwise if the user gives
    # something like stopping_rules = [SDDP.IterationLimit(100)], the vector
    # will be concretely typed and we can't add a TimeLimit.
    #stopping_rules = algo_params.stopping_rules
    #stopping_rules = convert(Vector{SDDP.AbstractStoppingRule}, stopping_rules)
    if length(algo_params.stopping_rules) == 0
        @warn(
            "You haven't specified a stopping rule! You can only terminate " *
            "the call to DynamicSDDiP.solve via a keyboard interrupt ([CTRL+C])."
        )
    end

    # Prepare options for logging
    #---------------------------------------------------------------------------
    options = DynamicSDDiP.Options(
        model,
        model.initial_root_state,
        time(),
        log,
        log_file_handle
    )

    ############################################################################
    # RE-INITIALIZE THE EXISTING VALUE FUNCTION AND PREPARE CUT SELECTION
    ############################################################################
    # Update the nodes with the selected cut type (SINGLE_CUT or MULTI_CUT)
    # and the cut deletion minimum.
    if algo_params.cut_selection_regime.cut_deletion_minimum < 0
        algo_params.cut_selection_regime.cut_deletion_minimum = typemax(Int)
    end

    # Change to nonconvex Bellman function to use non-convex cuts
    #---------------------------------------------------------------------------
    # fortunately, node.bellman_function requires no specific type
    for (key, node) in model.nodes
        node.bellman_function = DynamicSDDiP.initialize_bellman_function(bellman_function, model, node)
        node.bellman_function.cut_type = algo_params.cut_type
        node.bellman_function.global_theta.cut_oracle.deletion_minimum =
            algo_params.cut_selection_regime.cut_deletion_minimum
        for oracle in node.bellman_function.local_thetas
            oracle.cut_oracle.deletion_minimum = algo_params.cut_selection_regime.cut_deletion_minimum
        end
    end

    ############################################################################
    # MODEL START
    ############################################################################
    status = :not_solved
    try
        status = solve_DynamicSDDiP(parallel_scheme, model, options, algo_params, applied_solvers)
    catch ex
        if isa(ex, InterruptException)
            status = :interrupted
            interrupt(parallel_scheme)
        else
            close(log_file_handle)
            rethrow(ex)
        end
    finally
    end

    ############################################################################
    # lOG MODEL RESULTS
    ############################################################################
    results = DynamicSDDiP.Results(status, log)
    model.ext[:results] = results
    if print_level > 0
        print_helper(print_footer, log_file_handle, results)
        if print_level > 1
            print_helper(TimerOutputs.print_timer, log_file_handle, DynamicSDDiP_TIMER)
            print_helper(println, log_file_handle)
        end
    end
    close(log_file_handle)
    return
end

"""
Solves the `model` using DynamicSDDiP in a serial scheme.
"""

function solve_DynamicSDDiP(parallel_scheme::SDDP.Serial, model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options, algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers) where {T}

    ############################################################################
    # SET UP STATE VARIABLE INFORMATION
    ############################################################################
    for (node_index, children) in model.nodes
        node = model.nodes[node_index]

        # Set info for state.in as for previous stage's state.out
        #-----------------------------------------------------------------------
        if node_index > 1
            for (i, (name, state)) in enumerate(node.states)
                # Get correct state_info
                state_info = model.nodes[node_index-1].states[name].info.out

                if state_info.has_lb
                    JuMP.set_lower_bound(state.in, state_info.lower_bound)
                end
                if state_info.has_ub
                    JuMP.set_upper_bound(state.in, state_info.upper_bound)
                end
                if state_info.binary
                    JuMP.set_binary(state.in)
                elseif state_info.integer
                    JuMP.set_integer(state.in)
                end

                # Store info to reset it later
                state.info.in = state_info
            end
        end
    end

    @infiltrate algo_params.infiltrate_state == :all

    ############################################################################
    # CALL ACTUAL SOLUTION PROCEDURE
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "loop" begin
        status = master_loop(parallel_scheme, model, options, algo_params, applied_solvers)
    end
    return status

end


"""
Loop function of DynamicSDDiP.
"""

function master_loop(parallel_scheme::SDDP.Serial, model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options, algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    regularization_regime::DynamicSDDiP.Regularization) where {T}

    ############################################################################
    # INITIALIZE PARAMETERS REQUIRED FOR REFINEMENTS
    ############################################################################
    previous_solution = nothing
    previous_bound = nothing
    sigma_increased = false
    bound_check = false

    ############################################################################
    # INITIALIZE BEST KNOWN POINT AND OBJECTIVE VALUE
    ############################################################################
    model.ext[:best_objective] = model.objective_sense == JuMP.MOI.MIN_SENSE ? Inf : -Inf
    model.ext[:best_point] = Vector{Dict{Symbol,Float64}}()

    ############################################################################
    # ACTUAL LOOP
    ############################################################################
    while true
        # start an iteration
        TimerOutputs.@timeit DynamicSDDiP_TIMER "iterations" begin
            result = iteration(model, options, algo_params, applied_solvers, previous_solution, bound_check, sigma_increased)
        end

        # logging
        log_iteration(algo_params, options.log_file_handle, options.log)

        # initialize parameters
        previous_solution = result.current_sol
        previous_bound = result.lower_bound
        sigma_increased = false
        bound_check = true

        @infiltrate algo_params.infiltrate_state in [:all, :sigma]

        # check for convergence and if not achieved, update parameters
        convergence_results = convergence_handler(
            result,
            model,
            options,
            algo_params,
            applied_solvers,
            sigma_increased,
            bound_check,
            previous_solution,
            previous_bound,
            algo_params.regularization_regime
        )

        previous_solution = convergence_results.previous_solution
        previous_bound = convergence_results.previous_bound
        sigma_increased = convergence_results.sigma_increased
        bound_check = convergence_results.bound_check

    end

    return
end


"""
Convergence handler if regularization is used.
"""

function convergence_handler(result::DynamicSDDiP.IterationResult,
    model::SDDP.PolicyGraph{T}, options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    sigma_increased::Bool, bound_check::Bool,
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    previous_bound::Union{Float64,Nothing},
    regularization_regime::DynamicSDDiP.Regularization) where {T}

    ############################################################################
    # IF CONVERGENCE IS ACHIEVED, COMPARE TRUE AND REGULARIZED PROBLEM
    ############################################################################
    if result.has_converged
        TimerOutputs.@timeit DynamicSDDiP_TIMER "sigma_test" begin
            sigma_test_results = forward_sigma_test(model, options, algo_params, applied_solvers, result.scenario_path, sigma_increased)
        end
        sigma_increased = sigma_test_results.sigma_increased
        @infiltrate algo_params.infiltrate_state in [:all, :sigma]

        if sigma_increased
            # reset previous values, as model is changed and convergence not achieved
            previous_solution = nothing
            previous_bound = nothing
            # binary refinement only when no sigma refinement has been made
            bound_check = false
        else
            ####################################################################
            # THE ALGORITHM TERMINATES
            ####################################################################
            # return convergence status
            return result.status
        end

    ############################################################################
    # IF NO CONVERGENCE IS ACHIEVED, DO DIFFERENT CHECKS
    ############################################################################
    else
        # CHECK IF LOWER BOUND HAS IMPROVED
        ########################################################################
        # NOTE: If not, then the cut was (probably) not tight enough,
        # so the binary approximation should be refined in the next iteration.
        # As different trial states could yield the same lower bound, this
        # is only done if also the trial state does not change, though.
        if !isnothing(previous_bound)
            if !isapprox(previous_bound, result.lower_bound)
                bound_check = false
            end
        else
            bound_check = false
        end

        # CHECK IF SIGMA SHOULD BE INCREASED (DUE TO LB > UB)
        ########################################################################
        if result.upper_bound - result.lower_bound < - 1e-8 # NOTE
                regularization_regime.sigma = regularization_regime.sigma * regularization_regime.sigma_factor
                sigma_increased = true
                previous_solution = nothing
                previous_bound = nothing
                bound_check = false
        else
            sigma_increased = false
        end

    end

    return (
        sigma_increased = sigma_increased,
        previous_solution = previous_solution,
        previous_bound = previous_bound,
        bound_check = bound_check
    )
end


"""
Convergence handler if regularization is not used.
"""

function convergence_handler(result::DynamicSDDiP.IterationResult,
    model::SDDP.PolicyGraph{T}, options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    sigma_increased::Bool, bound_check::Bool,
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    previous_bound::Union{Float64,Nothing},
    regularization_regime::DynamicSDDiP.NoRegularization) where {T}

    ############################################################################
    # IF CONVERGENCE IS ACHIEVED, COMPARE TRUE AND REGULARIZED PROBLEM
    ############################################################################
    if result.has_converged
        @infiltrate algo_params.infiltrate_state in [:all, :sigma]
        ########################################################################
        # THE ALGORITHM TERMINATES
        ########################################################################
        return result.status

    ############################################################################
    # IF NO CONVERGENCE IS ACHIEVED, DO DIFFERENT CHECKS
    ############################################################################
    else
        # CHECK IF LOWER BOUND HAS IMPROVED
        ########################################################################
        # If not, then the cut was (probably) not tight enough,
        # so the binary approximation should maybe be refined in the next iteration.
        if !isnothing(previous_bound)
            if !isapprox(previous_bound, result.lower_bound)
                bound_check = false
            end
        else
            bound_check = false
        end

        # CHECK IF LB > UB
        ########################################################################
        if result.upper_bound - result.lower_bound < - 1e-8 # NOTE
            error("LB < UB for DynamicSDDiP. Terminating.")
        end

    end

    return (
        sigma_increased = sigma_increased,
        previous_solution = previous_solution,
        previous_bound = previous_bound,
        bound_check = bound_check
    )
end


"""
Executing an inner loop iteration for DynamicSDDiP.
"""
function iteration(
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    bound_check::Bool,
    sigma_increased::Bool
    ) where {T}

    ############################################################################
    # SET ITERATION COUNTER
    ############################################################################
    if haskey(model.ext, :iteration)
        model.ext[:iteration] += 1
    else
        model.ext[:iteration] = 1
    end

    ############################################################################
    # FORWARD PASS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "forward_pass" begin
        forward_trajectory = DynamicSDDiP.forward_pass(model, options, algo_params, applied_solvers, algo_params.forward_pass)
    end

    ############################################################################
    # BINARY REFINEMENT
    ############################################################################
    # If the forward pass solution and the lower bound did not change during
    # the last iteration, then increase the binary precision (for all stages)
    refinement_check = false
    binary_refinement = :none

    TimerOutputs.@timeit DynamicSDDiP_TIMER "bin_refinement" begin
        if !isnothing(previousSolution) && bound_check
                refinement_check = DynamicSDDiP.binary_refinement_check(
                    model,
                    previous_solution,
                    forward_trajectory.sampled_states,
                    refinement_check,
                    algo_params.state_approximation_regime
                )
        end
        if refinement_check
            binary_refinement = DynamicSDDiP.binary_refinement(
                model,
                algo_params.state_approximation_regime.binary_precision,
                binary_refinement
            )
        end
    end
    # bound_check = true
    @infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # BACKWARD PASS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "backward_pass" begin
        cuts = DynamicSDDiP.backward_pass(
            model,
            options,
            algo_params,
            applied_solvers,
            forward_trajectory.scenario_path,
            forward_trajectory.sampled_states,
            forward_trajectory.objective_states,
            forward_trajectory.belief_states,
        )
    end

    ############################################################################
    # CALCULATE LOWER BOUND
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "calculate_bound" begin
        first_stage_results = calculate_bound(model)
    end
    bound = first_stage_results.bound

    ############################################################################
    # CHECK IF BEST KNOWN SOLUTION HAS BEEN IMPROVED
    ############################################################################
    if model.objective_sense == JuMP.MOI.MIN_SENSE
        if forward_trajectory.cumulative_value < model.ext[:best_objective] || sigma_increased
            # udpate best upper bound
            model.ext[:best_objective] = forward_trajectory.cumulative_value
            # update best point so far
            model.ext[:best_point] = forward_trajectory.sampled_states
        end
    else
        if forward_trajectory.cumulative_value > model.ext[:best_objective] || sigma_increased
            # udpate best lower bound
            model.ext[:best_objective] = forward_trajectory.cumulative_value
            # update best point so far
            model.ext[:best_point] = forward_trajectory.sampled_states
        end
    end

    ############################################################################
    # PREPARE LOGGING
    ############################################################################
    subproblem_size = first_stage_results.problem_size[1]
    model.ext[:total_cuts] = 0
    model.ext[:active_cuts] = 0

    for (node_index, children) in model.nodes
        node = model.nodes[node_index]

        if length(node.children) == 0
            continue
        end
        model.ext[:total_cuts] += node.ext[:total_cuts]
        model.ext[:active_cuts] += node.ext[:active_cuts]
    end

    push!(
         options.log,
         Log(
             model.ext[:iteration],
             bound,
             model.ext[:best_objective],
             forward_trajectory.cumulative_value,
             forward_trajectory.sampled_states,
             time() - options.start_time,
             sigma_increased,
             binary_refinement,
             subproblem_size,
             algo_params,
             model.ext[:lag_iterations],
             model.ext[:lag_status],
             model.ext[:total_cuts],
             model.ext[:active_cuts],
         ),
     )

    ############################################################################
    # CHECK IF THE LOOP CONVERGED YET
    ############################################################################
    has_converged, status = convergence_test(model, options.log, algo_params.stopping_rules)

    @infiltrate algo_params.infiltrate_state in [:all]

    return DynamicSDDiP.IterationResult(
        bound,
        model.ext[:best_objective],
        model.ext[:best_point],
        forward_trajectory.scenario_path,
        has_converged,
        status,
        cuts,
    )
end
