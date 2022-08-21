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

const DynamicSDDiP_TIMER = TimerOutputs.TimerOutput()

"""
Solves the `model`. In contrast to SDDP.jl, all parameters configuring
the algorithm are given (and possibly pre-defined) in algo_params.
"""
function solve(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    problem_params::DynamicSDDiP.ProblemParams
    )

    ############################################################################
    # INITIALIZATION (SIMILAR TO SDDP.jl)
    ############################################################################
    # Reset the TimerOutput
    TimerOutputs.reset_timer!(DynamicSDDiP_TIMER)

    # Prepare logging
    log_file_handle = open(algo_params.log_file, "a")
    log = Log[]

    if algo_params.print_level > 0
        print_helper(print_banner, log_file_handle)
    end
    if algo_params.run_numerical_stability_report
        report = sprint(
            io -> SDDP.numerical_stability_report(
                io,
                model,
                print = algo_params.print_level > 0,
            ),
        )
        SDDP.print_helper(print, log_file_handle, report)
    end

    if algo_params.print_level > 1
        print_helper(print_parameters, log_file_handle, algo_params, applied_solvers, problem_params)
    end

    # Maybe add run_numerical_stability_report as in SDDP.jl later

    # Prepare stopping rules
    #---------------------------------------------------------------------------
    # Convert the vector to an AbstractStoppingRule. Otherwise if the user gives
    # something like stopping_rules = [SDDP.IterationLimit(100)], the vector
    # will be concretely typed and we can't add a TimeLimit.
    # stopping_rules = algo_params.stopping_rules
    # stopping_rules = convert(Vector{SDDP.AbstractStoppingRule}, stopping_rules)

    if length(algo_params.stopping_rules) == 0
        @warn(
            "You haven't specified a stopping rule! You can only terminate " *
            "the call to DynamicSDDiP.solve via a keyboard interrupt ([CTRL+C])."
        )
    end

    # Prepare binary_precision
    #---------------------------------------------------------------------------
    for cut_generation_regime in algo_params.cut_generation_regimes
        regime = cut_generation_regime.state_approximation_regime
        if isa(regime,DynamicSDDiP.BinaryApproximation) && isempty(regime.binary_precision)
            # If no binary_precision dict has been defined explicitly, it is
            # initialized as empty. Then, for each state take a default precision.
            for (name, state_comp) in model.nodes[1].states
                if JuMP.is_binary(state_comp.out) || JuMP.is_integer(state_comp.out)
                    regime.binary_precision[name] = 1.0
                else
                    ub = JuMP.upper_bound(state_comp.out)
                    lb = 0.0 # all states are assumed to satisfy non-negativity constraints
                    regime.binary_precision[name] = (ub-lb)/127.0
                end
            end
        end
    end

    # Prepare sigma
    #---------------------------------------------------------------------------
    regime = algo_params.regularization_regime
    if isa(regime,DynamicSDDiP.Regularization) && isempty(regime.sigma)
        for (node_index, _) in model.nodes
            if node_index == 1
                # first stage requires no regularization
                push!(regime.sigma, 0.0)
            else
                push!(regime.sigma, 1.0)
            end
        end
    end

    # Prepare options for logging
    #---------------------------------------------------------------------------
    options = DynamicSDDiP.Options(
        model,
        model.initial_root_state,
        algo_params.risk_measure,
        time(),
        log,
        log_file_handle
    )

    # Prepare dictionary to log status of Lagrangian dual solution
    #---------------------------------------------------------------------------
    model.ext[:lag_status_dict] = Dict([(:opt, 0), (:conv, 0), (:sub, 0), (:iter, 0), (:unbounded, 0), (:issues, 0), (:mn_opt, 0), (:mn_iter, 0), (:mn_issue, 0), (:subgr_stalling, 0)])

    ############################################################################
    # RE-INITIALIZE THE EXISTING VALUE FUNCTION AND PREPARE CUT SELECTION
    ############################################################################
    # Update the nodes with the selected cut type (SINGLE_CUT or MULTI_CUT)
    # and the cut deletion minimum.
    regime = algo_params.cut_selection_regime
    if regime == DynamicSDDiP.CutSelection && regime.cut_deletion_minimum < 0
        algo_params.cut_selection_regime.cut_deletion_minimum = typemax(Int)
    end

    # Change to nonconvex Bellman function to use non-convex cuts
    #---------------------------------------------------------------------------
    # fortunately, node.bellman_function requires no specific type
    for (key, node) in model.nodes
        # Set solver
        JuMP.set_optimizer(node.subproblem, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4,"TimeLimit"=>300
        ))
        JuMP.set_silent(node.subproblem)

        if key != model.root_node

            # Determine correct objective sense and bounds for bellman function
            if model.objective_sense == MOI.MIN_SENSE
                lower_bound = JuMP.lower_bound(node.bellman_function.global_theta.theta)
                upper_bound = Inf
            elseif model.objective_sense == MOI.MAX_SENSE
                upper_bound = JuMP.upper_bound(node.bellman_function.global_theta.theta)
                lower_bound = -Inf
            end

            # Initialize bellman function
            bellman_function = BellmanFunction(lower_bound = lower_bound, upper_bound = upper_bound)
            node.bellman_function = DynamicSDDiP.initialize_bellman_function(bellman_function, model, node)
            node.bellman_function.cut_type = algo_params.cut_type

            # Prepare multi-cut case
            if node.bellman_function.cut_type == SDDP.MULTI_CUT
                # Count number of local thetas (number of sets of dual variables later)
                counter = 0
                for child in node.children
                    if isapprox(child.probability, 0.0, atol = 1e-6)
                        continue
                    end
                    child_node = model[child.term]
                    for noise in SDDP.sample_backward_noise_terms(algo_params.backward_sampling_scheme, child_node)
                        counter = counter + 1
                    end
                end
                _add_locals_if_necessary(node, node.bellman_function, counter)
            end

            if isa(algo_params.cut_selection_regime, DynamicSDDiP.CutSelection)
                node.bellman_function.global_theta.deletion_minimum =
                    algo_params.cut_selection_regime.cut_deletion_minimum
                for oracle in node.bellman_function.local_thetas
                    oracle.deletion_minimum = algo_params.cut_selection_regime.cut_deletion_minimum
                end
            end

            # Set-up counter for Benders cuts (Chen & Luedtke approach)
            node.ext[:Benders_cuts_original] = Tuple{Int64, Symbol}[]
            node.ext[:Benders_cuts_binary] = Tuple{Int64, Symbol}[]

            # Identify constraints containing state variables and store them in node.ext[:state_constraints]
            identify_state_constraints!(node)

        end
    end

    ############################################################################
    # MODEL START
    ############################################################################
    status = :not_solved
    try
        status = solve_DynamicSDDiP(algo_params.parallel_scheme, model, options, algo_params, applied_solvers)
    catch ex
        if isa(ex, InterruptException)
            status = :interrupted
            interrupt(algo_params.parallel_scheme)
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
    if algo_params.print_level > 0
        print_helper(print_footer, log_file_handle, results)
        if algo_params.print_level > 1
            # print_helper(TimerOutputs.print_timer, log_file_handle, DynamicSDDiP_TIMER)
            TimerOutputs.print_timer(log_file_handle, DynamicSDDiP_TIMER, allocations=true)
            TimerOutputs.print_timer(stdout, DynamicSDDiP_TIMER, allocations=true)
        end

        print_helper(print_lag_status, log_file_handle, model.ext[:lag_status_dict])
    end
    close(log_file_handle)
    return
end

"""
Solves the `model` using DynamicSDDiP in a serial scheme.
"""

function solve_DynamicSDDiP(
    parallel_scheme::SDDP.Serial,
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers
    ) where {T}

    ############################################################################
    # SET UP STATE VARIABLE INFORMATION
    ############################################################################
    for (node_index, children) in model.nodes
        node = model.nodes[node_index]

        # Set info for state.in as for previous stage's state.out
        #-----------------------------------------------------------------------
        if node_index > 1
            for (i, (name, state)) in enumerate(node.states)
                state_out_previous_stage = model.nodes[node_index-1].states[name].out
                state_in = state.in

                set_up_state_in_info!(state_out_previous_stage, state_in)

            end
        end

        # Store info for all states (state.in, state.out) for later
        # This is required for variable fixing and unfixing
        #-----------------------------------------------------------------------
        node.ext[:state_info_storage] = Dict{Symbol,DynamicSDDiP.StateInfoStorage}()

        for (i, (name, state)) in enumerate(node.states)
                variable_info_in = get_variable_info(state.in)
                variable_info_out = get_variable_info(state.out)
                node.ext[:state_info_storage][name] = DynamicSDDiP.StateInfoStorage(variable_info_in, variable_info_out)
        end

    end

    Infiltrator.@infiltrate algo_params.infiltrate_state == :all

    ############################################################################
    # LOG ITERATION HEADER
    ############################################################################
    if algo_params.print_level > 0
        print_helper(io -> println(io, "Solver: ", parallel_scheme, "\n"), options.log_file_handle)
        print_helper(print_iteration_header, options.log_file_handle)
    end

    ############################################################################
    # CALL ACTUAL SOLUTION PROCEDURE
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "loop" begin
        status = master_loop(parallel_scheme, model, options, algo_params,
            applied_solvers)
    end
    return status

end


"""
Loop function of DynamicSDDiP.
"""

function master_loop(
    parallel_scheme::SDDP.Serial,
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers
    ) where {T}

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
        sigma_increased = false
        bound_check = true

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :sigma]

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

        if result.has_converged
            return result.status
        end

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

function convergence_handler(
    result::DynamicSDDiP.IterationResult,
    model::SDDP.PolicyGraph{T},
    options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    sigma_increased::Bool,
    bound_check::Bool,
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    previous_bound::Union{Float64,Nothing},
    regularization_regime::DynamicSDDiP.Regularization
    ) where {T}

    ############################################################################
    # IF CONVERGENCE IS ACHIEVED, COMPARE TRUE AND REGULARIZED PROBLEM
    ############################################################################
    if result.has_converged
        TimerOutputs.@timeit DynamicSDDiP_TIMER "sigma_test" begin
            sigma_test_results = forward_sigma_test(model, options, algo_params, applied_solvers, result.scenario_path, sigma_increased)
        end
        sigma_increased = sigma_test_results.sigma_increased
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :sigma]

        if sigma_increased
            # reset previous values, as model is changed and convergence not achieved
            previous_solution = nothing
            previous_bound = nothing
            # binary refinement only when no sigma refinement has been made
            bound_check = false
            # no convergence
            result.has_converged = false
        else
            ####################################################################
            # THE ALGORITHM WILL TERMINATE
            ####################################################################
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
        # Check if deterministic stopping is used
        deterministic_stopping = false
        for stopping_rule in algo_params.stopping_rules
            if isa(stopping_rule, DynamicSDDiP.DeterministicStopping)
                deterministic_stopping = true
            end
        end

        if result.upper_bound - result.lower_bound < - 1e-8 && deterministic_stopping
                regularization_regime.sigma = regularization_regime.sigma * regularization_regime.sigma_factor
                sigma_increased = true
                previous_solution = nothing
                previous_bound = nothing
                bound_check = false
        else
            sigma_increased = false
            # Update of previous_solution and previous_bound
            previous_solution = result.current_sol
            previous_bound = result.lower_bound
        end

    end

    return (
        sigma_increased = sigma_increased,
        previous_solution = previous_solution,
        previous_bound = previous_bound,
        bound_check = bound_check,
    )
end


"""
Convergence handler if regularization is not used.
"""

function convergence_handler(
    result::DynamicSDDiP.IterationResult,
    model::SDDP.PolicyGraph{T}, options::DynamicSDDiP.Options,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    sigma_increased::Bool, bound_check::Bool,
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    previous_bound::Union{Float64,Nothing},
    regularization_regime::DynamicSDDiP.NoRegularization
    ) where {T}

    ############################################################################
    # IF CONVERGENCE IS ACHIEVED, COMPARE TRUE AND REGULARIZED PROBLEM
    ############################################################################
    if result.has_converged
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :sigma]
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
        # Check if deterministic stopping is used
        deterministic_stopping = false
        for stopping_rule in algo_params.stopping_rules
            if isa(stopping_rule, DynamicSDDiP.DeterministicStopping)
                deterministic_stopping = true
            end
        end

        if result.upper_bound - result.lower_bound < - 1e-8 && deterministic_stopping
            error("LB < UB for DynamicSDDiP. Terminating.")
        else
            # Update previous_solution and previous_bound
            previous_solution = result.current_sol
            previous_bound = result.lower_bound
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
    solution_check = true
    binary_refinement = :none

    # TimerOutputs.@timeit DynamicSDDiP_TIMER "bin_refinement" begin
    #     binary_refinement = DynamicSDDiP.binary_refinement(
    #         model,
    #         previous_solution,
    #         forward_trajectory.sampled_states,
    #         algo_params,
    #         solution_check,
    #         binary_refinement,
    #         bound_check
    #     )
    # end

    # bound_check = true
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # BACKWARD PASS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "backward_pass" begin
        DynamicSDDiP.backward_pass(
            model,
            options,
            algo_params,
            applied_solvers,
            forward_trajectory.scenario_path,
            forward_trajectory.sampled_states,
            forward_trajectory.epi_states,
            # forward_trajectory.objective_states,
            # forward_trajectory.belief_states,
        )
    end

    ############################################################################
    # CALCULATE LOWER BOUND
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "calculate_bound" begin
        first_stage_results = calculate_bound(model)
    end
    bound = first_stage_results.bound

    #Infiltrator.@infiltrate

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
             model.ext[:agg_lag_iterations],
             model.ext[:lag_iterations],
             model.ext[:total_cuts],
             model.ext[:active_cuts],
             model.ext[:total_solves],
         ),
     )

    ############################################################################
    # CHECK IF THE LOOP CONVERGED YET
    ############################################################################
    has_converged, status = convergence_test(model, options.log, algo_params.stopping_rules)

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    return DynamicSDDiP.IterationResult(
        bound,
        model.ext[:best_objective],
        model.ext[:best_point],
        forward_trajectory.scenario_path,
        has_converged,
        status,
    )
end
