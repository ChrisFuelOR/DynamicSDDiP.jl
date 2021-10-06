# The functions
# > "solve",
# > "master_loop"
# > "inner_loop"
# are derived from similar named functions (solve, master_loop) in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Solves the `model`. In contrast to SDDP.jl, here in algo_params specific parameters
configuring the solution procedure are given to the function.
"""
function solve(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers;
    ############################################################################
    iteration_limit::Union{Int,Nothing} = nothing,
    time_limit::Union{Real,Nothing} = nothing,
    print_level::Int = 1,
    log_file::String = "DynamicSDDiP.log",
    log_frequency::Int = 1,
    run_numerical_stability_report::Bool = true,
    stopping_rules = SDDP.AbstractStoppingRule[],
    risk_measure = SDDP.Expectation(),
    sampling_scheme = SDDP.InSampleMonteCarlo(),
    cut_type = SDDP.SINGLE_CUT,
    cycle_discretization_delta::Float64 = 0.0,
    refine_at_similar_nodes::Bool = true,
    cut_deletion_minimum::Int = 1,
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme = SDDP.CompleteSampler(),
    dashboard::Bool = false,
    parallel_scheme::SDDP.AbstractParallelScheme = SDDP.Serial(),
    forward_pass::SDDP.AbstractForwardPass = SDDP.DefaultForwardPass(),
)

    # INITIALIZATION (AS IN SDDP)
    ############################################################################
    # Reset the TimerOutput and define logging-specific variables
    TimerOutputs.reset_timer!(DynamicSDDiP_TIMER)
    log_file_handle = open(log_file, "a")
    log = Log[]

    if print_level > 0
        print_helper(print_banner, log_file_handle)
    end

    if print_level > 1
        print_helper(print_parameters, log_file_handle, algo_params, applied_solvers)
    end

    # TODO: Add run_numerical_stability_report as in SDDP?

    # Convert the vector to an AbstractStoppingRule. Otherwise if the user gives
    # something like stopping_rules = [SDDP.IterationLimit(100)], the vector
    # will be concretely typed and we can't add a TimeLimit.
    stopping_rules = convert(Vector{SDDP.AbstractStoppingRule}, stopping_rules)
    # Add the limits as stopping rules. An IterationLimit or TimeLimit may
    # already exist in stopping_rules, but that doesn't matter.
    if iteration_limit !== nothing
        push!(stopping_rules, SDDP.IterationLimit(iteration_limit))
    end
    if time_limit !== nothing
        push!(stopping_rules, SDDP.TimeLimit(time_limit))
    end
    if length(stopping_rules) == 0
        @warn(
            "You haven't specified a stopping rule! You can only terminate " *
            "the call to DynamicSDDiP.solve via a keyboard interrupt ([CTRL+C])."
        )
    end

    # Update the nodes with the selected cut type (SINGLE_CUT or MULTI_CUT)
    # and the cut deletion minimum.
    if cut_deletion_minimum < 0
        cut_deletion_minimum = typemax(Int)
    end
    for (key, node) in model.nodes
        node.bellman_function.cut_type = cut_type
        node.bellman_function.global_theta.cut_oracle.deletion_minimum =
            cut_deletion_minimum
        for oracle in node.bellman_function.local_thetas
            oracle.cut_oracle.deletion_minimum = cut_deletion_minimum
        end
    end

    dashboard_callback = if dashboard
        launch_dashboard()
    else
        (::Any, ::Any) -> nothing
    end

    SDDP_options = DynamicSDDiP.Options(
        model,
        model.initial_root_state,
        sampling_scheme,
        backward_sampling_scheme,
        risk_measure,
        cycle_discretization_delta,
        refine_at_similar_nodes,
        stopping_rules,
        dashboard_callback,
        print_level,
        time(),
        log,
        log_file_handle,
        log_frequency,
        forward_pass,
    )

    # MODEL START
    ############################################################################
    status = :not_solved
    try
        status = solve_DynamicSDDiP(parallel_scheme, model, sddp_options, algo_params, applied_solvers)
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

    # SET UP SOME STUFF
    ############################################################################
    for (node_index, children) in model.nodes
        node = model.nodes[node_index]

        # Set info for x_in (taking bounds, binary, integer info from previous stage's x_out)
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

        # also re-initialize the existing value function such that nonlinear cuts are used
        # fortunately, node.bellman_function requires no specific type
        node.bellman_function = initialize_bellman_function_nonconvex(bellman_function, model, node)
        node.bellman_function.cut_type = cut_type
        node.bellman_function.global_theta.cut_oracle.deletion_minimum = deletion_minimum
        for oracle in node.bellman_function.local_thetas
            oracle.cut_oracle.deletion_minimum = deletion_minimum
        end
    end

    @infiltrate algo_params.infiltrate_state == :all

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
    applied_solvers::DynamicSDDiP.AppliedSolvers) where {T}

    # INITIALIZE PARAMETERS REQUIRED FOR REFINEMENTS
    ############################################################################
    previous_solution = nothing
    previous_bound = nothing
    sigma_increased = false
    bound_check = false

    # INITIALIZE BEST KNOWN POINT AND OBJECTIVE VALUE
    ############################################################################
    model.ext[:best_objective] = model.objective_sense == JuMP.MOI.MIN_SENSE ? Inf : -Inf
    model.ext[:best_point] = Vector{Dict{Symbol,Float64}}()

    # ACTUAL LOOP
    ############################################################################
    while true
        # start an iteration
        TimerOutputs.@timeit DynamicSDDiP_TIMER "iterations" begin
            result = iteration(model, options, algo_params, applied_solvers, previous_solution, bound_check, sigma_increased)
        end

        # set previous_solution and previous_bound using iteration results
        previous_solution = result.current_sol
        previous_bound = result.lower_bound

        # logging
        log_iteration(options, options.log)

        # initialize sigma_increased
        sigma_increased = false

        # initialize bound_check
        bound_check = true
        @infiltrate algo_params.infiltrate_state in [:all, :sigma]

        # IF CONVERGENCE IS ACHIEVED, CHECK IF ACTUAL OR ONLY REGULARIZED PROBLEM IS SOLVED
        ########################################################################
        if result.has_converged
            TimerOutputs.@timeit DynamicSDDiP_TIMER "sigma_test" begin
                sigma_test_results = forward_sigma_test(model, options, algo_params, applied_solvers, result.scenario_path, options.forward_pass, sigma_increased)
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
                # return convergence status
                return result.status
            end

        # IF NO CONVERGENCE IS ACHIEVED, DO DIFFERENT CHECKS
        ########################################################################
        else
            # CHECK IF LOWER BOUND HAS IMPROVED
            ############################################################################
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
            ############################################################################
            if result.upper_bound - result.lower_bound < - algo_params.opt_atol * 0.1
                algo_params.sigma = algo_params.sigma * algo_params.sigma_factor
                sigma_increased = true
                previous_solution = nothing
                previous_bound = nothing
                bound_check = false
            else
                sigma_increased = false
            end

        end
    end
end
