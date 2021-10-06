# The functions
# > "print_helper",
# > "print_banner",
# > "print_iteration_header",
# > "print_iteration",
# > "print_footer",
# > "log_iteration",
# > "write_log_to_csv",
# and structs
# > "Log"
# > "Options"
# > "Results"
# are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################
struct Log
    iteration::Int
    lower_bound::Float64
    best_upper_bound::Float64
    current_upper_bound::Float64
    current_state::Vector{Dict{Symbol,Float64}}
    time::Float64
    #total_solves::Int
    #sigma::Vector{Float64}
    #binaryPrecision::Dict{Symbol,Float64}
    sigma_increased::Union{Bool,Nothing}
    binary_refinement::Union{Symbol,Nothing}
    subproblem_size::Union{Dict{Symbol,Int64},Nothing}
    algo_params::DynamicSDDiP.AlgoParams
    lag_iterations::Union{Vector{Int},Nothing}
    lag_status::Union{Vector{Symbol},Nothing}
    total_cuts::Int
    active_cuts::Int
end


# Internal struct: storage for SDDP options and cached data. Users shouldn't
# interact with this directly.
struct Options{T}
    # The initial state to start from the root node.
    initial_state::Dict{Symbol,Float64}
    # The sampling scheme to use on the forward pass.
    sampling_scheme::SDDP.AbstractSamplingScheme
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme
    # Storage for the set of possible sampling states at each node. We only use
    # this if there is a cycle in the policy graph.
    starting_states::Dict{T,Vector{Dict{Symbol,Float64}}}
    # Risk measure to use at each node.
    risk_measures::Dict{T,SDDP.AbstractRiskMeasure}
    # The delta by which to check if a state is close to a previously sampled
    # state.
    cycle_discretization_delta::Float64
    # Flag to add cuts to similar nodes.
    refine_at_similar_nodes::Bool
    # The node transition matrix.
    Φ::Dict{Tuple{T,T},Float64}
    # A list of nodes that contain a subset of the children of node i.
    similar_children::Dict{T,Vector{T}}
    stopping_rules::Vector{SDDP.AbstractStoppingRule}
    dashboard_callback::Function
    print_level::Int
    start_time::Float64
    log::Vector{DynamicSDDiP.Log}
    log_file_handle
    log_frequency::Int
    forward_pass::SDDP.AbstractForwardPass

    # Internal function: users should never construct this themselves.
    function Options(
        model::SDDP.PolicyGraph{T},
        initial_state::Dict{Symbol,Float64},
        sampling_scheme::SDDP.AbstractSamplingScheme,
        backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
        risk_measures,
        cycle_discretization_delta::Float64,
        refine_at_similar_nodes::Bool,
        stopping_rules::Vector{SDDP.AbstractStoppingRule},
        dashboard_callback::Function,
        print_level::Int,
        start_time::Float64,
        log::Vector{DynamicSDDiP.Log},
        log_file_handle,
        log_frequency::Int,
        forward_pass::SDDP.AbstractForwardPass,
    ) where {T}
        return new{T}(
            initial_state,
            sampling_scheme,
            backward_sampling_scheme,
            SDDP.to_nodal_form(model, x -> Dict{Symbol,Float64}[]),
            SDDP.to_nodal_form(model, risk_measures),
            cycle_discretization_delta,
            refine_at_similar_nodes,
            SDDP.build_Φ(model),
            SDDP.get_same_children(model),
            stopping_rules,
            dashboard_callback,
            print_level,
            start_time,
            log,
            log_file_handle,
            log_frequency,
            forward_pass,
        )
    end
end

struct Results
    status::Symbol
    log::Vector{DynamicSDDiP.Log}
end

function print_helper(f, io, args...)
    f(stdout, args...)
    f(io, args...)
end

function print_banner(io)
    println(
        io,
        "--------------------------------------------------------------------------------",
    )
    println(io, "                      DynamicSDDiP.jl (c) Christian Füllner, 2021")
    println(io, "re-uses code from     SDDP.jl (c) Oscar Dowson, 2017-20")
    println(io)
    flush(io)
end

function print_parameters(io, algo_params::DynamicSDDiP.AlgoParams, applied_solvers::DynamicSDDiP.AppliedSolvers)

    # Printing the time
    println(io, Dates.now())

    # Printint the file name
    print(io, "calling ")
    print(io, @__DIR__)
    println(io)
    println(io)

    # Printing the parameters used
    println(io, Printf.@sprintf("opt_rtol: %1.4e", algo_params.opt_rtol))
    println(io, Printf.@sprintf("opt_atol: %1.4e", algo_params.opt_atol))
    println(io, Printf.@sprintf("iteration_limit: %5d", algo_params.iteration_limit))
    println(io, Printf.@sprintf("time_limit (sec): %6d", algo_params.time_limit))
    println(io, "------------------------------------------------------------------------")
    print(io, "Binary approximation used: ")
    println(io, algo_params.binary_approx)
    print(io, "Initial binary precision: ")
    println(io, algo_params.binary_precision)
    println(io, "------------------------------------------------------------------------")
    print(io, "Regularization used: ")
    println(io, algo_params.regularization)
    println(io, Printf.@sprintf("Initial sigma: %4.1e", algo_params.sigma))
    println(io, Printf.@sprintf("Sigma increase factor: %4.1e", algo_params.sigma:factor))
    println(io, "------------------------------------------------------------------------")
    println(io, Printf.@sprintf("Lagrangian rtol: %1.4e", algo_params.lagrangian_rtol))
    println(io, Printf.@sprintf("Lagrangian atol: %1.4e", algo_params.lagrangian_atol))
    println(io, Printf.@sprintf("iteration_limit: %5d", algo_params.iteration_limit))
    print(io, "Dual initialization: ")
    println(io, algo_params.dual_initialization_method)
    print(io, "Dual solution method: ")
    println(io, algo_params.dual_solution_method)
    print(io, "Dual status regime: ")
    println(io, algo_params.dual_status_regime)
    print(io, "Dual bound regime: ")
    println(io, algo_params.dual_bound_regime)
    print(io, "Cut type: ")
    println(io, algo_params.cut_type)
    print(io, "Magnanti and Wong approach used: ")
    println(io, algo_params.magnanti_wong)
    print(io, "Numerical focus used: ")
    println(io, algo_params.numerical_focus)
    println(io, "------------------------------------------------------------------------")
    if algo_params.dual_solution_method == :level_bundle
        println(io, Printf.@sprintf("Level parameter: %2.4e", algo_params.bundle_params.level_factor))
        println(io, Printf.@sprintf("Bundle alpha: %2.4e", algo_params.bundle_params.bundle_alpha))
        println(io, Printf.@sprintf("Bundle factor: %2.4e", algo_params.bundle_params.bundle_factor))
        println(io, "------------------------------------------------------------------------")
    end
    print(io, "Cut projection method: ")
    println(io, algo_params.cut_projection_method)
    print(io, "Cut selection used: ")
    println(io, algo_params.cut_selection)
    println(io, "------------------------------------------------------------------------")
    print(io, Printf.@sprintf("LP solver: %15s", applied_solvers.LP))
    print(io, Printf.@sprintf("MILP solver: %15s", applied_solvers.MILP))
    print(io, Printf.@sprintf("(MI)NLP solver: %15s", applied_solvers.NLP))
    print(io, Printf.@sprintf("Lagrange solver: %15s", applied_solvers.Lagrange))
    println(io, "------------------------------------------------------------------------")

    flush(io)
end


function print_iteration_header(io)
    println(
        io,
        " Inner_Iteration   Upper Bound    Best Upper Bound     Lower Bound      Gap       Time (s)         sigma_ref    bin_ref     tot_var     bin_var     int_var       con       cuts   active     Lag iterations & status     ",
    )
    flush(io)
end

function print_iteration(io, log::Log)
    print(io, lpad(Printf.@sprintf("%5d", log.iteration), 15))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.upper_bound), 13))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.best_upper_bound), 16))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.lower_bound), 13))
    print(io, "   ")

    gap = abs(log.best_upper_bound - log.lower_bound)/max(log.best_upper_bound, log.lower_bound)
    print(io, lpad(Printf.@sprintf("%3.4f", gap), 8))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.time), 13))
    print(io, "   ")
    if !isnothing(log.sigma_increased)
    	print(io, Printf.@sprintf("%9s", log.sigma_increased ? "true" : "false"))
    else
   	    print(io, lpad(Printf.@sprintf(""), 9))
    end
    print(io, "   ")
    if !isnothing(log.binary_refinement)
        print(io, Printf.@sprintf("%9s", log.binary_refinement))
    else
   	    print(io, lpad(Printf.@sprintf(""), 9))
    end
    print(io, "   ")
    if !isnothing(log.subproblem_size)
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:total_var]))
        print(io, "   ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:bin_var]))
        print(io, "   ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:int_var]))
        print(io, "   ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:total_con]))
    else
        print(io, lpad(Printf.@sprintf(""), 45))
    end
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%5d", log.total_cuts), 7))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%5d", log.active_cuts), 7))
    print(io, "     ")

    if !isnothing(log.lag_iterations)
        print(io, log.lag_iterations)
    else
        print(io, lpad(Printf.@sprintf(""), 19))
    end
    print(io, "   ")
    if !isnothing(log.lag_status)
        print(io, log.lag_status)
    else
        print(io, lpad(Printf.@sprintf(""), 19))
    end
    print(io, "   ")

    println(io)
    flush(io)
end


function print_footer(io, training_results)
    println(io, "\nTerminating DynamicSDDiP with status: $(training_results.status)")
    println(
        io,
        "------------------------------------------------------------------------------",
    )
    flush(io)
end

function log_iteration(options, log)
    options.dashboard_callback(log[end], false)
    if options.print_level > 0 && mod(length(log), options.log_frequency) == 0
        print_helper(print_iteration, options.log_file_handle, log[end])
    end
end


"""
    write_log_to_csv(model::PolicyGraph, filename::String)

Write the log of the most recent training to a csv for post-analysis.

Assumes that the model has been trained via [`DynamicSDDiP.solve`](@ref).
"""
function write_log_to_csv(model::SDDP.PolicyGraph, filename::String, algo_params::DynamicSDDiP.AlgoParams)
    # TODO
end


function print_lagrange_header(io)
    println(
        io,
        " Iteration     f_approx    best_actual     f_actual  ",
    )
    flush(io)
end

function print_lag_iteration(io, iter::Int, f_approx::Float64, best_actual::Float64, f_actual::Float64)
    print(io, lpad(Printf.@sprintf("%5d", iter), 15))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.10e", f_approx), 13))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.10e", best_actual), 13))
    print(io, "   ")
    print(io, lpad(Printf.@sprintf("%1.10e", f_actual), 13))

    println(io)
    flush(io)
end
