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
    sigma_increased::Union{Bool,Nothing}
    binary_refinement::Union{Symbol,Nothing}
    subproblem_size::Union{Dict{Symbol,Int64},Nothing}
    algo_params::DynamicSDDiP.AlgoParams
    agg_lag_iterations::Int64
    lag_iterations::Union{Vector{Int},Vector{Float64},Nothing}
    lag_status::Union{Vector{String},Nothing}
    total_cuts::Int
    active_cuts::Int
    total_solves::Int
end

# Internal struct: storage for SDDP options and cached data. Users shouldn't
# interact with this directly.
struct Options{T}
    # The initial state to start from the root node.
    initial_state::Dict{Symbol,Float64}
    # Storage for the set of possible sampling states at each node. We only use
    # this if there is a cycle in the policy graph.
    starting_states::Dict{T,Vector{Dict{Symbol,Float64}}}
    # Risk measure to use at each node.
    risk_measures::Dict{T,SDDP.AbstractRiskMeasure}
    # The node transition matrix.
    Φ::Dict{Tuple{T,T},Float64}
    # A list of nodes that contain a subset of the children of node i.
    similar_children::Dict{T,Vector{T}}
    start_time::Float64
    log::Vector{DynamicSDDiP.Log}
    log_file_handle::Any

    # Internal function: users should never construct this themselves.
    function Options(
        model::SDDP.PolicyGraph{T},
        initial_state::Dict{Symbol,Float64},
        risk_measures,
        start_time::Float64,
        log::Vector{DynamicSDDiP.Log},
        log_file_handle::Any
    ) where {T}
        return new{T}(
            initial_state,
            SDDP.to_nodal_form(model, x -> Dict{Symbol,Float64}[]),
            SDDP.to_nodal_form(model, risk_measures),
            SDDP.build_Φ(model),
            SDDP.get_same_children(model),
            start_time,
            log,
            log_file_handle
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
    println(io)
    println(io)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io, "DynamicSDDiP.jl (c) Christian Füllner, 2021")
    println(io, "re-uses code from SDDP.jl (c) Oscar Dowson, 2017-21")
    flush(io)
end

function print_parameters(io, algo_params::DynamicSDDiP.AlgoParams, applied_solvers::DynamicSDDiP.AppliedSolvers, problem_params::DynamicSDDiP.ProblemParams)

    # Printint the file name
    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "PATH")
    println(io, "calling ")
    println(io, @__DIR__)
    println(io, Base.source_path())

    # Printing the time
    println(io, "DATETIME")
    println(io, Dates.now())

    ############################################################################
    # PRINTING THE PARAMETERS USED
    ############################################################################
    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "STOPPING RULES")
    if isempty(algo_params.stopping_rules)
        println(io, "No stopping rule defined.")
    else
        for stopping_rule in algo_params.stopping_rules
            println(io, stopping_rule)
        end
    end

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "REGULARIZATION")
    println(io, algo_params.regularization_regime)

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "CUT GENERATION REGIMES")
    for regime in algo_params.cut_generation_regimes
        println(io, regime)
        println(io)
    end

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "CUT AGGREGATION")
    println(io, algo_params.cut_aggregation_regime)

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "CUT SELECTION")
    println(io, algo_params.cut_selection_regime)

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "APPLIED SOLVERS (LP, MILP, MIQCP, MINLP, NLP, Lagrange)")
    println(io, applied_solvers)
    println(io, algo_params.solver_approach)
    println(io, algo_params.numerical_focus)
    println(io, algo_params.silent)

    if !isnothing(algo_params.seed)
        println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
        println(io, "SAMPLING")
        println(io, "Used seed for sampling scenarios: ")
        print(io, rpad(Printf.@sprintf("%s", algo_params.seed), 10))
        println(io)
    end

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "RUN DESCRIPTION")
    println(io, algo_params.run_description)

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "PROBLEM DESCRIPTION")
    println(io, "Number of stages: ", problem_params.number_of_stages)
    println(io, "Number of realizations per stage: ", problem_params.number_of_realizations)

    if !isnothing(problem_params.tree_seed)
            println(io, "Seed for scenario tree sampling: ", problem_params.tree_seed)
    end

    flush(io)
end


function print_iteration_header(io)

    rule = "─"
    rule_length = 200

    #total_table_width = sum(textwidth.((sec_ncalls, time_headers, alloc_headers))) + 3
    printstyled(io, "", rule^rule_length, "\n"; bold=true)

    header = "It.#"
    print(io, rpad(Printf.@sprintf("%s", header), 6))
    print(io, "  ")
    header = "UB"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "Best UB"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "LB"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "Gap"
    print(io, lpad(Printf.@sprintf("%s", header), 8))
    print(io, "  ")
    header = "Time"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "It_Time"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "Refinements"
    print(io, lpad(Printf.@sprintf("%s", header), 20))
    print(io, "  ")
    header = "#Var."
    print(io, lpad(Printf.@sprintf("%s", header), 31))
    print(io, "  ")
    header = "#Constr."
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "#Cuts"
    print(io, lpad(Printf.@sprintf("%s", header), 16))
    print(io, "       ")
    header = "Lagrangian Dual"
    print(io, rpad(Printf.@sprintf("%s", header), 40))
    print(io, "  ")

    println(io)

    header = ""
    print(io, rpad(Printf.@sprintf("%s", header), 53))
    header = "[%]"
    print(io, lpad(Printf.@sprintf("%s", header), 8))
    print(io, "  ")
    header = "[s]"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "[s]"
    print(io, lpad(Printf.@sprintf("%s", header), 13))
    print(io, "  ")
    header = "σ"
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "Bin."
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "Total"
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "{0,1}"
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "ℤ"
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "Total"
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")
    header = "Total"
    print(io, lpad(Printf.@sprintf("%s", header), 7))
    print(io, "  ")
    header = "Active"
    print(io, lpad(Printf.@sprintf("%s", header), 7))
    print(io, "  ")
    header = "# It."
    print(io, lpad(Printf.@sprintf("%s", header), 9))
    print(io, "  ")

    println(io)

    printstyled(io, "", rule^rule_length, "\n"; bold=true)

    flush(io)
end

function print_iteration(io, log::Log, start_time::Float64)
    print(io, rpad(Printf.@sprintf("%-5d", log.iteration), 6))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.current_upper_bound), 13))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.best_upper_bound), 13))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.lower_bound), 13))
    print(io, "  ")

    gap = abs(log.best_upper_bound - log.lower_bound)/(max(log.best_upper_bound, log.lower_bound + 1e-10))

    print(io, lpad(Printf.@sprintf("%3.4f", gap), 8))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.time), 13))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%1.6e", log.time - start_time), 13))
    print(io, "  ")
    if !isnothing(log.sigma_increased)
    	print(io, Printf.@sprintf("%9s", log.sigma_increased ? "true" : "false"))
    else
   	    print(io, lpad(Printf.@sprintf(""), 9))
    end
    print(io, "  ")
    if !isnothing(log.binary_refinement)
        print(io, Printf.@sprintf("%9s", log.binary_refinement))
    else
   	    print(io, lpad(Printf.@sprintf(""), 9))
    end
    print(io, "  ")
    if !isnothing(log.subproblem_size)
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:total_var]))
        print(io, "  ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:bin_var]))
        print(io, "  ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:int_var]))
        print(io, "  ")
       	print(io, Printf.@sprintf("%9d", log.subproblem_size[:total_con]))
    else
        print(io, lpad(Printf.@sprintf(""), 45))
    end
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%5d", log.total_cuts), 7))
    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%5d", log.active_cuts), 7))

    print(io, "  ")
    print(io, lpad(Printf.@sprintf("%5d", log.agg_lag_iterations), 9))
    #print(io, "       ")

    # if !isnothing(log.lag_iterations)
    #     print(io, log.lag_iterations)
    # else
    #     print(io, lpad(Printf.@sprintf(""), 19))
    # end
    # print(io, "  ")
    # if !isnothing(log.lag_status)
    #     print(io, log.lag_status)
    # else
    #     print(io, lpad(Printf.@sprintf(""), 19))
    # end
    # print(io, "  ")

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

function log_iteration(algo_params::DynamicSDDiP.AlgoParams, log_file_handle::Any, log::Vector{DynamicSDDiP.Log})
    if algo_params.print_level > 0 && mod(length(log), algo_params.log_frequency) == 0
        # Get time() after last iteration to compute iteration specific time
        if lastindex(log) > 1
            start_time = log[end-1].time
        else
            start_time = 0.0
        end

        print_helper(print_iteration, log_file_handle, log[end], start_time)
    end
end


"""
    write_log_to_csv(model::PolicyGraph, filename::String)

Write the log of the most recent training to a csv for post-analysis.

Assumes that the model has been trained via [`DynamicSDDiP.solve`](@ref).
"""
function write_log_to_csv(model::SDDP.PolicyGraph, filename::String, algo_params::DynamicSDDiP.AlgoParams)
    # TO-DO
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

function print_simulation(io, algo_params::DynamicSDDiP.AlgoParams, μ::Float64, ci::Float64, lower_bound::Float64)

    println(io)
    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "SIMULATION RESULTS")
    println(io, algo_params.simulation_regime)
    println(io, "Lower bound: ", lower_bound)
    println(io, "Statistical upper bound (confidence interval): ", μ, " ± ", ci )
    println(io, "Pessimistic upper bound: ", μ + ci )
    flush(io)
end

function print_det_equiv(io, problem_params::DynamicSDDiP.ProblemParams, value::Float64)

    println(io)
    println(io)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io,"#########################################################################################################################################",)
    println(io, "DynamicSDDiP.jl (c) Christian Füllner, 2021")
    println(io, "re-uses code from SDDP.jl (c) Oscar Dowson, 2017-21")
    flush(io)

    # Printint the file name
    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "PATH")
    println(io, "calling ")
    println(io, @__DIR__)
    println(io, Base.source_path())

    # Printing the time
    println(io, "DATETIME")
    println(io, Dates.now())

    println(io, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
    println(io, "PROBLEM DESCRIPTION")
    println(io, "Number of stages: ", problem_params.number_of_stages)
    println(io, "Number of realizations per stage: ", problem_params.number_of_realizations)

    if !isnothing(problem_params.tree_seed)
            println(io, "Seed for scenario tree sampling: ", problem_params.tree_seed)
    end

    println(io, "SOLVED DETERMINISTIC EQUIVALENT")
    println(io, "Optimal value: ", value)

    flush(io)
end
