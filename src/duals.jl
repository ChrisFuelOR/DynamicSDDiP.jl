# The function
# > "get_dual_solution",
# is derived from the same function in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

#*******************************************************************************
# LP DUAL
#*******************************************************************************

"""
Solving the dual problem to obtain cut information - using LP relaxation
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.LinearDuality,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    subproblem = node.subproblem

    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "LP_relaxation" begin
        dual_results = solve_LP_relaxation(node, subproblem, algo_params, cut_generation_regime, applied_solvers, DynamicSDDiP.LPDuals())
    end
 
    dual_vars = dual_results.dual_vars
    dual_obj = dual_results.dual_obj

    """ Note that we always initialize π₀ as well, even if we do not change it
    from 1.0 in the classical framework."""
    dual_0_var = 1.0

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    store_dual_values!(node, dual_values, dual_vars, bin_state, cut_generation_regime.state_approximation_regime)

    return (
        dual_values=dual_values,
        dual_0_var=dual_0_var,
        bin_state=bin_state,
        intercept=dual_obj,
        iterations=1,
        lag_status=:B,
        add_cut_flag=add_cut_flag,
    )
end

function solve_LP_relaxation(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_initalization_regime::DynamicSDDiP.LPDuals,
)

    # Get number of states and create zero vector for duals
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
    dual_vars_initial = zeros(number_of_states)

    # Create LP Relaxation
    undo_relax = JuMP.relax_integrality(node.subproblem)

    # Define appropriate solver
    reset_solver!(subproblem, algo_params, applied_solvers, :LP_relax, algo_params.solver_approach)

    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_LP_relax" begin
        JuMP.optimize!(subproblem)
    end

    # Solve LP Relaxation
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    dual_obj = JuMP.objective_value(subproblem)

    # Get dual values (reduced costs) for binary states as initial solution
    get_and_set_dual_values!(node, dual_vars_initial, cut_generation_regime.state_approximation_regime)

    # Note: due to JuMP's dual convention, we need to flip the sign for
    # maximization problems.
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1
    dual_vars_initial = dual_sign * dual_vars_initial

    # Undo relaxation
    undo_relax()

    return(
        dual_obj = dual_obj,
        dual_vars = dual_vars_initial,
    )
end


#*******************************************************************************
# STRENGTHENED DUAL
#*******************************************************************************

"""
Solving the dual problem to obtain cut information - using LP relaxation
and strengthening by Lagrangian relaxation
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.StrengthenedDuality,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    subproblem = node.subproblem

    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "LP_relaxation" begin
        dual_results = solve_LP_relaxation(node, subproblem, algo_params, cut_generation_regime, applied_solvers, DynamicSDDiP.LPDuals())
    end

    dual_vars = dual_results.dual_vars
    dual_obj = dual_results.dual_obj

    """ Note that we always initialize π₀ as well, even if we do not change it
    from 1.0 in the classical framework."""
    dual_0_var = 1.0

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # GET STRENGTHENING INFORMATION
    ############################################################################
    # solve lagrangian relaxed problem for these dual values
    if node.has_integrality
        dual_obj = _getStrengtheningInformation(node, dual_vars, algo_params, cut_generation_regime, applied_solvers)
    end

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    store_dual_values!(node, dual_values, dual_vars, bin_state, cut_generation_regime.state_approximation_regime)

    return (
        dual_values=dual_values,
        dual_0_var=dual_0_var,
        bin_state=bin_state,
        intercept=dual_obj,
        iterations=1,
        lag_status=:SB,
        add_cut_flag=add_cut_flag,
    )
end


#*******************************************************************************
# LAGRANGIAN DUAL
#*******************************************************************************

"""
Solving the dual problem to obtain cut information - using Lagrangian dual
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.LagrangianDuality,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    subproblem = node.subproblem

    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)

    # storages for information on Lagrangian dual
    lag_obj = 0
    lag_iterations = 0
    lag_status = :none

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "dual_initialization" begin
        dual_vars = initialize_duals(node, subproblem, algo_params, cut_generation_regime, applied_solvers, duality_regime.dual_initialization_regime)
    end

    """ Note that we always initialize π₀ as well, even if we do not change it
    from 1.0 in the classical framework."""
    dual_0_var = 1.0

    ############################################################################
    # GET PRIMAL SOLUTION TO BOUND LAGRANGIAN DUAL, TRACK CONVERGENCE AND DEBUGGING
    ############################################################################
    # REGULARIZE PROBLEM IF REGULARIZATION IS USED
    node.ext[:regularization_data] = Dict{Symbol,Any}()
    regularize_bw!(node, node_index, subproblem, cut_generation_regime, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)
    #regularize_bw!(node, node_index, subproblem, cut_generation_regime, DynamicSDDiP.NoRegularization(), cut_generation_regime.state_approximation_regime)

    # RESET SOLVER (as it may have been changed in between for some reason)
    if isa(algo_params.regularization_regime, DynamicSDDiP.NoRegularization)
        reset_solver!(subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)
    else
        reset_solver!(subproblem, algo_params, applied_solvers, :reg, algo_params.solver_approach)
    end

    # SOLVE PRIMAL PROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    primal_obj = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    # DEREGULARIZE PROBLEM IF REQUIRED
    deregularize_bw!(node, subproblem, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)
    #deregularize_bw!(node, subproblem, DynamicSDDiP.NoRegularization(), cut_generation_regime.state_approximation_regime)

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # GET BOUNDS FOR LAGRANGIAN DUAL
    ############################################################################
    # primal_obj = 1.3735 #augmented case #1.3735, 1.3
    bound_results = get_dual_bounds(node, node_index, algo_params, primal_obj, duality_regime.dual_bound_regime)
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    try
        ########################################################################
        # CALL SOLUTION METHOD
        ########################################################################
        # Solve dual and return a dict with the multiplier of the copy constraints.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
            results = solve_lagrangian_dual(
                    node,
                    node_index,
                    primal_obj,
                    dual_vars,
                    bound_results,
                    algo_params,
                    cut_generation_regime,
                    applied_solvers,
                    duality_regime.dual_solution_regime
                    )
        end

        lag_obj = results.lag_obj
        lag_iterations = results.iterations
        lag_status = results.lag_status

        #println(primal_obj, ", ", lag_obj, ", ", lag_iterations, ", ", lag_status)

        subproblem.ext[:sddp_policy_graph].ext[:agg_lag_iterations] += results.iterations

        # Counter to compare only number of iterations for converged cases
        if lag_status in (:opt, :conv, :sub, :iter, :mn_opt, :mn_iter)
            subproblem.ext[:sddp_policy_graph].ext[:corr_lag_iterations] += results.iterations
            subproblem.ext[:sddp_policy_graph].ext[:corr_realizations] += 1
        end

        ########################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ########################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(subproblem.ext[:sddp_policy_graph], lag_status, duality_regime.dual_status_regime)

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    catch e
        #SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    store_dual_values!(node, dual_values, dual_vars, bin_state, cut_generation_regime.state_approximation_regime)

    return (
        dual_values=dual_values,
        dual_0_var=dual_0_var,
        bin_state=bin_state,
        intercept=lag_obj,
        iterations=lag_iterations,
        lag_status=lag_status,
        add_cut_flag=add_cut_flag,
    )
end


#*******************************************************************************
# LAGRANGIAN DUAL - UNIFIED SETTING - MULTI-CUT
#*******************************************************************************

"""
Solving the dual problem to obtain cut information - using Lagrangian dual
in the unified framework with multi-cuts
"""
function get_dual_solution(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    add_cut_flag::Bool,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    subproblem = node.subproblem

    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)

    # storages for information on Lagrangian dual
    lag_obj = 0
    lag_iterations = 0
    lag_status = :none

    ############################################################################
    # INITIALIZE DUALS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "dual_initialization" begin
        dual_vars = initialize_duals(node, subproblem, algo_params, cut_generation_regime, applied_solvers, duality_regime.dual_initialization_regime)
    end

    # Initialize π₀ as 1 (0 is not suitable for relatively complete recourse)
    dual_0_var = 1.0

    ############################################################################
    # GET TRUE PRIMAL OBJECTIVE VALUE
    ############################################################################
    # RESET SOLVER (as it may have been changed in between for some reason)
    reset_solver!(subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)

    # SOLVE PRIMAL SUBPROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    primal_original_obj = JuMP.objective_value(subproblem)

    ############################################################################
    # GET CORE POINT RESP. NORMALIZATION COEFFICIENTS
    ############################################################################
    normalization_regime = cut_generation_regime.duality_regime.normalization_regime

    """ If we do not use a normalization that requires a core point, we simply
    return nothing """
    normalization_coeff = get_normalization_coefficients(
        node,
        number_of_states,
        epi_state,
        primal_original_obj,
        algo_params,
        applied_solvers,
        cut_generation_regime.state_approximation_regime,
        normalization_regime,
    )
  
    ############################################################################
    # GET PRIMAL SOLUTION AND ADDRESS UNBOUNDEDNESS
    ############################################################################
    """
    For reverse polar cuts (linear normalization) cuts, we solve approximations of the 
    primal of the normalized Lagrangian dual problem first.
    
    First, their feasibility/infeasibility helps us to detect if a core point 
    candidate could be bad and lead to an unbounded dual problem. 
    Then, we may introduce artificial bounds to the dual problem.
    
    Second, the primal_obj value can be used as an upper bound for the normalized
    Lagrangian dual problem. For RP cuts such bound is required, as the outer
    problem may be unbounded.

    For deep cuts solving the primal projection problem is not implemented
    so far, so primal_obj is a user-specified bound (or nothing).

    """

    # Initialize values
    node.ext[:primal_data] = Dict{Symbol,Any}()
    primal_unified_obj = Inf
    dual_multiplier_bound = nothing

    # Call primal solution method to possibly detect unboundedness
    unbounded_result = detect_unboundedness(node, epi_state, normalization_coeff, cut_generation_regime.duality_regime, algo_params, applied_solvers, cut_generation_regime, normalization_regime)
    
    # Take a measure to address unboundedness
    if unbounded_result.unbounded_flag
        #println(normalization_coeff.ω, round(normalization_coeff.ω₀, digits=2), ", ", round(primal_original_obj, digits=2), ", ", normalization_coeff.ω₀ >= primal_original_obj, ", ", epi_state)

        if isa(normalization_regime.unbounded_regime, DynamicSDDiP.Unbounded_Opt_SB)
            # Get strengthened Benders cut
            dual_results = get_dual_solution(node, node_index, i, epi_state, add_cut_flag, algo_params, cut_generation_regime, applied_solvers, duality_regime)

            return (
                dual_values = dual_results.dual_values,
                dual_0_var = dual_results.dual_0_var,
                bin_state = dual_results.bin_state,
                intercept = dual_results.lag_obj,
                iterations = dual_results.lag_iterations,
                lag_status = :unbounded,
                add_cut_flag = dual_results.add_cut_flag,
            )

        elseif isa(normalization_regime.unbounded_regime, DynamicSDDiP.Unbounded_Opt_Bound)
            primal_unified_obj = normalization_regime.unbounded_regime.user_dual_objective_bound
            dual_multiplier_bound = normalization_regime.unbounded_regime.user_dual_multiplier_bound
            # Prepare bounds for Lagrangian dual
            bound_results = get_dual_bounds(node, node_index, algo_params, primal_unified_obj, duality_regime.dual_bound_regime)
        end

    else
        # Prepare bounds for Lagrangian dual
        primal_unified_obj = unbounded_result.primal_unified_obj
        dual_multiplier_bound = Inf
        bound_results = get_dual_bounds(node, node_index, algo_params, primal_unified_obj, duality_regime.dual_bound_regime)
    end
      
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    try
        ########################################################################
        # CALL SOLUTION METHOD
        ########################################################################
        # Solve dual and return a dict with the multiplier of the copy constraints.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
            results = solve_unified_lagrangian_dual(
                        node,
                        node_index,
                        i,
                        epi_state,
                        normalization_coeff,
                        primal_unified_obj,
                        dual_vars,
                        dual_0_var,
                        bound_results,
                        dual_multiplier_bound,
                        algo_params,
                        cut_generation_regime,
                        applied_solvers,
                        duality_regime.dual_solution_regime
                        )
        end

        lag_obj = results.lag_obj
        lag_iterations = results.iterations
        lag_status = results.lag_status
        dual_0_var = results.dual_0_var

        #println(primal_original_obj, ", ", primal_unified_obj, ", ", lag_obj, ", ", lag_iterations, ", ", lag_status, ", ", dual_0_var)    
        #println(lag_iterations, ", ", lag_status, ", ", dual_0_var, ", ", sum(abs.(dual_vars)))
        #println(node_index, " ,", i, " ,", primal_unified_obj, " ,", lag_obj, " ,", lag_status, " ,", lag_iterations)

        subproblem.ext[:sddp_policy_graph].ext[:agg_lag_iterations] += results.iterations

        # Re-set lag status if we had to introduce artificial bounds due to unboundedness
        if unbounded_result.unbounded_flag
            lag_status = :unbounded
        end

        # Counter to compare only number of iterations for converged cases
        if lag_status in (:opt, :conv, :sub, :iter, :mn_opt, :mn_iter)
            subproblem.ext[:sddp_policy_graph].ext[:corr_lag_iterations] += results.iterations
            subproblem.ext[:sddp_policy_graph].ext[:corr_realizations] += 1
        end

        ########################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ########################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(subproblem.ext[:sddp_policy_graph], lag_status, duality_regime.dual_status_regime)

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    catch e
        #SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    if isapprox(lag_obj, 0.0, atol=1e-8)

        add_cut_flag = false

    elseif !isnothing(normalization_coeff) && all(normalization_coeff.ω .== 0.0) && isapprox(normalization_coeff.ω₀, 0.0, atol=1e-8)
        """ If the linear pseudonorm is used, but all coefficients are zero,
        then the Lagrangian dual is not normalized, but unbounded. Analogously,
        a zero function is optimized over the reverse polar set, which can yield
        any point in this unbounded set. Therefore, we are not guaranteed to
        obtain any meaningful cut.

        Note that we artificially bound the dual objective, so that
        we obtain valid, but possibly very large cut coefficients.
        However, experiments show that the cuts still tend to be invalid due
        to numerical issues for these large coefficients. Therefore, in this
        case we do not construct a new cut at all (or at least restrict
        to redundant cut).
        """

        add_cut_flag = false

    elseif dual_0_var < 1e-6
        add_cut_flag = false

    else
        # We have to correct the intercept. We do this at this point, as (in
        # contrast to some other corrections) it is not required for Benders cuts.
        lag_obj = lag_obj + epi_state * dual_0_var
    end

    store_dual_values!(node, dual_values, dual_vars, bin_state, cut_generation_regime.state_approximation_regime)

    #println(dual_vars/dual_0_var)

    return (
        dual_values=dual_values,
        dual_0_var=dual_0_var,
        bin_state=bin_state,
        intercept=lag_obj,
        iterations=lag_iterations,
        lag_status=lag_status,
        add_cut_flag=add_cut_flag,
    )
end

function detect_unboundedness(
    node::SDDP.Node,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_Relint,DynamicSDDiP.Core_Conv}
)
    # Initialization
    unbounded_flag = false

    # Solve first auxiliary primal problem
    results_1 = solve_primal_projection_problem(node, epi_state, normalization_coeff, DynamicSDDiP.NoIntegerRelax(), normalization_regime.copy_regime, duality_regime, algo_params, applied_solvers, cut_generation_regime)

    # Solve second auxiliary primal problem
    node.ext[:primal_data] = Dict{Symbol,Any}()
    results_2 = solve_primal_projection_problem(node, epi_state, normalization_coeff, DynamicSDDiP.IntegerRelax(), normalization_regime.copy_regime, duality_regime, algo_params, applied_solvers, cut_generation_regime)

    # Unboundedness check
    if !results_1.unbounded_flag && !results_2.unbounded_flag
        # both problems feasible
        unbounded_flag = false
    elseif results_1.unbounded_flag && results_2.unbounded_flag
        # both problems infeasible
        unbounded_flag = true
    else
        # both problems feasible
        if normalization_regime.unbounded_regime.strict_proxy
            unbounded_flag = true
        else
            unbounded_flag = false
        end
    end

    return (unbounded_flag = unbounded_flag, primal_unified_obj = results_1.primal_unified_obj)
end

function detect_unboundedness(
    node::SDDP.Node,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    normalization_regime::Union{DynamicSDDiP.L₁_Deep,DynamicSDDiP.L₂_Deep,DynamicSDDiP.L∞_Deep,DynamicSDDiP.L₁∞_Deep,DynamicSDDiP.ChenLuedtke}
)

    return (unbounded_flag = false, primal_unified_obj = Inf)
end

function solve_primal_projection_problem(
    node::SDDP.Node,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
)

    # INITIALIZE VALUES
    unbounded_flag = false
    primal_unified_obj = Inf

    # PREPARE PRIMAL TO LAGRANGIAN DUAL PROBLEM (INCLUDES POSSIBLE REGULARIZATION)
    undo_relax = construct_unified_primal_problem!(node, node.index, node.subproblem, epi_state, normalization_coeff, integer_regime, copy_regime, duality_regime, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    # RESET SOLVER (as it may have been changed in between for some reason)
    reset_solver!(node.subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)

    # SOLVE THE PRIMAL PROBLEM
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(node.subproblem)
    end

    # CHECK TERMINATION STATUS TO GET BOUNDS FOR LAGRANGIAN DUAL
    if JuMP.termination_status(node.subproblem) == MOI.OPTIMAL
        primal_unified_obj = JuMP.objective_value(node.subproblem)
    else
        # Infeasibility is an indicator for unboundedness of the Lagrangian dual
        # We do not differentiate termination statuses here, because we need artificial bounds anyway
        unbounded_flag = true
    end

    # CHANGE THE PROBLEM TO THE PREVIOUS FORM AGAIN (INCLUDES POSSIBLE DEREGULARIZATION)
    deconstruct_unified_primal_problem!(node, node.subproblem, undo_relax, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    return (unbounded_flag = unbounded_flag, primal_unified_obj = primal_unified_obj)
end

