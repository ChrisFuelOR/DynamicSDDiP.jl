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
    set_solver!(subproblem, algo_params, applied_solvers, :LP_relax, algo_params.solver_approach)

    JuMP.optimize!(subproblem)

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
    DynamicSDDiP.set_solver!(subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)

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
        subproblem.ext[:sddp_policy_graph].ext[:agg_lag_iterations] += results.iterations
        lag_status = results.lag_status

        ########################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ########################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(lag_status, duality_regime.dual_status_regime)

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
    # GET PRIMAL SOLUTION, TRACK CONVERGENCE AND DEBUGGING
    ############################################################################
    """ Note that we also solve the primal problem if we generate cuts in
    the unified framework to allow for a possible tightness check.
    However, we cannot use the obtained value to bound the Lagrangian dual
    problem. """

    # REGULARIZE PROBLEM IF REGULARIZATION IS USED
    node.ext[:regularization_data] = Dict{Symbol,Any}()
    regularize_bw!(node, node_index, subproblem, cut_generation_regime, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    # RESET SOLVER (as it may have been changed in between for some reason)
    DynamicSDDiP.set_solver!(subproblem, algo_params, applied_solvers, :backward_pass, algo_params.solver_approach)

    # SOLVE PRIMAL PROBLEM (can be regularized or not)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_primal" begin
        JuMP.optimize!(subproblem)
    end

    # Try recovering from numerical issues
    if (JuMP.termination_status(subproblem) != MOI.OPTIMAL)
        elude_numerical_issues!(subproblem, algo_params)
    end

    #if node.subproblem.ext[:sddp_policy_graph].ext[:iteration] == 10
    #    Infiltrator.@infiltrate
    #end

    #@assert JuMP.termination_status(subproblem) == MOI.OPTIMAL
    primal_obj = JuMP.objective_value(subproblem)

    # DEREGULARIZE PROBLEM IF REQUIRED
    deregularize_bw!(node, subproblem, algo_params.regularization_regime, cut_generation_regime.state_approximation_regime)

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

    ############################################################################
    # GET BOUNDS FOR LAGRANGIAN DUAL
    ############################################################################
    bound_results = get_dual_bounds(node, node_index, algo_params, primal_obj, duality_regime.dual_bound_regime)
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    ############################################################################
    # GET CORE POINT RESP. NORMALIZATION COEFFICIENTS
    ############################################################################
    """ If we do not use a normalization that requires a core point, we simply
    return nothing """
    normalization_coeff = get_normalization_coefficients(
        node,
        number_of_states,
        epi_state,
        algo_params,
        cut_generation_regime.state_approximation_regime,
        duality_regime.normalization_regime,
        duality_regime.copy_regime
    )

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
                        primal_obj,
                        dual_vars,
                        dual_0_var,
                        bound_results,
                        algo_params,
                        cut_generation_regime,
                        applied_solvers,
                        duality_regime.dual_solution_regime
                        )
        end

        lag_obj = results.lag_obj
        lag_iterations = results.iterations

        #Infiltrator.@infiltrate

        subproblem.ext[:sddp_policy_graph].ext[:agg_lag_iterations] += results.iterations
        lag_status = results.lag_status
        dual_0_var = results.dual_0_var

        #if node.subproblem.ext[:sddp_policy_graph].ext[:iteration] == 4
        #    Infiltrator.@infiltrate
        #end

        ########################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ########################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(lag_status, duality_regime.dual_status_regime)

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

    # elseif isapprox(lag_obj, 0.0, atol=1e-8) && isapprox(dual_0_var, 0.0, atol=1e-8)
    #     """ The incumbent (state, epi_state) is contained in the epigraph (of the
    #     convex closure of the value function). Hence, we cannot obtain a cut
    #     to separate it from the epigraph.
    #     Furthermore, since dual_0_var is 0, the cut that we obtain is a (redundant)
    #     feasibility cut, e.g. state >= lower_bound(state). Since we only want
    #     to deal with optimality cuts in multistage problems, we do not want to
    #     use such cut at all. Therefore, by division by dual_0_var we construct
    #     classical Lagrangian optimality cuts. However, if dual_0_var is 0, this
    #     causes numerical issues.
    #     We avoid this by replacing the occuring redundant cut in this case by
    #     a different redundant cut."""
    #
    #     dual_0_var = 1.0
    #     dual_vars .= zeros(length(dual_vars))
    #     lag_obj = 0.0

    elseif !isnothing(normalization_coeff) && all(normalization_coeff.ω .== 0.0) && isapprox(normalization_coeff.ω₀, 0.0, atol=1e-8)
        """ If the linear pseudonorm is used, but all coefficients are zero,
        then the Lagrangian dual is not normalized, but unbounded. Analogously,
        a zero function is optimized over the reverse polar set, which can yield
        any point in this unbounded set. Therefore, we are not guaranteed to
        obtain any meaningful cut.

        Note that we artificially bound the dual objective with 1e-9, so that
        we obtain valid, but usually very large cut coefficients.
        However, experiments show that the cuts still tend to be invalid due
        to numerical issues for these large coefficients. Therefore, in this
        case we do not construct a new cut at all (or at least restrict
        to redundant cut).
        """

        add_cut_flag = false
        # dual_0_var = 1.0
        # dual_vars .= zeros(length(dual_vars))
        # lag_obj = 0.0

    else
        # We have to correct the intercept. We do this at this point, as (in
        # contrast to some other corrections) it is not required for Benders cuts.
        lag_obj = lag_obj + epi_state * dual_0_var
    end

    store_dual_values!(node, dual_values, dual_vars, bin_state, cut_generation_regime.state_approximation_regime)

    #Infiltrator.@infiltrate

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
# AUXILIARY METHODS
#*******************************************************************************

"""
Determining objective and/or variable bounds for the Lagrangian dual
if ValueBound is used.

Note that we always solve the primal problem, even if we do not use its
objective value as the objective bound, as this is more convenient for
debugging purposes, and does not take much time compared to solving the
Lagrangian dual.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.ValueBound,
    )

    return (
        obj_bound = primal_obj,
        dual_bound = Inf
    )

end

"""
Determining objective and/or variable bounds for the Lagrangian dual
if NormBound is used.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.NormBound,
    )

    dual_bound = Inf
    if isa(algo_params.regularization_regime, DynamicSDDiP.NoRegularization)
        # if no regularization is used, bounds should be Inf even if intended to use
        dual_bound = Inf
    else
        dual_bound = algo_params.regularization_regime.sigma[node_index]
    end

    return (
        obj_bound = Inf,
        dual_bound = dual_bound
    )

end

"""
Determining objective and/or variable bounds for the Lagrangian dual
if BothBound is used.
"""
function get_dual_bounds(
    node::SDDP.Node,
    node_index::Int64,
    algo_params::DynamicSDDiP.AlgoParams,
    primal_obj::Float64,
    dual_bound_regime::DynamicSDDiP.BothBounds,
    )

    dual_bound = Inf
    if isa(algo_params.regularization_regime, DynamicSDDiP.NoRegularization)
        # if no regularization is used, bounds should be Inf even if intended to use
        dual_bound = Inf
    else
        dual_bound = algo_params.regularization_regime.sigma[node_index]
    end

    return (
        obj_bound = primal_obj,
        dual_bound = dual_bound
    )

end


"""
Checking the status of the Lagrangian dual solution and throw an error if required
under rigorous regime.
"""
function lagrangian_status_check(
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Rigorous,
    )

    if lag_status == :conv
        error("Lagrangian dual converged to value < solver_obj.")
    elseif lag_status == :sub
        error("Lagrangian dual had subgradients zero without LB=UB.")
    elseif lag_status == :iter
        error("Solving Lagrangian dual exceeded iteration limit.")
    elseif lag_status == :unbounded
        error("Normalized Lagrangian dual unbounded and reached artificial bound.")
    elseif lag_status == :issues
        error("Lagrangian LB > UB due to numerical issues.")
    end

    return
end


"""
Trivial check of the status of the Lagrangian dual solution under lax regime.
"""
function lagrangian_status_check(
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Lax,
    )

    return
end


"""
Initializing duals with zero vector.
"""
function initialize_duals(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_initalization_regime::DynamicSDDiP.ZeroDuals,
)

    # Get number of states and create zero vector for duals
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
    dual_vars_initial = zeros(number_of_states)

    return dual_vars_initial

end


"""
Initializing duals by solving LP relaxation.
"""
function initialize_duals(
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
    undo_relax = JuMP.relax_integrality(subproblem);

    # Define appropriate solver
    set_solver!(subproblem, algo_params, applied_solvers, :LP_relax, algo_params.solver_approach)

    # Solve LP Relaxation
    JuMP.optimize!(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL
    # or MOI.FEASIBLE_POINT???

    # Get dual values (reduced costs) for binary states as initial solution
    get_and_set_dual_values!(node, dual_vars_initial, cut_generation_regime.state_approximation_regime)

    # Undo relaxation
    undo_relax()

    return dual_vars_initial

end

function get_and_set_dual_values!(
    node::SDDP.Node,
    dual_vars_initial::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
       variable_name = node.ext[:backward_data][:bin_states][name]
       reference_to_constr = JuMP.FixRef(variable_name)
       dual_vars_initial[i] = JuMP.dual(reference_to_constr)
    end

    return
end

function get_and_set_dual_values!(
    node::SDDP.Node,
    dual_vars_initial::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    for (i, name) in enumerate(keys(node.states))

        reference_to_constr = JuMP.FixRef(node.states[name].in)
        dual_vars_initial[i] = JuMP.dual(reference_to_constr)
    end

    return
end

function store_dual_values!(
    node::SDDP.Node,
    dual_values::Dict{Symbol, Float64},
    dual_vars::Vector{Float64},
    bin_state::Dict{Symbol, BinaryState},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    #old_rhs = node.ext[:backward_data][:old_rhs]

    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
        dual_values[name] = dual_vars[i]

        #value = old_rhs[i]
        value = JuMP.fix_value(node.ext[:backward_data][:bin_states][name])
        x_name = node.ext[:backward_data][:bin_x_names][name]
        k = node.ext[:backward_data][:bin_k][name]
        bin_state[name] = BinaryState(value, x_name, k)
    end

    return
end

function store_dual_values!(
    node::SDDP.Node,
    dual_values::Dict{Symbol, Float64},
    dual_vars::Vector{Float64},
    bin_state::Dict{Symbol, BinaryState},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    for (i, name) in enumerate(keys(node.states))
        dual_values[name] = dual_vars[i]
    end

    return
end

function get_number_of_states(
    node::SDDP.Node,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    return length(node.ext[:backward_data][:bin_states])
end

function get_number_of_states(
    node::SDDP.Node,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    return length(node.states)
end
