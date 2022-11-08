"""
Solving the aggregated Lagrangian relaxation problem, i.e. the inner problem of the
unified Lagrangian dual
"""

function _solve_Lagrangian_relaxation_aggregated!(
    node::SDDP.Node,
    π_k::Array{Float64,2},
    π0_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Array{Float64,2},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    w_k::Vector{Float64},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Initialize Lagrangian aggregated value
    L_k = 0
    w_k_agg = 0

    # Store the old objective
    old_obj = JuMP.objective_function(model)

    ############################################################################
    # Iterate over all noise terms and solve the single-scenario inner problem
    ############################################################################
    for j in 1:length(node.noise_terms)
        noise = node.noise_terms[j]

        # Parametrize the problem correctly
        parameterize(node, noise.term)

        # Set objective
        JuMP.set_objective_function(model, JuMP.@expression(model,  π0_k * noise.probability * w_expr + π_k[j,:]' * h_expr))

        # Solve the single Lagrangian problem
        Lag_results = _solve_Lagrangian_relaxation_single!(
            model,
            h_expr,
            h_k[j,:],
            w_expr,
            w_k[j],
            algo_params,
            applied_solvers,
            update_subgradients
            )

        # somehow with slicing an udpate directly in _solve_Lagrangian_relaxation_single does not seem to work
        L_k += Lag_results.L_k_j
        h_k[j,:] = Lag_results.h_k_j
        w_k[j] = Lag_results.w_k_j

        # Update w_k_agg
        w_k_agg += noise.probability * w_k[j]
    end

    # Reset old objective
    JuMP.set_objective_function(model, old_obj)

    return (L_k = L_k, w_k_agg = w_k_agg)
end

"""
Solving the Lagrangian relaxation problem, i.e. the inner problem of the
unified Lagrangian dual, for a single noise term
"""

function _solve_Lagrangian_relaxation_single!(
    model::JuMP.Model,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k_j::Vector{Float64},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    w_k_j::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    update_subgradients::Bool = true,
)

    # Optimization
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_Lag_inner" begin
        JuMP.optimize!(model)
    end
    # Try recovering from numerical issues
    if (JuMP.termination_status(model) != MOI.OPTIMAL)
        elude_numerical_issues!(model, algo_params)
    end

    # Store objective value for return
    L_k_j = JuMP.objective_value(model)

    # Update of subgradients
    if update_subgradients
        h_k_j .= JuMP.value.(h_expr)
        w_k_j = JuMP.value.(w_expr)
    end

    return (
    L_k_j = L_k_j,
    h_k_j = h_k_j,
    w_k_j = w_k_j
    )
end


"""
Kelley's method to solve an aggregated Lagrangian dual
"""
function solve_aggregated_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    outgoing_state::Dict{Symbol,Float64},
    scenario_path,
    i::Int64,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    primal_obj::Float64,
    π_k_agg::Vector{Float64},
    π0_k::Float64,
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.Kelley,
)

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_noise = length(SDDP.sample_backward_noise_terms(backward_sampling_scheme, node))
    number_of_states_per_noise = get_number_of_states(node, cut_generation_regime.state_approximation_regime)

    # The best estimate for the dual objective value (ignoring optimization sense)
    L_star = -Inf

    # The original value for x (former old_rhs); same for all noises
    x_in_value = zeros(number_of_states_per_noise)

    # Copy constraint multipliers aggregated
    #---------------------------------------------------------------------------
    # The current estimate for π (in our case determined in initialization)
    # π_k_agg

    # Copy constraint multipliers per noise term
    #---------------------------------------------------------------------------
    # The current estimate for π0 (in our case determined in initialization)
    # Initialized with same value (zeros) as π_k_agg
    π_k = zeros(number_of_noise, number_of_states_per_noise)
    # for j in size(π_k, 2)
    #  π_k[:,j] = π_k_agg
    # end

    # The best estimatce for π (former best_mult)
    π_star = zeros(number_of_noise, number_of_states_per_noise)

    # Copy constraint slacks and subgradients solve_aggregated_lagrangian
    #---------------------------------------------------------------------------
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states_per_noise)

    # Copy constraint slacks and subgradients per noise term
    #---------------------------------------------------------------------------
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_noise, number_of_states_per_noise)

    # Epi multiplier aggregated
    #---------------------------------------------------------------------------
    # The current estimate for π0 (in our case determined in initialization)
    # π0_k

    # The best estimate for π0 (former best_mult)
    π0_star = 1.0

    # Epi slacks and subgradients aggregated
    #---------------------------------------------------------------------------
    # The expression for ̄the epi slack
    w_expr = JuMP.AffExpr()
    # The current value of the epi slack (epi subgradient)
    w_k_agg = 0.0

    # Epi slacks and subgradients per noise term
    #---------------------------------------------------------------------------
    # The current value of the epi slack for a specific noise term
    w_k = zeros(number_of_noise)

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = cut_generation_regime.duality_regime.atol
    rtol = cut_generation_regime.duality_regime.rtol
    iteration_limit = cut_generation_regime.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.copy_regime)
    end
    node.ext[:backward_data][:old_rhs] = x_in_value

    ############################################################################
    # SET EPI SLACK
    ############################################################################
    w_expr = JuMP.@expression(node.subproblem, JuMP.objective_function(node.subproblem) - epi_state)

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_approx_model" begin
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        approx_model = JuMP.Model()
        JuMP.set_optimizer(approx_model, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4,"TimeLimit"=>300,"NumericFocus"=>algo_params.numerical_focus
        ))
        JuMP.set_silent(approx_model)

        approx_model.ext[:sddp_policy_graph] = node.subproblem.ext[:sddp_policy_graph]

        if isa(cut_generation_regime.duality_regime.normalization_regime, DynamicSDDiP.L₂_Deep)
            set_solver!(approx_model, algo_params, applied_solvers, :l₂, algo_params.solver_approach)
        else
            set_solver!(approx_model, algo_params, applied_solvers, :kelley, algo_params.solver_approach)
        end

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Norm Minimization part later
        JuMP.@variable(approx_model, π⁺[1:number_of_noise, 1:number_of_states_per_noise] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_noise, 1:number_of_states_per_noise] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        JuMP.@variable(approx_model, π₀ >= 0)

        # User-specific bounds
        ########################################################################
        # We cannot use the primal_obj as an obj_bound in the unified framework,
        # so we use an arbitrarily chosen upper bound or bound the dual multipliers.
        if !isnothing(cut_generation_regime.duality_regime.user_dual_objective_bound)
            set_objective_bound!(approx_model, s, cut_generation_regime.duality_regime.user_dual_objective_bound)
        end

        if !isnothing(cut_generation_regime.duality_regime.user_dual_multiplier_bound)
            bound = cut_generation_regime.duality_regime.user_dual_multiplier_bound
            JuMP.@constraint(approx_model, π⁺[1:number_of_noise, 1:number_of_states_per_noise] .<= bound)
            JuMP.@constraint(approx_model, π⁻[1:number_of_noise, 1:number_of_states_per_noise] .<= bound)
            JuMP.set_upper_bound(π₀, bound)
        end

        # Algorithm-specific bounds, e.g. due to regularization
        ########################################################################
        if !isinf(bound_results.dual_bound)
            set_multiplier_bounds!(node, approx_model, number_of_states_per_noise, number_of_noise, bound_results.dual_bound,
                algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
                cut_generation_regime.duality_regime)
        end

        # Space restriction by Chen & Luedtke
        ########################################################################
        # Add dual space restriction (span of earlier multipliers)
        dual_space_restriction_agg!(node, approx_model, i, number_of_states_per_noise, number_of_noise, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.dual_space_regime)

        # Normalization of Lagrangian dual
        ########################################################################
        # Add normalization constraint depending on abstract normalization regime
        add_normalization_constraint!(node, approx_model, number_of_states_per_noise, number_of_noise, normalization_coeff, cut_generation_regime.duality_regime.normalization_regime)

    end

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none
    feas_flag = false

    # set up optimal value of approx_model (former f_approx)
    t_k = Inf

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol_agg" begin
            relax_results = _solve_Lagrangian_relaxation_aggregated!(node, π_k, π0_k, h_expr, h_k, w_expr, w_k, algo_params, applied_solvers, true)
        end
        L_k = relax_results.L_k
        w_k_agg = relax_results.w_k_agg

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
            π0_star = π0_k
        end

        ########################################################################
        # ADD CUTTING PLANE
        ########################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "add_cut" begin
            JuMP.@constraint(approx_model, t <= s * (L_k + sum(h_k[j,:]' * (π[j,:] .- π_k[j,:]) for j in 1:number_of_noise) + w_k_agg' * (π₀ - π0_k)))
        end

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        TimerOutputs.@timeit DynamicSDDiP_TIMER "outer_sol" begin
            JuMP.optimize!(approx_model)

            # Try recovering from numerical issues
            if (JuMP.termination_status(approx_model) != MOI.OPTIMAL)
                #Infiltrator.@infiltrate
                #elude_numerical_issues!(approx_model, algo_params)
                feas_flag = true
                break
            end
        end

        t_k = JuMP.objective_value(approx_model)
        π_k .= JuMP.value.(π)
        π0_k = JuMP.value.(π₀)

        # if JuMP.value(π0) == 0
        #     π0_k = π0_k
        # else
        #     π0_k = JuMP.value.(π0)
        # end

        # Sometimes the solver (e.g. Gurobi) provides a float point approximation
        # of zero, which is slightly negative, e.g. 1.144917E-16, even though
        # π0_k >= 0 is enforced as a constraint.
        # In such case, the inner problem of the Lagrangian relaxation may
        # become unbounded. Therefore, we set π0_k manually to 0 then.
        if π0_k < 0
            π0_k = 0.0
        end

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        if L_star > t_k + atol/10.0
            #error("Could not solve for Lagrangian duals. LB > UB.")
            break
        end

    end

    ############################################################################
    # CONVERGENCE ANALYSIS
    ############################################################################
    obj_bound = 1e15
    if !isnothing(cut_generation_regime.duality_regime.user_dual_objective_bound)
        obj_bound = cut_generation_regime.duality_regime.user_dual_objective_bound
    end

    TimerOutputs.@timeit DynamicSDDiP_TIMER "convergence_check" begin
        if isapprox(L_star, t_k, atol = atol, rtol = rtol) && isapprox(L_star, obj_bound, rtol=1e-4)
            # UNBOUNDEDNESS DETECTED, which is only bounded by artificial
            # objective bound
            # Leads still to a valid cut, but no vertex of the reverse polar set
            lag_status = :unbounded
        elseif isapprox(L_star, t_k, atol = atol, rtol = rtol)
            # CONVERGENCE ACHIEVED
            # it does not make sense to compare with primal_obj here
            lag_status = :opt
        elseif all(h_k .== 0) && w_k == 0
            # NO OPTIMALITY ACHIEVED, BUT STILL ALL SUBGRADIENTS ARE ZERO
            # may occur due to numerical issues
            lag_status = :sub
        elseif iter == iteration_limit
            # TERMINATION DUE TO ITERATION LIMIT
            # stil leads to a valid cut
            lag_status = :iter
        elseif L_star > t_k + atol/10.0
            # NUMERICAL ISSUES, LOWER BOUND EXCEEDS UPPER BOUND
            lag_status = :bound_issues
        elseif feas_flag
            # NUMERICAL ISSUES, SOME SUBPROBLEM BECAME INFEASIBLE (OR UNBOUNDED)
            lag_status = :feas_issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    # if lag_status == :opt || lag_status == :unbounded
    # # In other cases we do not have an optimal solution from Kelley's method,
    # # so finding the minimal norm optimal solution does not make sense.
    #     TimerOutputs.@timeit DynamicSDDiP_TIMER "minimal_norm" begin
    #         mn_results = minimal_norm_choice_unified!(node, node_index, approx_model, π_k, π_star, π0_k, π0_star, t_k, h_expr, h_k, w_expr, w_k, s, L_star,
    #         iteration_limit, atol, rtol, cut_generation_regime.duality_regime.dual_choice_regime, iter, lag_status, algo_params, applied_solvers)
    #
    #         iter = mn_results.iter
    #         lag_status = mn_results.lag_status
    #     end
    # # elseif isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
    #     # println("Proceeding without minimal norm choice.")
    # end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, cut_generation_regime.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass, algo_params.solver_approach)

    ############################################################################
    # AGGREGATE DUAL COEFFICIENTS
    ############################################################################
    # Get aggregated dual_vars value by summing up the optimal solutions for each noise
    for i in size(π_star,2)
        for j in size(π_star,1)
            π_k_agg[i] += π_star[j,i]
        end
    end
    π_k_agg .= - π_k_agg
    π0_k = π0_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status, dual_0_var = π0_k)

end
