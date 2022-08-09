"""
Solving the Lagrangian relaxation problem, i.e. the inner problem of the
unified Lagrangian dual
"""

function _solve_unified_Lagrangian_relaxation!(
    node::SDDP.Node,
    π_k::Vector{Float64},
    π0_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    w_k::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Set the objective of the Lagrangian relaxation
    old_obj = JuMP.objective_function(model)
    JuMP.set_objective_function(model, JuMP.@expression(model, π0_k * w_expr + π_k' * h_expr))
    #Infiltrator.@infiltrate

    # Optimization
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_Lag_inner" begin
        JuMP.optimize!(model)
    end

    # Try recovering from numerical issues
    if (JuMP.termination_status(model) != MOI.OPTIMAL)
        elude_numerical_issues!(model, algo_params)
    end

    # Update the correct values
    L_k = JuMP.objective_value(model)
    if update_subgradients
        h_k .= JuMP.value.(h_expr)
        w_k = JuMP.value(w_expr)
    end

    # Reset old objective
    JuMP.set_objective_function(model, old_obj)

    return (L_k = L_k, w_k = w_k)
end


"""
Kelley's method to solve unified Lagrangian dual
"""
function solve_unified_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    primal_obj::Float64,
    π_k::Vector{Float64},
    π0_k::Float64,
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
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
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
    # The original value for x (former old_rhs)
    x_in_value = zeros(number_of_states)
    # The current estimate for π (in our case determined in initialization)
    # π_k
    # The best estimate for π (former best_mult)
    π_star = zeros(number_of_states)
    # The best estimate for the dual objective value (ignoring optimization sense)
    L_star = -Inf
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states)
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)
    #---------------------------------------------------------------------------
    # The current estimate for π0 (in our case determined in initialization)
    # π0_k
    # The best estimate for π0 (former best_mult)
    π0_star = 1.0
    # The expression for ̄the epi slack (will be determined later)
    w_expr = JuMP.AffExpr()
    # The current value of the epi slack
    w_k = 0.0

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
    # LOGGING OF LAGRANGIAN DUAL
    ############################################################################
    #lag_log_file_handle = open("C:/Users/cg4102/Documents/julia_logs/Lagrange.log", "a")
    #print_helper(print_lagrange_header, lag_log_file_handle)

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_approx_model" begin
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        JuMP.set_optimizer(approx_model, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),"MIPGap"=>1e-4
        ))
        JuMP.set_silent(node.subproblem)

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
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
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
            JuMP.@constraint(approx_model, π⁺[1:number_of_states] .<= bound)
            JuMP.@constraint(approx_model, π⁻[1:number_of_states] .<= bound)
            JuMP.set_upper_bound(π₀, bound)
        end

        # Algorithm-specific bounds, e.g. due to regularization
        ########################################################################
        set_multiplier_bounds!(node, approx_model, number_of_states, bound_results.dual_bound,
            algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
            cut_generation_regime.duality_regime)

        # Space restriction by Chen & Luedtke
        ########################################################################
        # Add dual space restriction (span of earlier multipliers)
        dual_space_restriction!(node, approx_model, i, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.dual_space_regime)

        # Normalization of Lagrangian dual
        ########################################################################
        # Add normalization constraint depending on abstract normalization regime
        add_normalization_constraint!(node, approx_model, number_of_states, normalization_coeff, cut_generation_regime.duality_regime.normalization_regime)

    end

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = Inf

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            relax_results = _solve_unified_Lagrangian_relaxation!(node, π_k, π0_k, h_expr, h_k, w_expr, w_k, algo_params, applied_solvers, true)
        end
        L_k = relax_results.L_k
        w_k = relax_results.w_k

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
            JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k) + w_k * (π₀ - π0_k)))
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
                elude_numerical_issues!(approx_model, algo_params)
            end
        end

        t_k = JuMP.objective_value(approx_model)
        π_k .= JuMP.value.(π)
        π0_k = JuMP.value.(π₀)

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
            lag_status = :issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    if lag_status == :opt || lag_status == :unbounded
    # In other cases we do not have an optimal solution from Kelley's method,
    # so finding the minimal norm optimal solution does not make sense.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "minimal_norm" begin
            mn_results = minimal_norm_choice_unified!(node, node_index, approx_model, π_k, π_star, π0_k, π0_star, t_k, h_expr, h_k, w_expr, w_k, s, L_star,
            iteration_limit, atol, rtol, cut_generation_regime.duality_regime.dual_choice_regime, iter, lag_status, algo_params, applied_solvers)

            iter = mn_results.iter
            lag_status = mn_results.lag_status
        end
    elseif isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
        println("Proceeding without minimal norm choice.")
    end

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
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)
    #if node.subproblem.ext[:sddp_policy_graph].ext[:iteration] == 11
    #    Infiltrator.@infiltrate
    #end
    # if all(π_star .== 0) && isapprox(π0_star, 0.0, atol=1e-8)
    #     println(node_index, ", ", node.subproblem.ext[:sddp_policy_graph].ext[:iteration], ", ", lag_status, ", ", epi_state,  ", ",primal_obj, ", ", L_star, ", ", t_k, ", ", π_star, ", ", π0_star)
    # end

    # Set dual_vars (here π_k) to the optimal solution
    π_k .= -π_star
    π0_k = π0_star

    #Infiltrator.@infiltrate

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status, dual_0_var = π0_k)

end


"""
Level Bundle method to solve unified Lagrangian dual
"""
function solve_unified_lagrangian_dual(
    node::SDDP.Node,
    node_index::Int64,
    i::Int64,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    primal_obj::Float64,
    π_k::Vector{Float64},
    π0_k::Float64,
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.LevelBundle
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_states = get_number_of_states(node, cut_generation_regime.state_approximation_regime)
    # The original value for x (former old_rhs)
    x_in_value = zeros(number_of_states)
    # The current estimate for π (in our case determined in initialization)
    # π_k
    # The best estimate for π (former best_mult)
    π_star = zeros(number_of_states)
    # The best estimate for the dual objective value (ignoring optimization sense)
    L_star = -Inf
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states)
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states)
    #---------------------------------------------------------------------------
    # The current estimate for π0 (in our case determined in initialization)
    # π0_k
    # The best estimate for π0 (former best_mult)
    π0_star = 1.0
    # The expression for ̄the epi slack (will be determined later)
    w_expr = JuMP.AffExpr()
    # The current value of the epi slack
    w_k = 0.0

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = cut_generation_regime.duality_regime.atol
    rtol = cut_generation_regime.duality_regime.rtol
    iteration_limit = cut_generation_regime.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)

    # Set bundle_parameters
    #---------------------------------------------------------------------------
    level_factor = dual_solution_regime.level_factor

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
        approx_model = JuMP.Model(Gurobi.Optimizer)
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
        JuMP.@variable(approx_model, π⁺[1:number_of_states] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states] >= 0)
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
            JuMP.@constraint(approx_model, π⁺[1:number_of_states] .<= bound)
            JuMP.@constraint(approx_model, π⁻[1:number_of_states] .<= bound)
            JuMP.set_upper_bound(π₀, bound)
        end

        # Algorithm-specific bounds, e.g. due to regularization
        ########################################################################
        set_multiplier_bounds!(node, approx_model, number_of_states, bound_results.dual_bound,
            algo_params.regularization_regime, cut_generation_regime.state_approximation_regime,
            cut_generation_regime.duality_regime)

        # Space restriction by Chen & Luedtke
        ########################################################################
        # Add dual space restriction (span of earlier multipliers)
        dual_space_restriction!(node, approx_model, i, cut_generation_regime.state_approximation_regime, cut_generation_regime.duality_regime.dual_space_regime)

        # Normalization of Lagrangian dual
        ########################################################################
        # Add normalization constraint depending on abstract normalization regime
        add_normalization_constraint!(node, approx_model, number_of_states, normalization_coeff, cut_generation_regime.duality_regime.normalization_regime)

    end

    ############################################################################
    # BUNDLE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = Inf

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            relax_results = _solve_unified_Lagrangian_relaxation!(node, π_k, π0_k, h_expr, h_k, w_expr, w_k, algo_params, applied_solvers, true)
        end
        L_k = relax_results.L_k
        w_k = relax_results.w_k

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]
        #Infiltrator.@infiltrate

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
            JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k) + w_k * (π₀ - π0_k)))
        end

        ########################################################################
        # RESET OBJECTIVE FOR APPROX_MODEL AFTER NONLINEAR MODEL
        ########################################################################
        JuMP.@objective(approx_model, Max, t)
        if isa(cut_generation_regime.duality_regime.normalization_regime, DynamicSDDiP.L₂_Deep)
            set_solver!(approx_model, algo_params, applied_solvers, :l₂, algo_params.solver_approach)
        else
            set_solver!(approx_model, algo_params, applied_solvers, :kelley, algo_params.solver_approach)
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
                elude_numerical_issues!(approx_model, algo_params)
            end
        end
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # COMPUTE GAP AND FORM A NEW LEVEL
        ########################################################################
        f_up = t_k # objective value of the approx_model
        f_down = L_star # best lower bound so far
        gap = f_up - f_down

        """
        We use a convex combination of f_up and f_down for the new level.
        level = (1-a) f_up + a f_down = f_up - a (f_up - f_down) = f_up - a gap

        For level_factor = 0 => level = f_up, i.e. t_k, that means, this
        reduces to the basic cutting-plane method.

        For level_factor = 1 => level = f_down, i.e. L_star. That means that
        we search for the closest multiplier, for which the level is larger than
        L_star. At least for L_star = L_k then the multiplier does not change
        anymore, because the condition is trivially satisfied for the current
        multiplier. Otherwise, this should yield one of the previous multipliers,
        so again no improvement.
        """

        # TODO: POSSIBLE IMPROVEMENTS/CHANGES
        """
        1.) STRUCTURE
        In the literature often a different order of the steps is proposed,
        so maybe we should change this.

        In particular, the gap is often determined before the approx_model
        is solved with the new multipliers. That means that the gap is determined
        using f_down and f_up = t_{k-1}. In this case, t_0 has to be determined
        appropriately, e.g. we can just choose primal_obj.

        2.) STOPPING
        In the literature it is often proposed to stop if the gap is sufficiently
        small. In our case, this doesn't make a difference, but since the gap
        can be determined differently as well, then our stopping criterion
        would change.

        3.) STABILITY CENTER
        We always choose the new multiplier obtained by solving the proximal
        problem as new stability center. It is also possible to change the
        stability center only if a sufficiently large improvement is achieved.

        4.) DETERMINING f_up
        Instead of solving approx_model, we could also determine f_up in the
        following way:
        > If the proximal problem is infeasible, set f_up = level.
        > If the proximal problem is feasible, do not change f_up.

        5.) DETERMINE gap
        Right now, the gap is determined as f_up - f_down with
        f_up = t_k
        f_down = L_star

        We could also choose f_up = primal_obj, since we already know that
        this is the best possible upper bound. Is this beneficial?
        """

        level = f_up - gap * level_factor
        # - atol/10.0 for numerical issues?
        JuMP.set_lower_bound(t, level)

        ########################################################################
        # DETERMINE NEXT ITERATION USING PROXIMAL PROBLEM
        ########################################################################
        # Objective function of approx model has to be adapted to new center
        JuMP.@objective(approx_model, Min, sum((π_k[i] - π[i])^2 for i in 1:number_of_states) + (π0_k - π₀)^2)
        set_solver!(approx_model, algo_params, applied_solvers, :level_bundle, algo_params.solver_approach)

        TimerOutputs.@timeit DynamicSDDiP_TIMER "bundle_sol" begin
            JuMP.optimize!(approx_model)
        end

        Infiltrator.@infiltrate

        """ This is removed, as sometimes the QCP is solved to optimality, but
        Gurobi is not able to get an optimality certificate
        (QCP status(13): Unable to satisfy optimality tolerances;
        a sub-optimal solution is available), even though other solvers prove
        that the solution is indeed optimal.
        Even if the solution is suboptimal, this should not be a problem, as long
        as it is feasible. Therefore, I think the algorithm should not
        stop here.
        """
        #@assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        # TODO: Maybe change this in the future again

        π_k .= JuMP.value.(π)
        π0_k = JuMP.value.(π₀)
        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        # Delete the level lower bound for the original approx_model again
        JuMP.delete_lower_bound(t)

        # Sometimes the solver (e.g. Gurobi) provides a float point approximation
        # of zero, which is slightly negative, e.g. 1.144917E-16, even though
        # π0_k >= 0 is enforced as a constraint.
        # In such case, the inner problem of the Lagrangian relaxation may
        # become unbounded. Therefore, we set π0_k manually to 0 then.
        if π0_k < 0
            π0_k = 0.0
        end

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        println(L_k, ", ", L_star, ", ", t_k)

        ########################################################################
        if L_star > t_k + atol/10.0
            # error("Could not solve for Lagrangian duals. LB > UB.")
            break
        end
    end

    println()

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
            # objective bound 1e-9
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
            lag_status = :issues
        end
    end

    ############################################################################
    # APPLY MINIMAL NORM CHOICE APPROACH IF INTENDED
    ############################################################################
    if lag_status == :opt || lag_status == :unbounded
        # In other cases we do not have an optimal solution from Kelley's method,
        # so finding the minimal norm optimal solution does not make sense.
        TimerOutputs.@timeit DynamicSDDiP_TIMER "minimal_norm" begin
            mn_results = minimal_norm_choice_unified!(node, node_index, approx_model, π_k, π_star, π0_k, π0_star, t_k, h_expr, h_k, w_expr, w_k, s, L_star,
            iteration_limit, atol, rtol, cut_generation_regime.duality_regime.dual_choice_regime, iter, lag_status, algo_params, applied_solvers)

            iter = mn_results.iter
            lag_status = mn_results.lag_status
        end
    elseif isa(cut_generation_regime.duality_regime.dual_choice_regime, DynamicSDDiP.MinimalNormChoice)
        println("Proceeding without minimal norm choice.")
    end

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
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    # Set dual_vars (here π_k) to the optimal solution
    π_k .= -π_star
    π0_k = π0_star

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status, dual_0_var = π0_k)

end

"""
Given the optimal dual objective value from the Kelley's method, try to find
the optimal dual multipliers with the smallest L1-norm.

Note that this is only done until the maximum number of iterations is
achieved in total.
"""

function minimal_norm_choice_unified!(
    node::SDDP.Node,
    node_index::Int64,
    approx_model::JuMP.Model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    π0_k::Float64,
    π0_star::Float64,
    t_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    w_k::Float64,
    s::Int,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.MinimalNormChoice,
    iter::Int,
    lag_status::Symbol,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    )

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    t = approx_model[:t]
    π = approx_model[:π]
    π₀ = approx_model[:π₀]

    # We actually would have to solve with the objective (π⁺ + π⁻)/π₀
    # To avoid a nonlinear problem, we fix the scaling factor π₀ to its optimal
    # value from the previous solution method.
    #JuMP.fix(π₀, π0_k, force=true)
    JuMP.fix(π₀, π0_star, force=true)
    #π0_k = π0_star

    # Reset objective
    JuMP.@objective(approx_model, Min, sum(π⁺) + sum(π⁻))
    JuMP.set_lower_bound(t, t_k)

    # The worst-case scenario in this for-loop is that we run through the
    # iterations without finding a new dual solution. However if that happens
    # we can just keep our current λ_star.
    for it in (iter+1):iteration_limit
        JuMP.optimize!(approx_model)

        try
            @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        catch err
            println("Proceeding without minimal norm choice.")
            return (iter=it, lag_status=:mn_issue)
        end

        π_k .= JuMP.value.(π)

        relax_results = _solve_unified_Lagrangian_relaxation!(node, π_k, π0_k, h_expr, h_k, w_expr, w_k, algo_params, applied_solvers, true)

        L_k = relax_results.L_k
        w_k = relax_results.w_k

        if isapprox(L_star, L_k, atol = atol, rtol = rtol)
            # At this point we tried the smallest ‖π‖ from the cutting plane
            # problem, and it returned the optimal dual objective value. No
            # other optimal dual vector can have a smaller norm.
            π_star .= π_k
            return (iter=it, lag_status=:mn_opt)
        end
        JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k) + w_k * (π₀ - π0_k)))
        # note that the last term is always zero, since π₀ is fixed

    end

    return (iter=iteration_limit, lag_status=:mn_iter)

end


function minimal_norm_choice_unified!(
    node::SDDP.Node,
    node_index::Int64,
    approx_model::JuMP.Model,
    π_k::Vector{Float64},
    π_star::Vector{Float64},
    π0_k::Float64,
    π0_star::Float64,
    t_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    w_k::Float64,
    s::Int,
    L_star::Float64,
    iteration_limit::Int,
    atol::Float64,
    rtol::Float64,
    dual_choice_regime::DynamicSDDiP.StandardChoice,
    iter::Int,
    lag_status::Symbol,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    )

    return (iter=iter, lag_status=lag_status)
end
