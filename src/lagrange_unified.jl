"""
Kelley's method to solve an aggregated Lagrangian dual
"""
function solve_aggregated_lagrangian(
    node::SDDP.Node,
    node_index::Int64,
    outgoing_state::Dict{Symbol,Float64},
    scenario_path,
    epi_state::Float64,
    primal_obj::Float64,
    π_k_agg::Vector{Float64},
    π0_k::Float64,
    bound_results::NamedTuple{(:obj_bound, :dual_bound),Tuple{Float64,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.Kelley
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################
    # A sign bit that is used to avoid if-statements in the models (see SDDP.jl)
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1

    # Storage for the cutting-plane method
    #---------------------------------------------------------------------------
    number_of_noise = length(SDDP.sample_backward_noise_terms(backward_sampling_scheme, node))
    number_of_states_per_noise = get_number_of_states(node, algo_params.state_approximation_regime)

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
    π_k = zeros(number_of_states_per_noise, number_of_noise)
    for j in size(π_k, 2)
        π_k[:,j] = π_k_agg
    end
    # The best estimate for π (former best_mult)
    π_star = zeros(number_of_states_per_noise, number_of_noise)

    # Copy constraint slacks and subgradients solve_aggregated_lagrangian
    #---------------------------------------------------------------------------
    # The expression for ̄x-z (former slacks)
    h_expr = Vector{JuMP.AffExpr}(undef, number_of_states_per_noise)

    # Copy constraint slacks and subgradients per noise term
    #---------------------------------------------------------------------------
    # The current value of ̄x-z (former subgradients)
    h_k = zeros(number_of_states_per_noise, number_of_noise)

    # Epi multiplier aggregated
    #---------------------------------------------------------------------------
    # The current estimate for π0 (in our case determined in initialization)
    # π0_k

    # The best estimate for π0 (former best_mult)
    π0_star = 0

    # Epi slacks and subgradients aggregated
    #---------------------------------------------------------------------------
    # The expression for ̄the epi slack
    w_expr = JuMP.AffExpr(undef)
    # The current value of the epi slack (epi subgradient)
    w_k_agg = 0

    # Epi slacks and subgradients per noise term
    #---------------------------------------------------------------------------
    # The current value of the epi slack for a specific noise term
    w_k = zeros(number_of_noise)

    # Set tolerances
    #---------------------------------------------------------------------------
    atol = algo_params.duality_regime.atol
    rtol = algo_params.duality_regime.rtol
    iteration_limit = algo_params.duality_regime.iteration_limit

    # Set solver for inner problem
    #---------------------------------------------------------------------------
    set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax)

    ############################################################################
    # RELAXING THE COPY CONSTRAINTS
    ############################################################################
    # Relax the copy constraints (same for all noise terms)
    TimerOutputs.@timeit DynamicSDDiP_TIMER "relax_copy" begin
        relax_copy_constraints!(node, x_in_value, h_expr, algo_params.state_approximation_regime)
    end
    node.ext[:backward_data][:old_rhs] = x_in_value
    # Set expression for w_expr
    w_expr = JuMP.@expression(node.subproblem, JuMP.objective_function(node.subproblem) - epi_state))

    ############################################################################
    # SET-UP THE APPROXIMATING CUTTING-PLANE MODEL
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "init_approx_model" begin
        # Approximation of Lagrangian dual by cutting planes
        # Optimizer is re-set anyway
        approx_model = JuMP.Model(GAMS.Optimizer)
        set_solver!(approx_model, algo_params, applied_solvers, :kelley)

        # Create the objective
        # Note that it is always formulated as a maximization problem, but that
        # s modifies the sense appropriately
        JuMP.@variable(approx_model, t)
        set_objective_bound!(approx_model, s, bound_results.obj_bound)
        JuMP.@objective(approx_model, Max, t)

        # Create the dual variables
        # Note that the real dual multipliers are split up into two non-negative
        # variables here, which is required for the Magnanti Wong part later
        JuMP.@variable(approx_model, π⁺[1:number_of_states_per_noise, number_of_noise] >= 0)
        JuMP.@variable(approx_model, π⁻[1:number_of_states_per_noise, number_of_noise] >= 0)
        JuMP.@expression(approx_model, π, π⁺ .- π⁻) # not required to be a constraint
        JuMP.@variable(approx_model, π0 >= 0)
        set_multiplier_bounds!(approx_model, number_of_states_per_noise, number_of_noise, bound_results.dual_bound)
    end

    ############################################################################
    # CUTTING-PLANE METHOD
    ############################################################################
    iter = 0
    lag_status = :none

    # set up optimal value of approx_model (former f_approx)
    t_k = 0 # why zero?
    #-inf is not possible, since then the while loop would not start at all

    while iter < iteration_limit && !isapprox(L_star, t_k, atol = atol, rtol = rtol)
        iter += 1

        ########################################################################
        # SOLVE LAGRANGIAN RELAXATION FOR GIVEN DUAL_VARS
        ########################################################################
        # Evaluate the inner problem and determine a subgradient
        TimerOutputs.@timeit DynamicSDDiP_TIMER "inner_sol" begin
            L_k = _solve_Lagrangian_relaxation!(
                node,
                π_k,
                π0_k,
                h_expr,
                h_k,
                w_expr,
                w_k_agg,
                w_k,
                true)
        end

        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        ########################################################################
        # UPDATE BEST FOUND SOLUTION SO FAR
        ########################################################################
        if s * L_k >= L_star
            L_star = s * L_k
            π_star .= π_k
            π0_star .= π0_k
        end

        ########################################################################
        # ADD CUTTING PLANE
        ########################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "add_cut" begin
            JuMP.@constraint(approx_model, t <= s * (L_k + sum(j in 1:number_of_noise, h_k[:,j]' * (π[:,j] .- π_k[:,j])) + w_k_agg' * (π0 - π0_k)))
        end

        ########################################################################
        # SOLVE APPROXIMATION MODEL
        ########################################################################
        # Get a bound from the approximate model
        TimerOutputs.@timeit DynamicSDDiP_TIMER "outer_sol" begin
            JuMP.optimize!(approx_model)
        end
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(approx_model)
        π_k .= JuMP.value.(π)
        π0_k .= JuMP.value.(π0)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        #print("UB: ", f_approx, ", LB: ", f_actual)
        #println()

        ########################################################################
        if L_star > t_k + atol/10.0
            error("Could not solve for Lagrangian duals. LB > UB.")
        end
    end

    ############################################################################
    # CONVERGENCE ANALYSIS
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "convergence_check" begin
        if isapprox(L_star, t_k, atol = atol, rtol = rtol)
            # CONVERGENCE ACHIEVED
            if isapprox(L_star, s * primal_obj, atol = atol, rtol = rtol)
                # CONVERGENCE TO TRUE OPTIMUM (APPROXIMATELY)
                lag_status = :opt
            else
                # CONVERGENCE TO A SMALLER VALUE THAN THE PRIMAL OBJECTIVE
                # sometimes this occurs due to numerical issues
                # still leads to a valid cut
                lag_status = :conv
            end
        elseif all(h_k .== 0)
            # NO OPTIMALITY ACHIEVED, BUT STILL ALL SUBGRADIENTS ARE ZERO
            # may occur due to numerical issues
            lag_status = :sub
        elseif iter == iteration_limit
            # TERMINATION DUE TO ITERATION LIMIT
            # stil leads to a valid cut
            lag_status = :iter
        end
    end

    ############################################################################
    # APPLY MAGNANTI AND WONG APPROACH IF INTENDED
    ############################################################################
    #TimerOutputs.@timeit DynamicSDDiP_TIMER "magnanti_wong" begin
    #    magnanti_wong!(node, approx_model, π_k, π_star, t_k, h_expr, h_k, s, L_star,
    #        iteration_limit, atol, rtol, algo_params.duality_regime.dual_choice_regime, iter)
    #end

    ############################################################################
    # RESTORE THE COPY CONSTRAINT x.in = value(x.in) (̄x = z)
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "restore_copy" begin
        restore_copy_constraints!(node, x_in_value, algo_params.state_approximation_regime)
    end

    ############################################################################
    # RESET SOLVER
    ############################################################################
    set_solver!(node.subproblem, algo_params, applied_solvers, :forward_pass)

    ############################################################################
    # LOGGING
    ############################################################################
    # print_helper(print_lag_iteration, lag_log_file_handle, iter, t_k, L_star, L_k)

    ############################################################################
    # AGGREGATE DUAL COEFFICIENTS
    ############################################################################
    # Get aggregated dual_vars value by summing up the optimal solutions for each noise
    for j in size(π_star, 2)
        π_k_agg[j] .+= π_star[:,j]
    end

    return (lag_obj = s * L_star, iterations = iter, lag_status = lag_status)

end


function set_multiplier_bounds!(approx_model::JuMP.Model, number_of_states_per_noise::Int, number_of_noise::Int, dual_bound::Float64)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]

    for i in 1:number_of_states_per_noise
        for j in 1:number_of_noise
            JuMP.set_upper_bound(π⁺[i,j], dual_bound)
            JuMP.set_upper_bound(π⁻[i,j], dual_bound)
        end
    end
end


"""
Solving the aggregated Lagrangian relaxation problem, i.e. the inner problem of the
Lagrangian dual
"""

function _solve_Lagrangian_relaxation_aggregated!(
    node::SDDP.Node,
    π_k::Array{Float64, 2},
    π0_k::Float64,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Array{Float64, 2},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    w_k_agg::Float64,
    w_k::Vector{Float64},
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Initialize Lagrangian aggregated value
    L_k = 0

    # Set the Lagrangian relaxation of the objective in the primal model
    old_obj = JuMP.objective_function(model)

    ############################################################################
    # Iterate over all noise terms and solve the single-scenario inner problem
    ############################################################################
    for j in 1:length(node.noise_terms)
        noise = node.noise_terms[j]

        # Parametrize the problem correctly
        parameterize(node, noise)

        # Set objective
        JuMP.set_objective_function(model, JuMP.@expression(model, - π0_k * noise.probability * w_expr - π_k[:,j]' * h_expr))

        # Solve the single Lagrangian problem
        L_k = += _solve_Lagrangian_relaxation_single!(
            model,
            h_expr,
            h_k[:,j],
            w_expr,
            w_k[j],
            update_subgradients
            )

        # Update w_k_agg
        w_k_agg += noise.probability * w_k[j]

    end

    # Reset old objective
    JuMP.set_objective_function(model, old_obj)

    return L_k
end


"""
Solving the Lagrangian relaxation problem, i.e. the inner problem of the
Lagrangian dual, for each noise term
"""

function _solve_Lagrangian_relaxation_single!(
    model::JuMP.Model,
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k_j::Vector{Float64},
    w_expr::JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    w_k_j::Float64,
    update_subgradients::Bool = true,
)

    # Optimization
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == MOI.OPTIMAL

    # Store objective value for return
    L_k_j = JuMP.objective_value(model)

    # Update of subgradients
    if update_subgradients
        h_k .= -JuMP.value.(h_expr)
        w_k_j .= -JuMP.value.(w_expr)
    end

    return L_k_j
end


# function magnanti_wong!(
#     node::SDDP.Node,
#     approx_model::JuMP.Model,
#     π_k::Vector{Float64},
#     π_star::Vector{Float64},
#     t_k::Float64,
#     h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
#     h_k::Vector{Float64},
#     s::Int,
#     L_star::Float64,
#     iteration_limit::Int,
#     atol::Float64,
#     rtol::Float64,
#     dual_choice_regime::DynamicSDDiP.MagnantiWongChoice,
#     iter::Int,
#     )
#
#     π⁺ = approx_model[:π⁺]
#     π⁻ = approx_model[:π⁻]
#     t = approx_model[:t]
#     π = approx_model[:π]
#
#     # Reset objective
#     JuMP.@objective(approx_model, Min, sum(π⁺) + sum(π⁻))
#     JuMP.set_lower_bound(t, t_k)
#
#     # The worst-case scenario in this for-loop is that we run through the
#     # iterations without finding a new dual solution. However if that happens
#     # we can just keep our current λ_star.
#     for _ in (iter+1):iteration_limit
#         JuMP.optimize!(approx_model)
#         @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
#         π_k .= value.(π)
#         L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)
#         if isapprox(L_star, L_k, atol = atol, rtol = rtol)
#             # At this point we tried the smallest ‖π‖ from the cutting plane
#             # problem, and it returned the optimal dual objective value. No
#             # other optimal dual vector can have a smaller norm.
#             π_star = π_k
#             return
#         end
#         JuMP.@constraint(approx_model, t <= s * (L_k + h_k' * (π .- π_k)))
#
#     end
#
#     return
# end
#
#
# function magnanti_wong!(
#     node::SDDP.Node,
#     approx_model::JuMP.Model,
#     π_k::Vector{Float64},
#     π_star::Vector{Float64},
#     t_k::Float64,
#     h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
#     h_k::Vector{Float64},
#     s::Int,
#     L_star::Float64,
#     iteration_limit::Int,
#     atol::Float64,
#     rtol::Float64,
#     dual_choice_regime::DynamicSDDiP.StandardChoice,
#     iter::Int,
#     )
#
#     return
# end
