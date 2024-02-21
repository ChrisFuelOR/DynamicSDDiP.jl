
function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

    for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(state)
        # Store expression for slack
        h_expr[i] = JuMP.@expression(node.subproblem, state - x_in_value[i])
        # Relax copy constraint (i.e. z does not have to take the value of ̄x anymore)
        JuMP.unfix(state)

        # Set bounds and integer constraints based on copy_regime
        follow_state_unfixing_binary!(state, copy_regime)
    end

    return
end


function relax_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

    for (i, (name, state)) in enumerate(node.states)
        # Store original value of ̄x, which z was fixed to
        x_in_value[i] = JuMP.fix_value(state.in)
        # Store expression for slack
        h_expr[i] = JuMP.@expression(node.subproblem, state.in - x_in_value[i])
        # Relax copy constraint (i.e. z does not have to take the value of ̄x anymore)
        JuMP.unfix(state.in)

        # Set bounds and integer constraints based on copy_regime
        variable_info = node.ext[:state_info_storage][name].in
        follow_state_unfixing!(state, variable_info, copy_regime)
    end

    return
end

function restore_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    for (i, (_, bin_state)) in enumerate(node.ext[:backward_data][:bin_states])
        prepare_state_fixing_binary!(node, bin_state)
        JuMP.fix(bin_state, x_in_value[i], force = true)
    end

    return
end

function restore_copy_constraints!(
    node::SDDP.Node,
    x_in_value::Vector{Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    for (i, (_, state)) in enumerate(node.states)
        prepare_state_fixing!(node, state)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end

    return
end

function set_objective_bound!(approx_model::JuMP.Model, s::Int, obj_bound::Float64)

    JuMP.set_upper_bound(approx_model[:t], s * obj_bound)

    return
end

function determine_weights!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation)

    weights = ones(number_of_states)

    for (i, (key, value)) in enumerate(node.ext[:backward_data][:bin_states])
        associated_original_state = node.ext[:backward_data][:bin_x_names][key]
    	beta = state_approximation_regime.binary_precision[associated_original_state]
    	associated_k = node.ext[:backward_data][:bin_k][key]
        #weights[i] = 2
        weights[i] = 2^(associated_k-1) * beta
    end

    return weights
end

function determine_weights!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation)

    weights = ones(number_of_states)
    return weights
end

function determine_weights!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, regularization_regime::Union{DynamicSDDiP.NoRegularization,DynamicSDDiP.Regularization},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation)

    weights = ones(number_of_states)
    return weights
end


function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.LagrangianDuality)

    weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, regularization_regime.norm_lifted, duality_regime)
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.LagrangianDuality)

    weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁(), duality_regime)

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.LagrangianDuality)

    weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, regularization_regime.norm_lifted, duality_regime)

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.LagrangianDuality)

    weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁(), duality_regime)

    return
end

function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states::Int,
    norm_lifted::DynamicSDDiP.L₁, duality_regime::DynamicSDDiP.LagrangianDuality)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]

    # This means that the supremum norm is bounded in the dual
    for i in 1:number_of_states
        JuMP.set_upper_bound(π⁺[i], dual_bound * weights[i])
        JuMP.set_upper_bound(π⁻[i], dual_bound * weights[i])
    end

    return
end

function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states::Int,
    norm_lifted::DynamicSDDiP.L∞, duality_regime::DynamicSDDiP.LagrangianDuality)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]

    # This means that the 1-norm is bounded in the dual
    JuMP.@constraint(approx_model, sum(weights[i] * (π⁺[i] + π⁻[i]) for i in 1:number_of_states) <= dual_bound)

    return
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
under rigorous regime. Moreover, we log the number of times a specific lag_status occured.
"""
function lagrangian_status_check(
    model::SDDP.PolicyGraph,
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Rigorous,
    )

    if lag_status == :opt
        model.ext[:lag_status_dict][:opt] += 1
    elseif lag_status == :conv
        error("Lagrangian dual converged to value < solver_obj.")
    elseif lag_status == :sub
        error("Lagrangian dual had subgradients zero without LB=UB.")
    elseif lag_status == :iter
        error("Solving Lagrangian dual exceeded iteration limit.")
    elseif lag_status == :unbounded
        error("Normalized Lagrangian dual unbounded and reached artificial bound.")
    elseif lag_status == :bound_issues
        error("Lagrangian LB > UB due to numerical issues.")
    elseif lag_status == :feas_issues
        error("Lagrangian subproblem became infeasible.")
    elseif lag_status == :mn_opt
        model.ext[:lag_status_dict][:mn_opt] += 1
    elseif lag_status == :mn_iter
        error("Solving Lagrangian dual with minimal norm choice exceeded iteration limit.")
    elseif lag_status == :mn_issue
        error("Numerical issue with minimal norm choice. Proceeded without this step.")
        # TODO: Should this be an error?
    elseif lag_status == :subgr_stalling
        model.ext[:lag_status_dict][:subgr_stalling] += 1
    end

    return
end


"""
Trivial check of the status of the Lagrangian dual solution under lax regime.
Moreover, we log the number of times a specific lag_status occured.
"""
function lagrangian_status_check(
    model::SDDP.PolicyGraph,
    lag_status::Symbol,
    dual_status_regime::DynamicSDDiP.Lax,
    )

    if lag_status == :opt
        model.ext[:lag_status_dict][:opt] += 1
    elseif lag_status == :conv
        model.ext[:lag_status_dict][:conv] += 1
    elseif lag_status == :sub
        model.ext[:lag_status_dict][:sub] += 1
    elseif lag_status == :iter
        model.ext[:lag_status_dict][:iter] += 1
    elseif lag_status == :unbounded
        model.ext[:lag_status_dict][:unbounded] += 1
        #println("unbounded")
    elseif lag_status == :bound_issues
        model.ext[:lag_status_dict][:bound_issues] += 1
        #println("bound_issues")
    elseif lag_status == :feas_issues
        model.ext[:lag_status_dict][:feas_issues] += 1
    elseif lag_status == :mn_opt
        model.ext[:lag_status_dict][:mn_opt] += 1
    elseif lag_status == :mn_iter
        model.ext[:lag_status_dict][:mn_iter] += 1
    elseif lag_status == :mn_issue
        model.ext[:lag_status_dict][:mn_issue] += 1
    elseif lag_status == :subgr_stalling
        model.ext[:lag_status_dict][:subgr_stalling] += 1
    end

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
    undo_relax = JuMP.relax_integrality(subproblem)

    # Define appropriate solver
    reset_solver!(subproblem, algo_params, applied_solvers, :LP_relax, algo_params.solver_approach)

    # Solve LP Relaxation
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_LP_relax" begin
        JuMP.optimize!(subproblem)
    end
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
