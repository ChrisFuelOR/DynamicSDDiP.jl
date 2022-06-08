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
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁, duality_regime)

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
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁, duality_regime)

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
