function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, regularization_regime.norm_lifted, duality_regime)

    return
end


function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁, DynamicSDDiP.LagrangianDuality())

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, regularization_regime.norm_lifted, duality_regime)

	return

end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁, DynamicSDDiP.LagrangianDuality())

    return
end


function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states::Int,
    norm_lifted::DynamicSDDiP.L₁, duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
	π₀ = approx_model[:π₀]

    # This means that the supremum norm is bounded in the dual
    for i in 1:number_of_states
        JuMP.@constraint(approx_model, π⁺[i] <= dual_bound * weights[i] * π₀)
        JuMP.@constraint(approx_model, π⁻[i] <= dual_bound * weights[i] * π₀)
    end

	return
end

function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states::Int,
    norm_lifted::DynamicSDDiP.L∞, duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
	π₀ = approx_model[:π₀]

	# This means that the 1-norm is bounded in the dual
	JuMP.@constraint(approx_model, sum(weights[i] * (π⁺[i] + π⁻[i]) for i in 1:number_of_states) <= dual_bound * π₀)

	return
end


function get_Benders_list!(
    node::SDDP.Node,
    K::Int64,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

    number_of_Benders_cuts = length(node.ext[:Benders_cuts_binary])

    if number_of_Benders_cuts == 0
            @error("Dual space restriction attempt, but no Benders cuts have been generated.")
    elseif  number_of_Benders_cuts <= K
            Benders_list = node.ext[:Benders_cuts_binary]
    else
            Benders_list = node.ext[:Benders_cuts_binary][number_of_Benders_cuts-K:number_of_Benders_cuts]
    end

    return (Benders_list = Benders_list, number_of_Benders_cuts = number_of_Benders_cuts)

end

function get_Benders_list!(
    node::SDDP.Node,
    K::Int64,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

    number_of_Benders_cuts = length(node.ext[:Benders_cuts_original])

    if number_of_Benders_cuts == 0
            @error("Dual space restriction attempt, but no Benders cuts have been generated.")
    elseif  number_of_Benders_cuts <= K
            Benders_list = node.ext[:Benders_cuts_original]
    else
            Benders_list = node.ext[:Benders_cuts_original][number_of_Benders_cuts-K:number_of_Benders_cuts]
    end

    return (Benders_list = Benders_list, number_of_Benders_cuts = number_of_Benders_cuts)

end

function dual_space_restriction!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    i::Int64,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.BendersSpanSpaceRestriction
)

    policy_graph = node.subproblem.ext[:sddp_policy_graph]
    ancestor_node = policy_graph.nodes[node.index-1]

    # Step 1: Get last K elements from Benders cuts
    ############################################################################
    K = dual_space_regime.K
    results = get_Benders_list!(ancestor_node, K, state_approximation_regime)
    Benders_list = results.Benders_list
    number_of_Benders_cuts = results.number_of_Benders_cuts

    # Step 2: Introduce a span variable
    ############################################################################
    span_variable = JuMP.@variable(approx_model, span_variable[1:number_of_Benders_cuts])

    # Step 3: Create span expression
    ############################################################################
    expr = JuMP.@expression(approx_model, 0)

    for k in 1:number_of_Benders_cuts
        coefficient = 0
        index = Benders_list[k][1]
        if Benders_list[k][2] == :single_cut
            coefficients = ancestor_node.bellman_function.global_theta.cuts[index].coefficients
        elseif Benders_list[k][2] == :multi_cut
            coefficients = ancestor_node.bellman_function.local_thetas[i].cuts[index].coefficients
        end
        coefficients_vector = zeros(length(coefficients))
        for (j, (key,value)) in enumerate(coefficients)
            coefficients_vector[j] = value
        end
        expr = JuMP.@expression(approx_model, expr .+ span_variable[k] * coefficients_vector)
    end

    # Step 4: Introduce a new constraint to restrict the dual space
    ############################################################################
    π = approx_model[:π]
    span_constraint = JuMP.@constraint(approx_model, π .== expr)

    return
end

function dual_space_restriction!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    i::Int64,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation, DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.NoDualSpaceRestriction
)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.L₁_Deep
)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    π₀ = approx_model[:π₀]

    JuMP.@constraint(approx_model, π₀ + sum(π⁺[i] + π⁻[i] for i in 1:number_of_states) <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.L₂_Deep
)

    π = approx_model[:π]
    π₀ = approx_model[:π₀]
    π_all = [π, π₀]

    #JuMP.@constraint(approx_model, [1, π_all] in SecondOrderCone())
    #JuMP.@NLconstraint(approx_model, sqrt(π₀^2 + sum(π[i]^2 for i in 1:number_of_states)) <= 1)
    JuMP.@constraint(approx_model, π₀^2 + sum(π[i]^2 for i in 1:number_of_states) <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.L∞_Deep
)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    π₀ = approx_model[:π₀]

    for i in 1:number_of_states
        JuMP.set_upper_bound(π⁺[i], 1)
        JuMP.set_upper_bound(π⁻[i], 1)
    end
    JuMP.set_upper_bound(π₀, 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.Fischetti
)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    π₀ = approx_model[:π₀]

    JuMP.@constraint(approx_model, π₀ + sum(π⁺[i] - π⁻[i] for i in 1:number_of_states) <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.Brandenberg
)

    # TODO

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
    normalization_regime::DynamicSDDiP.ChenLuedtke
)

    @assert :span_variable in keys(object_dictionary(approx_model))

    span_variable = approx_model[:span_variable]
    π₀ = approx_model[:π₀]

    abs_span_variable = JuMP.@variable(approx_model, abs_span_variable[1:length(span_variable)])
    JuMP.@constraint(approx_model, π₀ + sum(abs_span_variable[i] for i in 1:length(span_variable)) <= 1)
    JuMP.@constraint(approx_model, abs_span_1[i=1:length(span_variable)], -abs_span_variable[i] <= span_variable[i])
    JuMP.@constraint(approx_model, abs_span_2[i=1:length(span_variable)], abs_span_variable[i] >= span_variable[i])

    return
end
