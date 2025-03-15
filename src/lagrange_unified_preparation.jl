function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, regularization_regime.norm_lifted, duality_regime)

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states::Int, dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁(), DynamicSDDiP.LagrangianDuality())

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
	Benders_symbol::Symbol,
    K::Int64,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    )

	full_Benders_list = node.ext[:Benders_cuts_binary]
	specific_Benders_list = Tuple{Int64, Symbol}[]

	for k in length(full_Benders_list):-1:1
		if length(specific_Benders_list) == K
			break
		end

		if full_Benders_list[k][2] == Benders_symbol
			push!(specific_Benders_list, full_Benders_list[k])
		end
	end

	number_of_Benders_cuts = length(specific_Benders_list)

	@assert number_of_Benders_cuts <= K
	if number_of_Benders_cuts == 0
    	@error("Dual space restriction attempt, but no Benders cuts have been generated.")
	end

    return (Benders_list = specific_Benders_list, number_of_Benders_cuts = number_of_Benders_cuts)

end

function get_Benders_list!(
    node::SDDP.Node,
	Benders_symbol::Symbol,
    K::Int64,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    )

	full_Benders_list = node.ext[:Benders_cuts_original]
	specific_Benders_list = Tuple{Int64, Symbol}[]

	for k in length(full_Benders_list):-1:1
		if length(specific_Benders_list) == K
			break
		end

		if full_Benders_list[k][2] == Benders_symbol
			push!(specific_Benders_list, full_Benders_list[k])
		end
	end

	number_of_Benders_cuts = length(specific_Benders_list)

	@assert number_of_Benders_cuts <= K
	if number_of_Benders_cuts == 0
    	@error("Dual space restriction attempt, but no Benders cuts have been generated.")
	end

    return (Benders_list = specific_Benders_list, number_of_Benders_cuts = number_of_Benders_cuts)

end

function dual_space_restriction!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    i::Int64,
	cut_aggregation_regime::DynamicSDDiP.SingleCutRegime,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.BendersSpanSpaceRestriction
)

    policy_graph = node.subproblem.ext[:sddp_policy_graph]
    ancestor_node = policy_graph.nodes[node.index-1]

    # Step 1: Get last K elements from Benders cuts
    ############################################################################
    K = dual_space_regime.K
    results = get_Benders_list!(ancestor_node, :all, K, state_approximation_regime)
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
        coefficients = ancestor_node.bellman_function.global_theta.cuts[index].coefficients
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
	cut_aggregation_regime::DynamicSDDiP.MultiCutRegime,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.BendersSpanSpaceRestriction
)

    policy_graph = node.subproblem.ext[:sddp_policy_graph]
    ancestor_node = policy_graph.nodes[node.index-1]

    # Step 1: Get last K elements from Benders cuts
    ############################################################################
    K = dual_space_regime.K
    results = get_Benders_list!(ancestor_node, Symbol(i), K, state_approximation_regime)
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
        coefficients = ancestor_node.bellman_function.local_thetas[i].cuts[index].coefficients
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
	cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation, DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.NoDualSpaceRestriction
)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
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
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
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
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
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
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.L₁∞_Deep
)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    π₀ = approx_model[:π₀]

	l₁_norm_exp = JuMP.@expression(approx_model, π₀ + sum(π⁺[i] + π⁻[i] for i in 1:number_of_states))

	sup_norm_var = JuMP.@variable(approx_model, sup_norm_aux)
	for i in 1:number_of_states
		JuMP.@constraint(approx_model, sup_norm_var >= π⁺[i])
		JuMP.@constraint(approx_model, sup_norm_var >= π⁻[i])
	end
	JuMP.@constraint(approx_model, sup_norm_var >= π₀)

    JuMP.@constraint(approx_model, 0.5 * l₁_norm_exp + 0.5 * sup_norm_var <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.ChenLuedtke
)

    @assert :span_variable in keys(JuMP.object_dictionary(approx_model))

    span_variable = approx_model[:span_variable]
    π₀ = approx_model[:π₀]

    abs_span_variable = JuMP.@variable(approx_model, abs_span_variable[1:length(span_variable)])
    JuMP.@constraint(approx_model, π₀ + sum(abs_span_variable[i] for i in 1:length(span_variable)) <= 1)
    JuMP.@constraint(approx_model, abs_span_1[i=1:length(span_variable)], -abs_span_variable[i] <= span_variable[i])
    JuMP.@constraint(approx_model, abs_span_2[i=1:length(span_variable)], abs_span_variable[i] >= span_variable[i])

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Relint,DynamicSDDiP.Core_Conv},
)

	π⁺ = approx_model[:π⁺]
	π⁻ = approx_model[:π⁻]
	π₀ = approx_model[:π₀]

	ω₀ = normalization_coeff.ω₀
	ω = normalization_coeff.ω

	JuMP.@constraint(approx_model, ω₀ * π₀ + sum(ω[i] * (π⁺[i] - π⁻[i]) for i in 1:number_of_states) <= 1)

	return
end