function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
    number_of_states_per_noise::Int,
	number_of_noise::Int,
	dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states_per_noise, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states_per_noise, number_of_noise, regularization_regime.norm_lifted, duality_regime)

    return
end


function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
	number_of_states_per_noise::Int,
	number_of_noise::Int,
	dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states_per_noise, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states_per_noise, number_of_noise, DynamicSDDiP.L₁(), DynamicSDDiP.LagrangianDuality())

    return
end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
	number_of_states_per_noise::Int,
	number_of_noise::Int,
	dual_bound::Float64,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states_per_noise, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states_per_noise, number_of_noise, regularization_regime.norm_lifted, duality_regime)

	return

end

function set_multiplier_bounds!(node::SDDP.Node, approx_model::JuMP.Model,
	number_of_states_per_noise::Int,
	number_of_noise::Int,
	dual_bound::Float64,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

	weights = determine_weights!(node, approx_model, number_of_states_per_noise, regularization_regime, state_approximation_regime)
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states_per_noise, number_of_noise, DynamicSDDiP.L₁(), DynamicSDDiP.UnifiedLagrangianDuality())

    return
end


function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states_per_noise::Int,
	number_of_noise::Int, norm_lifted::DynamicSDDiP.L₁,
	duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

    π = approx_model[:π]
	π₀ = approx_model[:π₀]

	# TODO: REGULARIZED CASE HAS TO BE CHECKED!!!
	JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise])

    # This means that the supremum norm is bounded in the dual
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] <= dual_bound * weights[i] * π₀)
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], -abs_aux[i] <= sum(π[j,i] for j in 1:number_of_noise))
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] >= sum(π[j,i] for j in 1:number_of_noise))

	return
end

function add_norm_constraints!(node::SDDP.Node, approx_model::JuMP.Model,
    weights::Vector{Float64}, dual_bound::Float64, number_of_states_per_noise::Int,
	number_of_noise::Int, norm_lifted::DynamicSDDiP.L∞,
	duality_regime::DynamicSDDiP.UnifiedLagrangianDuality)

    π = approx_model[:π]
	π₀ = approx_model[:π₀]

	# TODO: REGULARIZED CASE HAS TO BE CHECKED!!!
	JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise])

	# This means that the 1-norm is bounded in the dual
	JuMP.@constraint(approx_model, sum(weights[i] * abs_aux[i] for i in 1:number_of_states_per_noise) <= dual_bound * π₀)
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], -abs_aux[i] <= sum(π[j,i] for j in 1:number_of_noise))
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] >= sum(π[j,i] for j in 1:number_of_noise))

	return
end


  function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
	number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.L₁_Deep
)

	π = approx_model[:π]
	π₀ = approx_model[:π₀]

	#JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise])
	#JuMP.@constraint(approx_model, sum(abs_aux[i] for i in 1:number_of_states_per_noise) <= 1)
	#JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], -abs_aux[i] <= sum(π[j,i] for j in 1:number_of_noise))
    #JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] >= sum(π[j,i] for j in 1:number_of_noise))

	JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise, j=1:number_of_noise])
	JuMP.@constraint(approx_model, sum(abs_aux[i,j] for i in 1:number_of_states_per_noise for j in 1:number_of_noise) <= 1)
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise, j in 1:number_of_noise], -abs_aux[i,j] <= π[j,i])
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise, j in 1:number_of_noise], abs_aux[i,j] >= π[j,i])


    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.L₂_Deep
)

    π = approx_model[:π]
    π₀ = approx_model[:π₀]
    π_all = [π, π₀]

    JuMP.@constraint(approx_model, π₀^2 + sum( (sum(π[j,i] for j in 1:number_of_noise))^2 for i in 1:number_of_states_per_noise) <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.L∞_Deep
)

    π = approx_model[:π]
    π₀ = approx_model[:π₀]

	JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise])
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] <= 1)
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], -abs_aux[i] <= sum(π[j,i] for j in 1:number_of_noise))
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] >= sum(π[j,i] for j in 1:number_of_noise))
	JuMP.set_upper_bound(π₀, 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.L₁∞_Deep
)

    π = approx_model[:π]
    π₀ = approx_model[:π₀]

	############################################################################
	# Modeling of absolute values
	JuMP.@variable(approx_model, abs_aux[i=1:number_of_states_per_noise])
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], -abs_aux[i] <= sum(π[j,i] for j in 1:number_of_noise))
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], abs_aux[i] >= sum(π[j,i] for j in 1:number_of_noise))

	############################################################################
	# Part for 1-norm
	l₁_norm_exp = JuMP.@expression(approx_model, π₀ + sum(abs_aux[i] for i in 1:number_of_states_per_noise))

	############################################################################
	# Part for supremum norm
	sup_norm_var = JuMP.@variable(approx_model, sup_norm_aux)

	JuMP.@constraint(approx_model, sup_norm_var >= π₀)
	JuMP.@constraint(approx_model, [i in 1:number_of_states_per_noise], sup_norm_var >= abs_aux[i])

	############################################################################
    JuMP.@constraint(approx_model, 0.5 * l₁_norm_exp + 0.5 * sup_norm_var <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.Fischetti
)

    π⁺ = approx_model[:π⁺]
    π⁻ = approx_model[:π⁻]
    π₀ = approx_model[:π₀]

    JuMP.@constraint(approx_model, -π₀ - sum(π⁺[j,i] - π⁻[j,i] for j in 1:number_of_noise for i in 1:number_of_states_per_noise) <= 1)

    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::DynamicSDDiP.ChenLuedtke
)

    @assert :span_variable in keys(JuMP.object_dictionary(approx_model))

	span_variable = approx_model[:span_variable]
	π₀ = approx_model[:π₀]

	@assert size(span_variable) == number_of_states_per_noise

    abs_span_variable = JuMP.@variable(approx_model, abs_span_variable[1:number_of_states_per_noise])
    JuMP.@constraint(approx_model, π₀ + sum(abs_span_variable[i] for i in 1:number_of_states_per_noise) <= 1)
    JuMP.@constraint(approx_model, abs_span_1[i=1:number_of_states_per_noise], -abs_span_variable[i] <= span_variable[i])
    JuMP.@constraint(approx_model, abs_span_2[i=1:number_of_states_per_noise], abs_span_variable[i] >= span_variable[i])
    return
end

function add_normalization_constraint!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    number_of_states_per_noise::Int, number_of_noise::Int,
	normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Optimal,DynamicSDDiP.Core_Relint},
)

	π = approx_model[:π]
	π₀ = approx_model[:π₀]

	ω₀ = normalization_coeff.ω₀
	ω = normalization_coeff.ω

	JuMP.@constraint(approx_model, ω₀ * π₀ + sum(ω[i] * sum(π[j,i] for j in 1:number_of_noise) for i in 1:number_of_states_per_noise) <= 1)

	return
end

"""
Get normalization coefficients for pseudo-deep cuts in the single-cut case
"""
function get_normalization_coefficients_agg(
    node::SDDP.Node,
    number_of_states_per_noise::Int, number_of_noise::Int,
	epi_state::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
	state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
	normalization_regime::AbstractNormalizationRegime,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

return

end


"""
Get normalization coefficients for pseudo-deep cuts in the single-cut case

Note that we have to use the single epi_state for the expected value function
for each of the realizations here.
"""
function get_normalization_coefficients_agg(
    node::SDDP.Node,
    number_of_states_per_noise::Int, number_of_noise::Int,
	epi_state::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Optimal, DynamicSDDiP.Core_Relint},
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Initialize direction
	ω₀_agg = 0
	ω_agg = zeros(number_of_states_per_noise)

    for noise in SDDP.sample_backward_noise_terms(backward_sampling_scheme, node)
		# Parameterize the subproblem for the current realization
		parameterize(node, noise.term)

		# Get core point and based on that normalization direction for the current realization
		current_direction = get_normalization_coefficients(
			node,
			number_of_states_per_noise,
			epi_state,
			algo_params,
			state_approximation_regime,
			normalization_regime,
			copy_regime,
		)

		# Update the aggregate direction using the probability
		ω_agg .= ω_agg + noise.probability * current_direction.ω
		ω₀_agg = ω₀_agg + noise.probability * current_direction.ω₀

    end

return (ω = ω_agg, ω₀ = ω₀_agg)

end

function dual_space_restriction_agg!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    i::Int64,
	number_of_states_per_noise::Int, number_of_noise::Int,
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
    number_of_Benders_cuts = results.number_of_Benders_cuts # K or less

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
    span_constraint = JuMP.@constraint(approx_model, [1:number_of_states_per_noise], sum(π[j,i] for j in number_of_noise) == expr[i])

    # Infiltrator.@infiltrate

    return
end

function dual_space_restriction_agg!(
    node::SDDP.Node,
    approx_model::JuMP.Model,
    i::Int64,
	number_of_states_per_noise::Int, number_of_noise::Int,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation, DynamicSDDiP.NoStateApproximation},
    dual_space_regime::DynamicSDDiP.NoDualSpaceRestriction
)

    return
end
