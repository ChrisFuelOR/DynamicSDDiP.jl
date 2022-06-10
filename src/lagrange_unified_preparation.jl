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
    add_norm_constraints!(node, approx_model, weights, dual_bound, number_of_states, DynamicSDDiP.L₁(), DynamicSDDiP.LagrangianDuality())

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
    normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Optimal},
)

	π⁺ = approx_model[:π⁺]
	π⁻ = approx_model[:π⁻]
	π₀ = approx_model[:π₀]

	ω₀ = normalization_coeff.ω₀
	ω = normalization_coeff.ω

	JuMP.@constraint(approx_model, ω₀ * π₀ + sum(ω[i] * (π⁺[i] - π⁻[i]) for i in 1:number_of_states) <= 1)

	return
end


"""
This method determines the midpoint of the state space in the backward pass
based on the StateApproximation approach. This midpoint is required for some
approaches to determine core points and normalizations in the unified cut
generation framework.
"""
function get_state_space_midpoint(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    state_space_midpoint = 0.5 * ones(number_of_states)

    return state_space_midpoint
end

function get_state_space_midpoint(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    state_space_midpoint = Vector{Float64}(undef, number_of_states)

	for (i, (name, state)) in enumerate(node.states)
		variable_info = node.ext[:state_info_storage][name].in

		if variable_info.has_ub  & variable_info.has_lb
            state_space_midpoint[i] = variable_info.lower_bound + 0.5 * (variable_info.upper_bound - variable_info.lower_bound)
        elseif variable_info.binary
            state_space_midpoint[i] = 0.5
        else
            Error("All states require bounds to determine state space midpoint.")
        end
   end

	return state_space_midpoint
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Midpoint.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.Core_Midpoint,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get the midpoint of the state space
	core_point_x = get_state_space_midpoint(node, number_of_states, state_approximation_regime)

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.Core_Midpoint,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get the midpoint of the state space
	core_point_x = get_state_space_midpoint(node, number_of_states, state_approximation_regime)

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Epsilon.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.Core_Epsilon,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get the current incumbent and perturb it with the predefined value
	# into the interior of [0,1]
	core_point_x = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
		incumbent = JuMP.fix_value(state)

		if incumbent == 1
			core_point_x[i] = incumbent - normalization_regime.perturb
		elseif incumbent == 0
			core_point_x[i] = incumbent + normalization_regime.perturb
		end
	end

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.Core_Epsilon,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	core_point_x = zeros(number_of_states)

	# First get the midpoint of the state space
	midpoint = get_state_space_midpoint(node, number_of_states, state_approximation_regime)

	# Get the incumbent. Each component is (by choice) perturbed towards the
	# midpoint
	incumbent = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.states)
		incumbent = JuMP.fix_value(state.in)

		if incumbent > midpoint[i]
			core_point_x[i] = incumbent - normalization_regime.perturb
		elseif incumbent < midpoint[i]
			core_point_x[i] = incumbent + normalization_regime.perturb
		else
			core_point_x[i] = incumbent
		end
	end

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Epsilon.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.Core_In_Out,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	core_point_x = zeros(number_of_states)

	# Get the current iteration
	policy_graph = node.subproblem.ext[:sddp_policy_graph]
	iteration = policy_graph.ext[:iteration]

	# Get the current incumbent.
	incumbent = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
		incumbent .= JuMP.fix_value(state)
	end

	# Set the core point to the current incumbent in the first iteration,
	# otherwise use a linear combination with the previous core point
	if iteration == 1
		core_point_x .= incumbent
	else
		core_point_x .= 0.5 * incumbent + 0.5 * node.ext[:last_core_point]
	end

	# Store core point for following iterations
	node.ext[:last_core_point] = core_point_x

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

	return (x = core_point_x, theta = core_point_theta)
end

function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.Core_In_Out,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	core_point_x = zeros(number_of_states)

	# Get the current iteration
	policy_graph = node.subproblem.ext[:sddp_policy_graph]
	iteration = policy_graph.ext[:iteration]

	# Get the current incumbent.
	incumbent = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.states)
		incumbent .= JuMP.fix_value(state.in)
	end

	# Set the core point to the current incumbent in the first iteration,
	# otherwise use a linear combination with the previous core point
	if iteration == 1
		core_point_x .= incumbent
	else
		core_point_x .= 0.5 * incumbent + 0.5 * node.ext[:last_core_point]
	end

	# Store core point for following iterations
	node.ext[:last_core_point] = core_point_x

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, state_approximation_regime)

	return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Optimal.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
    state_approximation_regime::Union{DynamicSDDiP.NoStateApproximation,DynamicSDDiP.BinaryApproximation},
	normalization_regime::DynamicSDDiP.Core_Optimal,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get core point and corresponding optimal value
	core_point = evaluate_approx_value_function_no_fix(node, number_of_states, state_approximation_regime, copy_regime)

	return (x = core_point.x, theta = core_point.theta)
end



"""
Use core point to obtain normalization coefficients
"""
function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Optimal},
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get core point
	core_point = get_core_point(node, number_of_states, state_approximation_regime, normalization_regime, copy_regime)

	# Get theta direction
	ω₀ = core_point.theta - epi_state

	# Get state direction
	ω = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
		incumbent = JuMP.fix_value(state)
		ω[i] = core_point.x[i] - incumbent
	end

	return (ω = ω, ω₀ = ω₀)
end

function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out,DynamicSDDiP.Core_Optimal},
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	# Get core point
	core_point = get_core_point(node, number_of_states, state_approximation_regime, normalization_regime, copy_regime)

	# Get theta direction
	ω₀ = core_point.theta - epi_state

	# Get state direction
	ω = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.states)
		incumbent = JuMP.fix_value(state.in)
		ω[i] = core_point.x[i] - incumbent
	end

	return (ω = ω, ω₀ = ω₀)

end


"""
Trivial case where no core point is required
"""
function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	return
end

function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
    )

	return
end

"""
Evaluate the approximate value function for a given core point
"""
function evaluate_approx_value_function(
	node::SDDP.Node,
	core_point_x::Vector{Float64},
	number_of_states::Int,
	state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	)

	subproblem = node.subproblem
	original_state_values = Vector{Float64}(undef, number_of_states)

	# Replace the currently fixed value with the one from the core point
	# Store the original state value
	for (i, (name, state_comp)) in enumerate(node.states)
        #prepare_state_fixing!(node, state_comp)
		original_state_values[i] = JuMP.fix_value(state_comp.in)
        JuMP.fix(state_comp.in, core_point_x[i], force=true)
    end

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    core_obj = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

	# Restore the original state value
	for (i, (name, state_comp)) in enumerate(node.states)
        #prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, original_state_values[i], force=true)
    end

	return core_obj
end

function evaluate_approx_value_function(
	node::SDDP.Node,
	core_point_x::Vector{Float64},
	number_of_states::Int,
	state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	)

	subproblem = node.subproblem
	original_state_values = Vector{Float64}(undef, number_of_states)

	# Replace the currently fixed value with the one from the core point
	# Store the original state value
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        #prepare_state_fixing!(node, state_comp)
		original_state_values[i] = JuMP.fix_value(state)
        JuMP.fix(state, core_point_x[i], force=true)
    end

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    core_obj = JuMP.objective_value(subproblem)
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

	# Restore the original state value
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        #prepare_state_fixing!(node, state)
        JuMP.fix(state, original_state_values[i], force=true)
    end

	return core_obj
end



"""
Trivial case where no core point is required
"""
function evaluate_approx_value_function_no_fix(
	node::SDDP.Node,
	number_of_states::Int,
	state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
	)

	subproblem = node.subproblem
	original_state_values = Vector{Float64}(undef, number_of_states)
	core_point_x = Vector{Float64}(undef, number_of_states)

	# Unfix the state variable and store the original value
	for (i, (name, state)) in enumerate(node.states)
        original_state_values[i] = JuMP.fix_value(state.in)
        JuMP.unfix(state.in)

        # Set bounds and integer constraints based on copy_regime
        variable_info = node.ext[:state_info_storage][name].in
        follow_state_unfixing!(state, variable_info, copy_regime)
    end

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

	@assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

	# Get the optimal solution
    core_obj = JuMP.objective_value(subproblem)
	for (i, (name, state)) in enumerate(node.states)
        core_point_x[i] = JuMP.value(state.in)
    end

	# Restore the original state value
	for (i, (name, state_comp)) in enumerate(node.states)
        prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, original_state_values[i], force=true)
    end

	return (x = core_point_x,  theta = core_obj)
end

function evaluate_approx_value_function_no_fix(
	node::SDDP.Node,
	number_of_states::Int,
	state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
	)

	subproblem = node.subproblem
	original_state_values = Vector{Float64}(undef, number_of_states)
	core_point_x = Vector{Float64}(undef, number_of_states)

	# Unfix the state variable and store the original value
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        original_state_values[i] = JuMP.fix_value(state)
        JuMP.unfix(state)

        # Set bounds and integer constraints based on copy_regime
        variable_info = node.ext[:state_info_storage][name].in
        follow_state_unfixing!(state, variable_info, copy_regime)
    end

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

	@assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

	# Get the optimal solution
    core_obj = JuMP.objective_value(subproblem)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        core_point_x[i] = JuMP.value(state)
    end

	# Restore the original state value
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        prepare_state_fixing!(node, state)
        JuMP.fix(state, original_state_values[i], force=true)
    end

	return (x = core_point_x,  theta = core_obj)
end
