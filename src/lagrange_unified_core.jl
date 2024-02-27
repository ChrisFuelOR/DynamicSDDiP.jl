function integer_relax(
    subproblem::JuMP.Model,
    integer_regime::DynamicSDDiP.NoIntegerRelax
)
    function do_nothing()
        return
    end

    return do_nothing
end

function integer_relax(
    subproblem::JuMP.Model,
    integer_regime::DynamicSDDiP.IntegerRelax
)
    return JuMP.relax_integrality(subproblem)
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


function update_core_point_candidate(
    core_point_candidate,
    epi_state::Float64,
    primal_obj::Float64,
    improvement_regime::DynamicSDDiP.NoImprovement
)

    return core_point_candidate
end

function update_core_point_candidate(
    core_point_candidate,
    epi_state::Float64,
    primal_obj::Float64,
    improvement_regime::DynamicSDDiP.EpiState
)

    if core_point_candidate.theta < epi_state
        return (x = core_point_candidate.x, theta = epi_state)
    else
        return core_point_candidate
    end
end

function update_core_point_candidate(
    core_point_candidate,
    epi_state::Float64,
    primal_obj::Float64,
    improvement_regime::DynamicSDDiP.PrimalObj
)

    if core_point_candidate.theta < primal_obj
        return (x = core_point_candidate.x, theta = primal_obj)
    else
        return core_point_candidate
    end
end

"""
Use core point to obtain normalization coefficients
"""
function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    primal_obj::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out, DynamicSDDiP.Core_Relint},
    )

	# Get core point
	core_point_candidate = get_core_point(node, number_of_states, algo_params, applied_solvers, state_approximation_regime, normalization_regime)

    # Update core point candidate based on improvement regime
    core_point_candidate = update_core_point_candidate(core_point_candidate, epi_state, primal_obj, normalization_regime.improvement_regime)

	# Get theta direction
	ω₀ = core_point_candidate.theta - epi_state

	""" For Core_Epsilon with parameter 0.0 we do not allow negative values, since
	they are related to numerical issues. """
	if isa(normalization_regime, DynamicSDDiP.Core_Epsilon) && isapprox(normalization_regime.perturb, 0.0) && ω₀ < 0
		ω₀ = 0.0
	end

	# Get state direction
	ω = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
		incumbent = JuMP.fix_value(state)
		ω[i] = core_point_candidate.x[i] - incumbent
	end

    # Normalization of the direction if intended
	if normalization_regime.normalize_direction
		norm_factor = sum(abs.(ω)) + ω₀
    	ω = ω ./ norm_factor
    	ω₀ = ω₀ / norm_factor
	end

	return (ω = ω, ω₀ = ω₀)
end

function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    primal_obj::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out, DynamicSDDiP.Core_Relint},
    )

	# Get core point
	core_point_candidate = get_core_point(node, number_of_states, algo_params, applied_solvers, state_approximation_regime, normalization_regime)
    Infiltrator.@infiltrate
    
    # Update core point candidate based on improvement regime
    core_point_candidate = update_core_point_candidate(core_point_candidate, epi_state, primal_obj, normalization_regime.improvement_regime)

	# Get theta direction
	ω₀ = core_point_candidate.theta - epi_state

	""" For Core_Epsilon with parameter 0.0 we do not allow negative values, since
	they are related to numerical issues. """
	if isa(normalization_regime, DynamicSDDiP.Core_Epsilon) && isapprox(normalization_regime.perturb, 0.0) && ω₀ < 0
		ω₀ = 0.0
	end

	# TODO: We could exclude very small coefficients in general based on AlgoParams.atol
	# in order to avoid numerical issues.

	# Get state direction
	ω = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.states)
		incumbent = JuMP.fix_value(state.in)
		ω[i] = core_point_candidate.x[i] - incumbent
	end

	# Normalization of the direction if intended
	if normalization_regime.normalize_direction
		norm_factor = sum(abs.(ω)) + ω₀
    	ω = ω ./ norm_factor
    	ω₀ = ω₀ / norm_factor
	end

	# println(core_point, ", ", ω, ", ", ω₀)
    Infiltrator.@infiltrate
	
	return (ω = ω, ω₀ = ω₀)

end


"""
Trivial case where no core point is required
"""
function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    primal_obj::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    )

	return
end

function get_normalization_coefficients(
    node::SDDP.Node,
    number_of_states::Int,
	epi_state::Float64,
    primal_obj::Float64,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
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
    algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out},
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

    # Possible integer relaxation
    undo_relax = integer_relax(subproblem, normalization_regime.integer_regime)

	# RESET SOLVER (as it may have been changed in between for some reason)
    reset_solver!(subproblem, algo_params, applied_solvers, :norm, algo_params.solver_approach)

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL
	core_obj = JuMP.objective_value(subproblem)

	# Restore the original state value
	for (i, (name, state_comp)) in enumerate(node.states)
        #prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, original_state_values[i], force=true)
    end

	# Remove possible integer relaxation
	undo_relax()

	return core_obj
end

function evaluate_approx_value_function(
	node::SDDP.Node,
	core_point_x::Vector{Float64},
	number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
	normalization_regime::Union{DynamicSDDiP.Core_Midpoint,DynamicSDDiP.Core_Epsilon,DynamicSDDiP.Core_In_Out},
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

	# Possible integer relaxation
    undo_relax = integer_relax(subproblem, normalization_regime.integer_regime)

	# RESET SOLVER (as it may have been changed in between for some reason)
    reset_solver!(subproblem, algo_params, applied_solvers, :norm, algo_params.solver_approach)

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

    # Maybe attempt numerical recovery as in SDDP
    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL
    core_obj = JuMP.objective_value(subproblem)

	# Restore the original state value
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
        #prepare_state_fixing!(node, state)
        JuMP.fix(state, original_state_values[i], force=true)
    end

    # Remove possible integer relaxation
	undo_relax()	

	return core_obj
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Midpoint.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
	normalization_regime::DynamicSDDiP.Core_Midpoint,
    )

	# Get the midpoint of the state space
	core_point_x = get_state_space_midpoint(node, number_of_states, state_approximation_regime)
    Infiltrator.@infiltrate

	# Get optimal value to core point
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, algo_params, applied_solvers, normalization_regime, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Epsilon.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.Core_Epsilon,
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
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, algo_params, applied_solvers, normalization_regime, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.Core_Epsilon,
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
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, algo_params, applied_solvers, normalization_regime, state_approximation_regime)

    return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Epsilon.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	normalization_regime::DynamicSDDiP.Core_In_Out,
    )

	core_point_x = zeros(number_of_states)

	# Get the current iteration
	policy_graph = node.subproblem.ext[:sddp_policy_graph]
	iteration = policy_graph.ext[:iteration]

	# Get the current incumbent.
	incumbent = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.ext[:backward_data][:bin_states])
		incumbent[i] = JuMP.fix_value(state)
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
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, algo_params, applied_solvers, normalization_regime, state_approximation_regime)

	return (x = core_point_x, theta = core_point_theta)
end

function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
	applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	normalization_regime::DynamicSDDiP.Core_In_Out,
    )

	core_point_x = zeros(number_of_states)

	# Get the current iteration
	policy_graph = node.subproblem.ext[:sddp_policy_graph]
	iteration = policy_graph.ext[:iteration]

	# Get the current incumbent.
	incumbent = zeros(number_of_states)
	for (i, (_, state)) in enumerate(node.states)
		incumbent[i] = JuMP.fix_value(state.in)
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
	core_point_theta = evaluate_approx_value_function(node, core_point_x, number_of_states, algo_params, applied_solvers, normalization_regime, state_approximation_regime)

	return (x = core_point_x, theta = core_point_theta)
end

"""
Compute a core point based on the used StateApproximationRegime and the
NormalizationRegime Core_Relint.
"""
function get_core_point(
    node::SDDP.Node,
    number_of_states::Int,
	algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
	normalization_regime::DynamicSDDiP.Core_Relint,
    )

	subproblem = node.subproblem

	# Storage for data
	node.ext[:relint_data] = Dict{Symbol,Any}()
	relint_data = node.ext[:relint_data]
    relint_data[:original_state_values] = Dict{Symbol,Float64}()
    relint_data[:old_objective_expression] = JuMP.GenericAffExpr
	relint_data[:old_constraints] = Dict{Symbol,Any}[]
    relint_data[:relint_variables] = JuMP.VariableRef[]
    relint_data[:relint_constraints] = JuMP.ConstraintRef[]

	core_point_x = Vector{Float64}(undef, number_of_states)

	############################################################################
	# CONSTRUCT FEASIBILITY PROBLEM
	############################################################################
	vars = construct_feasibility_problem!(node, state_approximation_regime, normalization_regime)
	theta = vars.theta
	scaling_var = vars.scaling_var
	slack_variables = vars.slack_variables

	############################################################################
	# COMPUTE THE CORE POINT
	############################################################################
	# RESET SOLVER (as it may have been changed in between for some reason)
    reset_solver!(subproblem, algo_params, applied_solvers, :norm, algo_params.solver_approach)

	# Solve the subproblem
	TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_core" begin
        JuMP.optimize!(subproblem)
    end

	@assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

	# Get the optimal solution
    core_obj = JuMP.value(theta) / JuMP.value(scaling_var)
	for (i, (name, state)) in enumerate(node.states)
        core_point_x[i] = JuMP.value(slack_variables[i]) / JuMP.value(scaling_var)
    end

	############################################################################
	# REGAIN THE ORIGINAL PROBLEM
	############################################################################
	regain_original_problem!(node, theta, scaling_var, vars.undo_relax, state_approximation_regime)

	return (x = core_point_x, theta = core_obj)
end


"""
Method to construct the feasibility model to compute a relative interior point
of the epigraph.
"""
function construct_feasibility_problem!(
    node::SDDP.Node,
    state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
    normalization_regime::DynamicSDDiP.Core_Relint,
    )

    subproblem = node.subproblem
    relint_data = node.ext[:relint_data]

    ############################################################################
    # ADD SCALING VARIABLE
    ############################################################################
    scaling_var = JuMP.@variable(subproblem, scaling_var >= 1)

    ############################################################################
    # UNFIX THE STATE VARIABLES AND ADD SLACKS FOR BOUNDS
    ############################################################################
    prepare_state_variables_for_feasibility_model!(node, state_approximation_regime, normalization_regime.copy_regime)

    ############################################################################
    # ADD SLACKS FOR CONSTRAINTS CONTAINING STATE VARIABLES
    ############################################################################
    for state_constraint in node.ext[:state_constraints]
        # Create storage for constraint data
        old_constraints_storage = Dict{Symbol,Any}()

        # Store the constraint name
        old_constraints_storage[:name] = JuMP.name(state_constraint)

        # Create a slack variable
        slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
        push!(relint_data[:relint_variables], slack_var)

        # Extract and store the expression of the constraint
        expr = JuMP.constraint_object(state_constraint).func
        old_constraints_storage[:expression] = expr

        # Extract and store the RHS of the constraint
        rhs = JuMP.normalized_rhs(state_constraint)
        old_constraints_storage[:rhs] = rhs

        # Extract type of constraint
        constraint_type = typeof(JuMP.constraint_object(state_constraint).set)
        old_constraints_storage[:type] = constraint_type

        # Add new constraint
        if constraint_type == MOI.GreaterThan{Float64}
            slack_con = JuMP.@constraint(subproblem, expr + slack_var >= rhs * scaling_var)
        elseif constraint_type == MOI.LessThan{Float64}
            slack_con = JuMP.@constraint(subproblem, expr + slack_var <= rhs * scaling_var)
        elseif constraint_type == MOI.EqualTo{Float64}
            slack_con = JuMP.@constraint(subproblem, expr + slack_var == rhs * scaling_var)
        end
        push!(relint_data[:relint_constraints], slack_con)

        # Delete old constraint
        JuMP.delete(subproblem, state_constraint)

        # push to old_constraints storage
        push!(relint_data[:old_constraints], old_constraints_storage)

        """ Note: Alternatively, we could just store the old rhs and then set
        it to a sufficiently large value to turn off the constraint. We could
        then introduce the new constraint in addition. After the core point
        computation we only would have to delete the new constraint and change
        the rhs of the old one. However, this approach does only work for
        inequality constraints, so we stick with the above approach, even if it
        may take considerable amounts of time."""

    end

    ############################################################################
    # MODIFY THE OBJECTIVE
    ############################################################################
    # Introduce an epigraph variable
    theta = JuMP.@variable(subproblem)

    # Extract and store the current objective
    old_objective = JuMP.objective_function(subproblem)
    relint_data[:old_objective_expression] = old_objective
    relint_data[:old_objective_sense] = JuMP.objective_sense(subproblem)

    # Add a slack variable
    slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
    push!(relint_data[:relint_variables], slack_var)

    # Add new epigraph constraint to the model
    epi_constraint = JuMP.@constraint(subproblem, -theta + old_objective + slack_var <= 0 * scaling_var)
    push!(relint_data[:relint_constraints], epi_constraint)

    # Add new objective
    slack_variables = relint_data[:relint_variables]
    JuMP.@objective(subproblem, Max, sum(slack_variables[i] for i in 1:length(slack_variables)))

    ############################################################################
    # INTEGER RELAXATION
    ############################################################################
    undo_relax = integer_relax(subproblem, normalization_regime.integer_regime)

    return (theta = theta, scaling_var = scaling_var, slack_variables = slack_variables, undo_relax = undo_relax)
end


"""
Method to prepare the state variables for the feasibility model based on the
chosen StateApproximation.
"""

function prepare_state_variables_for_feasibility_model!(
	node::SDDP.Node,
	state_approximation_regime::DynamicSDDiP.BinaryApproximation,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
	)

	relint_data = node.ext[:relint_data]

	for (i, (name, state)) in enumerate(node.ext[:backward_data][:bin_states])
		relint_data[:original_state_values][name] = JuMP.fix_value(state)
		JuMP.unfix(state)

		# Set bounds and integer constraints based on copy_regime
		follow_state_unfixing_binary_slack!(node, state, copy_regime)
	end

	return
end

function prepare_state_variables_for_feasibility_model!(
	node::SDDP.Node,
	state_approximation_regime::DynamicSDDiP.NoStateApproximation,
	copy_regime::DynamicSDDiP.AbstractCopyRegime,
	)

	relint_data = node.ext[:relint_data]

	for (i, (name, state_comp)) in enumerate(node.states)
		relint_data[:original_state_values][name] = JuMP.fix_value(state_comp.in)
		JuMP.unfix(state_comp.in)

		# Set bounds and integer constraints based on copy_regime
		variable_info = node.ext[:state_info_storage][name].in
		follow_state_unfixing_slack!(node, state_comp, variable_info, copy_regime)
	end

	return
end

"""
Method to regain the original subproblems after computing a relative interior
point of the epigraph.
"""
function regain_original_problem!(
	node::SDDP.Node,
	theta::JuMP.VariableRef,
	scaling_var::JuMP.VariableRef,
    undo_relax::Any,
	state_approximation_regime::Union{DynamicSDDiP.BinaryApproximation,DynamicSDDiP.NoStateApproximation},
	)

	subproblem = node.subproblem
	relint_data = node.ext[:relint_data]

	############################################################################
	# UNDO INTEGER_RELAX
	############################################################################
    undo_relax()

	############################################################################
	# FIX STATE VARIABLES AGAIN
	############################################################################
	# Restore the original state value
	prepare_state_variables_for_original_model!(node, state_approximation_regime)

	############################################################################
	# RECREATE ORIGINAL CONSTRAINTS
	############################################################################
	for old_constraints_storage in relint_data[:old_constraints]
		con_name = old_constraints_storage[:name]
		expr = old_constraints_storage[:expression]
		rhs = old_constraints_storage[:rhs]
		con_type = old_constraints_storage[:type]

		if con_type == MOI.GreaterThan{Float64}
			con = JuMP.@constraint(subproblem, expr >= rhs)
		elseif con_type == MOI.LessThan{Float64}
			con = JuMP.@constraint(subproblem, expr <= rhs)
		elseif con_type == MOI.EqualTo{Float64}
			con = JuMP.@constraint(subproblem, expr == rhs)
		end
		JuMP.set_name(con, con_name)
	end

	############################################################################
	# RECREATE ORIGINAL OBJECTIVE
	############################################################################
	if relint_data[:old_objective_sense] == MOI.MIN_SENSE
		JuMP.@objective(subproblem, Min, relint_data[:old_objective_expression])
	else
		JuMP.@objective(subproblem, Max, relint_data[:old_objective_expression])
	end

	############################################################################
	# REMOVE ALL NEW CONSTRAINTS AND VARIABLES AGAIN
	############################################################################
	JuMP.delete(subproblem, relint_data[:relint_variables])
	JuMP.delete(subproblem, scaling_var)
	JuMP.unregister(subproblem, :scaling_var)
	JuMP.delete(subproblem, theta)

	for con in relint_data[:relint_constraints]
		JuMP.delete(subproblem, con)
	end

	############################################################################
	# REMOVE ALL CACHED DATA
	############################################################################
	delete!(node.ext, :relint_data)

	############################################################################
	# RESTORE LIST OF CONSTRAINTS CONTAINING STATE VARIABLES
	############################################################################
	identify_state_constraints!(node)

	return

end

"""
Method to prepare the state variables for the original model based on the
chosen StateApproximation.
"""

function prepare_state_variables_for_original_model!(
	node::SDDP.Node,
	state_approximation_regime::DynamicSDDiP.BinaryApproximation
	)

	relint_data = node.ext[:relint_data]

	for (i, (name, state)) in enumerate(node.ext[:backward_data][:bin_states])
		prepare_state_fixing_binary!(node, state)
		JuMP.fix(state, relint_data[:original_state_values][name], force=true)
	end

	return
end

function prepare_state_variables_for_original_model!(
	node::SDDP.Node,
	state_approximation_regime::DynamicSDDiP.NoStateApproximation
	)

	relint_data = node.ext[:relint_data]

	for (i, (name, state_comp)) in enumerate(node.states)
        prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, relint_data[:original_state_values][name], force=true)
    end

	return
end

"""
Auxiliary method to relax state variable bounds for NoStateApproximation.
"""
function follow_state_unfixing_slack!(node::SDDP.Node, state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.StateSpaceCopy)

	follow_state_unfixing_slack!(node, state, variable_info, DynamicSDDiP.ConvexHullCopy())

    if variable_info.binary
        JuMP.set_binary(state.in)
    elseif variable_info.integer
        JuMP.set_integer(state.in)
    end

    return
end

function follow_state_unfixing_slack!(node::SDDP.Node, state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.ConvexHullCopy)

	subproblem = node.subproblem
	relint_data = node.ext[:relint_data]
	scaling_var = subproblem[:scaling_var]

	if variable_info.has_lb
		lb = variable_info.lower_bound
	elseif variable_info.binary
		lb = 0
	else
		lb = -1e9
	end

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, -state.in + slack_var <= -lb * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)


	if variable_info.has_ub
		ub = variable_info.upper_bound
	elseif variable_info.binary
		ub = 1
	else
		ub = 1e9
	end

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, state.in + slack_var <= ub * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

    return
end

function follow_state_unfixing_slack!(node::SDDP.Node, state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.NoBoundsCopy)

	subproblem = node.subproblem
	relint_data = node.ext[:relint_data]
	scaling_var = subproblem[:scaling_var]

	# To avoid unboundedness
	lb = -1e9
	ub = 1e9

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, -state.in + slack_var <= -lb * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, state.in + slack_var <= ub * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

    return
end


"""
Auxiliary method to relax state variable bounds for BinaryApproximation.
"""
function follow_state_unfixing_binary_slack!(node::SDDP.Node, state::JuMP.VariableRef, copy_regime::DynamicSDDiP.StateSpaceCopy)

	follow_state_unfixing_binary_slack!(node, state, DynamicSDDiP.ConvexHullCopy())

    JuMP.set_binary(state)

    return
end

function follow_state_unfixing_binary_slack!(node::SDDP.Node, state::JuMP.VariableRef, copy_regime::DynamicSDDiP.ConvexHullCopy)

	subproblem = node.subproblem
	relint_data = node.ext[:relint_data]
	scaling_var = subproblem[:scaling_var]

	lb = 0
	ub = 1

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, -state + slack_var <= -lb * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, state + slack_var <= ub * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

    return
end

function follow_state_unfixing_binary_slack!(node::SDDP.Node, state::JuMP.VariableRef, copy_regime::DynamicSDDiP.NoBoundsCopy)

	subproblem = node.subproblem
	relint_data = node.ext[:relint_data]
	scaling_var = subproblem[:scaling_var]

	# TODO: To avoid unboundedness
	lb = 0
	ub = 1

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, -state + slack_var <= -lb * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

	# Add slack variable
	slack_var = JuMP.@variable(subproblem, lower_bound = 0, upper_bound = 1)
	push!(relint_data[:relint_variables], slack_var)

	# Add constraint for state variable bounds
	slack_con = JuMP.@constraint(subproblem, state + slack_var <= ub * scaling_var)
	push!(relint_data[:relint_constraints], slack_con)

    return
end

"""
Modifying the (regularized) subproblem to obtain an (approximate) primal to the Lagrangian dual
if no regularization and no state approximation is used
"""
function construct_unified_primal_problem!(
    node::SDDP.Node,
	node_index::Int64,
    subproblem::JuMP.Model,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    duality_regime::DynamicSDDiP.UnifiedLagrangianDuality,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    primal_data = node.ext[:primal_data]
	primal_data[:fixed_state_value] = Dict{Symbol,Float64}()
    primal_data[:primal_variables] = JuMP.VariableRef[]
    primal_data[:primal_constraints] = JuMP.ConstraintRef[]
    old_obj = primal_data[:old_objective] = JuMP.objective_function(subproblem)
    primal_data[:slacks] = Any[]
	number_of_states = 0

    ω₀ = normalization_coeff.ω₀
    ω = normalization_coeff.ω

    # INTRODUCE A NEW VARIABLE ETA
    ############################################################################
    eta = JuMP.@variable(subproblem, eta >= 0)
    push!(primal_data[:primal_variables], eta)

	# change the objective
	JuMP.set_objective(subproblem, JuMP.objective_sense(subproblem), eta)

    # UNFIX THE STATE VARIABLES (RELAXATION)
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        primal_data[:fixed_state_value][name] = JuMP.fix_value(state_comp.in)
        push!(primal_data[:slacks], primal_data[:fixed_state_value][name] - state_comp.in)
        JuMP.unfix(state_comp.in)
        variable_info = node.ext[:state_info_storage][name].in
        follow_state_unfixing!(state_comp, variable_info, copy_regime)
		number_of_states = i
    end
    slack = primal_data[:slacks]

    # INTRODUCE NEW CONSTRAINTS
    ############################################################################
    # objective constraint
    const_obj = JuMP.@constraint(subproblem, ω₀ * eta >= old_obj - epi_state)
    push!(primal_data[:primal_constraints], const_obj)

    # state constraint
    const_state = JuMP.@constraint(subproblem, [i=1:number_of_states], -ω[i] * eta == slack[i])
    append!(primal_data[:primal_constraints], const_state)

	# RELAX INTEGER IF REQUIRED
    ############################################################################
    undo_relax = integer_relax(subproblem, integer_regime)

    return undo_relax
end

"""
Modifying the (regularized) subproblem to obtain an (approximate) primal to the Lagrangian dual
if no regularization, but binary state approximation is used
"""
function construct_unified_primal_problem!(
    node::SDDP.Node,
	node_index::Int64,
    subproblem::JuMP.Model,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    duality_regime::UnifiedLagrangianDuality,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

	bw_data = node.ext[:backward_data]
    binary_states = bw_data[:bin_states]

    primal_data = node.ext[:primal_data]
	primal_data[:fixed_state_value] = Dict{Symbol,Float64}()
    primal_data[:primal_variables] = JuMP.VariableRef[]
    primal_data[:primal_constraints] = JuMP.ConstraintRef[]
    old_obj = primal_data[:old_objective] = JuMP.objective_function(subproblem)
    primal_data[:slacks] = Any[]
	number_of_states = 0

    ω₀ = normalization_coeff.ω₀
    ω = normalization_coeff.ω

    # INTRODUCE A NEW VARIABLE ETA
    ############################################################################
    eta = JuMP.@variable(subproblem, eta >= 0)
    push!(primal_data[:primal_variables], eta)

	# change the objective
	JuMP.set_objective(subproblem, JuMP.objective_sense(subproblem), eta)

    # UNFIX THE STATE VARIABLES (RELAXATION)
    ############################################################################
    for (i, (name, state_comp)) in enumerate(binary_states)
        primal_data[:fixed_state_value][name] = JuMP.fix_value(state_comp)
        push!(primal_data[:slacks], primal_data[:fixed_state_value][name] - state_comp)
        JuMP.unfix(state_comp)
        follow_state_unfixing_binary!(state_comp, copy_regime)
        number_of_states = i
    end
    slack = primal_data[:slacks]

    # INTRODUCE NEW CONSTRAINTS
    ############################################################################
    # objective constraint
    const_obj = JuMP.@constraint(subproblem, ω₀ * eta >= old_obj - epi_state)
    push!(primal_data[:primal_constraints], const_obj)

    # state constraint
    const_state = JuMP.@constraint(subproblem, [i=1:number_of_states], -ω[i] * eta == slack[i])
    append!(primal_data[:primal_constraints], const_state)

    # RELAX INTEGER IF REQUIRED
    ############################################################################
    undo_relax = integer_relax(subproblem, integer_regime)

    return undo_relax
end

"""
Modifying the (regularized) subproblem to obtain an (approximate) primal to the Lagrangian dual
if regularization, but no state approximation is used
"""
function construct_unified_primal_problem!(
    node::SDDP.Node,
	node_index::Int64,
    subproblem::JuMP.Model,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    duality_regime::UnifiedLagrangianDuality,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    primal_data = node.ext[:primal_data]
	primal_data[:fixed_state_value] = Dict{Symbol,Float64}()
    primal_data[:primal_variables] = JuMP.VariableRef[]
    primal_data[:primal_constraints] = JuMP.ConstraintRef[]
	primal_data[:reg_variables] = JuMP.VariableRef[]
	primal_data[:reg_constraints] = JuMP.ConstraintRef[]
    old_obj = primal_data[:old_objective] = JuMP.objective_function(subproblem)
    primal_data[:slacks] = Any[]
	number_of_states = 0

    ω₀ = normalization_coeff.ω₀
    ω = normalization_coeff.ω

    # INTRODUCE A NEW VARIABLE ETA
    ############################################################################
    eta = JuMP.@variable(subproblem, eta >= 0)
    push!(primal_data[:primal_variables], eta)

	# change the objective
	JuMP.set_objective(subproblem, JuMP.objective_sense(subproblem), eta)

    # UNFIX THE STATE VARIABLES (RELAXATION)
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        primal_data[:fixed_state_value][name] = JuMP.fix_value(state_comp.in)
        push!(primal_data[:slacks], primal_data[:fixed_state_value][name] - state_comp.in)
        JuMP.unfix(state_comp.in)
        variable_info = node.ext[:state_info_storage][name].in
        follow_state_unfixing!(state_comp, variable_info, regularization_regime.copy_regime) # normalization_regime.copy_regime?
        number_of_states = i
    end
    slack = primal_data[:slacks]

    # INTRODUCE NEW CONSTRAINTS
    ############################################################################
    # Variable for objective
    v = JuMP.@variable(subproblem, base_name = "reg_v")
    push!(reg_data[:reg_variables], v)

    # Get sign for regularization term
    fact = (JuMP.objective_sense(subproblem) == JuMP.MOI.MIN_SENSE ? 1 : -1)

	# objective constraint
    const_obj = JuMP.@constraint(subproblem, ω₀ * eta >= old_obj + fact * regularization_regime.sigma[node_index] * v - epi_state)
    push!(primal_data[:primal_constraints], const_obj)

	# Variables for absolute values
    alpha = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "alpha")
    append!(primal_data[:reg_variables], alpha)

	# Variable for regularization (used instead of fixed value here)
	x_aux = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "x_aux")
    append!(primal_data[:primal_variables], x_aux)

    # Constraints for absolute values
    const_plus = JuMP.@constraint(subproblem, [i=1:number_of_states], -alpha[i] <= x_aux[i] - state_comp.in[i])
    const_minus = JuMP.@constraint(subproblem, [i=1:number_of_states], slack[i] <= x_aux[i] - state_comp.in[i])
    append!(primal_data[:reg_constraints], const_plus)
    append!(primal_data[:reg_constraints], const_minus)

	# Add norm specific constraint
	add_norm_constraint!(subproblem, v, alpha, primal_data, number_of_states, regularization_regime.norm)

    # state constraint
    const_state = JuMP.@constraint(subproblem, [i=1:number_of_states], -ω[i] * eta == primal_data[:fixed_state_value][name][i] - x_aux[i])
    append!(primal_data[:primal_constraints], const_state)

    # RELAX INTEGER IF REQUIRED
    ############################################################################
    undo_relax = integer_relax(subproblem, integer_regime)

    return undo_relax
end


"""
Modifying the (regularized) subproblem to obtain an (approximate) primal to the Lagrangian dual
if regularization and binary state approximation is used
"""
function construct_unified_primal_problem!(
    node::SDDP.Node,
	node_index::Int64,
    subproblem::JuMP.Model,
    epi_state::Float64,
    normalization_coeff::Union{Nothing,NamedTuple{(:ω, :ω₀),Tuple{Vector{Float64},Float64}}},
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    duality_regime::UnifiedLagrangianDuality,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

	bw_data = node.ext[:backward_data]
    binary_states = bw_data[:bin_states]

    primal_data = node.ext[:primal_data]
    primal_data[:fixed_state_value] = Dict{Symbol,Float64}()
    primal_data[:primal_variables] = JuMP.VariableRef[]
    primal_data[:primal_constraints] = JuMP.ConstraintRef[]
	primal_data[:reg_variables] = JuMP.VariableRef[]
	primal_data[:reg_constraints] = JuMP.ConstraintRef[]
    old_obj = primal_data[:old_objective] = JuMP.objective_function(subproblem)
    primal_data[:slacks] = Any[]
    primal_data[:weights] = Float64[]
    sigma_bin = regularization_regime.sigma[node_index]
	number_of_states = 0

    ω₀ = normalization_coeff.ω₀
    ω = normalization_coeff.ω

    # INTRODUCE A NEW VARIABLE ETA
    ############################################################################
    eta = JuMP.@variable(subproblem, eta >= 0)
    push!(primal_data[:primal_variables], eta)

	# change the objective
	JuMP.set_objective(subproblem, JuMP.objective_sense(subproblem), eta)

	# UNFIX THE STATE VARIABLES (RELAXATION)
    ############################################################################
    for (i, (name, state_comp)) in enumerate(binary_states)
        primal_data[:fixed_state_value][name] = JuMP.fix_value(state_comp)
        push!(primal_data[:slacks], primal_data[:fixed_state_value][name] - state_comp)
        JuMP.unfix(state_comp)
        follow_state_unfixing_binary!(state_comp, duality_regime.copy_regime)

		# determine and store the corresponding weight
		associated_original_state = node.ext[:backward_data][:bin_x_names][name]
		beta = state_approximation_regime.binary_precision[associated_original_state]
		associated_k = node.ext[:backward_data][:bin_k][name]
		push!(reg_data[:weights], 2^(associated_k-1) * beta)

        number_of_states = i
    end
    slack = primal_data[:slacks]

    # INTRODUCE NEW CONSTRAINTS
    ############################################################################
    # Variable for objective
    v = JuMP.@variable(subproblem, base_name = "reg_v")
    push!(reg_data[:reg_variables], v)

    # Get sign for regularization term
    fact = (JuMP.objective_sense(subproblem) == JuMP.MOI.MIN_SENSE ? 1 : -1)

	# objective constraint
    const_obj = JuMP.@constraint(subproblem, ω₀ * eta >= old_obj + fact * regularization_regime.sigma[node_index] * v - epi_state)
    push!(primal_data[:primal_constraints], const_obj)

	# Variables for absolute values
    alpha = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "alpha")
    append!(primal_data[:reg_variables], alpha)

	# Variable for regularization (used instead of fixed value here)
	x_aux = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "x_aux")
    append!(primal_data[:primal_variables], x_aux)

    # Constraints for absolute values
    const_plus = JuMP.@constraint(subproblem, [i=1:number_of_states], -alpha[i] <= x_aux[i] - state_comp.in[i])
    const_minus = JuMP.@constraint(subproblem, [i=1:number_of_states], slack[i] <= x_aux[i] - state_comp.in[i])
    append!(primal_data[:reg_constraints], const_plus)
    append!(primal_data[:reg_constraints], const_minus)

	# Add norm specific constraint
	add_norm_constraint_binary!(subproblem, v, alpha, primal_data, number_of_states, regularization_regime.norm_lifted)

    # state constraint
    const_state = JuMP.@constraint(subproblem, [i=1:number_of_states], -ω[i] * eta == primal_data[:fixed_state_value][name][i] - x_aux[i])
    append!(primal_data[:primal_constraints], const_state)

    # RELAX INTEGER IF REQUIRED
    ############################################################################
    undo_relax = integer_relax(subproblem, integer_regime)

    return undo_relax
end

"""
Re-modifying the (regularized) primal to the Lagrangian dual
if no regularization and binary state approximation is used
"""
function deconstruct_unified_primal_problem!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    undo_relax::Any,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    primal_data = node.ext[:primal_data]

    # UNDO INTEGER RELAX IF REQUIRED
    ############################################################################
    undo_relax()

    # FIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, primal_data[:fixed_state_value][name], force=true)
    end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, primal_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    JuMP.delete(subproblem, primal_data[:primal_variables])
	JuMP.unregister(subproblem, :eta)

    for constraint in primal_data[:primal_constraints]
        JuMP.delete(subproblem, constraint)
    end

    delete!(node.ext, :primal_data)

    return
end


"""
Re-modifying the (regularized) primal to the Lagrangian dual
if no regularization, but binary state approximation is used
"""
function deconstruct_unified_primal_problem!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    undo_relax::Any,
    regularization_regime::DynamicSDDiP.NoRegularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    primal_data = node.ext[:primal_data]
	bw_data = node.ext[:backward_data]

    # UNDO INTEGER RELAX IF REQUIRED
    ############################################################################
    undo_relax()

	# FIX THE STATE VARIABLES
	############################################################################
	for (i, (name, state_comp)) in enumerate(bw_data[:bin_states])
		prepare_state_fixing_binary!(node, state_comp)
		JuMP.fix(state_comp, reg_data[:fixed_state_value][name], force=true)
	end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, primal_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    JuMP.delete(subproblem, primal_data[:primal_variables])
	JuMP.unregister(subproblem, :eta)

    for constraint in primal_data[:primal_constraints]
        JuMP.delete(subproblem, constraint)
    end

    delete!(node.ext, :primal_data)

    return
end

"""
Re-modifying the (regularized) primal to the Lagrangian dual
if regularization, but no binary state approximation is used
"""
function deconstruct_unified_primal_problem!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    undo_relax::Any,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    primal_data = node.ext[:primal_data]

    # UNDO INTEGER RELAX IF REQUIRED
    ############################################################################
    undo_relax()

    # FIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, primal_data[:fixed_state_value][name], force=true)
    end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, primal_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    JuMP.delete(subproblem, primal_data[:primal_variables])
	JuMP.delete(subproblem, primal_data[:reg_variables])
	JuMP.unregister(subproblem, :eta)

    for constraint in primal_data[:primal_constraints]
        JuMP.delete(subproblem, constraint)
    end

	for constraint in primal_data[:reg_constraints]
        JuMP.delete(subproblem, constraint)
    end

    delete!(node.ext, :primal_data)

    return
end


"""
Re-modifying the (regularized) primal to the Lagrangian dual
if regularization and binary state approximation are used
"""
function deconstruct_unified_primal_problem!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    undo_relax::Any,
    regularization_regime::DynamicSDDiP.Regularization,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    primal_data = node.ext[:primal_data]
	bw_data = node.ext[:backward_data]

    # UNDO INTEGER RELAX IF REQUIRED
    ############################################################################
    undo_relax()

	# FIX THE STATE VARIABLES
	############################################################################
	for (i, (name, state_comp)) in enumerate(bw_data[:bin_states])
		prepare_state_fixing_binary!(node, state_comp)
		JuMP.fix(state_comp, reg_data[:fixed_state_value][name], force=true)
	end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, primal_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    JuMP.delete(subproblem, primal_data[:primal_variables])
	JuMP.delete(subproblem, primal_data[:reg_variables])
	JuMP.unregister(subproblem, :eta)

    for constraint in primal_data[:primal_constraints]
        JuMP.delete(subproblem, constraint)
    end

	for constraint in primal_data[:reg_constraints]
        JuMP.delete(subproblem, constraint)
    end

    delete!(node.ext, :primal_data)

    return
end
