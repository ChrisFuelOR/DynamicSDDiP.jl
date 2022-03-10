# The functions
# > "BellmanFunction",
# > "bellman_term",
# > "initialize_bellman_function"
# > "refine_bellman_function"
# > "_add_average_cut"
# > "_add_cut"
# > "_add_cut_constraints_to_models"
# and structs
# > "SampledState",
# > "LevelOneOracle",
# > "CutApproximation"
# > "BellmanFunction"
# are derived from similar named functions and structs in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

################################################################################
# DEFINING STRUCTURES TO STORE CUTS
################################################################################
mutable struct SampledState
    state::Dict{Symbol,Float64}
    dominating_cut::DynamicSDDiP.Cut
    best_objective::Float64
    #obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    #belief_y::Union{Nothing,Dict{T,Float64} where {T}}
end

mutable struct CutApproximation
    theta::JuMP.VariableRef
    states::Dict{Symbol,JuMP.VariableRef}
    # objective_states::Union{Nothing,NTuple{N,JuMP.VariableRef} where {N}}
    # belief_states::Union{Nothing,Dict{T,JuMP.VariableRef} where {T}}
    # Storage for cut selection
    cuts::Vector{DynamicSDDiP.Cut}
    sampled_states::Vector{DynamicSDDiP.SampledState}
    cuts_to_be_deleted::Vector{DynamicSDDiP.Cut}
    deletion_minimum::Int

    function CutApproximation(
        theta::JuMP.VariableRef,
        states::Dict{Symbol,JuMP.VariableRef},
        # objective_states,
        # belief_states,
        deletion_minimum::Int,
    )
        return new(
            theta,
            states,
            # objective_states,
            # belief_states,
            DynamicSDDiP.Cut[],
            DynamicSDDiP.SampledState[],
            DynamicSDDiP.Cut[],
            deletion_minimum
        )
    end
end

mutable struct BellmanFunction
    global_theta::CutApproximation
    local_thetas::Vector{CutApproximation}
    cut_type::SDDP.CutType
    # Cuts defining the dual representation of the risk measure.
    risk_set_cuts::Set{Vector{Float64}}
end

"""
Required to define Bellman functions.
"""
function BellmanFunction(;
    lower_bound = -Inf,
    upper_bound = Inf,
    deletion_minimum::Int = 1,
    cut_type::SDDP.CutType = SDDP.SINGLE_CUT,
)
    return SDDP.InstanceFactory{BellmanFunction}(
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        deletion_minimum = deletion_minimum,
        cut_type = cut_type,
    )
end

function bellman_term(bellman_function::DynamicSDDiP.BellmanFunction)
    return bellman_function.global_theta.theta
end

################################################################################
# INITIALIZE BELLMAN FUNCTION
################################################################################

"""
Initializing the bellman function for the subproblems.
"""
function initialize_bellman_function(
    factory::SDDP.InstanceFactory{BellmanFunction},
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
) where {T}
    lower_bound, upper_bound, deletion_minimum, cut_type = -Inf, Inf, 0, SDDP.SINGLE_CUT
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in BellmanFunction.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :deletion_minimum
            deletion_minimum = value
        elseif kw == :cut_type
            cut_type = value
        else
            error("Keyword $(kw) not recognised as argument to BellmanFunction.")
        end
    end
    if lower_bound == -Inf && upper_bound == Inf
        error("You must specify a finite bound on the cost-to-go term.")
    end
    if length(node.children) == 0
        lower_bound = upper_bound = 0.0
    end
    Θᴳ = JuMP.@variable(node.subproblem, base_name="Θᴳ")
    lower_bound > -Inf && JuMP.set_lower_bound(Θᴳ, lower_bound)
    upper_bound < Inf && JuMP.set_upper_bound(Θᴳ, upper_bound)
    # Initialize bounds for the objective states. If objective_state==nothing,
    # this check will be skipped by dispatch.
    # SDDP._add_initial_bounds(node.objective_state, Θᴳ)
    x′ = Dict(key => var.out for (key, var) in node.states)
    ## obj_μ = node.objective_state !== nothing ? node.objective_state.μ : nothing
    ## belief_μ = node.belief_state !== nothing ? node.belief_state.μ : nothing
    return BellmanFunction(
        ## CutApproximation(Θᴳ, x′, obj_μ, belief_μ, deletion_minimum),
        CutApproximation(Θᴳ, x′, deletion_minimum),
        CutApproximation[],
        cut_type,
        Set{Vector{Float64}}(),
    )
end

################################################################################
# REFINE BELLMAN FUNCTION
################################################################################

"""
Refining the bellman function of a node by constructing a new cut
"""
# Could also be shifted to SDDP.jl, since it overwrites an existing function,
# but with additional arguments. Therefore, both methods can be distinguished.
function refine_bellman_function(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    bellman_function::BellmanFunction,
    risk_measure::SDDP.AbstractRiskMeasure,
    trial_points::Dict{Symbol,Float64},
    anchor_points::Dict{Symbol,Float64},
    bin_states::Dict{Symbol,BinaryState},
    dual_variables::Vector{Dict{Symbol,Float64}},
    noise_supports::Vector,
    nominal_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
) where {T}

    ############################################################################
    # CHECK IF ALL ELEMENTS CONTAINING INFORMATION ON ALL BACKWARD OPENINGS
    # HAVE THE SAME/REQUIRED LENGTH
    ############################################################################
    @assert (length(dual_variables) == length(noise_supports)
                                    == length(nominal_probability)
                                    == length(objective_realizations)
                                    )

    ############################################################################
    # RISK-RELATED PREPARATIOn
    ############################################################################
    # Preliminaries that are common to all cut types.
    risk_adjusted_probability = similar(nominal_probability)
    offset = SDDP.adjust_probability(
        risk_measure,
        risk_adjusted_probability,
        nominal_probability,
        noise_supports,
        objective_realizations,
        model.objective_sense == MOI.MIN_SENSE,
    )

    ############################################################################
    # ADD A NEW CUT TO THE CURRENT APPROXIMATION
    ############################################################################
    """ Note that so far only single cuts are supported by DynamicSDDiP.
    For multi cuts some more implementation is required.
    """

    if bellman_function.cut_type == SDDP.SINGLE_CUT
        return _add_average_cut(
            node,
            node_index,
            trial_points,
            anchor_points,
            bin_states,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables,
            offset,
            algo_params,
            model.ext[:iteration],
            applied_solvers
        )
    else  # Add a multi-cut
        @assert bellman_function.cut_type == SDDP.MULTI_CUT
        # TODO: Not implemented so far, see SDDP.jl
    end
end


"""
Adding one more cut to the bellman function (taking expectations in stochastic case).
"""
function _add_average_cut(
    node::SDDP.Node,
    node_index::Int64,
    trial_points::Dict{Symbol,Float64},
    anchor_points::Dict{Symbol,Float64},
    bin_states::Dict{Symbol,BinaryState},
    risk_adjusted_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    dual_variables::Vector{Dict{Symbol,Float64}},
    offset::Float64,
    algo_params::DynamicSDDiP.AlgoParams,
    iteration::Int64,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
)

    # Some initializations
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(dual_variables)

    ############################################################################
    # EXPECTED INTERCEPT AND DUAL VARIABLES
    ############################################################################
    # Calculate the expected intercept and dual variables with respect to the
    # risk-adjusted probability distribution.
    πᵏ = set_up_dict_for_duals(bin_states, trial_points, algo_params.state_approximation_regime)
    θᵏ = offset

    for i in 1:length(objective_realizations)
        p = risk_adjusted_probability[i]
        θᵏ += p * objective_realizations[i]
        for (key, dual) in dual_variables[i]
            πᵏ[key] += p * dual
        end
    end

    ############################################################################
    # GET CORRECT SIGMA
    ############################################################################
    # As cuts are created for the value function of the following state,
    # we need the parameters for this stage.
    if isa(algo_params.regularization_regime, DynamicSDDiP.NoRegularization)
        sigma = nothing
    else
        sigma = algo_params.regularization_regime.sigma[node_index+1]
    end

    ############################################################################
    # ADD THE CUT USING THE NEW EXPECTED COEFFICIENTS
    ############################################################################
    _add_cut(
        node,
        node.bellman_function.global_theta,
        θᵏ,
        πᵏ,
        bin_states,
        anchor_points,
        trial_points,
        # obj_y,
        # belief_y,
        sigma,
        iteration,
        algo_params.infiltrate_state,
        algo_params,
        applied_solvers,
        algo_params.state_approximation_regime,
    )

    return (theta = θᵏ, pi = πᵏ, λ = bin_states)
end


"""
Adding one more cut based on dual information if BinaryApproximation is used.
"""
# Add the cut to the model and the convex approximation.
function _add_cut(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    θᵏ::Float64, # epigraph variable theta
    πᵏ::Dict{Symbol,Float64}, # dual multipliers (cut coefficients)
    λᵏ::Dict{Symbol,BinaryState}, # binary states (anchor point in binary space for cut using BinaryApproximation)
    xᵏ_b::Dict{Symbol,Float64}, # anchor point for cut using BinaryApproximation
    xᵏ::Dict{Symbol,Float64}, # trial point (anchor point for cut without BinaryApproximation), outgoing_state
    # obj_y::Union{Nothing,NTuple{N,Float64}},
    # belief_y::Union{Nothing,Dict{T,Float64}},
    sigma::Union{Nothing,Float64},
    iteration::Int64,
    infiltrate_state::Symbol,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
) where {N,T}

    ############################################################################
    # CORRECT THE INTERCEPT (WE USE A DIFFERENT CUT FORMULA)
    ############################################################################
    for (key, λ) in λᵏ
        θᵏ -= πᵏ[key] * λᵏ[key].value
    end
    @infiltrate infiltrate_state in [:bellman, :all]

    if isnothing(sigma)
        sigma_use = nothing
    else
        sigma_use = copy(sigma)
    end

    ############################################################################
    # CONSTRUCT NONLINEAR CUT STRUCT
    ############################################################################
    cut = DynamicSDDiP.NonlinearCut(
            θᵏ,
            πᵏ,
            xᵏ,
            xᵏ_b,
            λᵏ,
            copy(state_approximation_regime.binary_precision),
            sigma_use,
            JuMP.VariableRef[],
            JuMP.ConstraintRef[],
            # obj_y,
            # belief_y,
            1,
            iteration
            )

    ############################################################################
    # ADD CUT PROJECTION TO SUBPROBLEM (we are already at the previous stage)
    ############################################################################
    _add_cut_constraints_to_models(node, V, cut, algo_params, infiltrate_state)

    ############################################################################
    # UPDATE CUT SELECTION
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "cut_selection" begin
        DynamicSDDiP._cut_selection_update(node, V, cut, xᵏ_b, xᵏ,
            applied_solvers, algo_params, infiltrate_state,
            algo_params.cut_selection_regime)
    end

    return
end


"""
Adding the new cut to the subproblem in case of a NonlinearCut.
"""
function _add_cut_constraints_to_models(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.NonlinearCut,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    model = JuMP.owner_model(V.theta)
    @assert model == node.subproblem

    number_of_states = length(node.states)

    """
    Later on, in the cut formula, for each component of the dual vector
    (in the binary space), we want to multiply the correct cut_coefficient
    with a newly introduced (relaxed) binary variable.

    We use the vector all_lambda to store all these [0,1] variables.
    Additionally, we make sure that at the same index of the vector
    all_coefficients the corresponding cut_coefficient is stored.

    In some cases to represent the cut projection closure, we need similar
    vector storages all_eta and all_mu.

    All other constraints and variables are introduced per state component.

    K_tilde counts how many components of the duals (in the binary space)
    have been considered so far.
    """
    all_coefficients = Float64[]
    all_lambda = JuMP.VariableRef[]
    all_eta = JuMP.VariableRef[]
    all_mu = JuMP.VariableRef[]
    K_tilde = 0

    ############################################################################
    # ADD NEW VARIABLES AND CONSTRAINTS
    ############################################################################
    for (state_index, (state_name, state_comp)) in enumerate(node.states)
        """
        If a state is binary already, it may be possible to abstain from
        introducing the cut projection / KKT constraints. That's actually
        what I did in the NCNBD implementation. For example, it seems
        unnecessary to introduce a new (relaxed) binary variable then.
        However, I am not absolutely sure if we in fact do not require
        the other constraints. Therefore, in this version, I also included
        them for binary state variables. At least, if StrongDuality is used,
        they are required anyway.
        """

        variable_info = node.ext[:state_info_storage][state_name].out

        if (!isfinite(variable_info.upper_bound) || !variable_info.has_ub) && !variable_info.binary
            error("When using DynamicSDDiP, state variables require an upper bound.")
        end

        ####################################################################
        # DETERMINE K (number of [0,1] variables required) AND BETA
        ####################################################################
        if variable_info.binary
            beta = 1.0
            K = 1
        elseif variable_info.integer
            beta = 1.0
            K = SDDP._bitsrequired(variable_info.upper_bound)
        else
            beta = cut.binary_precision[state_name]
            K = SDDP._bitsrequired(round(Int, variable_info.upper_bound / beta))
        end

        ####################################################################
        # DEFINE CUT PROJECTION CLOSURE
        ####################################################################
        # For the current state component add all required variables
        # and constraints to represent the cut projection closure
        K_tilde = represent_cut_projection_closure!(
                    model,
                    node,
                    ########################################################
                    V.states[state_name], # state_comp
                    state_name,
                    state_index,
                    ########################################################
                    cut.coefficients,
                    cut.binary_state,
                    cut.sigma,
                    cut.cut_variables,
                    cut.cut_constraints,
                    cut.iteration,
                    beta,
                    ########################################################
                    all_coefficients,
                    all_lambda,
                    all_eta,
                    all_mu,
                    ########################################################
                    K,
                    K_tilde,
                    ########################################################
                    infiltrate_state,
                    ########################################################
                    algo_params.state_approximation_regime.cut_projection_regime
                    )

    end

    @infiltrate infiltrate_state in [:bellman, :all]

    ############################################################################
    # MAKE SOME VALIDITY CHECKS
    ############################################################################
    validity_checks!(cut, V, K_tilde, all_lambda, all_mu, all_eta,
        all_coefficients, number_of_states,
        algo_params.state_approximation_regime.cut_projection_regime)

    number_of_duals = size(all_lambda, 1)

    ############################################################################
    # ADD THE ORIGINAL CUT CONSTRAINT AS WELL
    ############################################################################
    expr = get_cut_expression(model, node, V, all_lambda, all_mu, all_eta,
        all_coefficients, number_of_states, number_of_duals,
        algo_params.state_approximation_regime.cut_projection_regime)

    constraint_ref = if JuMP.objective_sense(model) == MOI.MIN_SENSE
        JuMP.@constraint(model, expr >= cut.intercept)
    else
        JuMP.@constraint(model, expr <= cut.intercept)
    end
    push!(cut.cut_constraints, constraint_ref)

    ############################################################################
    # ADD SOS1 STRONG DUALITY CONSTRAINT
    ############################################################################
    #add_strong_duality_cut!(model, node, cut, V, all_lambda, all_mu, all_eta,
    #    all_coefficients, number_of_states, number_of_duals,
    #    algo_params.state_approximation_regime.cut_projection_regime)

    return

end


"""
Defining the constraints and variables corresponding to representing
the cut projection closure.
"""
function represent_cut_projection_closure!(
    model::JuMP.Model,
    node::SDDP.Node,
    ########################################################
    state_comp::JuMP.VariableRef,
    state_name::Symbol,
    state_index::Int64,
    ########################################################
    coefficients::Dict{Symbol,Float64},
    binary_state::Dict{Symbol,BinaryState},
    sigma::Union{Nothing,Float64},
    cut_variables::Vector{JuMP.VariableRef},
    cut_constraints::Vector{JuMP.ConstraintRef},
    iteration::Int64,
    beta::Float64,
    ########################################################
    all_coefficients::Vector{Float64},
    all_lambda::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    ########################################################
    K::Int64,
    K_tilde::Int64,
    ########################################################
    infiltrate_state::Symbol,
    ########################################################
    cut_projection_regime::Union{DynamicSDDiP.SOS1, DynamicSDDiP.BigM, DynamicSDDiP.KKT}
    )

    ############################################################################
    # STORE THE RELATED CUT COEFFICIENTS π
    ############################################################################
    related_coefficients = Vector{Float64}(undef, K)

    for (i, (name, value)) in enumerate(coefficients)
        if binary_state[name].x_name == state_name
            index = binary_state[name].k
            related_coefficients[index] = coefficients[name]
        end
    end
    append!(all_coefficients, related_coefficients)

    ############################################################################
    # ADD REQUIRED VARIABLES (CPC-CONSTRAINTS 4, 5)
    ############################################################################
    """
    Note that we use gamma instead of lambda here to avoid misconstructions
    compared to the binary approximation in the backward pass.
    """
    γ = JuMP.@variable(model, [k in 1:K], lower_bound=0, upper_bound=1, base_name = "γ_" * string(state_index) * "_it" * string(iteration))
    ν = JuMP.@variable(model, [k in 1:K], lower_bound=0, base_name = "ν_" * string(state_index) * "_it" * string(iteration))
    μ = JuMP.@variable(model, [k in 1:K], lower_bound=0, base_name = "μ_" * string(state_index) * "_it" * string(iteration))
    η = JuMP.@variable(model, base_name = "η_" * string(state_index) * "_it" * string(iteration))

    # Store those variables in the storing vectors
    append!(all_lambda, γ)
    append!(all_mu, μ)
    push!(all_eta, η)

    # Store those variables as cut_variables for the existing cut
    append!(cut_variables, γ)
    append!(cut_variables, ν)
    append!(cut_variables, μ)
    push!(cut_variables, η)

    ############################################################################
    # ADD BINARY EXPANSION CONSTRAINT (CPC-CONSTRAINT 3)
    ############################################################################
    binary_constraint = JuMP.@constraint(
        model,
        state_comp == SDDP.bincontract([γ[k] for k in 1:K], beta)
    )
    push!(cut_constraints, binary_constraint)

    ############################################################################
    # ADD DUAL FEASIBILITY CONSTRAINTS (CPC-CONSTRAINT 2)
    ############################################################################
    primal_feas_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        -related_coefficients[k] - ν[k] + μ[k] + 2^(k-1) * beta * η == 0
    )
    append!(cut_constraints, primal_feas_constraints)

    ############################################################################
    # ADD COMPLEMENTARITY CONSTRAINTS
    ############################################################################
    add_complementarity_constraints!(
        model,
        node,
        ########################################################################
        state_index,
        ########################################################################
        cut_variables,
        cut_constraints,
        sigma,
        iteration,
        beta,
        ########################################################################
        related_coefficients,
        ########################################################################
        γ,
        μ,
        ν,
        ########################################################################
        K,
        ########################################################################
        infiltrate_state,
        ########################################################################
        cut_projection_regime
    )

    ############################################################################
    # INCREASE K_tilde
    ############################################################################
    K_tilde += K

    return K_tilde

end

"""
Function representing complementarity constraints using projection_regime KKT.
"""
function add_complementarity_constraints!(
    model::JuMP.Model,
    node::SDDP.Node,
    ########################################################################
    state_index::Int64,
    ########################################################################
    cut_variables::Vector{JuMP.VariableRef},
    cut_constraints::Vector{JuMP.ConstraintRef},
    sigma::Union{Nothing,Float64},
    iteration::Int64,
    beta::Float64,
    ########################################################################
    related_coefficients::Vector{Float64},
    ########################################################################
    γ::Vector{JuMP.VariableRef},
    μ::Vector{JuMP.VariableRef},
    ν::Vector{JuMP.VariableRef},
    ########################################################################
    K::Int64,
    ########################################################################
    infiltrate_state::Symbol,
    ########################################################################
    cut_projection_regime::DynamicSDDiP.KKT,
)

    @infiltrate infiltrate_state in [:bellman, :all]

    ############################################################################
    # ADD COMPLEMENTARITY CONSTRAINTS
    ############################################################################
    complementarity_constraint_1 = JuMP.@constraint(
        model,
        [k=1:K],
        ν[k] * γ[k] == 0
    )
    append!(cut_constraints, complementarity_constraint_1)

    complementarity_constraint_2 = JuMP.@constraint(
        model,
        [k=1:K],
        μ[k] * (γ[k]-1) == 0
    )
    append!(cut_constraints, complementarity_constraint_2)

    return

end


"""
Function representing complementarity constraints using projection_regime BigM.
"""
function add_complementarity_constraints!(
    model::JuMP.Model,
    node::SDDP.Node,
    ########################################################################
    state_index::Int64,
    ########################################################################
    cut_variables::Vector{JuMP.VariableRef},
    cut_constraints::Vector{JuMP.ConstraintRef},
    sigma::Union{Nothing,Float64},
    iteration::Int64,
    beta::Float64,
    ########################################################################
    related_coefficients::Vector{Float64},
    ########################################################################
    γ::Vector{JuMP.VariableRef},
    μ::Vector{JuMP.VariableRef},
    ν::Vector{JuMP.VariableRef},
    ########################################################################
    K::Int64,
    ########################################################################
    infiltrate_state::Symbol,
    ########################################################################
    cut_projection_regime::DynamicSDDiP.BigM,
)

    @infiltrate infiltrate_state in [:bellman, :all]

    ############################################################################
    # ADD ADDITIONAL BINARY VARIABLES
    ############################################################################
    w = JuMP.@variable(model, [k in 1:K], binary=true, base_name = "w_" * string(state_index) * "_it" * string(iteration))
    u = JuMP.@variable(model, [k in 1:K], binary=true, base_name = "u_" * string(state_index) * "_it" * string(iteration))
    append!(cut_variables, w)
    append!(cut_variables, u)

    ############################################################################
    # DETERMINE BIG-M PARAMETER
    ############################################################################
    # TODO: Could be further improved later on
    bigM = get_bigM(node, sigma, beta, related_coefficients, K)

    ############################################################################
    # ADD BIG-M CONSTRAINTS
    ############################################################################
    bigM_11_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        γ[k] <= w[k]
    )
    append!(cut_constraints, bigM_11_constraints)

    bigM_12_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        ν[k] <= bigM * (1-w[k])
    )
    append!(cut_constraints, bigM_12_constraints)

    bigM_21_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        1 - γ[k] <= u[k]
    )
    append!(cut_constraints, bigM_21_constraints)

    bigM_22_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        μ[k] <= bigM * (1-u[k])
    )
    append!(cut_constraints, bigM_22_constraints)

    return

end


"""
Determine a reasonable bigM value based on the maximum upper bound of all state
components, sigma and beta. This could be improved later.
"""
function get_bigM(node::SDDP.Node, sigma::Union{Nothing,Float64}, beta::Float64, related_coefficients::Vector{Float64}, K::Int64)

    ############################################################################
    # DETERMINE U_MAX
    ############################################################################
    U_max = 0
    for (i, (name, state)) in enumerate(node.states)
        variable_info = node.ext[:state_info_storage][name].out

        upper_bound = -Inf
        if variable_info.binary
            upper_bound = 1
        elseif variable_info.has_ub
            upper_bound = variable_info.upper_bound
        end

        if upper_bound > U_max
            U_max = upper_bound
        end

    end

    ############################################################################
    # DETERMINE BIG-M
    ############################################################################
    bigM = 0
    if isnothing(sigma)
        # no regularization is used, so bigM is just bounded by an arbitrary value
        bigM = 1e4
    else
        bigM = 2 * sigma * U_max
    end

    return bigM
end


"""
Function representing complementarity constraints using projection_regime KKT.
"""
function add_complementarity_constraints!(
    model::JuMP.Model,
    node::SDDP.Node,
    ########################################################################
    state_index::Int64,
    ########################################################################
    cut_variables::Vector{JuMP.VariableRef},
    cut_constraints::Vector{JuMP.ConstraintRef},
    sigma::Union{Nothing,Float64},
    iteration::Int64,
    beta::Float64,
    ########################################################################
    related_coefficients::Vector{Float64},
    ########################################################################
    γ::Vector{JuMP.VariableRef},
    μ::Vector{JuMP.VariableRef},
    ν::Vector{JuMP.VariableRef},
    ########################################################################
    K::Int64,
    ########################################################################
    infiltrate_state::Symbol,
    ########################################################################
    cut_projection_regime::DynamicSDDiP.SOS1,
)

    @infiltrate infiltrate_state in [:bellman, :all]

    ############################################################################
    # AUXILIARY VARIABLE
    ############################################################################
    # First represent γ[k]-1 as a new variable
    ρ = JuMP.@variable(model, [k in 1:K], lower_bound=0, upper_bound=1, base_name = "ρ_" * string(state_index) * "_it" * string(iteration))
    append!(cut_variables, ρ)

    ρ_constraint = JuMP.@constraint(model, [k in 1:K], ρ[k] == 1 - γ[k])
    append!(cut_constraints, ρ_constraint)

    ############################################################################
    # ADD SOS1 CONSTRAINTS
    ############################################################################
    for k in 1:K
        SOS1_constraint_1 = JuMP.@constraint(model, [γ[k], ν[k]] in JuMP.SOS1())
        SOS1_constraint_2 = JuMP.@constraint(model, [μ[k], ρ[k]] in JuMP.SOS1())
        push!(cut_constraints, SOS1_constraint_1)
        push!(cut_constraints, SOS1_constraint_2)
    end

    return

end


"""
Defining the constraints and variables corresponding to representing
the cut projection closure if StrongDuality is exploited.
"""
function represent_cut_projection_closure!(
    model::JuMP.Model,
    node::SDDP.Node,
    ########################################################
    state_comp::JuMP.VariableRef,
    state_name::Symbol,
    state_index::Int64,
    ########################################################
    coefficients::Dict{Symbol,Float64},
    binary_state::Dict{Symbol,BinaryState},
    sigma::Union{Nothing,Float64},
    cut_variables::Vector{JuMP.VariableRef},
    cut_constraints::Vector{JuMP.ConstraintRef},
    iteration::Int64,
    beta::Float64,
    ########################################################
    all_coefficients::Vector{Float64},
    all_lambda::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    ########################################################
    K::Int64,
    K_tilde::Int64,
    ########################################################
    infiltrate_state::Symbol,
    ########################################################
    cut_projection_regime::DynamicSDDiP.StrongDuality
    )

    ############################################################################
    # STORE THE RELATED CUT COEFFICIENTS π
    ############################################################################
    related_coefficients = Vector{Float64}(undef, K)

    for (i, (name, value)) in enumerate(coefficients)
        if binary_state[name].x_name == state_name
            index = binary_state[name].k
            related_coefficients[index] = coefficients[name]
        end
    end
    append!(all_coefficients, related_coefficients)

    ############################################################################
    # ADD REQUIRED VARIABLES (CPC-CONSTRAINTS 4, 5)
    ############################################################################
    μ = JuMP.@variable(model, [k in 1:K], lower_bound=0, base_name = "μ_" * string(state_index) * "_it" * string(iteration))
    η = JuMP.@variable(model, base_name = "η_" * string(state_index) * "_it" * string(iteration))

    # Store those variables in the storing vectors
    append!(all_mu, μ)
    push!(all_eta, η)

    # Store those variables as cut_variables for the existing cut
    append!(cut_variables, μ)
    push!(cut_variables, η)

    ############################################################################
    # ADD DUAL FEASIBILITY CONSTRAINTS (CPC-CONSTRAINT 2d)
    ############################################################################
    primal_feas_constraints = JuMP.@constraint(
        model,
        [k=1:K],
        μ[k] + 2^(k-1) * beta * η >= related_coefficients[k]
    )
    append!(cut_constraints, primal_feas_constraints)

    ############################################################################
    # INCREASE K_tilde
    ############################################################################
    K_tilde += K

    return K_tilde

end


################################################################################
# AUXILIARY FUNCTIONS
################################################################################

function set_up_dict_for_duals(
    bin_states::Dict{Symbol,BinaryState},
    trial_points::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    return Dict(key => 0.0 for key in keys(bin_states))
end

function set_up_dict_for_duals(
    bin_states::Dict{Symbol,BinaryState},
    trial_points::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    return Dict(key => 0.0 for key in keys(trial_points))
end

function validity_checks!(
    cut::DynamicSDDiP.NonlinearCut,
    V::DynamicSDDiP.CutApproximation,
    K_tilde::Int64,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    cut_projection_regime::Union{DynamicSDDiP.SOS1,DynamicSDDiP.BigM,DynamicSDDiP.KKT},
    )

    @assert (K_tilde == size(collect(values(cut.coefficients)), 1)
                    == size(all_coefficients, 1)
                    == size(all_mu, 1)
                    == size(all_lambda, 1)
                    )
    @assert (number_of_states == size(all_eta, 1)
                              == length(V.states)
                              )

    return
end

function validity_checks!(
    cut::DynamicSDDiP.NonlinearCut,
    V::DynamicSDDiP.CutApproximation,
    K_tilde::Int64,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    cut_projection_regime::DynamicSDDiP.StrongDuality,
    )

    @assert (K_tilde == size(collect(values(cut.coefficients)), 1)
                    == size(all_coefficients, 1)
                    == size(all_mu, 1)
                    )

    @assert (number_of_states == size(all_eta, 1)
                              == length(V.states)
                              )

    return
end

function get_cut_expression(
    model::JuMP.Model,
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    number_of_duals::Int64,
    cut_projection_regime::Union{DynamicSDDiP.SOS1,DynamicSDDiP.BigM,DynamicSDDiP.KKT},
    )

    expr = JuMP.@expression(
        model,
        V.theta - sum(all_coefficients[j] * all_lambda[j]  for j in 1:number_of_duals)
    )

    return expr
end

function get_cut_expression(
    model::JuMP.Model,
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    number_of_duals::Int64,
    cut_projection_regime::DynamicSDDiP.StrongDuality,
    )

    expr = JuMP.@expression(
        model,
        V.theta - sum(all_mu[j]  for j in 1:size(all_mu, 1))
        - sum(x * all_eta[i]  for (i, (_,x)) in enumerate(V.states))
    )

    return expr
end

function add_strong_duality_cut!(
    model::JuMP.Model,
    node::SDDP.Node,
    cut::DynamicSDDiP.NonlinearCut,
    V::DynamicSDDiP.CutApproximation,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    number_of_duals::Int64,
    cut_projection_regime::DynamicSDDiP.SOS1,
    )

    strong_duality_expr = JuMP.@expression(
        model,
        sum(all_coefficients[j] * all_lambda[j]  for j in 1:number_of_duals)
        - sum(all_mu[j]  for j in 1:number_of_duals)
        - sum(node.ext[:state_info_storage][sym].out.lower_bound * all_eta[i]  for (i, (sym, x)) in enumerate(V.states))
    )

    constraint_ref = if JuMP.objective_sense(model) == MOI.MIN_SENSE
        JuMP.@constraint(model, strong_duality_expr >= 0)
    else
        JuMP.@constraint(model, strong_duality_expr <= 0)
    end
    push!(cut.cut_constraints, constraint_ref)

    return
end

function add_strong_duality_cut!(
    model::JuMP.Model,
    node::SDDP.Node,
    cut::DynamicSDDiP.NonlinearCut,
    V::DynamicSDDiP.CutApproximation,
    all_lambda::Vector{JuMP.VariableRef},
    all_mu::Vector{JuMP.VariableRef},
    all_eta::Vector{JuMP.VariableRef},
    all_coefficients::Vector{Float64},
    number_of_states::Int64,
    number_of_duals::Int64,
    cut_projection_regime::Union{DynamicSDDiP.BigM,DynamicSDDiP.KKT,DynamicSDDiP.StrongDuality},
    )

    return
end

################################################################################

"""
Adding one more cut based on dual information if NoStateApproximation is used.
"""
# Add the cut to the model and the convex approximation.
function _add_cut(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    θᵏ::Float64, # epigraph variable theta
    πᵏ::Dict{Symbol,Float64}, # dual multipliers (cut coefficients)
    λᵏ::Dict{Symbol,BinaryState}, # binary states (anchor point in binary space for cut using BinaryApproximation)
    xᵏ_b::Dict{Symbol,Float64}, # anchor point for cut using BinaryApproximation
    xᵏ::Dict{Symbol,Float64}, # trial point (anchor point for cut without BinaryApproximation), outgoing_state
    # obj_y::Union{Nothing,NTuple{N,Float64}},
    # belief_y::Union{Nothing,Dict{T,Float64}},
    sigma::Union{Nothing,Float64},
    iteration::Int64,
    infiltrate_state::Symbol,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
) where {N,T}

    ############################################################################
    # CORRECT THE INTERCEPT (WE USE A DIFFERENT CUT FORMULA)
    ############################################################################
    for (key, x) in xᵏ
        θᵏ -= πᵏ[key] * x
    end
    @infiltrate infiltrate_state in [:bellman, :all]

    if isnothing(sigma)
        sigma_use = nothing
    else
        sigma_use = copy(sigma)
    end

    ############################################################################
    # CONSTRUCT LINEAR CUT STRUCT
    ############################################################################
    cut = DynamicSDDiP.LinearCut(
            θᵏ,
            πᵏ,
            xᵏ,
            sigma_use,
            nothing,
            # obj_y,
            # belief_y,
            1,
            iteration
            )

    ############################################################################
    # ADD CUT TO SUBPROBLEM (we are already at the previous stage)
    ############################################################################
    _add_cut_constraints_to_models(node, V, cut, algo_params, infiltrate_state)

    ############################################################################
    # UPDATE CUT SELECTION
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "cut_selection" begin
        DynamicSDDiP._cut_selection_update(node, V, cut, xᵏ_b, xᵏ,
            applied_solvers, algo_params, infiltrate_state,
            algo_params.cut_selection_regime)
    end

    return
end


"""
Adding the new cut to the subproblem in case of a LinearCut.
"""
function _add_cut_constraints_to_models(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.LinearCut,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    )

    ############################################################################
    # SOME INITIALIZATIONS
    ############################################################################
    model = JuMP.owner_model(V.theta)
    @assert model == node.subproblem

    @infiltrate infiltrate_state in [:bellman, :all]

    ############################################################################
    # ADD THE LINEAR CUT CONSTRAINT
    ############################################################################
    expr = JuMP.@expression(
        model,
        V.theta - sum(cut.coefficients[i] * x for (i, x) in V.states)
    )

    cut.cut_constraint = if JuMP.objective_sense(model) == MOI.MIN_SENSE
        JuMP.@constraint(model, expr >= cut.intercept)
    else
        JuMP.@constraint(model, expr <= cut.intercept)
    end

    return

end
