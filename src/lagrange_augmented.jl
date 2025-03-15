"""
Considering an augmented Lagrangian relaxation problem, i.e. the inner problem of the
Lagrangian dual using the 1-norm
"""

function _augmented_Lagrangian_relaxation!(
    node::SDDP.Node,
    node_index::Int64,
    π_k::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    regularization_regime::DynamicSDDiP.NoRegularization,
    update_subgradients::Bool = true,
)
    L_k = _solve_Lagrangian_relaxation!(node, π_k, h_expr, h_k, true)

    return L_k
end

function _augmented_Lagrangian_relaxation!(
    node::SDDP.Node,
    node_index::Int64,
    π_k::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    regularization_regime::DynamicSDDiP.Regularization,
    update_subgradients::Bool = true,
)
    rho = regularization_regime.sigma[node_index]

    L_k = _solve_augmented_Lagrangian_relaxation!(node, π_k, h_expr, h_k, rho, true)

    return L_k
end

"""
Solving the augmented Lagrangian relaxation problem, i.e. the inner problem of the
Lagrangian dual using the 1-norm
"""

function _solve_augmented_Lagrangian_relaxation!(
    node::SDDP.Node,
    π_k::Vector{Float64},
    h_expr::Vector{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}},
    h_k::Vector{Float64},
    rho::Float64,
    update_subgradients::Bool = true,
)
    model = node.subproblem

    # Set the Lagrangian relaxation of the objective in the primal model
    old_obj = JuMP.objective_function(model)

    # Get number of states
    number_of_states = length(h_k)

    # Add new variable
    omega = JuMP.@variable(model, [1:number_of_states])

    # Adapt objective function
    #Infiltrator.@infiltrate
    JuMP.set_objective_function(model, JuMP.@expression(model, old_obj - π_k' * h_expr ))
    JuMP.set_objective_function(model, JuMP.@expression(model, old_obj - π_k' * h_expr + rho * sum(omega[i] for i in 1:number_of_states)))

    # Add norm constraints
    norm_1 = JuMP.@constraint(model, [i=1:number_of_states], -omega[i] <= h_expr[i])
    norm_2 = JuMP.@constraint(model, [i=1:number_of_states], omega[i] >= h_expr[i])

    # Optimization
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == MOI.OPTIMAL

    # Update the correct values
    L_k = JuMP.objective_value(model)
    if update_subgradients
        h_k .= -JuMP.value.(h_expr)
    end

    # Reset old objective
    JuMP.set_objective_function(model, old_obj)

    # Delete augmentation constraints and variables
    JuMP.delete(model, norm_1)
    JuMP.delete(model, norm_2)
    JuMP.delete(model, omega)

    return L_k
end
