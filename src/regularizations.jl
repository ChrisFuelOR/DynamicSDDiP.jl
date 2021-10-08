# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################


"""
Modifying the forward pass problem to include a regularization term
if regularization is used.
"""
function regularize_subproblem!(node::SDDP.Node, node_index::Int64,
    subproblem::JuMP.Model, regularization_regime::DynamicSDDiP.Regularization)

    # Note that the copy constraint is not modeled explicitly here. Instead,
    # the state variable is unfixed and takes the role of z in our paper.
    # It is then subtracted from the fixed value to obtain the so called slack.

    reg_data = node.ext[:regularization_data]
    reg_data[:fixed_state_value] = Dict{Symbol,Float64}()
    reg_data[:slacks] = Any[]
    reg_data[:reg_variables] = JuMP.VariableRef[]
    reg_data[:reg_constraints] = JuMP.ConstraintRef[]

    number_of_states = 0

    # UNFIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        reg_data[:fixed_state_value][name] = JuMP.fix_value(state_comp.in)
        push!(reg_data[:slacks], reg_data[:fixed_state_value][name] - state_comp.in)
        JuMP.unfix(state_comp.in)
        follow_state_unfixing!(state_comp)
        number_of_states = i
    end

    # STORE ORIGINAL OBJECTIVE FUNCTION
    ############################################################################
    old_obj = reg_data[:old_objective] = JuMP.objective_function(subproblem)

    # DEFINE NEW VARIABLES, CONSTRAINTS AND OBJECTIVE
    ############################################################################
    # These variables and constraints are used to define the norm of the slack as a MILP
    # Using the lifting approach without binary requirements
    slack = reg_data[:slacks]

    # Variable for objective
    v = JuMP.@variable(subproblem, base_name = "reg_v")
    push!(reg_data[:reg_variables], v)

    # Get sign for regularization term
    fact = (JuMP.objective_sense(subproblem) == JuMP.MOI.MIN_SENSE ? 1 : -1)

    # New objective
    new_obj = old_obj + fact * regularization_regime.sigma[node_index] * v
    JuMP.set_objective_function(subproblem, new_obj)

    # Variables
    alpha = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "alpha")
    append!(reg_data[:reg_variables], alpha)

    # Constraints
    const_plus = JuMP.@constraint(subproblem, [i=1:number_of_states], -alpha[i] <= slack[i])
    const_minus = JuMP.@constraint(subproblem, [i=1:number_of_states], slack[i] <= alpha[i])
    append!(reg_data[:reg_constraints], const_plus)
    append!(reg_data[:reg_constraints], const_minus)

    const_norm = JuMP.@constraint(subproblem, v >= sum(alpha[i] for i in 1:number_of_states))
    push!(reg_data[:reg_constraints], const_norm)

end

"""
Trivial modification of the forward pass problem if no regularization is used.
"""
function regularize_subproblem!(node::SDDP.Node, node_index::Int64,
    subproblem::JuMP.Model, regularization_regime::DynamicSDDiP.NoRegularization)
    return
end


"""
Modifying the forward pass problem to remove the regularization term
and regain the original model if regularization is used.
"""
function deregularize_subproblem!(node::SDDP.Node, subproblem::JuMP.Model, regularization_regime::DynamicSDDiP.Regularization)

    reg_data = node.ext[:regularization_data]

    # FIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(node.states)
        prepare_state_fixing!(node, state_comp)
        JuMP.fix(state_comp.in, reg_data[:fixed_state_value][name], force=true)
    end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, reg_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    delete(subproblem, reg_data[:reg_variables])

    for constraint in reg_data[:reg_constraints]
        delete(subproblem, constraint)
    end

    delete!(node.ext, :regularization_data)

end

"""
Trivial modification of the forward pass problem if no regularization is used.
"""
function deregularize_subproblem!(node::SDDP.Node, subproblem::JuMP.Model, regularization_regime::DynamicSDDiP.NoRegularization)
    return
end


"""
Introducing a regularizing term to the backward pass problem in binary space.
"""
function regularize_backward!(node::SDDP.Node, subproblem::JuMP.Model, sigma::Float64)

    bw_data = node.ext[:backward_data]
    binary_states = bw_data[:bin_states]

    number_of_states = size(collect(values(binary_states)), 1)

    reg_data = node.ext[:regularization_data]
    reg_data[:fixed_state_value] = Dict{Symbol,Float64}()
    reg_data[:slacks] = Any[]
    reg_data[:reg_variables] = JuMP.VariableRef[]
    reg_data[:reg_constraints] = JuMP.ConstraintRef[]

    # DETERMINE SIGMA TO BE USED IN BINARY SPACE
    ############################################################################
    Umax = 0
    for (i, (name, state_comp)) in enumerate(node.states)
        if state_comp.info.out.upper_bound > Umax
            Umax = state_comp.info.out.upper_bound
        end
    end
    # Here, not sigma, but a different regularization parameter is used
    sigma_bin = sigma * Umax

    # UNFIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(binary_states)
        reg_data[:fixed_state_value][name] = JuMP.fix_value(state_comp)
        push!(reg_data[:slacks], reg_data[:fixed_state_value][name] - state_comp)
        JuMP.unfix(state_comp)
        follow_state_unfixing_binary!(state_comp)
    end

    # STORE ORIGINAL OBJECTIVE FUNCTION
    ############################################################################
    old_obj = reg_data[:old_objective] = JuMP.objective_function(subproblem)

    # DEFINE NEW VARIABLES, CONSTRAINTS AND OBJECTIVE
    ############################################################################
    # These variables and constraints are used to define the norm of the slack as a MILP
    # Using the lifting approach without binary requirements
    slack = reg_data[:slacks]

    # Variable for objective
    v = JuMP.@variable(subproblem, base_name = "reg_v")
    push!(reg_data[:reg_variables], v)

    # Get sign for regularization term
    fact = (JuMP.objective_sense(subproblem) == JuMP.MOI.MIN_SENSE ? 1 : -1)

    # New objective
    new_obj = old_obj + fact * sigma_bin * v
    JuMP.set_objective_function(subproblem, new_obj)

    # Variables
    alpha = JuMP.@variable(subproblem, [i=1:number_of_states], base_name = "alpha")
    append!(reg_data[:reg_variables], alpha)

    # Constraints
    const_plus = JuMP.@constraint(subproblem, [i=1:number_of_states], -alpha[i] <= slack[i])
    const_minus = JuMP.@constraint(subproblem, [i=1:number_of_states], slack[i] <= alpha[i])
    append!(reg_data[:reg_constraints], const_plus)
    append!(reg_data[:reg_constraints], const_minus)

    const_norm = JuMP.@constraint(subproblem, v >= sum(alpha[i] for i in 1:number_of_states))
    push!(reg_data[:reg_constraints], const_norm)

end


"""
Regaining the unregularized problem in binary space.
"""
function deregularize_backward!(node::SDDP.Node, subproblem::JuMP.Model)

    reg_data = node.ext[:regularization_data]
    bw_data = node.ext[:backward_data]

    # FIX THE STATE VARIABLES
    ############################################################################
    for (i, (name, state_comp)) in enumerate(bw_data[:bin_states])
        prepare_state_fixing_binary!(node, state_comp)
        JuMP.fix(state_comp, reg_data[:fixed_state_value][name], force=true)
    end

    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    JuMP.set_objective_function(subproblem, reg_data[:old_objective])

    # DELETE ALL REGULARIZATION-BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    delete(subproblem, reg_data[:reg_variables])

    for constraint in reg_data[:reg_constraints]
        delete(subproblem, constraint)
    end

    delete!(node.ext, :regularization_data)

end
