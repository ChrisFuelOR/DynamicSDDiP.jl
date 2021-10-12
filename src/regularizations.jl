# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################


"""
Modifying the forward pass problem to include a regularization term.
"""
function regularize_subproblem!(node::SDDP.Node, subproblem::JuMP.Model, sigma::Float64)

    #NOTE: The copy constraint is not modeled explicitly here. Instead,
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
    new_obj = old_obj + fact * sigma * v
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
Modifying the forward pass problem to remove the regularization term
and regain the original model.
"""
function deregularize_subproblem!(node::SDDP.Node, subproblem::JuMP.Model)

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


"""
Check if an increase of the number of binary variables is required.
"""
function binary_refinement_check(
    model::SDDP.PolicyGraph{T},
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    sampled_states::Vector{Dict{Symbol,Float64}},
    solution_check::Bool,
    ) where {T}

    # Check if feasible solution has changed since last iteration
    # If not, then the cut was not tight enough, so binary approximation should be refined
    for i in 1:size(previous_solution,1)
        for (name, state_comp) in model.nodes[i].states
            current_sol = sampled_states[i][name]
            previous_sol = previous_solution[i][name]
            if ! isapprox(current_sol, previous_sol)
                solution_check = false
            end
        end
    end

    return solution_check

end


"""
Executing a binary refinement: Increasing the number of binary variables to
approximate the states
"""
function binary_refinement(
    model::SDDP.PolicyGraph{T},
    algo_params::NCNBD.AlgoParams,
    binary_refinement::Symbol,
    ) where {T}

    all_refined = Int[]

    # Consider stage 2 here (should be the same for all following stages)
    # Precision is only used (and increased) for continuous variables
    for (name, state_comp) in model.nodes[2].states
        if !state_comp.info.in.binary && !state_comp.info.in.integer
            current_prec = algo_params.binary_precision[name]
            ub = state_comp.info.in.upper_bound
            K = SDDP._bitsrequired(round(Int, ub / current_prec))
            new_prec = ub / sum(2^(k-1) for k in 1:K+1)

            # Only apply refinement if int64 is appropriate to represent this number
            if ub / new_prec > 2^63 - 1
                push!(all_refined, 0)
            else
                push!(all_refined, 1)
                algo_params.binary_precision[name] = new_prec
            end

        else
            push!(all_refined, 2)
        end
    end

    # Check refinement status
    if 0 in all_refined
        if 1 in all_refined
            binary_refinement = :partial
        else
            binary_refinement = :impossible
        end
    else
        binary_refinement = :all
    end

    return binary_refinement

end
