# The functions
# > "set_incoming_state!",
# > "setup_state",
# > "get_outgoing_state",
# and structs
# > "State",
# > "StateInfo"
# are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
This struct is used to store the bound and integer information of a state
variable (either .in or .out).
This is required to be able to retrieve this information after fixing/unfixing
the states for regularization or relaxation purposes.

In NCNBD.jl this was done by introducing a struct NCNBD.State with a field
stateInfo, however only for the linearized model. Overwriting the SDDP.State in
general without rewriting a lot of SDDP code does not seem possible.
"""
mutable struct VariableInfo
    has_lb::Bool
    lower_bound::Float64
    has_ub::Bool
    upper_bound::Float64
    binary::Bool
    integer::Bool

    function VariableInfo(
        has_lb = false,
        lower_bound = NaN,
        has_ub = false,
        upper_bound = NaN,
        binary = false,
        integer = false,
    )
        return new(
            has_lb, lower_bound, has_ub, upper_bound, binary, integer
        )
    end
end

"""
This struct is used to store the bound and integer information of a state
variable (both .in or .out) by using VariableInfo structs.
This is required to be able to retrieve this information after fixing/unfixing
the states for regularization or relaxation purposes.
"""
mutable struct StateInfoStorage
    in::DynamicSDDiP.VariableInfo
    out::DynamicSDDiP.VariableInfo
end

"""
Function to set up the storage in VariableInfo.
"""
function get_variable_info(state::JuMP.VariableRef)

    variable_info = DynamicSDDiP.VariableInfo()

    if JuMP.has_lower_bound(state)
        variable_info.has_lb = true
        variable_info.lower_bound = JuMP.lower_bound(state)
    end
    if JuMP.has_upper_bound(state)
        variable_info.has_ub = true
        variable_info.upper_bound = JuMP.upper_bound(state)
    end
    if JuMP.is_fixed(state)
        variable_info.fixed = true
        variable_info.fixed_value = JuMP.fix_value(state)
    end
    if JuMP.is_binary(state)
        variable_info.binary = true
    end
    if JuMP.is_integer(state)
        variable_info.integer = true
    end

    return variable_info
end


################################################################################

"""
Set the incoming state variable of the node to the values contained in the
state dict. This basically means that the variable is fixed.

This requires to delete the binary or integer type (with function
prepare_state_fixing!()), cause otherwise problems occur when fixing the
variables. The bounds are not removed explicitly, though, since this can
be automatically done via force=true in the fix command.
"""
function set_incoming_state!(node::SDDP.Node, state::Dict{Symbol,Float64})
    for (state_name, value) in state
        prepare_state_fixing!(node, state_name)
        JuMP.fix(node.states[state_name].in, value, force=true)
    end
    return
end

"""
Preparation of fixing a state variable with three different types of argument.

Note that this differentiation is required since we use the term "state" for
different things. One time for the dict of SDDP.States in our model, each
basically containing two variable references (.in and .out).
Another time for a current state (incoming_state, outgoing_state), which is
a specific allocation of values to one of these variables (e.g. representing
the current trial point). This is a Dict{Symbol,Float64} then. The third case
is used for binary approximations of the state space.

The first method with the symbol is used for setting the incoming state, since
there we loop over all states getting the state_name (::Symbol) and the
value (::Float64) from the current state dict (which is the incoming_state_value
dict, actually).

The second method with the SDDP.State is used for resetting the model after
a regularization or a (Lagrangian) relaxation. In this case, we reset the
integer and binary type of the relaxed variables and the bounds and then fix
them again to their originally fixed values.

The third one is used for the binary case.
"""
function prepare_state_fixing!(node::SDDP.Node, state_name::Symbol)
    if JuMP.is_binary(node.states[state_name].in)
        JuMP.unset_binary(node.states[state_name].in)
    elseif JuMP.is_integer(node.states[state_name].in)
        JuMP.unset_integer(node.states[state_name].in)
    end
    return
end

function prepare_state_fixing!(node::SDDP.Node, state::SDDP.State)
    if JuMP.is_binary(state.in)
        JuMP.unset_binary(state.in)
    elseif JuMP.is_integer(state.in)
        JuMP.unset_integer(state.in)
    end
end

function prepare_state_fixing_binary!(node::SDDP.Node, state::JuMP.VariableRef)
    if JuMP.is_binary(state)
        JuMP.unset_binary(state)
    elseif JuMP.is_integer(state)
        JuMP.unset_integer(state)
    end
     return
end

################################################################################

"""
Get the outgoing state which can be used on the following stage
to set the incoming state.

Requires node.subproblem to have been solved with PrimalStatus == FeasiblePoint.

I'm not sure if I actually need this function outside of NCNBD or if I could
simply use the SDDP one.
"""
function get_outgoing_state(node::SDDP.Node)
    values = Dict{Symbol,Float64}()
    for (name, state_comp) in node.states
        # To fix some cases of numerical infeasiblities, if the outgoing value
        # is outside its bounds, project the value back onto the bounds. There
        # is a pretty large (×5) penalty associated with this check because it
        # typically requires a call to the solver. It is worth reducing
        # infeasibilities though.
        outgoing_value = JuMP.value(state_comp.out)
        if JuMP.has_upper_bound(state_comp.out)
            current_bound = JuMP.upper_bound(state_comp.out)
            if current_bound < outgoing_value
                outgoing_value = current_bound
            end
        end
        if JuMP.has_lower_bound(state_comp.out)
            current_bound = JuMP.lower_bound(state_comp.out)
            if current_bound > outgoing_value
                outgoing_value = current_bound
            end
        end
        values[name] = outgoing_value
    end
    return values
end

################################################################################

"""
This method is the counterpart to prepare_state_fixing!().

It makes sure that if a state variable is unfixed (e.g. during the regularization
or binarization process), the bounds and integer type associated with this
state originally are reintroduced if required.

In the forward pass, we can decide on this using the copy_regime in the
regularization_regime. In the backward_pass, the original state variables are
just relaxed without any bounds, integrality constraints etc. Those constraints
are only required for the (relaxed) binary variables.
"""
function follow_state_unfixing!(state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.StateSpaceCopy)

    if variable_info.has_lb
        JuMP.set_lower_bound(state.in, variable_info.lower_bound)
    else
        # avoid unboundedness
        JuMP.set_lower_bound(state.in, -1e9)
    end
    if variable_info.has_ub
        JuMP.set_upper_bound(state.in, variable_info.upper_bound)
    else
        # avoid unboundedness
        JuMP.set_upper_bound(state.in, 1e9)
    end
    if variable_info.binary
        JuMP.set_binary(state.in)
    elseif variable_info.integer
        JuMP.set_integer(state.in)
    end

    return
end

function follow_state_unfixing!(state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.ConvexHullCopy)

    if variable_info.has_lb
        JuMP.set_lower_bound(state.in, variable_info.lower_bound)
    elseif variable_info.binary
        JuMP.set_lower_bound(state.in, 0)
    else
        # avoid unboundedness
        JuMP.set_lower_bound(state.in, -1e9)
    end
    if variable_info.has_ub
        JuMP.set_upper_bound(state.in, variable_info.upper_bound)
    elseif variable_info.binary
        JuMP.set_upper_bound(state.in, 1.0)
    else
        # avoid unboundedness
        JuMP.set_upper_bound(state.in, 1e9)
    end

    return
end

function follow_state_unfixing!(state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.NoBoundsCopy)

    # avoid unboundedness
    JuMP.set_lower_bound(state.in, -1e9)
    JuMP.set_upper_bound(state.in, 1e9)

    return
end

"""
This method is the counterpart to prepare_state_fixing!(), but for variables
in the lifted binary space.

It makes sure that if a state variable is unfixed (e.g. during the regularization
or binarization process), the bounds and integer type associated with this
state originally are reintroduced if required.

We can decide using the copy_regime parameter in the duality_regime, whether
all constraints, only bounds or no constraints are imposed.
"""

function follow_state_unfixing_binary!(state::JuMP.VariableRef, copy_regime::DynamicSDDiP.StateSpaceCopy)

    JuMP.set_binary(state)

    return
end

function follow_state_unfixing_binary!(state::JuMP.VariableRef, copy_regime::DynamicSDDiP.ConvexHullCopy)

    JuMP.set_lower_bound(state, 0)
    JuMP.set_upper_bound(state, 1)

    return
end

function follow_state_unfixing_binary!(state::JuMP.VariableRef, copy_regime::DynamicSDDiP.NoBoundsCopy)

    # Aim: Impose no constraints for the state variable.
    # However, we have to set some bounds to ensure feasibility.
    JuMP.set_lower_bound(state, 0)
    JuMP.set_upper_bound(state, 1e9)

    return
end



################################################################################

"""
Function to store the bound and integer information of the .out state of the
previous stage in the .in state of the current stage.
This is required as we - in contrast to SDDP.jl - need also bounds for the
incoming states when they are relaxed.
"""
function set_up_state_in_info!(state_out_previous_stage::JuMP.VariableRef, state_in::JuMP.VariableRef)

    if JuMP.has_lower_bound(state_out_previous_stage)
        lower_bound = JuMP.lower_bound(state_out_previous_stage)
        JuMP.set_lower_bound(state_in, lower_bound)
    end
    if JuMP.has_upper_bound(state_out_previous_stage)
        upper_bound = JuMP.upper_bound(state_out_previous_stage)
        JuMP.set_upper_bound(state_in, upper_bound)
    end
    if JuMP.is_binary(state_out_previous_stage)
        JuMP.set_binary(state_in)
    elseif JuMP.is_integer(state_out_previous_stage)
        JuMP.set_integer(state_in)
    end

    return
end


"""
Function to identify which constraints contain the in-component of the state variables.
Note that we consider mixed-integer linear programs, so this only includes affine
constraints so far.
"""
function identify_state_constraints!(node::SDDP.Node)

    node.ext[:state_constraints] = JuMP.ConstraintRef[]

    for constraint in JuMP.all_constraints(node.subproblem, JuMP.AffExpr, MOI.GreaterThan{Float64})
        check_constraint_for_state(node, constraint, node.ext[:state_constraints])
    end

    for constraint in JuMP.all_constraints(node.subproblem, JuMP.AffExpr, MOI.LessThan{Float64})
        check_constraint_for_state(node, constraint, node.ext[:state_constraints])
    end

    for constraint in JuMP.all_constraints(node.subproblem, JuMP.AffExpr, MOI.EqualTo{Float64})
        check_constraint_for_state(node, constraint, node.ext[:state_constraints])
    end

end


function check_constraint_for_state(
    node::SDDP.Node,
    constraint::JuMP.ConstraintRef,
    storage_vector::Vector{JuMP.ConstraintRef},
    )

    for (name, state) in node.states
        if JuMP.normalized_coefficient(constraint, state.in) != 0.0
            push!(storage_vector, constraint)
            break
        end
    end

end


"""
Functions to binarize the state space (statically!) after a given number
of iterations
"""
function apply_late_binarization_nodes!(model::SDDP.PolicyGraph{T}, algo_params::DynamicSDDiP.AlgoParams,) where {T}

    # Go through all nodes of the problem and apply binarization
    for (node_index, children) in model.nodes
        node = model.nodes[node_index]

        # Apply binarization to the subproblem of this node
        apply_late_binarization_node!(node, node_index, algo_params.late_binarization_regime.K)

        # Update the Bellman function state attribute
        new_states = Dict(key => var.out for (key, var) in node.states)
        if algo_params.cut_type == SDDP.SINGLE_CUT
            node.bellman_function.global_theta.states = new_states
        else
            for theta in node.bellman_function.local_thetas
                theta.states = new_states
            end
        end

    end

    return
end


function apply_late_binarization_node!(node::SDDP.Node, node_index::Int64, K::Int64)

    model = SDDP.get_policy_graph(node.subproblem)

    # Specific step for first node
    ############################################################################
    if node_index == 1
        # Fix original state variable to initial_root_state
        for (i, (name, state)) in enumerate(model.nodes[1].states)
            JuMP.fix(state.in, model.initial_root_state[name])
        end

        # Re-set initial root_state
        model.initial_root_state = Dict{Symbol,Float64}()
    end

    # Step I for all nodes: Copy the original states
    ############################################################################
    """ We need this because later on we cannot iterate over node.states
    anymore, since we have already added the new binary states."""
    original_states = copy(node.states)

    # Step II for all nodes: Get the number of required binary variables
    ############################################################################
    K_total = 0
    for (i, (name, state)) in enumerate(original_states)
        # get variable_info (NOTE: This is problematic, since in principle the in-variable can have different properties)
        variable_info = node.ext[:state_info_storage][name].out
        K_total += calculate_K(variable_info, K)
    end

    # Step III for all nodes: Add the new state variables
    ############################################################################
    """ Note that we have to do this all at once instead of state by state,
    since we cannot use anonymous variables for SDDP.State, but also cannot
    register the same name again and again in subproblem """
    binary_var = JuMP.@variable(node.subproblem, binary_var[i in 1:K_total], SDDP.State, Bin, initial_value = 0, base_name = "stat_bin_")

    # Step IV for all nodes: Add the binary approximation constraints
    ############################################################################
    # states to delete
    states_to_delete = Symbol[]
    K_so_far = 0

    # iterate over all original states
    for (i, (name, state)) in enumerate(original_states)
        # get variable_info (NOTE: This is problematic, since in principle the in-variable can have different properties)
        variable_info = node.ext[:state_info_storage][name].out
        # used K for this state
        K_current = calculate_K(variable_info, K)

        # apply binarization to state
        K_so_far = apply_late_binarization_state!(node, node_index, state, variable_info, binary_var, K_current, K_so_far)

        # STORE STATE TO BE DELETED
        # (should not be done here, as it changes node.states while iterating over it)
        ########################################################################
        push!(states_to_delete, name)
    end

    @assert(K_so_far == K_total)

    # Delete states from node.states
    ############################################################################
    for name in states_to_delete
        delete!(node.states, name)
        delete!(node.ext[:state_info_storage], name)
    end

    # For remaining (i.e. new) states, add variable info to dictionary
    ############################################################################
    for (i, (name, state)) in enumerate(node.states)
        variable_info_in = get_variable_info(state.in)
        variable_info_out = get_variable_info(state.out)
        node.ext[:state_info_storage][name] = DynamicSDDiP.StateInfoStorage(variable_info_in, variable_info_out)
    end

    return
end


function calculate_K(
    variable_info::DynamicSDDiP.VariableInfo,
    K::Int64,
    )

    if variable_info.binary
        """ Note that we still introduce a new binary variable for completeness """
        return 1
    else
        if !isfinite(variable_info.upper_bound) || !variable_info.has_ub
            error("When using a binary expansion, state variables require an upper bound.")
        end

        if variable_info.integer
            # no matter which value K takes we use an exact representation here
            return SDDP._bitsrequired(variable_info.upper_bound)

        else #continuous case
            # Get binary precision beta based on K
            return K
        end
    end
end


function apply_late_binarization_state!(
    node::SDDP.Node,
    node_index::Int64,
    state_variable::SDDP.State,
    variable_info::DynamicSDDiP.VariableInfo,
    binary_var::Vector{SDDP.State{JuMP.VariableRef}},
    K_current::Int64,
    K_so_far::Int64,
    )

    subproblem = node.subproblem

    """ Note that the new state variables are automatically added as states
    in node.states and that the initial_root_state is automatically stored."""

    if variable_info.binary
        JuMP.@constraint(subproblem, state_variable.out == binary_var[K_so_far+1].out)
        # NOTE: This may cause issues if the state's initial value is greater than 0
        JuMP.@constraint(subproblem, state_variable.in == binary_var[K_so_far+1].in)

    else
        if variable_info.integer
            JuMP.@constraint(subproblem, state_variable.out == SDDP.bincontract([binary_var[i].out for i in K_so_far+1:K_so_far+K_current]))
            # NOTE: This may cause issues if the state's initial value is greater than 0
            JuMP.@constraint(subproblem, state_variable.in == SDDP.bincontract([binary_var[i].in for i in K_so_far+1:K_so_far+K_current]))

        else #continuous case
            # Get binary precision beta based on K
            #beta = variable_info.upper_bound / (2^K_current - 1)
            beta = 1.0
            JuMP.@constraint(subproblem, state_variable.out == SDDP.bincontract([binary_var[i].out for i in K_so_far+1:K_so_far+K_current], beta))
            # NOTE: This may cause issues if the state's initial value is greater than 0
            JuMP.@constraint(subproblem, state_variable.in == SDDP.bincontract([binary_var[i].in for i in K_so_far+1:K_so_far+K_current], beta))

        end
    end

    # We should also unfix the original states (they are somehow still fixed from the last iteration)
    JuMP.unfix(state_variable.in)
    follow_state_unfixing!(state_variable, variable_info, DynamicSDDiP.StateSpaceCopy())

    return K_so_far + K_current
end
