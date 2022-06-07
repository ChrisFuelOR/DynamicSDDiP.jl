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
    end
    if variable_info.has_ub
        JuMP.set_upper_bound(state.in, variable_info.upper_bound)
    end

    return
end

function follow_state_unfixing!(state::SDDP.State, variable_info::DynamicSDDiP.VariableInfo, copy_regime::DynamicSDDiP.NoBoundsCopy)

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
