# This file redefines some functions from the 'JuMP.jl' package by
# Copyright (c) 2017: Iain Dunning, Joey Huchette, Miles Lubin, and contributors
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2017: Iain Dunning, Joey Huchette, Miles Lubin, and contributors

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

function JuMP.build_variable(
    _error::Function,
    info::JuMP.VariableInfo,
    ::Type{SDDP.State};
    initial_value = NaN,
    kwargs...,
)

    if isnan(initial_value)
        _error(
            "When creating a state variable, you must set the " *
            "`initial_value` keyword to the value of the state variable at" *
            " the root node.",
        )
    end
    return SDDP.StateInfo(
        JuMP.VariableInfo(
            false,
            NaN, # lower bound
            false,
            NaN, # upper bound
            false,
            NaN,  # fixed value
            false,
            NaN,  # start value
            false,
            false, # binary and integer
        ),
        info,
        initial_value,
        kwargs,
    )
end

function JuMP.add_variable(
    subproblem::JuMP.Model,
    state_info::SDDP.StateInfo,
    name::String,
)
    state = SDDP.State(
        JuMP.add_variable(
            subproblem,
            JuMP.ScalarVariable(state_info.in),
            name * "_in",
        ),
        JuMP.add_variable(
            subproblem,
            JuMP.ScalarVariable(state_info.out),
            name * "_out",
        ),
    )
    node = SDDP.get_node(subproblem)
    sym_name = Symbol(name)
    @assert !haskey(node.states, sym_name)  # JuMP prevents duplicate names.
    node.states[sym_name] = state
    graph = SDDP.get_policy_graph(subproblem)
    graph.initial_root_state[sym_name] = state_info.initial_value
    return state
end


JuMP.variable_type(model::JuMP.Model, ::Type{State}) = SDDP.State

function JuMP.value(state::State{JuMP.VariableRef})
    return State(JuMP.value(state.in), JuMP.value(state.out))
end
