# The functions
# > "set_objective",
# > "parameterize",
# are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

# Internal function: set the objective of node to the stage objective, plus the
# cost/value-to-go term.
function set_objective(subproblem::JuMP.Model)
    node = SDDP.get_node(subproblem)
    objective_state_component = SDDP.get_objective_state_component(node)
    belief_state_component = SDDP.get_belief_state_component(node)
    if objective_state_component != JuMP.AffExpr(0.0) ||
       belief_state_component != JuMP.AffExpr(0.0)
        node.stage_objective_set = false
    end
    if !node.stage_objective_set
        JuMP.set_objective(
            subproblem,
            JuMP.objective_sense(subproblem),
            @expression(
                subproblem,
                node.stage_objective +
                objective_state_component +
                belief_state_component +
                bellman_term(node.bellman_function)
            )
        )
    end
    node.stage_objective_set = true
    return
end

function parameterize(node::SDDP.Node, noise)
    node.parameterize(noise)
    # set objective function and Bellman function for MILP
    set_objective(node.subproblem)
    return
end
