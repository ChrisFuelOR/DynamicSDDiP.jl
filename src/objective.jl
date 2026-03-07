# The functions
# > "set_objective",
# > "parameterize",
# are derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2026 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2026 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

# Internal function: set the objective of node to the stage objective, plus the
# cost/value-to-go term.
function set_objective(subproblem::JuMP.Model)
    node = SDDP.get_node(subproblem)
    if !node.stage_objective_set
        JuMP.set_objective(
            subproblem,
            JuMP.objective_sense(subproblem),
            JuMP.@expression(
                subproblem,
                node.stage_objective +
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
