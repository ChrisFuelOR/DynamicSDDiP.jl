# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function SDDP.bellman_term(bellman_function::DynamicSDDiP.BellmanFunction)
    return bellman_function.global_theta.theta
end