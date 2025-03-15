# The function
# > "convergence_test"
# is derived from similar named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

# ======================= Deterministic Stopping Rule ====================== #

stopping_rule_status(::DeterministicStopping) = :DeterministicStopping

function convergence_test(graph::SDDP.PolicyGraph, log::Vector{Log}, rule::DeterministicStopping)

    bool_rtol = abs(log[end].best_upper_bound - log[end].lower_bound)/abs(max(log[end].best_upper_bound, log[end].lower_bound)) <= rule.rtol
    bool_atol = log[end].best_upper_bound - log[end].lower_bound <= rule.atol
    bool_neg = log[end].best_upper_bound - log[end].lower_bound >= -1e-4

    return (bool_rtol || bool_atol) && bool_neg
end

# ======================= Iteration Limit Stopping Rule ====================== #
stopping_rule_status(::SDDP.IterationLimit) = :iteration_limit

function convergence_test(graph::SDDP.PolicyGraph, log::Vector{Log}, rule::SDDP.IterationLimit)
    return log[end].iteration >= rule.limit
end

# ========================= Time Limit Stopping Rule ========================= #
stopping_rule_status(::SDDP.TimeLimit) = :time_limit

function convergence_test(graph::SDDP.PolicyGraph, log::Vector{Log}, rule::SDDP.TimeLimit)
    return log[end].time >= rule.limit
end

# ========================= Bound Stalling Stopping Rule ========================= #
stopping_rule_status(::SDDP.BoundStalling) = :bound_stalling

function convergence_test(
    graph::SDDP.PolicyGraph{T},
    log::Vector{Log},
    rule::SDDP.BoundStalling,
) where {T}
    if length(log) < rule.num_previous_iterations + 1
        return false
    end
    # No change in the bound. There are three possibilities:
    #  1) we haven't added enough cuts
    #  2) the problem was deterministic or myopic
    #  3) there were existing cuts
    existing_solves = log[1].total_solves > log[end].total_solves / length(log)
    if !existing_solves && isapprox(log[1].lower_bound, log[end].lower_bound; atol = 1e-6)
        return all(l -> isapprox(l.lower_bound, l.best_upper_bound; atol = 1e-6), log)
    end
    for i in 1:rule.num_previous_iterations
        if abs(log[end-i].lower_bound - log[end-i+1].lower_bound) > rule.tolerance
            return false
        end
    end
    return true
end

# ============================================================================ #
function convergence_test(
    graph::SDDP.PolicyGraph,
    log::Vector{Log},
    stopping_rules::Vector{SDDP.AbstractStoppingRule},
)
    for stopping_rule in stopping_rules
        if convergence_test(graph, log, stopping_rule)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end
