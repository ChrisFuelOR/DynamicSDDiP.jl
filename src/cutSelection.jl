# The functions
# > "_cut_selection_update"
# > "_eval_height"
# are derived from equally named functions in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Trivial cut selection function if no cut selection is used for nonlinear cuts.
"""
# Internal function: update the Level-One datastructures inside `bellman_function`.
function _cut_selection_update(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.NonlinearCut,
    anchor_state::Dict{Symbol,Float64},
    trial_state::Dict{Symbol,Float64},
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    cut_selection_regime::DynamicSDDiP.NoCutSelection
)

    ############################################################################
    # ADD CUTS AND STATES TO THE ORACLE
    ############################################################################
    sampled_state_anchor = DynamicSDDiP.SampledState(anchor_state, cut, NaN)
    sampled_state_anchor.best_objective = _eval_height(node, cut, sampled_state_anchor, applied_solvers, algo_params)
    sampled_state_trial = DynamicSDDiP.SampledState(trial_state, cut, NaN)
    sampled_state_trial.best_objective = _eval_height(node, cut, sampled_state_trial, applied_solvers, algo_params)

    push!(V.sampled_states, sampled_state_anchor)
    push!(V.sampled_states, sampled_state_trial)

    push!(V.cuts, cut)

    ############################################################################
    # DETERMINE NUMBER OF CUTS FOR LOGGING
    ############################################################################
    count_cuts(node, V)

    return
end

"""
Trivial cut selection function if no cut selection is used for linear cuts.
"""
# Internal function: update the Level-One datastructures inside `bellman_function`.
function _cut_selection_update(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.LinearCut,
    anchor_state::Dict{Symbol,Float64},
    trial_state::Dict{Symbol,Float64},
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    cut_selection_regime::DynamicSDDiP.NoCutSelection
)

    ############################################################################
    # ADD CUTS AND STATES TO THE ORACLE
    ############################################################################
    sampled_state_trial = DynamicSDDiP.SampledState(trial_state, cut, NaN)
    sampled_state_trial.best_objective = _eval_height(node, cut, sampled_state_trial, applied_solvers, algo_params)

    push!(V.sampled_states, sampled_state_trial)

    push!(V.cuts, cut)

    ############################################################################
    # DETERMINE NUMBER OF CUTS FOR LOGGING
    ############################################################################
    count_cuts(node, V)

    return
end


"""
Simple cut selection feature for nonlinear cuts.
Note that we only compare nonlinear cuts with each other.
"""
function _cut_selection_update(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.NonlinearCut,
    anchor_state::Dict{Symbol,Float64},
    trial_state::Dict{Symbol,Float64},
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    cut_selection_regime::DynamicSDDiP.CutSelection
)

    ############################################################################
    # GET MODEL INFORMATION
    ############################################################################
    model = JuMP.owner_model(V.theta)
    is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE

    ############################################################################
    # GET TRIAL STATE AND ANCHOR STATE
    ############################################################################
    """
    By considering both type of states, we have way more states than cuts,
    so we may not eliminiate that many cuts. On the other hand,
    only considering the trial points is not sufficient, because it may lead to
    creating the same cuts over and over again if the binary approximation
    is not refined. Only considering the anchor points is also not sufficient,
    since we are actually interested in good approximations in the trial points.
    """
    sampled_state_anchor = DynamicSDDiP.SampledState(anchor_state, cut, NaN)
    sampled_state_anchor.best_objective = _eval_height(node, cut, sampled_state_anchor, applied_solvers, algo_params)
    sampled_state_trial = DynamicSDDiP.SampledState(trial_state, cut, NaN)
    sampled_state_trial.best_objective = _eval_height(node, cut, sampled_state_trial, applied_solvers, algo_params)

    ############################################################################
    # LOOP THROUGH PREVIOUSLY VISITED STATES (ANCHOR OR TRIAL STATES)
    ############################################################################
    """
    We loop through previously sampled states and compare the height of the
    most recent cut with the current best.
    If the new cut yields an improvement, store this one instead.
    """
    for old_state in V.sampled_states
        height = _eval_height(node, cut, old_state, applied_solvers, algo_params)
        if SDDP._dominates(height, old_state.best_objective, is_minimization)
            old_state.dominating_cut.non_dominated_count -= 1
            cut.non_dominated_count += 1
            old_state.dominating_cut = cut
            old_state.best_objective = height
        end
    end
    push!(V.sampled_states, sampled_state_anchor)
    push!(V.sampled_states, sampled_state_trial)

    ############################################################################
    # LOOP THROUGH PREVIOUSLY VISITED CUTS
    ############################################################################
    """
    Now we loop through previously discovered cuts and compare their height
    at the existing states. If a cut is an improvement, add it to a queue to be added.
    """
    for old_cut in V.cuts
        if isa(old_cut,DynamicSDDiP.LinearCut)
            # We only care about other nonlinear cuts.
            continue
        end

        if !isempty(old_cut.cut_constraints)
            # We only care about cuts not currently in the model.
            continue
        end

        # For anchor state
        height = _eval_height(node, old_cut, sampled_state_anchor, applied_solvers, algo_params)
        if SDDP._dominates(height, sampled_state_anchor.best_objective, is_minimization)
            sampled_state_anchor.dominating_cut.non_dominated_count -= 1
            old_cut.non_dominated_count += 1
            sampled_state_anchor.dominating_cut = old_cut
            sampled_state_anchor.best_objective = height
            _add_cut_constraints_to_models(node, V, old_cut, algo_params, cut_generation_regime, infiltrate_state)
        end

        # For trial state
        height = _eval_height(node, old_cut, sampled_state_trial, applied_solvers, algo_params)
        if SDDP._dominates(height, sampled_state_trial.best_objective, is_minimization)
            sampled_state_trial.dominating_cut.non_dominated_count -= 1
            old_cut.non_dominated_count += 1
            sampled_state_trial.dominating_cut = old_cut
            sampled_state_trial.best_objective = height
            _add_cut_constraints_to_models(node, V, old_cut, algo_params, cut_generation_regime, infiltrate_state)
        end

    end
    push!(V.cuts, cut)

    ############################################################################
    # DETERMINE CUTS TO BE DELETED
    ############################################################################
    for cut in V.cuts
        if cut.non_dominated_count < 1
            if !isempty(cut.cut_constraints)
                push!(V.cuts_to_be_deleted, cut)
            end
        end
    end

    ############################################################################
    # DELETE CUTS
    ############################################################################
    if length(V.cuts_to_be_deleted) >= V.deletion_minimum
        for cut in V.cuts_to_be_deleted
            for constraint_ref in cut.cut_constraints
                JuMP.delete(model, constraint_ref)
            end
            for variable_ref in cut.cut_variables
                JuMP.delete(model, variable_ref)
            end
            cut.cut_variables = JuMP.VariableRef[]
            cut.cut_constraints = JuMP.ConstraintRef[]

            cut.non_dominated_count = 0
        end
    end
    empty!(V.cuts_to_be_deleted)

    ############################################################################
    # DETERMINE NUMBER OF CUTS FOR LOGGING
    ############################################################################
    count_cuts(node, V)

    return
end


"""
Simple cut selection feature for linear cuts.
Note that we only compare linear cuts with each other.
"""
# Internal function: update the Level-One datastructures inside `bellman_function`.
function _cut_selection_update(
    node::SDDP.Node,
    V::DynamicSDDiP.CutApproximation,
    cut::DynamicSDDiP.LinearCut,
    anchor_state::Dict{Symbol,Float64},
    trial_state::Dict{Symbol,Float64},
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams,
    infiltrate_state::Symbol,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    cut_selection_regime::DynamicSDDiP.CutSelection
)

    ############################################################################
    # GET MODEL INFORMATION
    ############################################################################
    model = JuMP.owner_model(V.theta)
    is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE

    ############################################################################
    # GET TRIAL STATE
    ############################################################################
    sampled_state_trial = DynamicSDDiP.SampledState(trial_state, cut, NaN)
    sampled_state_trial.best_objective = _eval_height(node, cut, sampled_state_trial, applied_solvers, algo_params)

    ############################################################################
    # LOOP THROUGH PREVIOUSLY VISITED STATES (ANCHOR OR TRIAL STATES)
    ############################################################################
    """
    We loop through previously sampled states and compare the height of the
    most recent cut with the current best.
    If the new cut yields an improvement, store this one instead.
    """
    for old_state in V.sampled_states
        height = _eval_height(node, cut, old_state, applied_solvers, algo_params)
        if SDDP._dominates(height, old_state.best_objective, is_minimization)
            old_state.dominating_cut.non_dominated_count -= 1
            cut.non_dominated_count += 1
            old_state.dominating_cut = cut
            old_state.best_objective = height
        end
    end
    push!(V.sampled_states, sampled_state_trial)

    ############################################################################
    # LOOP THROUGH PREVIOUSLY VISITED CUTS
    ############################################################################
    """
    Now we loop through previously discovered cuts and compare their height
    at the existing states. If a cut is an improvement, add it to a queue to be added.
    """
    for old_cut in V.cuts
        if isa(old_cut,DynamicSDDiP.NonlinearCut)
            # We only care about other linear cuts.
            continue
        end

        if !isnothing(old_cut.cut_constraint)
            # We only care about cuts not currently in the model.
            continue
        end

        # For trial state
        height = _eval_height(node, old_cut, sampled_state_trial, applied_solvers, algo_params)
        if SDDP._dominates(height, sampled_state_trial.best_objective, is_minimization)
            sampled_state_trial.dominating_cut.non_dominated_count -= 1
            old_cut.non_dominated_count += 1
            sampled_state_trial.dominating_cut = old_cut
            sampled_state_trial.best_objective = height
            _add_cut_constraints_to_models(node, V, old_cut, algo_params, cut_generation_regime, infiltrate_state)
        end
    end
    push!(V.cuts, cut)

    ############################################################################
    # DETERMINE CUTS TO BE DELETED
    ############################################################################
    for cut in V.cuts
        if cut.non_dominated_count < 1
            if cut.cut_constraint !== nothing
                push!(V.cuts_to_be_deleted, cut)
            end
        end
    end

    ############################################################################
    # DELETE CUTS
    ############################################################################
    if length(V.cuts_to_be_deleted) >= V.deletion_minimum
        for cut in V.cuts_to_be_deleted
            JuMP.delete(model, cut.cut_constraint)
            cut.cut_constraint = nothing
            cut.non_dominated_count = 0
        end
    end
    empty!(V.cuts_to_be_deleted)

    ############################################################################
    # DETERMINE NUMBER OF CUTS FOR LOGGING
    ############################################################################
    count_cuts(node, V)

end


function count_cuts(node::SDDP.Node, V::DynamicSDDiP.CutApproximation)
    node.ext[:total_cuts] += size(V.cuts, 1)
    counter = 0

    for cut in V.cuts
        counter = count_cut(cut, counter)
    end
    node.ext[:active_cuts] += counter

    return

end

function count_cut(cut::DynamicSDDiP.NonlinearCut, counter::Int64)

    if !isempty(cut.cut_constraints)
        counter += 1
    end

    return counter
end

function count_cut(cut::DynamicSDDiP.LinearCut, counter::Int64)

    if cut.cut_constraint !== nothing
        counter += 1
    end

    return counter
end


"""
Evaluating some nonlinear cut at a specific point.
"""
# Internal function: calculate the height of `cut` evaluated at `state`.
function _eval_height(node::SDDP.Node, cut::DynamicSDDiP.NonlinearCut,
    sampled_state::DynamicSDDiP.SampledState, applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams)

    ############################################################################
    # CREATE CUT-PROJECTION CLOSURE LP TO EVALUATE CUT
    ############################################################################
    # Create a new JuMP model to evaluate the height of a non-convex cut
    model = JuMP.Model()
    model.ext[:sddp_policy_graph] = node.subproblem.ext[:sddp_policy_graph]
    set_solver!(model, algo_params, applied_solvers, :cut_selection, algo_params.solver_approach)

    # Storages for coefficients and binary states
    binary_state_storage = JuMP.VariableRef[]
    all_coefficients = Float64[]
    binary_variables_so_far = 0

    for (i, (state_name, value)) in enumerate(sampled_state.state)
        variable_info = node.ext[:state_info_storage][state_name].out

        if variable_info.binary
            ####################################################################
            # BINARY CASE
            ####################################################################
            # introduce one binary variable to the model
            binary_var = JuMP.@variable(model)
            push!(binary_state_storage, binary_var)
            binary_variables_so_far += 1

            # introduce binary expansion constraint to the model
            binary_constraint = JuMP.@constraint(model, binary_var == value)

            # determine the correct cut coefficient
            related_coefficient = 0.0
            for (bin_name, value) in cut.coefficients
                if cut.binary_state[bin_name].x_name == state_name
                    related_coefficient = cut.coefficients[bin_name]
                end
            end
            push!(all_coefficients, related_coefficient)

        else
            ####################################################################
            # NON-BINARY CASE
            ####################################################################
            if !isfinite(variable_info.upper_bound) || !variable_info.has_ub
                error("When using SDDiP, state variables require an upper bound.")
            end

            # Get K and beta
            if variable_info.integer
                beta = 1
                K = SDDP._bitsrequired(variable_info.upper_bound)
            else
                beta = cut.binary_precision[state_name]
                K = SDDP._bitsrequired(round(Int, variable_info.upper_bound / beta))
            end

            # introduce binary variables to the model
            binary_var = JuMP.@variable(model, [k in 1:K], lower_bound=0, upper_bound=1)
            append!(binary_state_storage, binary_var)
            binary_variables_so_far += K

            # introduce binary expansion constraint to the model
            binary_constraint = JuMP.@constraint(model, SDDP.bincontract([binary_var[k] for k in 1:K], beta) == value)

            # determine the correct cut coefficient
            related_coefficients = Vector{Float64}(undef, K)
            for (bin_name, value) in cut.coefficients
                if cut.binary_state[bin_name].x_name == state_name
                    index = cut.binary_state[bin_name].k
                    related_coefficients[index] = cut.coefficients[bin_name]
                end
            end
            append!(all_coefficients, related_coefficients)

        end
    end

    ############################################################################
    # SOME VALIDITY CHECKS
    ############################################################################
    @assert(size(all_coefficients, 1) == size(binary_state_storage, 1)
                == binary_variables_so_far
                == size(collect(values(cut.coefficients)),1)
                )

    ############################################################################
    # ADD OBJECTIVE TO THE MODEL
    ############################################################################
    objective_sense_stage = JuMP.objective_sense(node.subproblem)
    eval_sense = (
        objective_sense_stage == JuMP.MOI.MIN_SENSE ? JuMP.MOI.MAX_SENSE : JuMP.MOI.MIN_SENSE
    )

    JuMP.@objective(
        model, eval_sense,
        cut.intercept + sum(all_coefficients[j] * binary_state_storage[j] for j in 1:binary_variables_so_far)
    )

    ############################################################################
    # SOLVE MODEL AND RETURN SOLUTION
    ############################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "solver_call_cut_selection" begin
        JuMP.optimize!(model)
    end
    height = JuMP.objective_value(model)

    # redo scaling by π₀ cause otherwise, cuts are not comparable to each other
    height = height / cut.scaling_coeff

    return height

end


"""
Evaluating some linear cut at a specific point.
"""
# Internal function: calculate the height of `cut` evaluated at `state`.
function _eval_height(node::SDDP.Node, cut::DynamicSDDiP.LinearCut,
    sampled_state::DynamicSDDiP.SampledState, applied_solvers::DynamicSDDiP.AppliedSolvers,
    algo_params::DynamicSDDiP.AlgoParams)

    # Start with intercept and increase by evaluating the states
    height = cut.intercept
    for (key, value) in cut.coefficients
        height += value * sampled_state.state[key]
    end

    # redo scaling by π₀ cause otherwise, cuts are not comparable to each other
    height = height / cut.scaling_coeff

    return height

end
