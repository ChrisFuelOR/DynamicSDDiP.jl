# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Process for binary refinements.
"""
function binary_refinement(
    model::SDDP.PolicyGraph{T},
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    sampled_states::Vector{Dict{Symbol,Float64}},
    algo_params::DynamicSDDiP.AlgoParams,
    solution_check::Bool,
    binary_refinement_status::Symbol,
    bound_check::Bool,
    ) where {T}

    # (1) BINARY APPROXIMATION USED?
    ############################################################################
    # Check if at least one cut_generation_regime uses BinaryApproximation
    refinement_check_required = false

    for cut_generation_regime in algo_params.cut_generation_regimes
        if isa(cut_generation_regime.state_approximation_regime, DynamicSDDiP.BinaryApproximation)
            refinement_check_required = true
        end
    end

    # (2) BINARY REFINEMENT REQUIRED?
    ############################################################################
    # Check if the binary precision has to be refined.
    # This is the case if the forward pass solution and the lower bound did not change during
    # the last iteration.
    # This is independent of the individual cut_generation_regime.

    if refinement_check_required
        if !isnothing(previous_solution)
                solution_check = DynamicSDDiP.binary_refinement_check(
                    model,
                    previous_solution,
                    sampled_states,
                    solution_check,
                    DynamicSDDiP.BinaryApproximation()
                )
        end

        # (3) BINARY REFINEMENT
        ########################################################################
        # This is done for all cut_generation_regimes with BinaryApproximation
        # and for all stages.
        if solution_check || bound_check
            for cut_generation_regime in algo_params.cut_generation_regime
                if isa(cut_generation_regime.state_approximation_regime)
                    binary_refinement = DynamicSDDiP.binary_refinement(
                        model,
                        cut_generation_regime.state_approximation_regime.binary_precision,
                        binary_refinement_status,
                    )
                end
            end
        end
    end

    return binary_refinement_status

    # TODO binary_refinement symbol should exist separately for all regimes
end


"""
Check if an increase of the number of binary variables is required
if binary approximation is used.
"""
function binary_refinement_check(
    model::SDDP.PolicyGraph{T},
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    sampled_states::Vector{Dict{Symbol,Float64}},
    solution_check::Bool,
    #bound_check::Bool,
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
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
Check if an increase of the number of binary variables is required
if no binary approximation is used.
"""
function binary_refinement_check(
    model::SDDP.PolicyGraph{T},
    previous_solution::Union{Vector{Dict{Symbol,Float64}},Nothing},
    sampled_states::Vector{Dict{Symbol,Float64}},
    refinement_check::Bool,
    #bound_check::Bool,
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
    ) where {T}

    return false
end

"""
Executing a binary refinement: Increasing the number of binary variables to
approximate the states
"""
function binary_refinement_execution(
    model::SDDP.PolicyGraph{T},
    binary_precision::Dict{Symbol, Float64},
    binary_refinement::Symbol,
    ) where {T}

    all_refined = Int[]

    # Consider stage 2 here (should be the same for all following stages)
    # Precision is only used (and increased) for continuous variables
    for (name, state_comp) in model.nodes[2].states

        variable_info = model.nodes[2].ext[:state_info_storage][name].in

        if !variable_info.binary && !variable_info.integer
            current_prec = binary_precision[name]
            ub = variable_info.upper_bound
            K = SDDP._bitsrequired(round(Int, ub / current_prec))
            new_prec = ub / sum(2^(k-1) for k in 1:K+1)

            # Only apply refinement if int64 is appropriate to represent this number
            if ub / new_prec > 2^63 - 1
                push!(all_refined, 0)
            else
                push!(all_refined, 1)
                binary_precision[name] = new_prec
            end

        else
            push!(all_refined, 2)
        end
    end

    # Re-set the list of stored Benders coefficients for all nodes
    for (key, node) in model.nodes
        if key != model.root_node
            node.ext[:Benders_cuts_binary] = Tuple{Int64, Symbol}[]
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
