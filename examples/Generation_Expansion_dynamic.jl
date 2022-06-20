# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

import JuMP
import SDDP
import DynamicSDDiP
using Revise
import GAMS
import Infiltrator
import Random


function model_config()

    # Stopping rules to be used
    stopping_rules = [SDDP.IterationLimit(20)]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax()
    dual_choice_regime = DynamicSDDiP.StandardChoice()
    #normalization_regime = DynamicSDDiP.Core_Epsilon(perturb=0.0)
    normalization_regime = DynamicSDDiP.Core_Relint()
    #dual_subproblemace_regime = DynamicSDDiP.BendersSpanSpaceRestriction(10, :multi_cut)
    dual_space_regime = DynamicSDDiP.NoDualSpaceRestriction()
    copy_regime = DynamicSDDiP.ConvexHullCopy()

    duality_regime_1 = DynamicSDDiP.UnifiedLagrangianDuality(
        atol = 1e-4,
        rtol = 1e-4,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
        normalization_regime = normalization_regime,
        dual_space_regime = dual_space_regime,
        copy_regime = copy_regime,
    )

    duality_regime_2 = DynamicSDDiP.LagrangianDuality(
        atol = 1e-4,
        rtol = 1e-4,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
        copy_regime = copy_regime,
    )

    duality_regime_3 = DynamicSDDiP.LinearDuality()
    #duality_regime = DynamicSDDiP.StrengthenedDuality()

    # State approximation and cut projection configuration
    state_approximation_regime_1 = DynamicSDDiP.NoStateApproximation()

    # Cut generation regimes
    cut_generation_regime_1 = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime_1,
        duality_regime = duality_regime_1,
    )

    cut_generation_regimes = [cut_generation_regime_1]

    # Regularization configuration
    regularization_regime = DynamicSDDiP.NoRegularization()

    # Cut aggregation regime
    cut_aggregation_regime = DynamicSDDiP.MultiCutRegime()

    cut_type = SDDP.SINGLE_CUT
    if isa(cut_aggregation_regime, DynamicSDDiP.MultiCutRegime)
        cut_type = SDDP.MULTI_CUT
    end

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/Generation_Expansion_dynamic.log"

    # Suppress solver output
    silent = true

    # Infiltration for debugging
    infiltrate_state = :none

    # Solver approach
    solver_approach = DynamicSDDiP.Direct_Solver()

    # Define solvers to be used
    applied_solvers = DynamicSDDiP.AppliedSolvers(
        LP = "Gurobi",
        MILP = "Gurobi",
        MIQCP = "Gurobi",
        MINLP = "Gurobi",
        NLP = "Gurobi",
        Lagrange = "Gurobi",
    )

    # Definition of algo_params
    algo_params = DynamicSDDiP.AlgoParams(
        stopping_rules = stopping_rules,
        regularization_regime = regularization_regime,
        cut_aggregation_regime = cut_aggregation_regime,
        cut_selection_regime = cut_selection_regime,
        cut_generation_regimes = cut_generation_regimes,
        cut_type = cut_type,
        log_file = log_file,
        silent = silent,
        infiltrate_state = infiltrate_state,
        solver_approach = solver_approach,
        numerical_focus = false,
    )

    # Start model with used configuration
    model_starter(algo_params, applied_solvers)
end


function model_starter(
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_definition()

    ## For repeatability
    Random.seed!(11119)

    ############################################################################
    # SOLVE MODEL
    ############################################################################
    DynamicSDDiP.solve(model, algo_params, applied_solvers)
end


function model_definition()

    build_cost = 1
    use_cost = 4 / 10000.0
    num_units = 5
    capacities = ones(num_units)
    demand_vals =
        0.5 * [
            5 5 5 5 5 5 5 5
            4 3 1 3 0 9 8 17
            0 9 4 2 19 19 13 7
            25 11 4 14 4 6 15 12
            6 7 5 3 8 4 17 13
        ]

    # 5 5 5 5 5 5 5 5
    # 4 3 1 3 0 9 8 17
    # 0 9 4 2 19 19 13 7
    # 25 11 4 14 4 6 15 12
    # 6 7 5 3 8 4 17 13

    # Cost of unmet demand
    penalty = 50
    # Discounting rate
    rho = 0.99

    number_of_stages = 5

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min,
    ) do subproblem, stage

        # DEFINE STAGE-t MODEL
        ########################################################################
        # State variables
        JuMP.@variable(
            subproblem,
            0 <= invested[1:num_units] <= 1,
            SDDP.State,
            Int, #TODO: Why not Bin
            initial_value = 0
        )

        # Local variables
        JuMP.@variables(subproblem, begin
            generation >= 0
            unmet >= 0
            demand
        end)

        # Constraints
        JuMP.@constraints(
            subproblem,
            begin
                # Can't un-invest
                investment[i in 1:num_units], invested[i].out >= invested[i].in
                # Generation capacity
                sum(capacities[i] * invested[i].out for i = 1:num_units) >= generation
                # Meet demand or pay a penalty
                unmet >= demand - sum(generation)
                # For fewer iterations order the units to break symmetry, units are identical (tougher numerically)
                [j in 1:(num_units-1)], invested[j].out <= invested[j+1].out
            end
        )

        # Demand is uncertain
        SDDP.parameterize(ω -> JuMP.fix(demand, ω), subproblem, demand_vals[stage, :])

        # Objective function
        JuMP.@expression(
            subproblem,
            investment_cost,
            build_cost * sum(invested[i].out - invested[i].in for i = 1:num_units)
        )
        SDDP.@stageobjective(
            subproblem,
            (investment_cost + generation * use_cost) * rho^(stage - 1) + penalty * unmet
        )
    end

    return model
end

model_config()
