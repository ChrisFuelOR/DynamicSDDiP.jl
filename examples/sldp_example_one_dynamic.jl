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
    stopping_rules = [SDDP.IterationLimit(12)]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Rigorous()
    dual_choice_regime = DynamicSDDiP.StandardChoice()
    #normalization_regime = DynamicSDDiP.Core_Epsilon(perturb=0.0)
    normalization_regime = DynamicSDDiP.L₁_Deep()
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
    duality_regime_4 = DynamicSDDiP.StrengthenedDuality()

    # State approximation and cut projection configuration
    state_approximation_regime_1 = DynamicSDDiP.NoStateApproximation()

    # Cut generation regimes
    cut_generation_regime = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime_1,
        duality_regime = duality_regime_1,
    )

    cut_generation_regimes = [cut_generation_regime]

    # Regularization configuration
    regularization_regime = DynamicSDDiP.NoRegularization()

    # Cut aggregation regime
    cut_aggregation_regime = DynamicSDDiP.MultiCutRegime()

    cut_type = SDDP.SINGLE_CUT
    if isa(cut_aggregation_regime, DynamicSDDiP.MultiCutRegime)
        cut_type = SDDP.MULTI_CUT
    end

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.CutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/sldp_example_one_dynamic.log"

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
        run_numerical_stability_report = false,
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
    #Random.seed!(11111) # issues
    Random.seed!(11112) # issues
    #Random.seed!(11113)

    ############################################################################
    # SOLVE MODEL
    ############################################################################
    DynamicSDDiP.solve(model, algo_params, applied_solvers)
end


function model_definition()

    number_of_stages = 8

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min,
    ) do subproblem, t

        JuMP.@variable(subproblem, x, SDDP.State, initial_value = 2.0)
        JuMP.@variables(subproblem, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)
        SDDP.@stageobjective(subproblem, 0.9^(t - 1) * (x⁺ + x⁻))
        JuMP.@constraints(subproblem, begin
            x.out == x.in + 2 * u - 1 + ω
            x⁺ >= x.out
            x⁻ >= -x.out
        end)
        points = [
            -0.3089653673606697,
            -0.2718277412744214,
            -0.09611178608243474,
            0.24645863921577763,
            0.5204224537256875,
        ]
        JuMP.set_silent(subproblem)
        SDDP.parameterize(φ -> JuMP.fix(ω, φ), subproblem, [points; -points])

        return
    end

    return model
end

model_config()
