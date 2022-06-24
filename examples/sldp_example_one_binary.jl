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
    stopping_rules = [SDDP.TimeLimit(3600)]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Rigorous()
    dual_choice_regime = DynamicSDDiP.StandardChoice()
    #normalization_regime = DynamicSDDiP.Core_Epsilon(perturb=1e-2)
    normalization_regime = DynamicSDDiP.Core_Optimal()
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
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/sldp_example_one_binary.log"

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

    ############################################################################
    # SOLVE MODEL
    ############################################################################
    DynamicSDDiP.solve(model, algo_params, applied_solvers)
end


function model_definition()

    number_of_stages = 8
    K = 8
    beta = 20/(2^K -1)
    #beta = 0.1 #0.01
    #K = 8 #11
    #slack_penalize = 100

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # Auxiliary binary state variables
        JuMP.@variable(subproblem, λ[1:K], SDDP.State, Bin, initial_value = 0)

        # Add former state variable as ordinary variable
        JuMP.@variable(subproblem, x_out, lower_bound = -10.0, upper_bound = 10.0)
        JuMP.@variable(subproblem, x_in, lower_bound = -10.0, upper_bound = 10.0)

        # Add slack variables to ensure complete continuous recourse
        JuMP.@variable(subproblem, slack⁺ >= 0)
        JuMP.@variable(subproblem, slack⁻ >= 0)

        # Add binary approximation constraints
        JuMP.@constraint(subproblem, x_out == -10 + sum(2^(k-1) * beta * λ[k].out for k in 1:K))
        if t==1
            JuMP.@constraint(subproblem, x_in == 2.0)
        else
            JuMP.@constraint(subproblem, x_in == -10 + sum(2^(k-1) * beta * λ[k].in for k in 1:K))
        end

        # Keep the rest of the constraints as it is
        JuMP.@variables(subproblem, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)

        # Add slack penalization factor (see SLDP.jl)
        slack_penalize = 1.0 + (1.0 - 0.9^(number_of_stages+1-t))/(1 - 0.9)*0.9^(t-1)

        SDDP.@stageobjective(subproblem, 0.9^(t - 1) * (x⁺ + x⁻) + slack_penalize * (slack⁺ + slack⁻))
        JuMP.@constraints(subproblem, begin
            x_out == x_in + 2 * u - 1 + ω + slack⁺ - slack⁻
            x⁺ >= x_out
            x⁻ >= -x_out
        end)
        points = [
            -0.3089653673606697,
            -0.2718277412744214,
            -0.09611178608243474,
            0.24645863921577763,
            0.5204224537256875,
        ]
        SDDP.parameterize(φ -> JuMP.fix(ω, φ), subproblem, [points; -points])

        JuMP.set_silent(subproblem)

        return
    end

    return model
end


Random.seed!(11111)
model_config()

# Random.seed!(11112)
# model_config()
#
# Random.seed!(11113)
# model_config()
#
# Random.seed!(11114)
# model_config()
#
# Random.seed!(11115)
# model_config()
