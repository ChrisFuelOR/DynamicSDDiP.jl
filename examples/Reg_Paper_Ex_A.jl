# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

using JuMP
using SDDP
using DynamicSDDiP
using Revise
#using Gurobi
using GAMS
#using SCIP
using Infiltrator


function model_config()

    # Stopping rules to be used
    stopping_rules = [DynamicSDDiP.DeterministicStopping()]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax()
    dual_choice_regime = DynamicSDDiP.StandardChoice()
    copy_regime = DynamicSDDiP.ConvexHullCopy()

    duality_regime = DynamicSDDiP.LagrangianDuality(
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

    # State approximation and cut projection configuration
    #cut_projection_regime = DynamicSDDiP.BigM()
    #binary_precision = Dict{Symbol, Float64}()

    # State approximation and cut projection configuration
    state_approximation_regime = DynamicSDDiP.NoStateApproximation()

    # Cut generation regimes
    cut_generation_regime = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime,
        duality_regime = duality_regime,
    )

    cut_generation_regimes = [cut_generation_regime]

    # Regularization configuration
    #regularization_regime = DynamicSDDiP.NoRegularization()
    regularization_regime = DynamicSDDiP.Regularization(sigma=[0.0,2.0], sigma_factor=5.0)

    # Cut aggregation regime
    cut_aggregation_regime = DynamicSDDiP.SingleCutRegime()

    cut_type = SDDP.SINGLE_CUT
    if isa(cut_aggregation_regime,DynamicSDDiP.MultiCutRegime)
        cut_type = SDDP.MULTI_CUT
    end

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/Reg_Paper_Ex_A.log"

    # Suppress solver output
    silent = false

    # Infiltration for debugging
    infiltrate_state = :none

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
    )

    # Define solvers to be used
    applied_solvers = DynamicSDDiP.AppliedSolvers(
        LP = "Gurobi",
        MILP = "Gurobi",
        MIQCP = "Gurobi",
        MINLP = "SCIP",
        NLP = "Gurobi",
        Lagrange = "Gurobi",
    )

    # Start model with used configuration
    model_starter(
        algo_params,
        applied_solvers,
    )
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

    number_of_stages = 2

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min
    ) do subproblem, t

        ########################################################################
        # DEFINE STAGE-t MODEL
        ########################################################################
        # State variables
        JuMP.@variable(subproblem, 0.0 <= x <= 1.0, SDDP.State, Bin, initial_value = 0)

        if t == 1

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out >= x.in)
            #JuMP.@constraint(subproblem, x.out <= 0.5)

            # Stage objective
            SDDP.@stageobjective(subproblem, -x.out)

        else
            # Local variables
            JuMP.@variable(subproblem, 0.0 <= y[i=1:2])
            JuMP.set_integer(y[1])
            JuMP.set_upper_bound(y[1], 2)
            JuMP.set_upper_bound(y[2], 3)

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out == 0)
            JuMP.@constraint(subproblem, con, 2*y[1] + y[2] >= 3*x.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, y[1] + y[2])

        end

    end

    return model
end

model_config()
