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
    # iteration limit 50? time limit 36000?

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Rigorous()
    dual_choice_regime = DynamicSDDiP.MagnantiWongChoice()
    duality_regime = DynamicSDDiP.LagrangianDuality(
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
    )

    # State approximation and cut projection configuration
    cut_projection_regime = DynamicSDDiP.BigM()
    binary_precision = Dict{Symbol, Float64}()

    state_approximation_regime = DynamicSDDiP.BinaryApproximation(
                                    binary_precision = binary_precision,
                                    cut_projection_regime = cut_projection_regime)

    # Regularization configuration
    regularization_regime = DynamicSDDiP.Regularization(sigma = [0.0, 10.0], sigma_factor = 5.0)

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/UC_2_Leo.log"

    # Suppress solver output
    silent = true

    # Infiltration for debugging
    infiltrate_state = :none

    # Definition of algo_params
    algo_params = DynamicSDDiP.AlgoParams(
        stopping_rules = stopping_rules,
        state_approximation_regime = state_approximation_regime,
        regularization_regime = regularization_regime,
        duality_regime = duality_regime,
        cut_selection_regime = cut_selection_regime,
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

    DynamicSDDiP.solve(model, algo_params, applied_solvers)
end


function model_definition()

    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min
    ) do subproblem, t

        ########################################################################
        # ALL STAGES
        ########################################################################
        # State variables
        JuMP.@variable(
            subproblem,
            0.0 <= commit <= 1.0,
            SDDP.State,
            Bin,
            initial_value = 0
        )

        JuMP.@variable(
            subproblem,
            0.0 <= gen <= 5000,
            SDDP.State,
            initial_value = 0
        )

        ########################################################################
        # STAGE 1
        ########################################################################
        if t == 1

            SDDP.@stageobjective(subproblem, 2*gen.out)

            JuMP.@constraint(subproblem, gen.out == 116.5)
            JuMP.@constraint(subproblem, commit.out == 1.0)
        end

        ########################################################################
        # STAGE 2
        ########################################################################
        if t == 2

            # demand slack
            JuMP.@variable(subproblem, demand_slack >= 0.0)
            JuMP.@variable(subproblem, neg_demand_slack >= 0.0)
            demand_penalty = 10000.0

            # generation bounds
            JuMP.@constraint(subproblem, genmin, gen.out >= commit.out * 0)
            JuMP.@constraint(subproblem, genmax, gen.out <= commit.out * 5000.0)

            # ramping
            JuMP.@constraint(subproblem, rampup, gen.out - gen.in <= 1000.0)
            JuMP.@constraint(subproblem, rampdown, gen.in - gen.out <= 1000.0)

            # load balance
            demand_1 = 130.5073368017375    # 126.81106813153158   # 130.5073368017375
            demand_2 = 69.45250729123958    # 68.97237105734341    # 69.45250729123958
            JuMP.@constraint(subproblem, load, gen.out + demand_slack - neg_demand_slack == demand_1 + demand_2 )

            demand_slack = subproblem[:demand_slack]
            neg_demand_slack = subproblem[:neg_demand_slack]
            SDDP.@stageobjective(subproblem, 2*gen.out +  demand_slack * demand_penalty + neg_demand_slack * demand_penalty)
        end

    end

    return model
end

model_config()
