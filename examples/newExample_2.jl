# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

using JuMP
using SDDP
using DynamicSDDiP
using Revise
using Gurobi
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
    regularization_regime = DynamicSDDiP.Regularization(sigma = [0.0, 1.0], sigma_factor = 5.0)
    # regularization_regime = DynamicSDDiP.NoRegularization()

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/newExample_2.log"

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

    ############################################################################
    # DEFINE BINARY APPROXIMATION IF INTENDED
    ############################################################################
    # for (name, state_comp) in model.nodes[1].ext[:lin_states]
    #     ub = JuMP.upper_bound(state_comp.out)
    #
    #     string_name = string(name)
    #     if occursin("gen", string_name)
    #         binaryPrecision[name] = binaryPrecisionFactor * ub
    #     else
    #         binaryPrecision[name] = 1
    #     end
    # end

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
        JuMP.@variable(subproblem, 0.0 <= x <= 1.0, SDDP.State, initial_value = 0)

        if t == 1

            # Constraints
            x = subproblem[:x]
            JuMP.@variable(subproblem, v)
            JuMP.@constraint(subproblem, v >= 0.7 - 6/5*x.out)
            JuMP.@constraint(subproblem, v >= -1.1 + 2.4*x.out)

            # Stage objective
            SDDP.@stageobjective(subproblem, v)

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
