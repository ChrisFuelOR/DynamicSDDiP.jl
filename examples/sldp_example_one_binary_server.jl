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
using Printf


function model_config(
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    seed::Int,
    )

    Random.seed!(seed)

    # Stopping rules to be used
    stopping_rules = [SDDP.TimeLimit(time_limit)]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax()
    dual_choice_regime = DynamicSDDiP.MinimalNormChoice()
    #dual_subproblemace_regime = DynamicSDDiP.BendersSpanSpaceRestriction(10, :multi_cut)
    dual_space_regime = DynamicSDDiP.NoDualSpaceRestriction()
    copy_regime = DynamicSDDiP.ConvexHullCopy()

    if duality_regime_sym == :uni_lag
        duality_regime = DynamicSDDiP.UnifiedLagrangianDuality(
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
            user_dual_objective_bound = 1e4,
        )
    elseif duality_regime_sym == :lag
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
    end

    # State approximation and cut projection configuration
    state_approximation_regime = DynamicSDDiP.NoStateApproximation()

    # Cut generation regimes
    cut_generation_regime = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime,
        duality_regime = duality_regime,
    )

    cut_generation_regimes = [cut_generation_regime]

    # Regularization configuration
    regularization_regime = DynamicSDDiP.NoRegularization()

    # Cut aggregation regime
    cut_type = SDDP.SINGLE_CUT
    if isa(cut_aggregation_regime, DynamicSDDiP.MultiCutRegime)
        cut_type = SDDP.MULTI_CUT
    end

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
        seed = seed,
        run_description = "SLDP Example with finer binary approximation and minimal norm choice"
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
    K = 11
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
        # points = [
        #     0.74142424122442444,
        #     0.59237849150344212,
        #     0.15298371036066613,
        #     -0.24645863921577763,
        #     -0.36119334780580125,
        # ]
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

function model_starter_single(
    num::Int,
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    seed::Int
    )

    try
        model_config(duality_regime_sym, normalization_regime, cut_aggregation_regime, cut_selection_regime, log_file, time_limit, seed)
    catch e
        @printf "Case %d terminated with error" num
        println()
        #throw(error(e))
        showerror(stdout, e, catch_backtrace())
        println()
        println("#############################################################")
        println()
    end
end


function model_starter_server()

    model_starter_single(0,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 200, 11111)

    # model_starter_single(1,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11111)
    # model_starter_single(2,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11111)
    # model_starter_single(3,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11111)
    # model_starter_single(4,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11111)
    # model_starter_single(5,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11111)
    model_starter_single(6,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11111)
    model_starter_single(7,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11111)
    model_starter_single(8,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11111)

    # model_starter_single(9,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11112)
    # model_starter_single(10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11112)
    # model_starter_single(11,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11112)
    # model_starter_single(12,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11112)
    # model_starter_single(13,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11112)
    # model_starter_single(14,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11112)
    # model_starter_single(15,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11112)
    # model_starter_single(16,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11112)
    #
    # model_starter_single(17,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11113)
    # model_starter_single(18,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11113)
    # model_starter_single(19,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11113)
    # model_starter_single(20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11113)
    # model_starter_single(21,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11113)
    # model_starter_single(22,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11113)
    # model_starter_single(23,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11113)
    # model_starter_single(24,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11113)

    model_starter_single(25,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11111)
    # model_starter_single(26,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11112)
    # model_starter_single(27,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11113)

    # model_starter_single(28,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11111)
    # model_starter_single(29,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11112)
    # model_starter_single(30,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11113)

    # model_starter_single(31,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11113)
    # model_starter_single(32,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11113)

    # model_starter_single(6,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11111)
    # model_starter_single(7,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11111)
    #
    # model_starter_single(14,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11112)
    # model_starter_single(15,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11112)
    #
    # model_starter_single(22,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11113)
    # model_starter_single(23,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11113)
    #
    # model_starter_single(25,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11111)
    # model_starter_single(26,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11112)
    # model_starter_single(27,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11113)
    #
    # model_starter_single(28,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11111)
    # model_starter_single(29,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11112)
    # model_starter_single(30,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11113)

    # model_starter_single(31,:uni_lag, DynamicSDDiP.Fischetti(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11113)


end

model_starter_server()
