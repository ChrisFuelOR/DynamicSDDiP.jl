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
import Distributions
using Printf


struct GeneratorType
    build_cost::Float64
    fuel_price::Float64
    heat_rate::Float64
    effiency::Float64
    om_cost::Float64
    av_capacity::Float64
    maximum_number::Int
end



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

    # Generator definition (there are six types: Base load, CC, CT, Nuclear, Wind, IGCC)
    generators = [
        Generator(1446000.0, 3.37, 8844.0, 0.4, 4.7, 1130.0, 4), # Base load
        Generator(795000.0, 9.11, 7196.0, 0.56, 2.11, 390.0, 10), # CC
        Generator(575000.0, 9.11, 10842.0, 0.4, 3.66, 380.0, 10), # CT
        Generator(1613000.0, 0.00093, 10400.0, 0.45, 0.51, 1180.0, 1), # Nuclear
        Generator(1650000.0, 0.0, 1.0, 1.0, 5.0, 175.0, 45), # Wind
        Generator(1671000.0, 3.37, 8613.0, 0.48, 2.98, 560.0, 4), #IGCC
    ]

    num_gen_types = length(generators)

    # Time horizon
    number_of_stages = 8

    # Annual interest rate
    r = 0.08

    # Demand penalty
    demand_penalty = 100000

    # Define number of binary expansion for each generator type
    # (we do not need beta, since the states are pure integer)
    K = [3, 4, 4, 1, 6, 3]

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # BINARY STATE VARIABLES
        ########################################################################
        # Note that we introduce a different number of binary variables for each generator type
        JuMP.@variable(
            subproblem,
            λ[i in 1:num_gen_types, j in 1:K[i]],
            SDDP.State,
            Bin,
            initial_value = 0
        )

        # ORIGINAL STATE VARIABLES
        ########################################################################
        # Number of generators built up to (and including) stage t
        JuMP.@variable(
            subproblem,
            0 <= gen_built_agg_out[i in 1:num_gen_types] <= generators[i][:maximum_number],
            #SDDP.State,
            Int,
            initial_value = 0
        )

        JuMP.@variable(
            subproblem,
            0 <= gen_built_agg_in[i in 1:num_gen_types] <= generators[i][:maximum_number],
            #SDDP.State,
            Int,
            initial_value = 0
        )

        # BINARY EXPANSION CONSTRAINTS
        ########################################################################
        JuMP.@constraint(subproblem, [i in 1:num_gen_types], gen_built_agg_out == sum(2^(k-1) * λ[i,k].out for k in 1:K[i]))
        if t==1
            JuMP.@constraint(subproblem, gen_built_agg_in == 0.0)
        else
            JuMP.@constraint(subproblem, [i in 1:num_gen_types], gen_built_agg_in == sum(2^(k-1) * λ[i,k].in for k in 1:K[i]))
        end

        # LOCAL VARIABLES
        ########################################################################
        # Number of generators built in stage t
        JuMP.@variable(
            subproblem,
            0 <= gen_built[i in 1:num_gen_types],
            Int
        )

        # Power generation in stage t
        JuMP.@variable(
            subproblem,
            0 <= production[i in 1:num_gen_types]
        )

        # Slack variable for unmet demand (ensure complete continuous recourse)
        JuMP.@variable(
            subproblem,
            0 <= unmet_demand
        )

        # RANDOM VARIABLES
        ########################################################################
        # RHS uncertainty of demand
        JuMP.@variable(
            subproblem,
            demand
        )

        # CONSTRAINTS
        ########################################################################
        # Generation capacity constraint
        JuMP.@constraint(
            subproblem,
            capacity[i in 1:num_gen_types],
            gen_built_agg[i].out * generators[i][:av_capacity] >= production[i]
        )

        # # Limitation on total number of generators
        # JuMP.@constraint(
        #     subproblem,
        #     gen_limit[i in 1:num_gen_types],
        #     gen_built_agg[i].out <= generators[:maximum_number]
        # )

        # Demand satisfaction
        JuMP.@constraint(
            subproblem,
            load,
            sum(production[i].out for i in 1:num_gen_types) + unmet_demand == demand
        )

        # State equation
        JuMP.@constraint(
            subproblem,
            state_equation[i in 1:num_gen_types],
            gen_built_agg[i].out == gen_built_agg[i].in + gen_built[i]
        )

        # STAGE OBJECTIVE
        ########################################################################
        # First define the fuel costs
        fuel_costs = Array{Float64,1}(undef, num_gen_types)

        for i in num_gen_types
            fuel_costs[i] = generators[i][:fuel_price] * generators[i][:heat_rate] * 1/generators[i][:efficiency]
        end

        SDDP.@stageobjective(
            subproblem,
            (sum((fuel_costs[i] + generators[i][:om_cost]) * production[i].out + generators[i][:build_cost] * gen_built[i] for i in 1:num_gen_types)
            + unmet_demand * demand_penalty) / (1+r)^(t-1)
        )



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
    # model_starter_single(6,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11111)
    # model_starter_single(7,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11111)
    # model_starter_single(8,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11111)

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

    # model_starter_single(25,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11111)
    # model_starter_single(26,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11112)
    # model_starter_single(27,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11113)

    model_starter_single(28,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11111)
    # model_starter_single(29,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11112)
    # model_starter_single(30,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11113)

    # model_starter_single(31,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11113)
    # model_starter_single(32,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11113)

end

model_starter_server()
