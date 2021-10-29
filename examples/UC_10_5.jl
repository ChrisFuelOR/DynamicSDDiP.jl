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


struct Generator
    comm_ini::Int
    gen_ini::Float64
    pmax::Float64
    pmin::Float64
    fuel_cost::Float64
    om_cost::Float64
    su_cost::Float64
    sd_cost::Float64
    ramp_up::Float64
    ramp_dw::Float64
end

struct Storage
    level_max::Float64
    level_ini::Float64
    level_end::Float64
    gen_max::Float64
    pump_max::Float64
    gen_eff::Float64
    pump_eff::Float64
end


function model_config()

    # Stopping rules to be used
    stopping_rules = [DynamicSDDiP.DeterministicStopping(atol=5e-4,rtol=5e-4), SDDP.TimeLimit(7200)]
    # iteration limit 50? time limit 36000?

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.LevelBundle()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax()
    dual_choice_regime = DynamicSDDiP.MagnantiWongChoice()
    duality_regime = DynamicSDDiP.LagrangianDuality(
        atol = 1e-4,
        rtol = 1e-4,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
    )

    # State approximation and cut projection configuration
    cut_projection_regime = DynamicSDDiP.SOS1()
    binary_precision = Dict{Symbol, Float64}()

    state_approximation_regime = DynamicSDDiP.BinaryApproximation(
                                    binary_precision = binary_precision,
                                    cut_projection_regime = cut_projection_regime)

    # Regularization configuration
    regularization_regime = DynamicSDDiP.Regularization(sigma = [0.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0], sigma_factor = 2.0)

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.CutSelection()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/UC_10_5.log"

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

    generators = [
        Generator(1, 1.06, 1.19, 0.37, 52.0, 0.0, 177.68, 17.0, 0.31, 0.36),
        Generator(0, 0.0, 1.13, 0.48, 54.0, 0.0, 171.60, 17.0, 0.28, 0.27),
        Generator(0, 0.0, 1.02, 0.47, 49.4, 0.0, 168.04, 17.0, 0.22, 0.275),
        Generator(1, 2.2, 2.82, 0.85, 61.6, 0.0, 486.81, 49.0, 0.9, 0.79),
        Generator(0, 0.0, 3.23, 0.84, 54.9, 0.0, 503.34, 50.0, 1.01, 1.00),
    ]

    demand_penalty = 5e2
    demand = [4.27 4.01 3.69 3.66 3.72 4.01 4.79 5.85 6.84 7.14]

    storages = [
        Storage(1.2, 0.5, 0.7, 0.45, 0.4, 0.9, 0.85),
        Storage(0.8, 0.3, 0.25, 0.35, 0.3, 0.92, 0.87),
    ]

    inflow = [0.2 0.3 0.4; 0.1 0.05 0.1]

    number_of_generators = size(generators,1)
    number_of_storages = 0
    number_of_stages = 10

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
        JuMP.@variable(
            subproblem,
            0.0 <= commit[i = 1:number_of_generators] <= 1.0,
            SDDP.State,
            Bin,
            initial_value = generators[i].comm_ini
        )

        JuMP.@variable(
            subproblem,
            0.0 <= gen[i = 1:number_of_generators] <= generators[i].pmax,
            SDDP.State,
            initial_value = generators[i].gen_ini
        )

        # start-up variables
       JuMP.@variable(subproblem, up[i=1:number_of_generators], Bin)
       JuMP.@variable(subproblem, down[i=1:number_of_generators], Bin)

       # demand slack
       JuMP.@variable(subproblem, demand_slack >= 0.0)
       JuMP.@variable(subproblem, neg_demand_slack >= 0.0)

       # cost variables
       JuMP.@variable(subproblem, startup_costs[i=1:number_of_generators] >= 0.0)
       JuMP.@variable(subproblem, shutdown_costs[i=1:number_of_generators] >= 0.0)
       JuMP.@variable(subproblem, fuel_costs[i=1:number_of_generators] >= 0.0)
       JuMP.@variable(subproblem, om_costs[i=1:number_of_generators] >= 0.0)

       # generation bounds
       JuMP.@constraint(subproblem, genmin[i=1:number_of_generators], gen[i].out >= commit[i].out * generators[i].pmin)
       JuMP.@constraint(subproblem, genmax[i=1:number_of_generators], gen[i].out <= commit[i].out * generators[i].pmax)

       # ramping
       # we do not need a case distinction as we defined initial_values
       JuMP.@constraint(subproblem, rampup[i=1:number_of_generators], gen[i].out - gen[i].in <= generators[i].ramp_up * commit[i].in + generators[i].pmin * (1-commit[i].in))
       JuMP.@constraint(subproblem, rampdown[i=1:number_of_generators], gen[i].in - gen[i].out <= generators[i].ramp_dw * commit[i].out + generators[i].pmin * (1-commit[i].out))

       # start-up and shut-down
       # we do not need a case distinction as we defined initial_values
       JuMP.@constraint(subproblem, startup[i=1:number_of_generators], up[i] >= commit[i].out - commit[i].in)
       JuMP.@constraint(subproblem, shutdown[i=1:number_of_generators], down[i] >= commit[i].in - commit[i].out)

       # additional storage state
        JuMP.@variable(
            subproblem,
            0.0 <= storage_level[j = 1:number_of_storages] <= storages[j].level_max,
            SDDP.State,
            initial_value = storages[j].level_ini,
        )

        # additional storage generation
        JuMP.@variable(
            subproblem,
            0.0 <= storage_gen[j = 1:number_of_storages] <= storages[j].gen_max,
        )

        # additional storage pumping
        JuMP.@variable(
            subproblem,
            0.0 <= storage_pump[j = 1:number_of_storages] <= storages[j].pump_max,
        )

        # additional storage level balance
        JuMP.@constraint(
            subproblem,
            level_balance[j=1:number_of_storages], storage_level[j].out == storage_level[j].in + storage_pump[j] * storages[j].pump_eff - storage_gen[j] / storages[j].gen_eff + inflow[j,t]
        )

        # additional storage end level
        if t == number_of_stages
            JuMP.@constraint(subproblem, storage_end[j=1:number_of_storages], storage_level[j].out >= storages[j].level_end)
        end

        # load balance
        JuMP.@constraint(subproblem, load, sum(gen[i].out for i in 1:number_of_generators) + demand_slack - neg_demand_slack + sum(storage_gen[j] - storage_pump[j] for j in 1:number_of_storages) == demand[t] )

        # costs
        JuMP.@constraint(subproblem, startupcost[i=1:number_of_generators], number_of_stages/24 * generators[i].su_cost * up[i] == startup_costs[i])
        JuMP.@constraint(subproblem, shutdowncost[i=1:number_of_generators], generators[i].sd_cost * down[i] == shutdown_costs[i])
        JuMP.@constraint(subproblem, fuelcost[i=1:number_of_generators], generators[i].fuel_cost * gen[i].out == fuel_costs[i])
        JuMP.@constraint(subproblem, omcost[i=1:number_of_generators], generators[i].om_cost * gen[i].out == om_costs[i])

        # define stage objective
        su_costs = subproblem[:startup_costs]
        sd_costs = subproblem[:shutdown_costs]
        f_costs = subproblem[:fuel_costs]
        om_costs = subproblem[:om_costs]
        demand_slack = subproblem[:demand_slack]
        neg_demand_slack = subproblem[:neg_demand_slack]
        SDDP.@stageobjective(subproblem,
                        sum(su_costs[i] + sd_costs[i] + f_costs[i] + om_costs[i] for i in 1:number_of_generators)
                        + demand_slack * demand_penalty + neg_demand_slack * demand_penalty)

    end

    return model
end

model_config()
