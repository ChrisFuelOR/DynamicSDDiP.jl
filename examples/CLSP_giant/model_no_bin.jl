import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

include("scenario_tree.jl")


function model_no_bin_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)

    """
    This model is the same as analyzed in the paper on Stochastic Lipschitz
    Dynamic Programming by Ahmed, Cabral and Costa (2022).
    The data was provided by Filipe Cabral.

    This is the model version without binary approximation of the state variables.
    """

    # Further parameters
    number_of_products = 20

    # Time parameters
    setup_time = [15 10 15 10 10 20 20 15 10 15 15 18 20 25 20 15 10 10 20 15]
    production_time = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]

    # Cost parameters
    setup_cost = [60 120 80 100 70 85 130 90 50 75 70 100 80 55 60 120 105 85 90 65]
    inventory_cost = [2 3 1 2 2 3 4 1 1 2 3 1 2 4 1 1 1 2 3 2]
    lostsale_cost = [20 30 10 15 10 20 10 30 40 25 10 20 15 15 10 20 35 40 40 30]

    # Average demand
    demand_avg = [40 35 0 20 30 45 25 0 10 20 15 0 30 20 20 50 35 30 0 10;
                  45 35 60 30 25 45 60 50 40 65 30 30 55 40 50 60 45 40 20 20;
                  65 55 45 35 60 70 55 35 50 40 40 50 55 35 55 75 40 50 30 15;
                  35 35 80 30 55 50 45 45 60 40 30 35 50 20 45 60 25 30 25 20;
                  #40 35 0 20 30 45 25 0 10 20;
                  45 35 60 30 25 45 60 50 40 65 30 30 55 40 50 60 45 40 20 20;
                  65 55 45 35 60 70 55 35 50 40 40 50 55 35 55 75 40 50 30 15;
                  35 35 80 30 55 50 45 45 60 40 30 35 50 20 45 60 25 30 25 20;
                  #40 35 0 20 30 45 25 0 10 20;
                  45 35 60 30 25 45 60 50 40 65 30 30 55 40 50 60 45 40 20 20;
                  65 55 45 35 60 70 55 35 50 40 40 50 55 35 55 75 40 50 30 15;
                  35 35 80 30 55 50 45 45 60 40 30 35 50 20 45 60 25 30 25 20;
                  #40 35 0 20 30 45 25 0 10 20;
                  45 35 60 30 25 45 60 50 40 65 30 30 55 40 50 60 45 40 20 20;
                  65 55 45 35 60 70 55 35 50 40 40 50 55 35 55 75 40 50 30 15;
                  35 35 80 30 55 50 45 45 60 40 30 35 50 20 45 60 25 30 25 20;
                  #40 35 0 20 30 45 25 0 10 20;
                  45 35 60 30 25 45 60 50 40 65 30 30 55 40 50 60 45 40 20 20;
                  65 55 45 35 60 70 55 35 50 40 40 50 55 35 55 75 40 50 30 15;
                  35 35 80 30 55 50 45 45 60 40 30 35 50 20 45 60 25 30 25 20
                  ]

    # Production time capacity
    capacity = 175/3*20

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # Add former state variable as ordinary variable (inventory)
        JuMP.@variable(subproblem, inventory[1:number_of_products], SDDP.State, lower_bound = 0.0, upper_bound = 600.0, initial_value = 0)

        # Production variables
        JuMP.@variable(subproblem, production[1:number_of_products] >= 0)

        # Lost sales variables
        JuMP.@variable(subproblem, lostsales[1:number_of_products] >= 0)

        # Set-up variables
        JuMP.@variable(subproblem, setup[1:number_of_products], Bin)

        # Random variable for demand
        JuMP.@variable(subproblem, demand[1:number_of_products])

        # Production only if machine is set-up (using Big-M formulation with M = 175)
        JuMP.@constraint(subproblem, [i=1:number_of_products], production[i] <= 175 * setup[i])

        # Time constraint
        JuMP.@constraint(subproblem, sum( production_time[i] * production[i] + setup_time[i] * setup[i] for i in 1:number_of_products) <= capacity)

        # Inventory balance equation
        JuMP.@constraint(subproblem, [i=1:number_of_products], inventory[i].out - lostsales[i] == inventory[i].in + production[i] - demand[i])

        # Stage objective
        SDDP.@stageobjective(subproblem, sum( setup_cost[i] * setup[i] + inventory_cost[i] * inventory[i].out + lostsale_cost[i] * lostsales[i] for i in 1:number_of_products))

        # PARAMETERIZE THE RANDOM VARIABLES
        ########################################################################
        # Get the support and probability for the current stage
        support = scenario_tree.support_array[t]
        probability = scenario_tree.probabilities_array[t]

        # Parameterize the demand
        SDDP.parameterize(subproblem, support, probability) do ω
            JuMP.fix(demand[1], ω.xi1 * demand_avg[t,1])
            JuMP.fix(demand[2], ω.xi2 * demand_avg[t,2])
            JuMP.fix(demand[3], ω.xi3 * demand_avg[t,3])
            JuMP.fix(demand[4], ω.xi4 * demand_avg[t,4])
            JuMP.fix(demand[5], ω.xi5 * demand_avg[t,5])
            JuMP.fix(demand[6], ω.xi6 * demand_avg[t,6])
            JuMP.fix(demand[7], ω.xi7 * demand_avg[t,7])
            JuMP.fix(demand[8], ω.xi8 * demand_avg[t,8])
            JuMP.fix(demand[9], ω.xi9 * demand_avg[t,9])
            JuMP.fix(demand[10], ω.xi10 * demand_avg[t,10])
            JuMP.fix(demand[11], ω.xi1 * demand_avg[t,11])
            JuMP.fix(demand[12], ω.xi2 * demand_avg[t,12])
            JuMP.fix(demand[13], ω.xi3 * demand_avg[t,13])
            JuMP.fix(demand[14], ω.xi4 * demand_avg[t,14])
            JuMP.fix(demand[15], ω.xi5 * demand_avg[t,15])
            JuMP.fix(demand[16], ω.xi6 * demand_avg[t,16])
            JuMP.fix(demand[17], ω.xi7 * demand_avg[t,17])
            JuMP.fix(demand[18], ω.xi8 * demand_avg[t,18])
            JuMP.fix(demand[19], ω.xi9 * demand_avg[t,19])
            JuMP.fix(demand[20], ω.xi10 * demand_avg[t,20])
        end

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_no_bin_set_up(
    number_of_stages::Int,
    number_of_realizations::Int;
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
    tree_seed::Int = 12345
)

    ############################################################################
    # DEFINE PROBLEM PARAMS
    ############################################################################
    problem_params = DynamicSDDiP.ProblemParams(number_of_stages, number_of_realizations, tree_seed = tree_seed)

    ############################################################################
    # GET FINITE SCENARIO TREE FOR MODEL
    ############################################################################
    scenario_tree = get_recombining_scenario_tree(algo_params, problem_params)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_no_bin_definition(problem_params, scenario_tree)

    return (model = model, problem_params = problem_params)
end
