import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

include("scenario_tree.jl")


function model_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)

    """
    This model is the same as analyzed in the paper on Stochastic Lipschitz
    Dynamic Programming by Ahmed, Cabral and Costa (2022).
    The data was provided by Filipe Cabral.
    """

    # Define number of binary expansion variables
    K = 10
    beta = 1.0 #600/(2^K-1)

    # Further parameters
    number_of_products = 3

    # Time parameters
    setup_time = [15 10 15]
    production_time = [1 1 1]

    # Cost parameters
    setup_cost = [60 120 80]
    inventory_cost = [2 3 1]
    lostsale_cost = [20 30 10]

    # Average demand
    demand_avg = [40 35 0;
                  45 35 60;
                  65 55 45;
                  35 35 80;
                  #40 35 0;
                  45 35 60;
                  65 55 45;
                  35 35 80;
                  #40 35 0;
                  45 35 60;
                  65 55 45;
                  35 35 80]

    # Production time capacity
    capacity = 175

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # Auxiliary binary state variables
        JuMP.@variable(subproblem, λ[1:number_of_products, 1:K], SDDP.State, Bin, initial_value = 0)

        # Add former state variable as ordinary variable (inventory)
        JuMP.@variable(subproblem, inventory_out[1:number_of_products], lower_bound = 0.0, upper_bound = 600.0)
        JuMP.@variable(subproblem, inventory_in[1:number_of_products], lower_bound = 0.0, upper_bound = 600.0)

        # Add binary approximation constraints
        JuMP.@constraint(subproblem, [i=1:number_of_products], inventory_out[i] == sum(2^(k-1) * beta * λ[i,k].out for k in 1:K))
        if t==1
            JuMP.@constraint(subproblem, [i=1:number_of_products], inventory_in[i] == 0.0)
        else
            JuMP.@constraint(subproblem, [i=1:number_of_products], inventory_in[i] == sum(2^(k-1) * beta * λ[i,k].in for k in 1:K))
        end

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
        JuMP.@constraint(subproblem, [i=1:number_of_products], inventory_out[i] - lostsales[i] == inventory_in[i] + production[i] - demand[i])

        # Stage objective
        SDDP.@stageobjective(subproblem, sum( setup_cost[i] * setup[i] + inventory_cost[i] * inventory_out[i] + lostsale_cost[i] * lostsales[i] for i in 1:number_of_products))

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
        end

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_set_up(
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
    model = model_definition(problem_params, scenario_tree)

    return (model = model, problem_params = problem_params)
end
