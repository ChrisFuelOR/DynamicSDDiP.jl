import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

include("scenario_tree.jl")

struct GeneratorType
    build_cost::Float64
    fuel_price::Float64
    heat_rate::Float64
    efficiency::Float64
    om_cost::Float64
    inst_capacity::Float64
    av_capacity::Float64
    maximum_number::Int
end

function model_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)

    """
    This model is based on the papers by Jin et al. (2011) and Zou et al. (2019).
    We consider six different types of generators:
    base load (coal), CC, CT, nuclear, wind and IGCC.

    Note that we define all costs as specific costs per hour (in order to have
    smaller coefficients in the model), so the actual total cost is obtained by
    multiplication by 8760. Moreover, all costs are expressed in 1000 dollars per
    GWh for better scaling of coefficients.

    This means that also all other parameters and variables, such as demand,
    generator production and capacity, build/fuel/O&M cost are expressed with
    respect to 1000 dollars or GWh, respectively.
    """

    generators = [
        GeneratorType(1446000.0/8760.0, 3.37, 8.844, 0.4, 4.7, 1.200, 1.13, 4), # Base load (Coal)
        GeneratorType(795000.0/8760.0, 9.11, 7.196, 0.56, 2.11, 0.4, 0.39, 10), # CC
        GeneratorType(575000.0/8760.0, 9.11, 10.842, 0.4, 3.66, 0.4, 0.38, 10), # CT
        GeneratorType(1613000.0/8760.0, 0.00093, 10.4, 0.45, 0.51, 1.2, 1.18, 1), # Nuclear
        GeneratorType(1650000.0/8760.0, 0.0, 1.0, 1.0, 5.0, 0.5, 0.175, 45), # Wind
        GeneratorType(1671000.0/8760.0, 3.37, 8.613, 0.48, 2.98, 0.6, 0.56, 4), #IGCC
    ]

    num_gen_types = length(generators)

    # Annual interest rate
    r = 0.08

    # Demand penalty
    demand_penalty = 1000

    # Define number of binary expansion for each generator type
    # (we do not need beta, since the states are pure integer)
    K = [3, 4, 4, 1, 6, 3]

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
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
            0 <= gen_built_agg_out[i in 1:num_gen_types] <= generators[i].maximum_number,
            #SDDP.State,
            Int
            #initial_value = 0
        )

        JuMP.@variable(
            subproblem,
            0 <= gen_built_agg_in[i in 1:num_gen_types] <= generators[i].maximum_number,
            #SDDP.State,
            Int
            #initial_value = 0
        )

        # BINARY EXPANSION CONSTRAINTS
        ########################################################################
        JuMP.@constraint(subproblem, [i in 1:num_gen_types], gen_built_agg_out[i] == sum(2^(k-1) * λ[i,k].out for k in 1:K[i]))
        if t==1
            JuMP.@constraint(subproblem, [i in 1:num_gen_types], gen_built_agg_in[i] == 0.0)
        else
            JuMP.@constraint(subproblem, [i in 1:num_gen_types], gen_built_agg_in[i] == sum(2^(k-1) * λ[i,k].in for k in 1:K[i]))
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

        # fuel cost variable
        JuMP.@variable(
            subproblem,
            fuel_cost[i in 1:num_gen_types]
        )

        # CONSTRAINTS
        ########################################################################
        # Generation capacity constraint
        JuMP.@constraint(
            subproblem,
            capacity[i in 1:num_gen_types],
            gen_built_agg_out[i] * generators[i].av_capacity >= production[i]
        )

        # Demand satisfaction
        JuMP.@constraint(
            subproblem,
            load,
            sum(production[i] for i in 1:num_gen_types) + unmet_demand == demand
        )

        # State equation
        JuMP.@constraint(
            subproblem,
            state_equation[i in 1:num_gen_types],
            gen_built_agg_out[i] == gen_built_agg_in[i] + gen_built[i]
        )

        # PARAMETERIZE THE RANDOM VARIABLES
        ########################################################################
        # Determine the original fuel costs for all (especially non-gas) generators
        JuMP.@variable(subproblem, fuel_costs[i in 1:num_gen_types])
        JuMP.@constraint(subproblem, fuel_cost_con[i in 1:num_gen_types], fuel_costs[i] == (generators[i].fuel_price * generators[i].heat_rate * 1/generators[i].efficiency + generators[i].om_cost) * production[i] )

        # Get the support and probability for the current stage
        support = scenario_tree.support_array[t]
        probability = scenario_tree.probabilities_array[t]

        # Parameterize the model using the uncertain demand and the natural gas prices
        SDDP.parameterize(subproblem, support, probability) do ω
               JuMP.fix(demand, ω.dem)
               JuMP.set_normalized_coefficient(fuel_cost_con[2], production[2], -ω.gas / 1.028 * generators[2].heat_rate * 1/generators[2].efficiency + generators[2].om_cost)
               JuMP.set_normalized_coefficient(fuel_cost_con[3], production[3], -ω.gas / 1.028 * generators[3].heat_rate * 1/generators[3].efficiency + generators[3].om_cost)
        end

        # STAGE OBJECTIVE
        ########################################################################
        SDDP.@stageobjective(
            subproblem,
            (sum(fuel_costs[i] + generators[i].build_cost * gen_built[i] * generators[i].inst_capacity for i in 1:num_gen_types)
            + unmet_demand * demand_penalty) / (1+r)^(t-1)
        )

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
