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
    The data was provided with the SLDP.jl package.
    """

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # State variables
        JuMP.@variable(subproblem, x, SDDP.State, lower_bound = -10.0, upper_bound = 10.0, initial_value = 0.0)

        # Add slack variables to ensure complete continuous recourse
        JuMP.@variable(subproblem, slack⁺ >= 0)
        JuMP.@variable(subproblem, slack⁻ >= 0)

        # Keep the rest of the constraints as it is
        JuMP.@variables(subproblem, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)

        # Add slack penalization factor (see SLDP.jl)
        slack_penalize = 1.0 + (1.0 - 0.9^(problem_params.number_of_stages+1-t))/(1 - 0.9)*0.9^(t-1)

        SDDP.@stageobjective(subproblem, 0.9^(t - 1) * (x⁺ + x⁻) + slack_penalize * (slack⁺ + slack⁻))
        JuMP.@constraints(subproblem, begin
            x.out == x.in + 2 * u - 1 + ω + slack⁺ - slack⁻
            x⁺ >= x.out
            x⁻ >= -x.out
        end)

        # Parameterize the random variable
        # Get the support and probability for the current stage
        support = scenario_tree.support_array[t]
        probability = scenario_tree.probabilities_array[t]
        if t == 1
            support = [0.0]
            probability = [1.0]
        end

        # Parameterize the model using the uncertain demand and the natural gas prices
        SDDP.parameterize(φ -> JuMP.fix(ω, φ), subproblem, support, probability)

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
    model = model_definition(problem_params, scenario_tree)

    return (model = model, problem_params = problem_params)
end
