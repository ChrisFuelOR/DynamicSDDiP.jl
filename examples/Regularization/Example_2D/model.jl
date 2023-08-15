import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

function model_definition(problem_params::DynamicSDDiP.ProblemParams)

    """
    This model is the same as analyzed in the paper on Stochastic Lipschitz
    Dynamic Programming by Ahmed, Cabral and Costa (2022).
    The data was provided by Filipe Cabral.
    """

    number_of_stages = 2

    model = SDDP.LinearPolicyGraph(
        stages = number_of_stages,
        lower_bound = -100.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min
    ) do subproblem, t

        ########################################################################
        # DEFINE STAGE-t MODEL
        ########################################################################
        # State variables
        JuMP.@variable(subproblem, 0.0 <= x[i=1:2] <= 5.0, SDDP.State, Int, initial_value = 0)

        if t == 1

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x[2].out <= 9/2 - x[1].out)

            # Stage objective
            SDDP.@stageobjective(subproblem, -1.5 * x[1].out - 4 * x[2].out)

            #JuMP.@constraint(subproblem, x[1].out == 2.0)
            #JuMP.@constraint(subproblem, x[2].out == 1.0)

        else
            # Local variables
            JuMP.@variable(subproblem, y[i=1:4], Bin)

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, con1, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= 10 - 1/3 * x[1].in - 2/3 * x[2].in)
            JuMP.@constraint(subproblem, con2, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= 10 - 2/3 * x[1].in - 1/3 * x[2].in)

            # Important: We also impose the local state constraint on x.in
            # redundant for fixed value, but not redundant in case of copy constraints
            # JuMP.@constraint(subproblem, x[2].in <= 9/2 - x[1].in)
            #JuMP.@constraint(subproblem, x[2].in <= 4 - x[1].in) # convex hull

            # Stage objective
            #SDDP.@stageobjective(subproblem, -16 * y[1] - 19 * y[2] - 23 * y[3] - 28 * y[4])
            SDDP.@stageobjective(subproblem, -1.6 * y[1] - 2.0 * y[2] - 2.8 * y[3] - 4.4 * y[4])

        end

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return

    end

    return model
end

# support = scenario_tree.support_array[t]
# probability = scenario_tree.probabilities_array[t]
#
# # Parameterize the demand
# SDDP.parameterize(subproblem, support, probability) do ω
#        JuMP.fix(demand[1], ω.xi1 * demand_avg[t,1])
#        JuMP.fix(demand[2], ω.xi2 * demand_avg[t,2])
#        JuMP.fix(demand[3], ω.xi3 * demand_avg[t,3])
# end


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
    # DEFINE MODEL
    ############################################################################
    model = model_definition(problem_params)

    return (model = model, problem_params = problem_params)
end
