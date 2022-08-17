import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

function model_definition(problem_params::DynamicSDDiP.ProblemParams)

    """
    This model is the same as Example 2 (illustrative example) in the SDDiP
    paper by Zou, Ahmed and Sun from 2019.
    """

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # State variables
        JuMP.@variable(subproblem, 0.0 <= x <= 1.0, SDDP.State, Bin, initial_value = 0)

        if t == 1
            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out >= x.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, -x.out)

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

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_set_up(
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
)

    ############################################################################
    # DEFINE PROBLEM PARAMS
    ############################################################################
    problem_params = DynamicSDDiP.ProblemParams(2, 1, tree_seed = tree_seed)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_definition(problem_params)

    return (model = model, problem_params = problem_params)
end
