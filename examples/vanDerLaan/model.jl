import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise


function model_definition(problem_params::DynamicSDDiP.ProblemParams)

    """
    This model is the same as Example in the paper on Scaled Cuts
    by van der Laan and Romeijnders.
    """

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        ########################################################################
        # DEFINE STAGE-t MODEL
        ########################################################################
        # State variables
        JuMP.@variable(subproblem, 0.0 <= x <= 4.0, SDDP.State, initial_value = 0)

        if t == 1

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out == 0)

            # Stage objective
            SDDP.@stageobjective(subproblem, x.out)

        else
            # Local variables
            JuMP.@variable(subproblem, y >= 0)
            JuMP.set_integer(y)

            # Uncertainty variable
            JuMP.@variable(subproblem, rhs)

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out == 0)
            JuMP.@constraint(subproblem, con, y >= rhs - x.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, 2*y)

            # Parameterize the RHS
            SDDP.parameterize(subproblem, [2.5, 3], [0.5, 0.5]) do ω
                   JuMP.fix(rhs, ω)
            end

        end

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
    # DEFINE MODEL
    ############################################################################
    model = model_definition(problem_params)

    return (model = model, problem_params = problem_params)
end
