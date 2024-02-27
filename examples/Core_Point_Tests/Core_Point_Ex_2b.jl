import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
using Printf
using Gurobi
import Random


function algo_config(
    core_point_approach::Symbol,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    improvement_regime::DynamicSDDiP.AbstractImprovementRegime,
    unbounded_regime::DynamicSDDiP.AbstractUnboundedRegime,
    )

    # Stopping rules to be used
    stopping_rules = [SDDP.IterationLimit(1), DynamicSDDiP.DeterministicStopping()]
    
    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax() #Rigorous
    dual_choice_regime = DynamicSDDiP.StandardChoice()
    dual_space_regime = DynamicSDDiP.NoDualSpaceRestriction()

    # Normalization regime
    if core_point_approach == :midpoint
        dual_normalization_regime = DynamicSDDiP.Core_Midpoint(copy_regime=copy_regime, integer_regime=integer_regime, unbounded_regime=unbounded_regime, improvement_regime=improvement_regime, normalize_direction=false)
    elseif core_point_approach == :relint
        dual_normalization_regime = DynamicSDDiP.Core_Relint(copy_regime=copy_regime, integer_regime=integer_regime, unbounded_regime=unbounded_regime, improvement_regime=improvement_regime, normalize_direction=false)
    elseif core_point_approach == :in_out
        dual_normalization_regime = DynamicSDDiP.Core_In_Out(copy_regime=copy_regime, integer_regime=integer_regime, unbounded_regime=unbounded_regime, improvement_regime=improvement_regime, normalize_direction=false)
    elseif core_point_approach == :eps
        dual_normalization_regime = DynamicSDDiP.Core_Epsilon(perturb=0.01, copy_regime=copy_regime, integer_regime=integer_regime, unbounded_regime=unbounded_regime, improvement_regime=improvement_regime, normalize_direction=false)
    end

    duality_regime = DynamicSDDiP.UnifiedLagrangianDuality(
        atol = 1e-4,
        rtol = 1e-4,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
        copy_regime = copy_regime,
        dual_space_regime = dual_space_regime,
        normalization_regime = dual_normalization_regime,
    )

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
    cut_type = SDDP.MULTI_CUT
    cut_aggregation_regime = DynamicSDDiP.MultiCutRegime()

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # Simulation
    simulation_regime = DynamicSDDiP.NoSimulation()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/Core_Point_1.log"

    # Suppress solver output
    silent = true

    # Infiltration for debugging
    infiltrate_state = :none

    # Solver approach
    solver_approach = DynamicSDDiP.Direct_Solver()

    # Define solvers to be used
    applied_solvers = DynamicSDDiP.AppliedSolvers()

    # Definition of algo_params
    algo_params = DynamicSDDiP.AlgoParams(
        stopping_rules = stopping_rules,
        regularization_regime = regularization_regime,
        cut_aggregation_regime = cut_aggregation_regime,
        cut_selection_regime = cut_selection_regime,
        cut_generation_regimes = cut_generation_regimes,
        simulation_regime = simulation_regime,
        late_binarization_regime = DynamicSDDiP.NoLateBinarization(),
        cut_type = cut_type,
        log_file = log_file,
        silent = silent,
        infiltrate_state = infiltrate_state,
        solver_approach = solver_approach,
        numerical_focus = false,
        run_numerical_stability_report = false,
        seed = 12345,
        run_description = ""
    )

    # Return algo_params and applied_solvers
    return(algo_params = algo_params, applied_solvers = applied_solvers)
end


function model_definition(fixed_x::Float64, problem_params::DynamicSDDiP.ProblemParams)

    """
    2nd illustrative example from the paper on Lipschitz regularization (Füllner, Sun & Rebennack).
    """

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # State variables
        JuMP.@variable(subproblem, 0.0 <= x <= 7.0, SDDP.State, Int, initial_value = 0)

        if t == 1

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, x.out == fixed_x)

            # Stage objective
            SDDP.@stageobjective(subproblem, 1)

        else
            # Local variables
            JuMP.@variable(subproblem, 0.0 <= y <= 4.0, Int)
            JuMP.@variable(subproblem, u, Bin)
            JuMP.@variable(subproblem, z, Int)

            # Constraints
            x = subproblem[:x]
            JuMP.@constraint(subproblem, con1, y >= 1/4 * z)
            JuMP.@constraint(subproblem, con2, y >= 1/2 * z - 1)
            JuMP.@constraint(subproblem, con3, y >= z - 4)
            JuMP.@constraint(subproblem, con4, u >= 1/8 * z + 1/8)
            JuMP.@constraint(subproblem, eq, z == x.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, y + u)

            Infiltrator.@infiltrate
        end
            
    end

    return model
end


function model_set_up(
    number_of_stages::Int,
    fixed_x::Float64;
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
)

    ############################################################################
    # DEFINE PROBLEM PARAMS
    ############################################################################
    problem_params = DynamicSDDiP.ProblemParams(number_of_stages, 1, tree_seed = 12345)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_definition(fixed_x, problem_params)

    return (model = model, problem_params = problem_params)
end


function model_starter(
    fixed_x::Float64,
    core_point_approach::Symbol,
    copy_regime::DynamicSDDiP.AbstractCopyRegime,
    integer_regime::DynamicSDDiP.AbstractIntegerRegime,
    improvement_regime::DynamicSDDiP.AbstractImprovementRegime,
    unbounded_regime::DynamicSDDiP.AbstractUnboundedRegime,
    )
    try
        ########################################################################
        # DEFINE ALGO PARAMS
        ########################################################################
        algo_config_output = algo_config(core_point_approach, copy_regime, integer_regime, improvement_regime, unbounded_regime)
        algo_params = algo_config_output.algo_params
        applied_solvers = algo_config_output.applied_solvers

        ########################################################################
        # DEFINE MODEL
        ########################################################################
        model_output = model_set_up(2, fixed_x, algo_params=algo_params, applied_solvers=applied_solvers)
        model = model_output.model
        problem_params = model_output.problem_params

        ########################################################################
        # SOLVE (TRAIN) MODEL
        ########################################################################
        Random.seed!(12345)
        DynamicSDDiP.solve(model, algo_params, applied_solvers, problem_params)

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

"""
model_starter(fixed_x, core_point_approach, copy_regime, integer_regime, improvement_regime, unbounded_regime)

Note that norm_reg_lifted will always be transformed to a weighted norm as stated in the paper. Otherwise, this should be changed manually. 
TODO: Maybe change this in the future.

model_starter(3.0, :midpoint, DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.NoIntegerRelax(), DynamicSDDiP.NoImprovement(), DynamicSDDiP.Unbounded_Opt_Bound(strict=true, user_dual_multiplier_bound=10.0))



"""

model_starter(0.0, :relint, DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.NoIntegerRelax(), DynamicSDDiP.NoImprovement(), DynamicSDDiP.Unbounded_Opt_Bound(strict=true, user_dual_multiplier_bound=3.0))
