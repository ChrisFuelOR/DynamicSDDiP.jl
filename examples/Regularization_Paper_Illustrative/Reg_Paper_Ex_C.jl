import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
using Printf
using Gurobi
import Random


function algo_config(
    reg::Bool, 
    sigma::Float64, 
    norm_reg::Union{Nothing, DynamicSDDiP.AbstractNorm}, 
    norm_reg_lifted::Union{Nothing, DynamicSDDiP.AbstractNorm}, 
    copy_regime_reg::DynamicSDDiP.AbstractCopyRegime,
    copy_regime_bw::DynamicSDDiP.AbstractCopyRegime,
    binarization::Bool,
    K::Int64,
    augmented::Bool
    )

    # Stopping rules to be used
    #stopping_rules = [SDDP.TimeLimit(time_limit), SDDP.IterationLimit(1000), SDDP.BoundStalling(20,1e-4), DynamicSDDiP.DeterministicStopping()]
    stopping_rules = [SDDP.IterationLimit(1), DynamicSDDiP.DeterministicStopping()]
    
    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.LevelBundle()
    dual_bound_regime = DynamicSDDiP.ValueBound()
    dual_status_regime = DynamicSDDiP.Lax() #Rigorous
    dual_choice_regime = DynamicSDDiP.MinimalNormChoice()

    duality_regime = DynamicSDDiP.LagrangianDuality(
        atol = 1e-4,
        rtol = 1e-4,
        iteration_limit = 1000,
        dual_initialization_regime = dual_initialization_regime,
        dual_bound_regime = dual_bound_regime,
        dual_solution_regime = dual_solution_regime,
        dual_choice_regime = dual_choice_regime,
        dual_status_regime = dual_status_regime,
        copy_regime = copy_regime_bw,
        augmented = augmented,
    )

    # State approximation and cut projection configuration
    if binarization
        cut_projection_regime = DynamicSDDiP.BigM()
        binary_precision = Dict{Symbol, Float64}()
        binary_precision[:b] = 2/(2^K-1)
        state_approximation_regime = DynamicSDDiP.BinaryApproximation(binary_precision = binary_precision, cut_projection_regime = cut_projection_regime)


    else
        state_approximation_regime = DynamicSDDiP.NoStateApproximation()
    end

    # Cut generation regimes
    cut_generation_regime = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime,
        duality_regime = duality_regime,
    )

    cut_generation_regimes = [cut_generation_regime]

    # Regularization configuration
    if reg
        regularization_regime = DynamicSDDiP.Regularization(sigma=[0.0,sigma], sigma_factor=5.0, norm = norm_reg, norm_lifted = norm_reg_lifted, copy_regime = copy_regime_reg)
    else
        regularization_regime = DynamicSDDiP.NoRegularization()
    end

    # Cut aggregation regime
    cut_type = SDDP.SINGLE_CUT
    cut_aggregation_regime = DynamicSDDiP.SingleCutRegime()

    # Cut selection configuration
    cut_selection_regime = DynamicSDDiP.NoCutSelection()

    # Simulation
    simulation_regime = DynamicSDDiP.NoSimulation()

    # File for logging
    log_file = "C:/Users/cg4102/Documents/julia_logs/Reg_Paper_Ex_C.log"

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
        JuMP.@variable(subproblem, 0.0 <= b <= 2.0, SDDP.State, initial_value = 0)

        if t == 1

            # Constraints
            b = subproblem[:b]
            JuMP.@constraint(subproblem, b.out == fixed_x + b.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, 1)

        else
            # Local variables
            JuMP.@variable(subproblem, 0.0 <= x[i=1:4])
            JuMP.set_integer(x[1])
            JuMP.set_integer(x[2])

            # Constraints
            b = subproblem[:b]
            JuMP.@constraint(subproblem, b.out == 0)
            JuMP.@constraint(subproblem, con, 1.25*x[1] - x[2] + 0.5*x[3] + 1/3*x[4] == b.in)

            # Stage objective
            SDDP.@stageobjective(subproblem, x[1] - 0.75*x[2] + 0.75*x[3] + 2.5*x[4])

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
    reg::Bool, 
    sigma::Float64, 
    norm_reg::Union{Nothing, DynamicSDDiP.AbstractNorm}, 
    norm_reg_lifted::Union{Nothing, DynamicSDDiP.AbstractNorm}, 
    copy_regime_reg::DynamicSDDiP.AbstractCopyRegime,
    copy_regime_bw::DynamicSDDiP.AbstractCopyRegime,
    binarization::Bool,
    K::Int64,
    augmented::Bool
    )
    try
        ########################################################################
        # DEFINE ALGO PARAMS
        ########################################################################
        algo_config_output = algo_config(reg, sigma, norm_reg, norm_reg_lifted, copy_regime_reg, copy_regime_bw, binarization, K, augmented)
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
model_starter(fixed_x, reg, sigma/rho, norm_reg, norm_reg_lifted, copy_regime_reg, copy_regime_bw, binarization, K, augmented)

Note that norm_reg_lifted will always be transformed to a weighted norm as stated in the paper. Otherwise, this should be changed manually. 
TODO: Maybe change this in the future.

LINEAR CUT - TESTS
model_starter(6/5, false, Inf, nothing, nothing, DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.StateSpaceCopy(), false, 0, false)                             CORRECT
model_starter(6/5, true, 1.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.StateSpaceCopy(), false, 0, false)          CORRECT
model_starter(6/5, true, 1/2, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.StateSpaceCopy(), false, 0, false)          CORRECT

CPC - TEST FOR BINARY REFINEMENTS:    
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 2, false)           CORRECT [v4]      
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 3, false)           CORRECT [v13]
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 4, false)           CORRECT [v6]

CPC - TEST FOR SIGMA REFINEMENTS:    
model_starter(6/5, true, 4/5, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 4, false)           CORRECT [v14]
model_starter(6/5, true, 1.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 4, false)           CORRECT [v5]
model_starter(6/5, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 4, false)           CORRECT [v7]
model_starter(6/5, true, 10.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 4, false)          [v18]

CPC - TEST FOR SUPREMUM NORM IN LIFTED SPACE: 
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L∞(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 3, false)           CORRECT [v15]
model_starter(6/5, true, 1.0, DynamicSDDiP.L₁(), DynamicSDDiP.L∞(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 3, false)           CORRECT

CPC - TEST FOR Z = {0,1} IN LIFTED SPACE: 
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.StateSpaceCopy(), true, 3, false)           DIFFERENT [v1]

CPC - TEST FOR BOUNDED SLOPE: 
model_starter(1.249, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)         DIFFERENT [v16_new]
model_starter(1.249, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)         DIFFERENT [V17_new]
=> manually change to non-weighted norm in lifted space (Reg + Bw), TODO: automize this

CPC - TEST FOR ALD CUT COMPARISON:
model_starter(1.249, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)         CORRECT [v16*]
model_starter(1.249, true, 10.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)        CORRECT [v16*]
model_starter(1.249, true, 20.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)        CORRECT [v16*]

model_starter(1.249, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)         CORRECT [v16*]
model_starter(1.249, true, 10.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)        CORRECT [v16*]
model_starter(1.249, true, 20.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)        CORRECT [v16*]

model_starter(6/5, true, 4/5, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)           CORRECT [v14*]
model_starter(6/5, true, 1.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)           CORRECT [v5*]
model_starter(6/5, true, 2.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)           CORRECT [v4*]
model_starter(6/5, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)           CORRECT [v7*]
model_starter(6/5, true, 10.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), false, 0, true)          [v18*]

# All ALD runs with only ValueBound.
# Second batch of tests with backward regularization (sigma = rho used in augmented Lagrangian dual), but with using the primal_obj bound from the unregularized case.  
# Note that for the augmented Lagrangian dual, we always use the L₁-norm.
"""

model_starter(1.249, true, 5.0, DynamicSDDiP.L₁(), DynamicSDDiP.L₁(), DynamicSDDiP.StateSpaceCopy(), DynamicSDDiP.ConvexHullCopy(), true, 8, false)   # ReLU comparison

