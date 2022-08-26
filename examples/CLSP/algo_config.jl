import SDDP
import DynamicSDDiP
import GAMS
import Gurobi
import Infiltrator
using Revise

function algo_config(
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    forward_seed::Int,
    )

    # Stopping rules to be used
    #stopping_rules = [SDDP.TimeLimit(time_limit), SDDP.IterationLimit(1000), SDDP.BoundStalling(20,1e-4)]
    stopping_rules = [SDDP.TimeLimit(time_limit), SDDP.BoundStalling(20,1e-4)]
    #stopping_rules = [ SDDP.IterationLimit(20), SDDP.BoundStalling(20,1e-4)]

    # Duality / Cut computation configuration
    dual_initialization_regime = DynamicSDDiP.ZeroDuals()
    dual_solution_regime = DynamicSDDiP.Kelley()
    dual_bound_regime = DynamicSDDiP.BothBounds()
    dual_status_regime = DynamicSDDiP.Lax()

    dual_choice_regime = DynamicSDDiP.StandardChoice()
    if isa(normalization_regime, DynamicSDDiP.L∞_Deep)
        dual_choice_regime = DynamicSDDiP.MinimalNormChoice()
    end

    #dual_space_regime = DynamicSDDiP.BendersSpanSpaceRestriction(10, :multi_cut)
    dual_space_regime = DynamicSDDiP.NoDualSpaceRestriction()
    copy_regime = DynamicSDDiP.ConvexHullCopy()

    if duality_regime_sym == :uni_lag
        duality_regime = DynamicSDDiP.UnifiedLagrangianDuality(
            atol = 1e-4,
            rtol = 1e-4,
            iteration_limit = 1000,
            dual_initialization_regime = dual_initialization_regime,
            dual_bound_regime = dual_bound_regime,
            dual_solution_regime = dual_solution_regime,
            dual_choice_regime = dual_choice_regime,
            dual_status_regime = dual_status_regime,
            normalization_regime = normalization_regime,
            dual_space_regime = dual_space_regime,
            copy_regime = copy_regime,
            #user_dual_multiplier_bound = 10.0, # 10.0
            user_dual_objective_bound = 1e4,
        )
    elseif duality_regime_sym == :lag
        duality_regime = DynamicSDDiP.LagrangianDuality(
            atol = 1e-4,
            rtol = 1e-4,
            iteration_limit = 1000,
            dual_initialization_regime = dual_initialization_regime,
            dual_bound_regime = dual_bound_regime,
            dual_solution_regime = dual_solution_regime,
            dual_choice_regime = dual_choice_regime,
            dual_status_regime = dual_status_regime,
            copy_regime = copy_regime,
        )
    elseif duality_regime_sym == :SB
        duality_regime = DynamicSDDiP.StrengthenedDuality()
    elseif duality_regime_sym == :B
        duality_regime = DynamicSDDiP.LinearDuality()
    end

    # State approximation and cut projection configuration
    state_approximation_regime = DynamicSDDiP.NoStateApproximation()

    # Cut generation regimes
    cut_generation_regime_2 = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime,
        duality_regime = duality_regime,
        #cut_away_approach = true,
        #iteration_to_start = 31,
    )

    cut_generation_regime_1 = DynamicSDDiP.CutGenerationRegime(
        state_approximation_regime = state_approximation_regime,
        duality_regime = DynamicSDDiP.StrengthenedDuality(),
    )

    cut_generation_regimes = [cut_generation_regime_2]

    # Regularization configuration
    regularization_regime = DynamicSDDiP.NoRegularization()

    # Cut aggregation regime
    cut_type = SDDP.SINGLE_CUT
    if isa(cut_aggregation_regime, DynamicSDDiP.MultiCutRegime)
        cut_type = SDDP.MULTI_CUT
    end

    # Simulation regime
    #simulation_regime = DynamicSDDiP.Simulation(sampling_scheme=DynamicSDDiP.OutOfSampleMonteCarlo(number_of_realizations=10,simulation_seed=232323),number_of_replications=1000)
    simulation_regime = DynamicSDDiP.Simulation(sampling_scheme=DynamicSDDiP.InSampleMonteCarlo(),number_of_replications=1000)
    #simulation_regime = DynamicSDDiP.NoSimulation()

    # Suppress solver output
    silent = true

    # Infiltration for debugging
    infiltrate_state = :none

    # Solver approach
    solver_approach = DynamicSDDiP.Direct_Solver()

    # Define solvers to be used
    applied_solvers = DynamicSDDiP.AppliedSolvers(
        LP = "Gurobi",
        MILP = "Gurobi",
        MIQCP = "Gurobi",
        MINLP = "Gurobi",
        NLP = "Gurobi",
        Lagrange = "Gurobi",
    )

    # Definition of algo_params
    algo_params = DynamicSDDiP.AlgoParams(
        stopping_rules = stopping_rules,
        regularization_regime = regularization_regime,
        cut_aggregation_regime = cut_aggregation_regime,
        cut_selection_regime = cut_selection_regime,
        cut_generation_regimes = cut_generation_regimes,
        simulation_regime = simulation_regime,
        cut_type = cut_type,
        log_file = log_file,
        silent = silent,
        infiltrate_state = infiltrate_state,
        solver_approach = solver_approach,
        numerical_focus = false,
        run_numerical_stability_report = false,
        seed = forward_seed,
        run_description = ""
    )

    # Return algo_params and applied_solvers
    return(algo_params = algo_params, applied_solvers = applied_solvers)
end
