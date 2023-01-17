import DynamicSDDiP
import Infiltrator
import MathOptInterface
using Revise
using Printf

include("algo_config.jl")
include("scenario_tree.jl")
include("model.jl")
include("model_no_bin.jl")
include("simulation.jl")

function model_starter(
    num::Int,
    number_of_stages::Int,
    number_of_realizations::Int,
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    forward_seed::Int,
    tree_seed::Int
)

    try
        ########################################################################
        # DEFINE ALGO PARAMS
        ########################################################################
        algo_config_output = algo_config(duality_regime_sym, normalization_regime, cut_aggregation_regime, cut_selection_regime, log_file, time_limit, forward_seed)
        algo_params = algo_config_output.algo_params
        applied_solvers = algo_config_output.applied_solvers

        ########################################################################
        # DEFINE MODEL
        ########################################################################
        model_output = model_no_bin_set_up(number_of_stages, number_of_realizations, algo_params=algo_params, applied_solvers=applied_solvers, tree_seed=tree_seed)
        model = model_output.model
        problem_params = model_output.problem_params

        ########################################################################
        # SOLVE (TRAIN) MODEL
        ########################################################################
        Random.seed!(forward_seed)
        DynamicSDDiP.solve(model, algo_params, applied_solvers, problem_params)

        ########################################################################
        # SIMULATE MODEL
        ########################################################################
        simulate(model, algo_params, problem_params, algo_params.simulation_regime)

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


function det_equiv_no_bin_starter(
    num::Int,
    number_of_stages::Int,
    number_of_realizations::Int,
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    forward_seed::Int,
    tree_seed::Int
)

    try
        ############################################################################
        # DEFINE ALGO PARAMS
        ############################################################################
        algo_config_output = algo_config(duality_regime_sym, normalization_regime, cut_aggregation_regime, cut_selection_regime, log_file, time_limit, forward_seed)
        algo_params = algo_config_output.algo_params
        applied_solvers = algo_config_output.applied_solvers

        ############################################################################
        # DEFINE MODEL
        ############################################################################
        model_output = model_no_bin_set_up(number_of_stages, number_of_realizations, algo_params=algo_params, applied_solvers=applied_solvers, tree_seed=tree_seed)
        model = model_output.model
        problem_params = model_output.problem_params

        ############################################################################
        # SOLVE (TRAIN) MODEL
        ############################################################################
        det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer, time_limit = 7200.0)
        JuMP.set_objective_sense(det_equiv, MathOptInterface.MIN_SENSE)
        JuMP.optimize!(det_equiv)
        print(JuMP.objective_value(det_equiv))

        ############################################################################
        # LOGGING
        ############################################################################
        log_file_handle = open(algo_params_SB_tol.log_file, "a")
        DynamicSDDiP.print_helper(DynamicSDDiP.print_det_equiv, log_file_handle, problem_params, JuMP.objective_value(det_equiv))
        close(log_file_handle)

    catch e
        @printf "Case %d (deterministic equivalent) terminated with error" num
        println()
        #throw(error(e))
        showerror(stdout, e, catch_backtrace())
        println()
        println("#############################################################")
        println()
    end

end

function model_starter_runs()

    """
    Specification of model runs that should be run one after the other.
    """

    #det_equiv_no_bin_starter(0,14,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_SB_tol.log", 7200, 11111, 12345)

    #Dummy run
    #model_starter(0,14,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_SB_tol.log", 100, 11111, 12345)

    model_starter(1,14,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Single_SB_tol.log", 10800, 11111, 12345)
    model_starter(2,14,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Multi_SB_tol.log", 10800, 11111, 12345)
    model_starter(3,14,10,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Mid_SB_tol.log", 10800, 11111, 12345)
    model_starter(4,14,10,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_InOut_SB_tol.log", 10800, 11111, 12345)
    model_starter(5,14,10,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Relint_SB_tol.log", 10800, 11111, 12345)
    model_starter(6,14,10,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Eps_SB_tol.log", 10800, 11111, 12345)
    model_starter(7,14,10,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_Lsup_SB_tol.log", 10800, 11111, 12345)
    model_starter(8,14,10,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_L1_SB_tol.log", 10800, 11111, 12345)
    model_starter(9,14,10,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_L1sup_SB_tol.log", 10800, 11111, 12345)

    #model_starter(10,14,10,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_B_tol.log", 10800, 11111, 12345)
    #model_starter(11,14,10,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_SB_tol.log", 10800, 11111, 12345)
    #model_starter(12,14,10,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_B_Multi_tol.log", 10800, 11111, 12345)
    #model_starter(13,14,10,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_SB_Multi_tol.log", 10800, 11111, 12345)

    #model_starter(14,5,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_L1_SB_tol.log", 10800, 11111, 12345)
    #model_starter(15,5,50,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_Span_SB_tol.log", 10800, 11111, 12345)
    #model_starter(16,5,50,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_Mid_SB_tol.log", 10800, 11111, 12345)
    #model_starter(17,5,50,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_Eps_SB_tol.log", 10800, 11111, 12345)
    #model_starter(18,5,50,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_Relint_SB_tol.log", 10800, 11111, 12345)
    #model_starter(19,5,50,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/AirRev_O_14_10_CL_Lsup_SB_tol.log", 10800, 11111, 12345)

end

#model_starter_runs()