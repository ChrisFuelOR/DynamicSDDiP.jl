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
        model_output = model_set_up(number_of_stages, number_of_realizations, algo_params=algo_params, applied_solvers=applied_solvers, tree_seed=tree_seed)
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

        ########################################################################
        # SIMULATE MODEL USING FULL SCENARIO TREE
        ########################################################################
        simulate(model, algo_params, problem_params, DynamicSDDiP.HistoricalSample())

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


function det_equiv_starter(
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
        model_output = model_set_up(number_of_stages, number_of_realizations, algo_params=algo_params, applied_solvers=applied_solvers, tree_seed=tree_seed)
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
        log_file_handle = open(algo_params.log_file, "a")
        DynamicSDDiP.print_helper(DynamicSDDiP.print_det_equiv, log_file_handle, problem_params, JuMP.objective_value(det_equiv), JuMP.objective_bound(det_equiv))
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

    #det_equiv_starter(0,2,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,2,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_Single.log", 3600, 11111, 12345)
    model_starter(3,2,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_Mid.log", 3600, 11111, 12345)
    model_starter(4,2,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_L1.log", 3600, 11111, 12345)
    model_starter(5,2,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_B.log", 3600, 11111, 12345)
    model_starter(6,2,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_2_20_SB.log", 3600, 11111, 12345)

    det_equiv_starter(0,3,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,3,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_Single.log", 3600, 11111, 12345)
    model_starter(3,3,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_Mid.log", 3600, 11111, 12345)
    model_starter(4,3,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_L1.log", 3600, 11111, 12345)
    model_starter(5,3,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_B.log", 3600, 11111, 12345)
    model_starter(6,3,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_3_20_SB.log", 3600, 11111, 12345)

    det_equiv_starter(0,4,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,4,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_Single.log", 3600, 11111, 12345)
    model_starter(3,4,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_Mid.log", 3600, 11111, 12345)
    model_starter(4,4,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_L1.log", 3600, 11111, 12345)
    model_starter(5,4,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_B.log", 3600, 11111, 12345)
    model_starter(6,4,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_4_20_SB.log", 3600, 11111, 12345)


    #det_equiv_starter(0,10,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_binary_single.log", 7200, 11111, 12345)

    # Dummy run
    # model_starter(0,10,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_single_10.log", 600, 11111, 12345)

    #model_starter(9,10,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_single_10.log", 18000, 11111, 12345)
    #model_starter(0,10,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_multi_10.log", 18000, 11111, 12345)
    #model_starter(1,10,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_Lsup_10.log", 18000, 11111, 12345)
    #model_starter(3,10,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_mid_10.log", 18000, 11111, 12345)
    #model_starter(2,10,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_binary_in_out_10.log", 18000, 11111, 12345)
    #model_starter(4,10,20,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_relint_10.log", 18000, 11111, 12345)
    #model_starter(5,10,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_eps_10.log", 18000, 11111, 12345)
    #model_starter(6,10,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_L1_10.log", 18000, 11111, 12345)
    #model_starter(7,10,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_L1sup_10.log", 18000, 11111, 12345)
    #model_starter(8,10,20,:uni_lag, DynamicSDDiP.Core_Optimal(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_binary_opt_10.log", 18000, 11111, 12345)

    #model_starter(10,10,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_Benders_10.log", 18000, 11111, 12345)
    #model_starter(11,10,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_Benders_10.log", 18000, 11111, 12345)
    #model_starter(12,10,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_Benders_10.log", 18000, 11111, 12345)
    #model_starter(13,10,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_Benders_10.log", 18000, 11111, 12345)

    # model_starter(14,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_L1_v2.log", 18000, 11111, 12345)
    # model_starter(15,16,20,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_Span_v2.log", 18000, 11111, 12345)
    # model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_Mid_v2.log", 18000, 11111, 12345)
    # model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_Eps_v2.log", 18000, 11111, 12345)
    # model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_Relint_v2.log", 18000, 11111, 12345)
    # model_starter(15,16,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_16_20_CL_Lsup_v2.log", 18000, 11111, 12345)
    #

end

#model_starter_runs()
