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
        det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer, time_limit = 7500.0)
        JuMP.set_objective_sense(det_equiv, MathOptInterface.MIN_SENSE)
        JuMP.set_optimizer_attribute(det_equiv, "TimeLimit", 7200.0)
        JuMP.optimize!(det_equiv)
        print(JuMP.objective_value(det_equiv))

        ############################################################################
        # LOGGING
        ############################################################################
        log_file_handle = open(algo_params.log_file, "a")
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

    #det_equiv_no_bin_starter(0,10,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_single_16_new_no_bin.log", 2400, 11111, 12345)

    det_equiv_no_bin_starter(0,2,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,2,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_Single.log", 3600, 11111, 12345)
    model_starter(3,2,50,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_Mid.log", 3600, 11111, 12345)
    model_starter(4,2,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_L1.log", 3600, 11111, 12345)
    model_starter(5,2,50,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_B.log", 3600, 11111, 12345)
    model_starter(6,2,50,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_2_50_SB.log", 3600, 11111, 12345)

    det_equiv_no_bin_starter(0,3,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,3,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_Single.log", 3600, 11111, 12345)
    model_starter(3,3,50,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_Mid.log", 3600, 11111, 12345)
    model_starter(4,3,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_L1.log", 3600, 11111, 12345)
    model_starter(5,3,50,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_B.log", 3600, 11111, 12345)
    model_starter(6,3,50,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_3_50_SB.log", 3600, 11111, 12345)

    det_equiv_no_bin_starter(0,4,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_det_equiv.log", 3600, 11111, 12345)
    model_starter(1,4,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_Single.log", 3600, 11111, 12345)
    model_starter(3,4,50,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_Mid.log", 3600, 11111, 12345)
    model_starter(4,4,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_L1.log", 3600, 11111, 12345)
    model_starter(5,4,50,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_B.log", 3600, 11111, 12345)
    model_starter(6,4,50,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_O_4_50_SB.log", 3600, 11111, 12345)

    # Dummy run
    # model_starter(0,16,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/clsp_large_binary_single_16_new_no_bin.log", 600, 11111, 12345)
    #
    #model_starter(9,16,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_Single_2h.log", 7200, 11111, 12345)
    #model_starter(0,16,50,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_Multi_2h.log", 7200, 11111, 12345)
    #model_starter(3,16,50,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_Mid_SB.log", 7200, 11111, 12345)
    #model_starter(2,16,50,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Large_O_16_50_InOut_MNC.log", 3600, 11111, 12345)
    #model_starter(4,16,50,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_Relint_SB.log", 7200, 11111, 12345)
    #model_starter(5,16,50,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_Eps_SB.log", 7200, 11111, 12345)
    #model_starter(6,16,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_L1_SB.log", 7200, 11111, 12345)
    #model_starter(7,16,50,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_L1sup_SB.log", 7200, 11111, 12345)
    #model_starter(1,16,50,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_16_50_prim_Lsup_SB.log", 7200, 11111, 12345)

    #model_starter(10,16,50,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_B_CS.log", 3600, 11111, 12345)
    #model_starter(11,16,50,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_CS.log", 3600, 11111, 12345)
    #model_starter(12,16,50,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_B_Multi_CS.log", 3600, 11111, 12345)
    #model_starter(13,16,50,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Multi_CS.log", 3600, 11111, 12345)

    #model_starter(14,16,50,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_L1_CL_CS.log", 3600, 11111, 12345)
    #model_starter(15,16,50,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Span_CL_CS.log", 3600, 11111, 12345)
    #model_starter(15,16,50,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Mid_CL_CS.log", 3600, 11111, 12345)
    #model_starter(15,16,50,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Eps_CL_CS.log", 3600, 11111, 12345)
    #model_starter(15,16,50,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Relint_CL_CS.log", 3600, 11111, 12345)
    #model_starter(15,16,50,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_Real_LB_O_16_50_CL_Lsup_CL_CS.log", 3600, 11111, 12345)


end

#model_starter_runs()
