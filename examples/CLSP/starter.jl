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
        #simulate(model, algo_params, problem_params, DynamicSDDiP.HistoricalSample())

        ########################################################################
        # DETERMINISTIC EQUIVALENT (INCLUDING CUTS!!!)
        ########################################################################
        #det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer, time_limit = 7200.0)
        #JuMP.set_objective_sense(det_equiv, MathOptInterface.MIN_SENSE)
        #JuMP.optimize!(det_equiv)
        #print(JuMP.objective_value(det_equiv))

        #log_file_handle = open(algo_params_CL3.log_file, "a")
        #DynamicSDDiP.print_helper(DynamicSDDiP.print_det_equiv, log_file_handle, problem_params, JuMP.objective_value(det_equiv))
        #close(log_file_handle)

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
        JuMP.set_optimizer_attribute(det_equiv, "TimeLimit", 7200.0)
        JuMP.optimize!(det_equiv)
        print(JuMP.objective_value(det_equiv))

        ############################################################################
        # LOGGING
        ############################################################################
        log_file_handle = open(algo_params_CL3.log_file, "a")
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

    #model_starter(14,2,10,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_2_1_Eps.log", 1800, 11111, 12345)

    model_starter(13,4,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_B_Multi.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_SB_Multi.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_Lsup.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_L1sup.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_4_20_Eps.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_4_20_Relint.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_4_20_InOut.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_4_20_Mid.log", 10800, 11111, 12345)

    #model_starter(2,6,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_6_20_Multi.log", 14400, 11111, 12345)
    #model_starter(4,6,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_6_20_Mid.log", 14400, 11111, 12345)
    #model_starter(4,6,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_6_20_InOut.log", 14400, 11111, 12345)
    #model_starter(5,6,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_6_20_Relint.log", 14400, 11111, 12345)
    #model_starter(6,6,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_6_20_Eps.log", 14400, 11111, 12345)
    #model_starter(7,6,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_6_20_Lsup.log", 14400, 11111, 12345)
    #model_starter(9,6,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_6_20_L1sup.log", 14400, 11111, 12345)
    #model_starter(11,6,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_6_20_B_Multi.log", 14400, 11111, 12345)
    #model_starter(12,6,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_6_20_SB_Multi.log", 14400, 11111, 12345)

    #model_starter(2,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_SSC_16_20_Single_SB.log", 28800, 11111, 12345)
    #model_starter(2,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_SSC_16_20_Multi_SB.log", 28800, 11111, 12345)
    #model_starter(4,10,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_Mid.log", 18000, 11111, 12345)
    #model_starter(4,10,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_InOut.log", 18000, 11111, 12345)
    #model_starter(5,10,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_Relint.log", 18000, 11111, 12345)
    #model_starter(10,10,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_Eps_SB.log", 18000, 11111, 12345)
    #model_starter(9,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_SSC_16_20_L1_SB.log", 28800, 11111, 12345)
    #model_starter(7,16,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_SSC_16_20_Lsup_SB.log", 28800, 11111, 12345)
    #model_starter(9,16,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_SSC_16_20_L1sup_SB.log", 28800, 11111, 12345)
    #model_starter(11,10,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_B_Multi.log", 18000, 11111, 12345)
    #model_starter(12,10,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_SB_Multi.log", 18000, 11111, 12345)
    #model_starter(11,10,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_B.log", 18000, 11111, 12345)
    #model_starter(12,10,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_SB.log", 18000, 11111, 12345)

    #model_starter(2,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_Single_SB.log", 28800, 11111, 12345)
    #model_starter(2,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_Multi_SB.log", 28800, 11111, 12345)
    #model_starter(4,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_Mid_SB.log", 28800, 11111, 12345)
    #model_starter(4,16,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_InOut.log", 28800, 11111, 12345)
    #model_starter(5,16,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_Relint_SB.log", 28800, 11111, 12345)
    #model_starter(10,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_Eps_SB.log", 28800, 11111, 12345)
    #model_starter(9,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_L1_SB.log", 28800, 11111, 12345)
    #model_starter(7,16,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_Lsup_SB.log", 28800, 11111, 12345)
    #model_starter(9,16,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_L1sup_SB.log", 28800, 11111, 12345)
    #model_starter(11,16,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_B_Multi.log", 28800, 11111, 12345)
    #model_starter(12,16,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_SB_Multi.log", 28800, 11111, 12345)
    #model_starter(11,16,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_B.log", 28800, 11111, 12345)
    #model_starter(12,16,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_SB.log", 28800, 11111, 12345)

    #model_starter(9,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_L1_Lsup.log", 28800, 11111, 12345)

    #model_starter(14,10,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_CL_L1.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_CL_Span.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.Core_Midpoint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_10_20_CL_Mid.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_10_20_CL_Eps.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_10_20_CL_Relint.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_10_20_CL_Lsup.log", 18000, 11111, 12345)
    #model_starter(15,10,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_10_20_CL_Lsup_no_SB.log", 18000, 11111, 12345)

    #model_starter(14,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_CL_L1.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_CL_Span.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_CL_Mid.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_CL_Eps.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_CL_Relint.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_16_20_CL_Lsup.log", 28800, 11111, 12345)
    #model_starter(15,16,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_16_20_CL_Lsup_no_SB.log", 28800, 11111, 12345)



end

#model_starter_runs()
