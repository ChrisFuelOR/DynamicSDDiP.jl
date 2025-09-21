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
        #simulate(model, algo_params, problem_params, algo_params.simulation_regime)

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

    #model_starter(1,16,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_B_2025.log", 14400, 11111, 12345)
    #model_starter(2,16,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_SB_2025.log", 14400, 11111, 12345)      
    #model_starter(3,16,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_B_Multi_2025.log", 14400, 11111, 12345)
    #model_starter(4,16,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_SB_Multi_2025.log", 14400, 11111, 12345)      
    #model_starter(5,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Single_2025.log", 14400, 11111, 12345)
    #model_starter(6,16,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Multi_2025.log", 14400, 11111, 12345)
  
    #model_starter(7,16,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_L1_2025.log", 14400, 11111, 12345)
    model_starter(10,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,copy_regime=DynamicSDDiP.ConvexHullCopy(),integer_regime=DynamicSDDiP.NoIntegerRelax(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Eps_2025.log", 14400, 11111, 12345)
    model_starter(11,16,20,:uni_lag, DynamicSDDiP.Core_Relint(copy_regime=DynamicSDDiP.ConvexHullCopy(),integer_regime=DynamicSDDiP.NoIntegerRelax(),improvement_regime=DynamicSDDiP.PrimalObj(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Relint_2025.log", 14400, 11111, 12345)
    model_starter(12,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(copy_regime=DynamicSDDiP.ConvexHullCopy(),integer_regime=DynamicSDDiP.NoIntegerRelax(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Mid_2025.log", 14400, 11111, 12345)
    
    
    #model_starter(8,100,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_Lsup.log", 14400, 11111, 12345)
    #model_starter(9,100,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_L1sup.log", 14400, 11111, 12345)
    #model_starter(10,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_regime=DynamicSDDiP.IntegerRelax(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_Eps_N_sbopt.log", 14400, 11111, 12345)
    #model_starter(11,16,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_Relint_N_sbopt.log", 14400, 11111, 12345)
    #model_starter(12,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_regime=DynamicSDDiP.IntegerRelax(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/Test.log", 14400, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Span_CL.log", 14400, 11111, 12345)
    
     #model_starter(13,4,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.5,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/Test.log", 14400, 11111, 12345)
     #model_starter(13,4,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.75,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/Test.log", 14400, 11111, 12345)
    #model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.9,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_Conv_90_N.log", 14400, 11111, 12345)
    #model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.99,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_NS_Conv_99_N.log", 14400, 11111, 12345)
    #model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Relint(copy_regime=DynamicSDDiP.StateSpaceCopy(), integer_regime=DynamicSDDiP.IntegerRelax()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_Core_Relint_SB.log", 14400, 11111, 12345)
 

    #model_starter(14,2,10,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=false,normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_2_1_Eps.log", 1800, 11111, 12345)

    #model_starter(13,4,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_B_Multi.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_SB_Multi.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_Lsup.log", 10800, 11111, 12345)
    #model_starter(14,4,20,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR_4_20_L1sup.log", 10800, 11111, 12345)
    #model_starter(14,16,20,:uni_lag, DynamicSDDiP.Core_Epsilon(copy_regime=DynamicSDDiP.ConvexHullCopy(),perturb=1e-2,integer_regime=DynamicSDDiP.NoIntegerRelax(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_Eps_sbopt_SB.log", 28800, 11111, 12345)
    #model_starter(14,16,20,:uni_lag, DynamicSDDiP.Core_Relint(copy_regime=DynamicSDDiP.ConvexHullCopy(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_Relint_sbopt_SB.log", 28800, 11111, 12345)
    #model_starter(14,6,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=false,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR3_4_20_InOut.log", 10800, 11111, 12345)
    #model_starter(14,16,20,:uni_lag, DynamicSDDiP.Core_Midpoint(copy_regime=DynamicSDDiP.ConvexHullCopy(),integer_regime=DynamicSDDiP.NoIntegerRelax(),normalize_direction=false,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_RR2_16_20_Mid_sbopt_SB.log", 28800, 11111, 12345)

    #model_starter(13,16,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.5,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Conv_50_SB.log", 28800, 11111, 12345)
    #model_starter(13,16,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.75,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_16_20_Conv_75_SB.log", 28800, 11111, 12345)
    #model_starter(13,10,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.9,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_10_20_Conv_90_SB.log", 18000, 11111, 12345)
    #model_starter(13,10,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.99,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CLSP_10_20_Conv_99_SB.log", 18000, 11111, 12345)   
end

#model_starter_runs()
