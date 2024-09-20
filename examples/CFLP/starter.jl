import DynamicSDDiP
import Infiltrator
import MathOptInterface
import Random
using Revise
using Printf

include("algo_config.jl")
include("scenario_tree.jl")
include("model.jl")
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

    #det_equiv_starter(0,5,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.CutSelection(), "C:/Users/cg4162/Documents/julia_logs/CFLP_5_10_det_equiv.log", 3600, 11111, 12345)

    #model_starter(1,100,20,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_B.log", 3600, 11111, 12345)
    #model_starter(2,5,10,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_SB.log", 3600, 11111, 12345)      
    #model_starter(3,5,10,:B, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_B_Multi.log", 3600, 11111, 12345)
    #model_starter(4,5,10,:SB, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_SB_Multi.log", 3600, 11111, 12345)      
    #model_starter(5,100,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_Single.log", 3600, 11111, 12345)
    #model_starter(6,5,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_Multi_SB.log", 3600, 11111, 12345)
  
    #model_starter(7,5,10,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_L1_SB.log", 3600, 11111, 12345)
    #model_starter(8,5,10,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_Lsup_SB.log", 3600, 11111, 12345)
    #model_starter(9,5,10,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_5_10_L1sup_SB.log", 3600, 11111, 12345)
    model_starter(10,100,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_regime=DynamicSDDiP.IntegerRelax(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Eps_N_sbopt.log", 14400, 11111, 12345)
    model_starter(11,100,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Relint_N_sbopt.log", 14400, 11111, 12345)
    #model_starter(12,100,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=true,normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_InOut.log", 7200, 11111, 12345)
    #model_starter(11,100,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_Relint_SB_NoNorm.log", 14400, 11111, 12345)
    #model_starter(11,100,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_Relint_SB_NoNorm.log", 14400, 11111, 12345)
    #model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=true,normalize_direction=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_Mid_SB_NoNorm.log", 14400, 11111, 12345)
    model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_regime=DynamicSDDiP.IntegerRelax(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Mid_N_sbopt.log", 14400, 11111, 12345)
    
    model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.5,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Conv_50_N_sbopt.log", 14400, 11111, 12345)
    model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.75,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Conv_75_N_sbopt.log", 14400, 11111, 12345)
    model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.9,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Conv_90_N_sbopt.log", 14400, 11111, 12345)
    model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Conv(lambda=0.99,copy_regime=DynamicSDDiP.StateSpaceCopy(),normalize_direction=true,unbounded_regime=DynamicSDDiP.Unbounded_Opt_SB()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_S_Conv_99_N_sbopt.log", 14400, 11111, 12345)
    #model_starter(13,100,20,:uni_lag, DynamicSDDiP.Core_Relint(copy_regime=DynamicSDDiP.StateSpaceCopy(), integer_regime=DynamicSDDiP.IntegerRelax()), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_Core_Relint_SB.log", 14400, 11111, 12345)
 
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_L1.log", 7200, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.ChenLuedtke(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_Span.log", 7200, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.Core_Midpoint(integer_relax=true,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_Mid.log", 7200, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2,integer_relax=true,normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_Epsilon.log", 7200, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.Core_Relint(normalize_direction=false), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_Relint.log", 7200, 11111, 12345)
    #model_starter(15,100,20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/CFLP_New_CL_Lsup.log", 7200, 11111, 12345)
    

end

#model_starter_runs()
