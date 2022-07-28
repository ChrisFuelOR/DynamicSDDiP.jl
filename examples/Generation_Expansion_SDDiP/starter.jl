import DynamicSDDiP
import Infiltrator
import MathOptInterface
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
    det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer)
    JuMP.set_objective_sense(det_equiv, MathOptInterface.MIN_SENSE)
    JuMP.optimize!(det_equiv)
    print(JuMP.objective_value(det_equiv))

    ############################################################################
    # LOGGING
    ############################################################################
    log_file_handle = open(algo_params.log_file, "a")
    DynamicSDDiP.print_helper(DynamicSDDiP.print_det_equiv, log_file_handle, problem_params, JuMP.objective_value(det_equiv))
    close(log_file_handle)

end


function model_starter_runs()

    """
    Specification of model runs that should be run one after the other.
    """

    #det_equiv_starter(0,3,20,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/GenExp.log", 50, 11111, 12345)

    model_starter(0,3,20,:uni_lag, DynamicSDDiP.Core_In_Out(integer_relax=true), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/GenExp.log", 600, 11111, 12345)

    # model_starter(1,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11111)
    # model_starter(2,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11111)
    # model_starter(3,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11111)
    # model_starter(4,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11111)
    # model_starter(5,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11111)
    # model_starter(6,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11111)
    # model_starter(7,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11111)
    # model_starter(8,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11111)

    # model_starter(9,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11112)
    # model_starter(10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11112)
    # model_starter(11,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11112)
    # model_starter(12,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11112)
    # model_starter(13,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11112)
    # model_starter(14,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11112)
    # model_starter(15,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11112)
    # model_starter(16,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11112)
    #
    # model_starter(17,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 3600, 11113)
    # model_starter(18,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11113)
    # model_starter(19,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11113)
    # model_starter(20,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11113)
    # model_starter(21,:uni_lag, DynamicSDDiP.L₁∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1sup.log", 3600, 11113)
    # model_starter(22,:uni_lag, DynamicSDDiP.Core_In_Out(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_in_out.log", 3600, 11113)
    # model_starter(23,:uni_lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_mid.log", 3600, 11113)
    # model_starter(24,:uni_lag, DynamicSDDiP.Core_Epsilon(perturb=1e-2), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_eps.log", 3600, 11113)

    # model_starter(25,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11111)
    # model_starter(26,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11112)
    # model_starter(27,:uni_lag, DynamicSDDiP.Core_Relint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_relint.log", 3600, 11113)

    # model_starter(28,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11111)
    # model_starter(29,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11112)
    # model_starter(30,:uni_lag, DynamicSDDiP.Core_Optimal(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_opt.log", 3600, 11113)

    # model_starter(31,:uni_lag, DynamicSDDiP.L∞_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_Lsup.log", 3600, 11113)
    # model_starter(32,:uni_lag, DynamicSDDiP.L₁_Deep(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_L1.log", 3600, 11113)

end

#model_starter_runs()
