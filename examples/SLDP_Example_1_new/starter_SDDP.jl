import DynamicSDDiP
import SDDP
import Gurobi
import Infiltrator
import MathOptInterface
using Revise
using Printf

include("scenario_tree.jl")
include("model.jl")
include("simulation.jl")

const GRB_ENV_2 = Gurobi.Env()

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

        algo_params.stopping_rules = [SDDP.IterationLimit(20)]

        if duality_regime_sym == :B
            duality_handler = SDDP.ContinuousConicDuality()
        elseif duality_regime_sym == :SB
            duality_handler = SDDP.StrengthenedConicDuality()
        elseif duality_regime_sym == :lag
            duality_handler = SDDP.LagrangianDuality(atol=1e-4,rtol=1e-4, iteration_limit=1000)
        end

        ########################################################################
        # DEFINE MODEL
        ########################################################################
        model_output = model_set_up(number_of_stages, number_of_realizations, algo_params=algo_params, applied_solvers=applied_solvers, tree_seed=tree_seed)
        model = model_output.model
        problem_params = model_output.problem_params

        for (node_index, node) in model.nodes
            node.optimizer = () -> Gurobi.Optimizer(GRB_ENV_2)
        end

        ########################################################################
        # SOLVE (TRAIN) MODEL
        ########################################################################
        Random.seed!(forward_seed)

        SDDP.train(
            model,
            time_limit = time_limit,
            print_level = 2,
            log_file = "SDDP.log",
            log_frequency = 1,
            stopping_rules = algo_params.stopping_rules,
            cut_type = algo_params.cut_type,
            cut_deletion_minimum = 1000, # can we turn off cut selection?
            duality_handler = duality_handler
        )

        # Which solver? Silent?

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


function model_starter_runs()

    """
    Specification of model runs that should be run one after the other.
    """

    # Dummy run
    # model_starter(0,8,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 100, 11111, 12345)

    model_starter(9,8,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.SingleCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_single.log", 1500, 11111, 12345)
    # model_starter(10,8,10,:lag, DynamicSDDiP.Core_Midpoint(), DynamicSDDiP.MultiCutRegime(), DynamicSDDiP.NoCutSelection(), "C:/Users/cg4102/Documents/julia_logs/sldp_binary_multi.log", 3600, 11111, 12345)

end

#model_starter_runs()
