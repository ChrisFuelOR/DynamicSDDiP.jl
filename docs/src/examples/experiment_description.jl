# # Computational experiments

# To test our proposed variant of SDDiP we conducted experiments on two different test problems:

#  * a capacitated lot-sizing problem (CLSP)
#  * a capacitated facility location problem (CFLP)

# For CLSP we performed tests for two different problem sizes, with 3 state variables (folder `CLSP`) and with 10 state variables (folder `CLSP_Large`). 

# For each example, the code consists of the following files:

#  * A **starter file** (e.g. `starter.jl`): This file is used to run one or several experiments for the test problem.
#  * A **config file** (`algo_config.jl`) setting up the parameters for SDDiP.
#  * A **model file** (e.g. `model.jl`) in which the multistage optimization problem is defined.
#  * A **scenario file** (e.g. `scenario_tree.jl`) preparing the scenarios for the model.
#  * A **simulation file** (`simulation.jl`) required for simulations after SDDiP has terminated.

# Note that for CLSP we ran tests with a static binary approximation of the state space (CLSP-Bin; files `model.jl` and `starter.jl`) and with the original state space (CLSP; `model_no_bin.jl`, `starter_no_bin.jl`). Therefore, the corresponding folders contain different model and starter files.


# ## The starter file

# In the starter file, the user can start a specific test run by calling the `model_starter` function (or `det_equiv_starter` if the determinist equivalent should be solved). 

# The function comes with several arguments that can be specified by the user and (for the most part) are passed to the `algo_params` struct later. They allow to define several model runs (e.g. for different cut or normalization techniques) at once and run them after each other.

# For each run, the following steps are executed:

#  * The algorithm is configured by calling the `algo_config` function from `algo_config.jl`. This yields a struct of type `AlgoParams` that is passed to the algorithm.
#  * The `model_set_up` function from `model.jl` is called to create the multistage optimization problem. It also contains the data (except for the scenario data) for the model.
#  * A seed for the forward pass sampling is specified.
#  * The `solve` function from `algorithm.jl` is called.
#  * The policy obtained with SDDiP is simulated on in-sample or out-of-sample data.


# ## The config file

# Most of the model parameter choices are hard-coded in function `algo_config`, as they do not change much between different runs of our experiments. Others which change more frequently between runs are passed to this function as arguments from the starter file.

using DynamicSDDiP
function algo_config(
    duality_regime_sym::Symbol,
    normalization_regime::DynamicSDDiP.AbstractNormalizationRegime,
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime,
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime,
    log_file::String,
    time_limit::Int,
    forward_seed::Int,
    )

# We performed a multitude of experiments. However, many of the parameters were chosen the same for all runs reported in the paper, and only a few crucial parameters different between runs. We give a detailed overview on the parameters for each run in [TODO](TODO). Also note that the parameters are logged and part of the log-file of each individual run.

# For explanations of the parameters that can be set in `algo_config.jl` see [Setting algorithmic parameters](params.md).


# !!! note "Remarks"
#     The code also contains the folder `Regularization_Paper_Illustrative` which we used for some toy examples in the [preprint](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/)].