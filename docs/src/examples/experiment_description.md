```@meta
EditURL = "experiment_description.jl"
```

# Computational experiments

To test our proposed variant of SDDiP we conducted experiments on two different test problems:

 * a capacitated lot-sizing problem (CLSP)
 * a capacitated facility location problem (CFLP)

For CLSP we performed tests for two different problem sizes, with 3 state variables (folder `CLSP`) and with 10 state variables (folder `CLSP_Large`).

For each example, the code consists of the following files:

 * A **starter file** (e.g. `starter.jl`): This file is used to run one or several experiments for the test problem.
 * A **config file** (`algo_config.jl`) setting up the parameters for SDDiP.
 * A **model file** (e.g. `model.jl`) in which the multistage optimization problem is defined.
 * A **scenario file** (e.g. `scenario_tree.jl`) preparing the scenarios for the model.
 * A **simulation file** (`simulation.jl`) required for simulations after SDDiP has terminated.

Note that for CLSP we ran tests with a static binary approximation of the state space (CLSP-Bin; files `model.jl` and `starter.jl`) and with the original state space (CLSP; `model_no_bin.jl`, `starter_no_bin.jl`). Therefore, the corresponding folders contain different model and starter files.

## The starter file

In the starter file, the user can start a specific test run by calling the `model_starter` function (or `det_equiv_starter` if the deterministic equivalent should be solved).

The function comes with several arguments that can be specified by the user and (for the most part) are passed to the `algo_params` struct later. They allow to define several model runs (e.g. for different cut or normalization techniques) at once and run them after each other.

For each run, the following steps are executed:

 * The algorithm is configured by calling the `algo_config` function from `algo_config.jl`. This yields a struct of type `AlgoParams` that is passed to the algorithm.
 * The `model_set_up` function from `model.jl` is called to create the multistage optimization problem. It also contains the data (except for the scenario data) for the model.
 * A seed for the forward pass sampling is specified.
 * The `solve` function from `algorithm.jl` is called.
 * The policy obtained with SDDiP is simulated on in-sample or out-of-sample data.

## The config file

Most of the model parameter choices are hard-coded in function `algo_config`, as they do not change much between different runs of our experiments. Others which change more frequently between runs are passed to this function as arguments from the starter file.

````@example experiment_description
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
end
````

Note that all the parameters for a given are also stored in the log-file of this run, see [Understanding the logging output](logging.md). We also describe the runs in a bit more detail below.

For explanations of the parameters that can be set in `algo_config.jl` see [Setting algorithmic parameters](../params.md).

## The results

The full results of our experiments are provided in folder `results` of this repository. For each run presented in the paper, there is a dedicated logging file which contains the parameter configuration of the run as well as the results. To understand the logging files, please see the description in [Understanding the logging output](logging.md).

## Reproducing the results

### CFLP

To reproduce the results from our experiments for CFLP, use the `starter.jl` file in folder `CFLP` with the pre-defined runs. The `algo_config.jl` file also contains the exact configuration to reproduce our results.

For runs in combination with SB cuts, make sure to change the `algo_config.jl` file by uncommenting `iteration_to_start = 21` for `cut_generation_regime_2` and by replacing `cut_generation_regimes = [cut_generation_regime_2]` with `cut_generation_regimes = [cut_generation_regime_1, cut_generation_regime_2]`.

If you want to use different logging files for different runs, the file name arguments in `starter.jl` should be adjusted.

### CLSP with Binarization

To reproduce the results from our experiments for CLSP with binarization, use the `starter.jl` file in folder `CLSP` with the pre-defined runs. The `algo_config.jl` file also contains the exact configuration to reproduce our results.

For 4 or 10 stages, make sure to adjust the variables `stages` and `time_limit` accordingly.

For runs in combination with SB cuts, make sure to change the `algo_config.jl` file by uncommenting `iteration_to_start = 21` for `cut_generation_regime_2` and by replacing `cut_generation_regimes = [cut_generation_regime_2]` with `cut_generation_regimes = [cut_generation_regime_1, cut_generation_regime_2]`.

If you want to use different logging files for different runs, the file name arguments in `starter.jl` should be adjusted.

### CLSP without Binarization

To reproduce the results from our experiments for CLSP without binarization, use the `starter_no_bin.jl` file in folder `CLSP_Large` with the pre-defined runs. The `algo_config.jl` file also contains the exact configuration to reproduce our results.

For runs in combination with SB cuts, make sure to change the `algo_config.jl` file by uncommenting `iteration_to_start = 21` for `cut_generation_regime_2` and by replacing `cut_generation_regimes = [cut_generation_regime_2]` with `cut_generation_regimes = [cut_generation_regime_1, cut_generation_regime_2]`.

If you want to use different logging files for different runs, the file name arguments in `starter_no_bin.jl` should be adjusted.

!!! note "Remarks"
    Note that you may not be able to reproduce the exact same results from our paper, as the output may differ if DynamicSDDiP is ran on a different machine. However, if the machine has similar specifications as the one we used for our experiments (for technical details, see the paper), the results should not differ by too much].

!!! note "Remarks"
    The code also contains the folder `Regularization_Paper_Illustrative` which we used for some toy examples in [our second preprint](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/)].

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

