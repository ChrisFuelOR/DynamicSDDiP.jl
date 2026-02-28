# # Implementation of the loglinearSDDP algorithm

# The implementation of loglinearSDDP is contained in the `src` folder and consists of the following files:

#  * `algorithm.jl`:
#      * Contains the main functions from (an earlier version of) SDDP.jl such as `train`, `backward_pass`, `solve_subproblem` etc. that are required for SDDP to work.
#      * Includes a few adjustments to our specific version of SDDP.
#      * In `parameterize` it is made sure that the cut intercepts are re-evaluated for a given scenario.
#      * Before calling the `train` method, the `train_loglinear` function is called to initialize the history of the stochastic process, to compute the cut exponents $\Theta$ (see [cut representation](cuts.md)), adjust the value functions (Bellman functions) to our needs etc.
#  * `ar_preparations.jl`:              
#      * Contains function `initialize_process_state` which uses the history from the `AutoregressiveProcess` struct to set up the initial state of the stochastic process. This state is then used as an input for following stages.
#      * If not enough history is provided by the user, based on the maximum dimension $L$ and the lag order $p$ some default history will be created by the code.
#  * `bellman.jl` and `bellman_redefine.jl`:
#      * Contains the functionality of the value functions (Bellman functions) that are approximated within SDDP, mostly borrowed from SDDP.jl
#      * Includes a few adjustments for our case, e.g. making sure that the code works with our extended `Cut` structs and our specific cut formulas (see [cut representation](cuts.md)).
#  * `cut_computations.jl`
#      * Contains most of the computations that are related to our special version of SDDP
#      * `compute_cut_exponents` computes the exponents $\Theta$ that are used in the cut formulas (see [cut representation](cuts.md)). They only have to be computed once in advance before the iteration loop in SDDP is started. Note that, compared to the paper, the indexing is a bit different. Precisely, $\theta(t,\tau,\ell,m,k)$ from the paper (with $k$ being a stage) translates to `cut_exponents_stage[t][τ,ℓ,m,κ]` (with κ being the lag and κ $= t-k$).
#      * `evaluate_cut_intercepts` updates the existing cut intercepts to a given scenario.
#      * `evaluate_stochastic_cut_intercept_tight` evaluates the cut intercept at the incumbent (including the current process state) where it is constructed.
#      * `update_process_state` updates the state of the stochastic process for a specific realization that was sampled. The current process state is always stored in `node.ext[:process_state]`.
#      * `compute_scenario_factors` computes the scenario-specific factors $\prod_{k=t-p}^{t-1} \prod_{m=1}^{L_k} \xi_{km}^{\Theta(t,\tau,\ell,m,k)}$ that are required in the cut intercept formula (see [cut representation](cuts.md)). They are the same for all cuts.
#      * `adapt_intercepts` iterates over all existing cuts to adjust the intercepts. To this end, first the final intercept value is computed using the `scenario_factors` and then `cut_intercept_variable` is fixed to this value.
#  * `duals.jl`:
#      * Contains the backward pass functionality of SDDP to compute optimal dual multipliers / cut coefficients. Compared to SDDP.jl this is enhanced a lot to be tailored to our special version of SDDP. In particular, the cut intercept factors $\alpha$ that are required in the cut formulas (see [cut representation](cuts.md)) are computed.
#      * `get_alpha` controls the process of computing these factors $\alpha$.
#      * `compute_alpha_t` computes the value of $\alpha^{(t)}$.
#      * `compute_alpha_tau` computes the value $\alpha^{(\tau)}$. If `simplified` is set to `True` in the `AutoregressiveProcess` struct, then a simplified version of this computation can be applied to accelerate the solution process 
#      * `get_existing_cut_factors` computes the first factor $\sum_{r \in R_{t+1}} \rho^*_{rtj} \alpha_{r,t+1,\ell}^{(\tau)}$ required for the computation of $\alpha$. It requires to iterate over all previously generated cuts and to multiply the corresponding dual multiplier with a cut intercept factor of that cut.
#  * `logging.jl`: Controls (most of) the logging of SDDP.
#  * `sampling_schemes.jl`:
#      * Contains functionality for the correct sampling of scenarios within SDDP.
#      * Function `sample_scenario` is taken from SDDP.jl, but adjusted to make sure that it fits log-linear AR processes. It computes a new realization of $\xi_t$ as required in the forward pass. To to dis, it is first sampled only from the stagewise independent process $\eta_t$ using functionality from SDDP.jl. Using the current state of the process and the process formula, then the new value of $\xi_t$ is computed (see formula (15) in the paper).
#      * `update_process_state` uses the new realization to update the state of the process for the following stage.
#      * `sample_backward_noise_terms` is the analogue to `sample_scenario` but for the backward pass of SDDP.
#  * `simulate.jl`: Takes the simulation functionality from SDDP.jl and adjusts it to fit our purposes, e.g. by using our tailored `sample_scenario` function.
#  * `typedefs.jl`: Contains the definition of the `ProblemParams`, `AlgoParams`, `AutoregressiveProcess` and `AutoregressiveProcessStage` structs discussed above.

# !!! note "Remarks"
#     1. We have used packages Tullio.jl and LoopVectorization.jl (with macro `@turbo`) to speed up the computations within `cut_computations.jl` and `duals.jl` as much and keep the overhead to standard SDDP as small as possible.
#     2. We have also implemented some Gurobi-specific methods, e.g. to accelerate the variable fixing or dual muliplier evaluation. These variants can only be used if the correct Gurobi-internal indices for the cut intercept variable, the coupling constraints or the cut constraints are provided using parameters `gurobi_fix_start` `gurobi_coupling_index_start` and  `gurobi_cut_index_start`. For our [hydrothermal scheduling problem](example/experiment_description.md) these values have been identified in advance.