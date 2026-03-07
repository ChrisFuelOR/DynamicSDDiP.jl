```@meta
EditURL = "algorithm.jl"
```

# Implementation of the dynamic SDDiP algorithm

The implementation of DynamicSDDiP is contained in the `src` folder and consists of a multitude of files.

## General functions

The first group of files is not related to steps of SDDiP but provides some general and auxiliary functions:

 * `DynamicSDDiP.jl`: General definition of the module.
 * `typedefs.jl`: Defines all kinds of structs that are either required for the configuration of parameters for the algorithm or required within the algorithm, e.g. to store cut information.
 * `solverHandling.jl`: Contains functionality to configure the solvers for the subproblems occuring in SDDiP, including setting some solver options. The solvers and (some of) their options can be specified by the user in the `AppliedSolvers` struct and are then processed in this file.
 * `logging.jl`: Defines the logging process for our version of SDDiP. Mostly borrowed from SDDP.jl, but adjusted to the specifics of our requirements. In particular, much more information is logged per iteration, see [Interpreting the logging results](examples/logging.md).
 * `stopping.jl`: Contains functionality for the stopping behavior of SDDiP. Mostly borrowed from SDDP.jl, but adapted to fit our logging process. Additionally, a deterministic stopping criterion can be used if the problem is not stochastic.
 * `objective.jl`: Sets the objective function for a subproblem. Mostly borrowed from SDDP.jl, but adapted to our setting.
 * `state.jl`: Contains functionality to deal with state variables. The main functionality is borrowed from SDDP.jl, but we had to make some adjustments. In particular, when setting up subproblems with Lipschitz regularization, state binarization and setting up auxiliary problems we need to make sure that we properly store state variable bounds, integer or binary requirements etc. in order to be able to restore them later.

## State space handling and regularization

The second group of files provides some functionality for state binarization (see [Binarization and non-convex cuts](binarization.md)) and Lipschitz regularization (see [Setting algorithmic parameters](params.md)).

 * `binarization.jl`: Contains functionality for applying a temporary state binarization in the backward pass subproblems.
     * The code is written in such a way that the binary approximation is applied before a subproblem (or its dual) is solved and removed again afterwards. This procedure requires to cache some information on the original state variables in order to be able to recreate them later.
 * `regularizations.jl`: Contains functionality for applying a Lipschitz regularization to the subproblems.
     * The code is written in such a way that a regularization is applied before a subproblem is solved and removed again afterwards. The reason is that in the backward pass a different version of the subproblem has to be considered, which is easier to set up from the initial problem. This procedure requires to cache some information on the state variables and the original objective in order to recreate the original subproblem afterwards.
     * The Lipschitz regularization initially leads to a non-linear objective function. The subproblem is then linearized by introducing additional linear constraints (with the specific steps depending on the chosen regularization norm).
     * When a regularization is applied in the backward pass, a different set of functions is used, as it has to be accounted for a potential state binarization.
 * `sigmaTest.jl`: If Lipschitz regularization is applied in the forward pass, the function in this file performs a forward pass without regularization to see if the sigma parameter is sufficiently large to ensure equivalence between the regularized and the non-regularized problem. This file is not relevant for the experiments in our paper.

## Algorithmic backbone

The third group of files provides the main algorithmic backbone of SDDiP and mostly re-uses code from SDDP.jl.

 * `algorithmMain.jl`: Contains the main functions to start the SDDiP solution process.
     * `solve`: Function that is called to start the solution process.
         * Initializes the logging, the stopping rules, the binary approximation of the state variables (if `state_approximation_regime == BinaryApproximation`) and the Lipschitz regularization (if `regularization_regime == Regularization`).
         * Note that it is required to re-initialize the Bellman functions after the standard SDDP.jl set-up in order to use our extended Bellman function and cut specification, which also allows for non-convex cuts.
         * Start the finalization of the log-files after SDDiP has terminated.
     * `solve_DynamicSDDiP`
         * Makes sure to store some state variable information in a way that is required for the algorithm, especially if a temporary binarization is used in the cut generation process.
         * Starts the logging process.
         * Calls master_loop to start the actual solution process.
     * `master_loop`
         * Actual loop that calls the iteration function for each iteration.
         * Calls the iteration logger.
         * Calls `convergence_handler`, a function that is required for the case where non-convex cuts are generated, and thus a Lipschitz regularization and a dynamic state binarization are used. The function checks if convergence according to the stopping criterion has been achieved and if the parameters `sigma` or `binary_precision` should be updated. For the experiments in this paper, this part is irrelevant.
     * `iteration`
         * Calls the `forward_pass` function.
         * Performs a refinement of the state binarization if required.
         * Calls the `backward_pass` function.
         * Calls the `calculate_bound` function to solve the first-stage problem and compute a lower bound (for minimization; otherwise an upper bound).
         * Updates the best found solution.
         * Prepares and triggers iteration logging.
 * `forward_pass.jl`: Contains the main forward pass functionality, i.e. sampling a scenario and solving the forward pass subproblems. This is mostly borrowed from SDDP.jl. However, there are some notable changes:
     * Values $\theta_n^i$ of the cut approximation are stored in `epi_state`.
     * Applies a Lipschitz regularization to the forward pass subproblem if required. Removes it afterwards.
     * No coverage of objective or belief states.
 * `backward_pass.jl`: Contains the main backward pass functionality, i.e. considering all realizations, calling functions from `dual.jl` to solve the dual subproblems, updating the approximations of the value functions by calling functions `bellman.jl` and solving the first-stage problem in `calculate_bound`. This is mostly borrowed from SDDP.jl. However, there are some notable changes:
     * Passes the `epi_state` information to the subproblem solution.
     * No coverage of objective or belief states.
     * Contains a loop to iterate over a list of `cut_generation_regime` that has been specified in the AlgoParams. If certain criteria are met, cuts are generated using one regime after the other. For more information, see [Setting algorithmic parameters](params.md).
     * If a temporary state binarization is used, the `anchor_state` (which may deviate from the incumbent) is computed here.
     * If the dual restriction approach from [Chen and Luedtke's paper](https://pubsonline.informs.org/doi/10.1287/ijoc.2022.1185) is applied, the dual space will be restricted to a span of coefficients from previous Benders cuts. Using function `update_Benders_cut_list` it is ensured that these coefficients are stored when (strengthened) Benders cuts are generated.

## Dual problems and cut generation

The fourth group of files provides the functionality for solving dual subproblems and generating different types of cuts.

 * `dual.jl`: Contains functionality for solving the dual subproblems. Especially, for each type of `duality_regime` (type of cuts) there exists a variant of function `get_dual_solution.jl`.
     * Includes initializing the dual multipliers according to `dual_initialization_regime`.
     * Includes applying a potential Lipschitz regularization in the backward pass according to `state_approximation_regime` and `regularization_regime`.
     * Includes solving the primal problem to get a bound on the dual objective value. * #      * Includes setting user-defined dual bounds from `dual_bound_regime`.
     * Includes calling a function for solving a Lagrangian dual problem if required. Those functions are included in `lagrange.jl` (for standard Lagrangian cuts) or `lagrange_unified.jl` (for our new Lagrangian cuts).
     * Includes calling a function to determine coefficients for the linear normalization function if Lagrangian LN cuts should be generated.
     * Includes storing the (optimal) dual multipliers in the right format afterwards.

!!! note "Remark"
    As described in our paper, a badly chosen core point candidate may lead to an unbounded normalized Lagrangian dual problem, and a failure of generating a cut. Therefore, we try to identify potential unboundedness in advance. We do so by checking for infeasibility of a related primal problem (for theoretical details, we refer to our paper). This check is done in function `detect_unboundedness`. If potential unboundedness is detected, we may either generate a strengthened Benders cut instead of an LN Lagrangian cut, or introduce artificial bounds to the dual problem.

!!! note "Remark"
    Note that in certain cases ($\pi_{n0} \approx 0$, (u_n, u_{n0}) \approx 0) we do not generate cuts to avoid numerical issues or adding redundant cuts.

 * `lagrange.jl`: Contains functionality to solve the standard Lagrangian dual problem from SDDiP.
     * Includes solving the inner relaxation and the outer problem using Kelley's cutting-plane method or a level bundle method.
     * Includes the option to solve an augmented problem.
     * Includes functionality for `MinimalNormChoice`.
     * Includes option to use suboptimal multipliers as well for approximations in the outer problem.
 * `lagrange_preparation.jl`: Contains some auxiliary functions required for `lagrange.jl`.
     * Includes relaxing and restoring constraints.
     * Includes setting multiplier bounds.
     * Includes regularization preparation. Note that here a weighted variant of `norm_lifted` is used.
     * Includes initializing dual multipliers.
     * Includes storing optimal dual multipliers and status check for the solution.
 * `lagrange_unified.jl`: Same as `lagrange.jl` but catered to normalized Lagrangian dual problems. This means that an additional multiplier $\pi_{n0}$ is considered, normalizations are used and dual space restriction are possible.
 * `lagrange_unified_preparation.jl`: Same as `lagrange_preparation.jl` but catered to normalized Lagrangian dual problems. In particular, contains functions to execute the normalization or dual space restriction.
 * `lagrange_unified_core.jl`: Contains functionality for the computation of core point candidates when generating LN Lagrangian cuts.
     * Includes functions `get_core_point` and `get_normalization_coefficients` for different type of heuristics.
     * Includes functionality to solve primal auxiliary problems for unboundedness detection.

!!! note "Remark"
    Note that function `get_core_point` for `Core_Conv` contains some code that is catered to solving the CLSP problem in our computational experiments. If other problems are solved, this code should be commented out.

## Updating the value function approximations

The fifth group of files provides the functionality for updating the approximations of the value functions by adding new cut constraints.

 * `bellman_jl`: Contains all kinds of functions to model cuts given the dual multipliers obtained in `dual.jl`. The main structure of the code follows that in SDDP.jl. However, the code contains many adjustments and extensions for handling our new types of Lagrangian cuts and non-convex cuts (see [Binarization and non-convex cuts](binarization.md)).
     * Extension of general SDDP.jl functions to the case of `NonlinearCut` (e.g. using `anchor_points`, `sigma`). In particular, the function `represent_cut_projection_closure` is implemented for `SOS1` and `BigM` as `cut_projection_regime`. It models the (linearized) KKT constraints representing a non-convex cut.
     * Extension of general SDDP.jl functions to the case of deep or LN Lagrangian cuts (e.g. using `epi_states` and `dual_0_var`).
     * Some added functionality to handle numerical issues and cut redundance (e.g. `add_cut_flag`, `check_for_cut_away`).
     * No coverage of objective or belief states.

 * `cut_selection.jl`: Contains functionality for cut selection schemes. Mostly borrowed from SDDP.jl but slightly adjusted so that also non-convex cuts can be used.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

