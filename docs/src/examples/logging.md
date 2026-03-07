```@meta
EditURL = "logging.jl"
```

# Understanding the logging results

For each run in our experiments, several outputs are logged in a dedicated log-file.
It consists of the following elements:

 * A section printing the main model run parameters and some general information
     * the date and the time of the run
     * the main algorithmic parameters stored in `AlgoParams` (`stopping rules`, `regularization_regime`, vector of `cut_generation_regime` including the `dual_regime` with all parameters and the `state_approximation_regime`, `cut_aggregation_regime`, `cut_selection_regime`, late binarization information (no longer part of the code))
     * the applied solvers and solver parameters
     * the seed for the simulation
     * an optional run description by the user
     * the main problem parameters stored in `ProblemParams` (number of stages, number of realizations per stage, seed for th forward pass)
 * A section logging information from the SDDiP iterations. Each row contains
     * the iteration number
     * the upper bound estimator for the respective iteration (this is a very bad estimator and only relevant for deterministic problems)
     * the max upper bound estimator found so far (this is a very bad estimator and only relevant for deterministic problems)
     * the deterministic lower bound
     * the optimality gap (this is only relevant for deterministic problems)
     * the deterministic lower bound
     * the total time in seconds
     * the time in seconds for the respective iteration
     * some information on refinements of the regularization parameter or the state binarization that took place (if `Regularization` and `BinaryApproximation` are used at all)
     * the total number of variables in each subproblem
     * the number of binary variables in each subproblem
     * the total number of original constraints in each subproblem
     * the total number of cuts in the subproblems (aggregated over all stages)
     * the active number of cuts in the subproblems (aggregated over all stages; relevant if `CutSelection` is used)
     * some information on the total number of iterations required to solve all the Lagrangian dual problems (if any) in the respective iteration
 * A table summarizing timing and memory allocation information for different steps of the algorithm
 * A summary of the final statuses of solving the Lagrangian dual problems (if any)
     * opt: the problem was solved to optimality
     * iter: the solution stopped at the iteration limit
     * conv: the solution method converged to a suboptimal solution
     * unbounded: the problem was detected as unbounded (relevant for linear normalization); a SB cut is generated instead
     * bound_issues or feas_issues: the solution method showed some numerical issues; a cut is generated using the best dual multipliers at the time when the issues are detected
     * sub, subgr_stalling: were relevant in an earlier implementation of a subgradient method; not relevant anymore
     * mn_opt: a second optimization step was used (``minimal norm choice''); this step terminated with optimality
     * mn_iter: a second optimization step was used (``minimal norm choice''); this step terminated with the iteration limit
     * mn_issue: a second optimization step was used (``minimal norm choice''); this step terminated with numerical issues
 * A section containing the results of an in-sample or out-of-sample simulation conducted after SDDiP has terminated. In each case the information contains
     * the deterministic lower bound
     * the simulated upper bound and a confidence interval

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

