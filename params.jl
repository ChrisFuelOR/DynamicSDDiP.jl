# # Setting algorithmic parameters

# When using DynamicSDDiP.jl, the user has the freedom to specify several algorithmic parameters. To hand these parameters to the algorithm, they should be stored in a struct of type `AlgoParams`, which is defined in `typedefs.jl`.

using DynamicSDDiP
using SDDP

mutable struct AlgoParams
    stopping_rules::Vector{SDDP.AbstractStoppingRule}
    regularization_regime::DynamicSDDiP.AbstractRegularizationRegime
    cut_aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime
    cut_selection_regime::DynamicSDDiP.AbstractCutSelectionRegime
    cut_generation_regimes::Vector{DynamicSDDiP.CutGenerationRegime}
    simulation_regime::DynamicSDDiP.AbstractSimulationRegime
    risk_measure::SDDP.AbstractRiskMeasure
    forward_pass::SDDP.AbstractForwardPass
    sampling_scheme::SDDP.AbstractSamplingScheme
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme
    parallel_scheme::SDDP.AbstractParallelScheme
    cut_type::SDDP.CutType
    refine_at_similar_nodes::Bool
    cycle_discretization_delta::Float64
    print_level::Int
    log_frequency::Int
    log_file::String
    run_numerical_stability_report::Bool
    numerical_focus::Bool
    silent::Bool
    infiltrate_state::Symbol
    seed::Union{Nothing,Int}
    run_description::String
    solver_approach::Union{DynamicSDDiP.GAMS_Solver,DynamicSDDiP.Direct_Solver}
end

# Note that all the parameters starting from `risk_measures` to `cycle_discretization_delta` are standard parameters of the `train` function in SDDP.jl which are required in our code as we are re-using some functionality from SDDP.jl.
# We have not changed them within our experiments for this paper, and they also should not be changed from their default setting for DynamicSDDiP.jl to work.

# We now discuss the parameters that actually can and should be adjusted by the user.
# To get an idea on how this is done, see the files for our computational experiments in the `examples` folder. For each run (or batch of runs) there exists an `algo_config.jl` file with the sole purpose to configure the parameters. 

# ## Stopping rules

# As for SDDP.jl, we can define a list of stopping rules for Dynamic SDDiP. In our experiments, we always used a time limit (`SDDP.TimeLimit`) (which differed for different test cases) and additionally stopped SDDP when the lower bounds did not improve by more than 1e-4 for 20 iterations (`SDDP.BoundStalling`). The second stopping criterion was mostly relevant when using only Benders or strengthened Benders cuts. 

using SDDP
time_limit = 10800
stopping_rules = [SDDP.TimeLimit(time_limit), SDDP.BoundStalling(20,1e-4)]

# Alternatively, it is also possible to define an iteration limit (`SDDP.IterationLimit`) or use a deterministic stopping criterion when no uncertainty is prevalent (`DeterministicStopping`; defined in `typedefs.jl`).

using SDDP
mutable struct DeterministicStopping <: SDDP.AbstractStoppingRule
    rtol::Float64
    atol::Float64
    function DeterministicStopping(;
        rtol=1e-8,
        atol=1e-8
    )
        return new(rtol, atol)
    end
end


# ## Cut generation approach

# A key aspect of DynamicSDDiP is that we can specify exactly how cuts should be generated in the backward pass of SDDiP. This can be controlled using two important concepts.

# First, we can specify which type of dual problem should be solved. There exist four main structs (`LinearDuality`, `StrengthenedDuality`, `LagrangianDuality`, `UnifiedLagrangianDuality`), each a subtype of `AbstractDualityRegime`. They determine the general solution approach and offer a variety of parameters for further configuration. 

# For instance, if we want to generate Benders cuts, we have to use

using DynamicSDDiP
duality_regime = DynamicSDDiP.LinearDuality()

# Second, we can define if a binary approximation of the state space should be applied in the cut generation process. This can be done by setting a `state_approximation_regime` variable to an object of either `NoStateApproximation` or `BinaryApproximation`.

# In our experiments for this paper, we did not use this second feature, thus always setting

using DynamicSDDiP
state_approximation_regime = DynamicSDDiP.NoStateApproximation()

# in the files `algo_config.jl`.

# !!! note "Remark"
#     This selection only refers to a temporary binary approximation of the state space within the cut generation process that will be applied automatically. In the original state space, this will lead to non-convex approximations of the value functions (**non-convex cuts**). 
#
#     Users can still apply a static and permanent binary approximation of the state space by adjusting their problem formulation accordingly.

# We discuss both concepts in more detail in [Configuring the cut generation](cut_generation.md) and [Binarization and non-convex cuts](binarization.md).

# With the `duality_regime` and the `state_approximation_regime` set, we can define an overall `cut_generation_regime`. For this purpose, there exists the struct `CutGenerationRegime`. Each object of this type has to contain a specific `duality_regime` and `state_approximation_regime`.

using DynamicSDDiP
mutable struct CutGenerationRegime
    state_approximation_regime::DynamicSDDiP.AbstractStateApproximationRegime
    duality_regime::DynamicSDDiP.AbstractDualityRegime
    iteration_to_start::Int64
    iteration_to_stop::Union{Int64,Float64}
    cut_away_approach::Bool
    cut_away_tol::Float64
end

# It is also possible to use several cut generation regimes in the same run, for instance if (strengthened) Benders cuts and Lagrangian cuts should be used within SDDiP. As you can see above, the `AlgoParams` struct always expects a vector of `CutGenerationRegime`.

# If not specified otherwise, all types of cuts are generated in each iteration. However, this can be fine-tuned. 
# The parameters `iteration_to_start` and `iteration_to_stop` allow to restrict the cut generation regime to a subset of iterations. 
# If `Cut_away_approach` is set to `true`, cuts of a regime will only be added to the subproblems if they lead to an improvement, i.e. cut away the current incumbent $(x_{a(n)}^i, \theta_n^i)$ by at least `cut_away_tol`. This is supposed to prevent adding redundant cuts.

# As an example, we can define to generate both strengthened Benders cuts and Lagrangian cuts, but the latter only starting from iteration 20.

using DynamicSDDiP
state_approximation_regime = DynamicSDDiP.NoStateApproximation()

cut_generation_regime_1 = DynamicSDDiP.CutGenerationRegime(
    state_approximation_regime = state_approximation_regime,
    duality_regime = DynamicSDDiP.StrengthenedDuality(),
)

cut_generation_regime_2 = DynamicSDDiP.CutGenerationRegime(
    state_approximation_regime = state_approximation_regime,
    duality_regime = DynamicSDDiP.LagrangianDuality(),
    iteration_to_start = 20,
)

cut_generation_regimes = [cut_generation_regime_1, cut_generation_regime_2]


# ## Lipschitz regularization

# Our code allows to apply a Lipschitz regularization. This can be done by setting a `regularization_regime` variable to an object of `Regularization`. Otherwise, `NoRegularization` has to be chosen.

using DynamicSDDiP

mutable struct Regularization <: DynamicSDDiP.AbstractRegularizationRegime
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    norm :: DynamicSDDiP.AbstractNorm
    norm_lifted :: DynamicSDDiP.AbstractNorm
    copy_regime :: DynamicSDDiP.AbstractCopyRegime
end

mutable struct NoRegularization <: DynamicSDDiP.AbstractRegularizationRegime end

# If `Regularization` is used, this means that in the forward pass of SDDiP a Lipschitz regularization with Lipschitz constant $\sigma_n$ (defined by parameter `sigma`) and the norm defined in `norm` is applied for each subproblem.
# The copy constraint $z_n = x_{a(n)}^i$ is removed and the expression $\sigma_n \lVert z_n - x_{a(n)}^i \rVert$ is added to the objective function. In addition, with parameter `copy_regime` we can specify if the variables $z_n$ should satisfy certain constraints in the forward pass problem. For details, we refer to [Handling Copy Constraints](copy_constraints.md).

# In the backward pass subproblems a Lipschitz regularization with parameter `sigma` is applied as well. However, here we can use a different norm which is specified by `norm_lifted`.

# !!! note "remark"
#     Note that if `sigma` is not chosen sufficiently large, it is not guaranteed that the original MS-MILP is solved when the regularization is applied. Therefore, we can specify a factor by which `sigma` is increased whenever the algorithm gets stuck without converging (`sigma_factor`). This has to be carefully checked (see file `sigmaTest.jl`) and is particularly challenging for stochastic problems. So far, we have only used it for multistage deterministic problems.

# !!! note "remark" 
#     Whereas a regularization can always be applied, this is mostly relevant for the generation of non-convex cuts, and thus in combination with a temporary state binarization (see above). This also explains why we can specify different norms for the forward and backward pass: In the backward pass we may work in a lifted space that requires a different norm (also note that our code automatically uses a weighted variant of `norm_lifted`). For more details, we refer to our [preprint on this topic](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/).

# For the experiments in this paper, the regularization_regime was always set to DynamicSDDiP.NoRegularization().

regularization_regime = DynamicSDDiP.NoRegularization()


# ## Cut aggregation

# We can decide whether cuts should be considered separately for each scenario (multi-cut approach) or aggregated such that there is only one set of cuts per stage (single-cut approach). To control this we can either choose `SingleCutRegime` or `MultiCutRegime`. Internally, these parameters are automatically related to the SDDP.jl parameters `SDDP.SINGLE_CUT` and `SDDP.MULTI_CUT`.

using DynamicSDDiP
mutable struct SingleCutRegime <: DynamicSDDiP.AbstractCutAggregationRegime end
mutable struct MultiCutRegime <: DynamicSDDiP.AbstractCutAggregationRegime end

# !!! note "remark" 
#     The new types of Lagrangian cuts introduced in our paper require to use a multi-cut approach. Therefore, for the experiments in our paper, only comparative runs using (strengthened) Benders cuts or standard Lagrangian cuts were conducted with `SingleCutRegime`.


# ## Cut selection

# To reduce the size of the subproblems, it is possible to apply a basic cut selection scheme which is borrowed from SDDP.jl. It tries to keep the number of active cuts at a minimum by removing cuts that are dominated at previously visited incumbents. This can be controlled by either choosing `CutSelection` or `NoCutSelection`. 

using DynamicSDDiP
mutable struct CutSelection <: DynamicSDDiP.AbstractCutSelectionRegime
    cut_deletion_minimum::Int
end

mutable struct NoCutSelection <: DynamicSDDiP.AbstractCutSelectionRegime end

# The parameter `cut_deletion_minimum` is taken from SDDP.jl and specifies a minimum number of cuts to cache before deleting cuts from the subproblem.

# In the experiments for which we report results in this paper, we did not apply cut selection. The reason is that we did not observe improvements in preliminary tests.


# ## Simulation

# After SDDiP is finished, it is typical to run a simulation to evaluate the computed policy. This can be done using in-sample data (`SDDP.InSampleMonteCarlo`) or out-of-sample data (`SDDP.OutOfSampleMonteCarlo`). 

# This step is important if reasonable statistical upper bounds shall be reported, as computing them is not possible from within the algorithm if only one scenario is sampled per forward pass (which is the default for SDDP.jl and also for DynamicSDDiP.jl).

# For the experiments in this paper, we always applied an in-sample simulation with 1,000 scenarios.


# ## Further parameters

# The final parameters concern some further details of the code, but not the algorithm itself. `Silent` is a boolean that defines if the solver output is printed to the console or not. `infiltrate_state` is only relevant for debugging purposes. `Solver_approach` should always be set to `Direct_Solver()` unless the GAMS.jl package is used. As for SDDP.jl, `log_file` allows to specify the path to a file in which the main run parameters and the output are logged. 

using DynamicSDDiP
silent = true
infiltrate_state = :none
solver_approach = DynamicSDDiP.Direct_Solver()


# All previously defined parameters are stored in a struct of type AlgoParams and then passed to the SDDiP algorithm in the run-file.


# ## Problem parameters

# In addition to the algorithmic parameters, some parameters of the model can be stored in struct `ProblemParams`. 

struct ProblemParams
    number_of_stages::Int
    number_of_realizations::Int
    tree_seed::Union{Nothing,Int}
end

# Here, `number_of_stages` refers to the number $T$ of stages in the problem. `number_of_realizations` refers to the finite number of realizations of the stagewise independent uncertainty. `tree_seed` refers to the seed that is used to generate the finitely many realizations per stage that represent the uncertainty in the model (i.e. the finite scenario tree). This should not be confused with the seed in `AlgoParams` that refers to how we sample from these realizations in the forward pass of SDDiP.


# ## Solvers

# The struct `AppliedSolvers` allows to specify details of the solvers that are used for the subproblems within SDDiP, including some solver options. In particular, it is possible to define different solvers for different types of problems (LPs, MILPs, MINLPs etc) but also for different phases of the algorithm (subproblems vs. Lagrangian relaxation). 

struct AppliedSolvers
    general_solver :: Bool
    solver_subproblem :: Any
    solver_Lagrange_relax :: Any
    solver_Lagrange_approx :: Any
    solver_Lagrange_Bundle :: Any
    solver_Lagrange_subgradient :: Any
    solver_cut_selection :: Any
    solver_LP_relax :: Any
    solver_reg :: Any
    solver_norm :: Any
    solver_tol :: Float64
    solver_time :: Int64
end

# Note that we can also specify the tolerance `solver_tol` for solving subproblems. This is the same for all problems, though. Moreover, (for now only for Gurobi) we can specify a time limit `solver_time` in seconds.

# In our experiments for this paper, we used Gurobi for all occuring subproblems.

# !!! note "Remark"
#     When `GAMS.jl` should be used, we have to use `GAMS_Solver` instead of `Direct_Solver` for the `AbstractSolverApproach`. Then, we can still define the actual solvers to be used within GAMS in structs of type `AppliedSolvers`. 
