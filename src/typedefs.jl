# The structs
# > "NonlinearCut",
# > "IterationResult"
# > "BackwardPassItems"
# are derived from similar named structs (Cut, IterationResult, BackwardPassItems)
# in the 'SDDP.jl' package by
# Oscar Dowson and released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson <o.dowson@gmail.com>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

import JuMP
import Revise

################################################################################
# BINARY STATE
################################################################################
"""
Struct to store information on the [0,1] (or binary) variables created
in the backward pass in case of BinaryApproximation.

value:  the value of the original state (which has been unfixed)
x_name: the name of the original state, the BinaryState is associated with
k:      the number of components of the [0,1] variable
"""
struct BinaryState
    value::Float64
    x_name::Symbol
    k::Int64
end

################################################################################
# STOPPING RULES
################################################################################
"""
    DeterministicStopping()

Terminate the algorithm once optimality is reached.
"""
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

################################################################################
# INITIALIZIING DUAL MULTIPLIERS
################################################################################
abstract type AbstractDualInitializationRegime end

"""
ZeroDuals means that the dual multipliers are initialized as zero.
LPDuals means that the LP relaxation is solved and the corresponding optimal
    dual multipliers are used as an initial solution.
Default is ZeroDuals.
"""

mutable struct ZeroDuals <: AbstractDualInitializationRegime end
mutable struct LPDuals <: AbstractDualInitializationRegime end

################################################################################
# SOLUTION METHOD FOR LAGRANGIAN DUAL
################################################################################
abstract type AbstractDualSolutionRegime end

"""
Kelley means that a classical cutting-plane method is used to solve the dual
    problem.
LevelBundle means that a level bundle method with specified parameters is
    used to solve the dual problem.
Subgradient means that a basic subgradient method (as in Bertsekas "Nonlinear
    Programming") is used to solve the dual problem.
Default is Kelley.
"""

mutable struct Kelley <: AbstractDualSolutionRegime end

mutable struct LevelBundle <: AbstractDualSolutionRegime
    level_factor::Float64
    switch_to_kelley::Bool
    function LevelBundle(;
        level_factor = 0.5,
        switch_to_kelley = true,
        )
        return new(level_factor, switch_to_kelley)
    end
end

mutable struct Subgradient <: AbstractDualSolutionRegime
    beta_up::Float64
    beta_down::Float64
    gamma::Float64
    wait::Int
    max_times_unchanged::Int
    function Subgradient(;
        beta_up = 1.5,
        beta_down = 0.95,
        gamma = 2,
        wait = 10,
        max_times_unchanged = 10,
        )
        return new(beta_up, beta_down, gamma, wait, max_times_unchanged)
    end
end



################################################################################
# BOUNDS IN LAGRANGIAL DUAL
################################################################################
abstract type AbstractDualBoundRegime end

"""
ValueBound means that the optimal value of the Lagrangian dual is bounded from
    the start using the optimal value of the primal problem (as in SDDP.jl).
    However, the dual multipliers are not bounded.
NormBound means that the dual multipliers in the Lagrangian dual are bounded
    in their norm. This makes most sense when using a regularization and choosing
    a dual bound related to the regularization parameter sigma.
BothBounds means that ValueBound and NormBound are both used.
Default is BothBounds.
"""

mutable struct ValueBound <: AbstractDualBoundRegime end
mutable struct NormBound <: AbstractDualBoundRegime end
mutable struct BothBounds <: AbstractDualBoundRegime end

################################################################################
# HOW TO HANDLE DUAL SOLUTIONS
################################################################################
abstract type AbstractDualStatusRegime end

"""
Rigorous means that the Lagrangian cuts are only used if the Lagrangian
    duals are solved as attempted (i.e. to approximate optimality). Otherwise,
    the algorithm terminates with an error.
Lax means that even if the Lagrangian dual is not solved as attempted
    (e.g. due to reaching the iteration limit without convergence,
    due to subgradients being zero despite LB!=UB, due to numerical issues,
    due to bound stalling), the current values of the multipliers are used
    to define a valid cut for the value function, hoping that this cut
    will suffice to improve the current incumbent.
Default is Rigorous.
"""

mutable struct Rigorous <: AbstractDualStatusRegime end
mutable struct Lax <: AbstractDualStatusRegime end

################################################################################
# HOW TO HANDLE DEGENERATE SOLUTIONS OF THE DUAL
################################################################################
abstract type AbstractDualChoiceRegime end

"""
StandardChoice means that the Lagrangian multipliers are used as determined
    by the cutting-plane or level bundle method.
MinimalNormChoice means that a second step is added to minimize the norm
    among all dual optimal solutions.
Default is MinimalNormChoice.
"""

mutable struct StandardChoice <: AbstractDualChoiceRegime end
mutable struct MinimalNormChoice <: AbstractDualChoiceRegime end

################################################################################
# CUT PROJECTION APPROACH
################################################################################
abstract type AbstractCutProjectionRegime end

"""
BigM means that the complementarity constraints of the KKT conditions of
    the cut projection closure are reformulated using a big-M approach.
SOS1 means that the complemnentarity constraints of the KKT conditions of
    the cut projection closure are reformulated using SOS1 constraints.
KKT means that the complementarity constraints of the KKT conditions are used
    in their original bilinear form.
StrongDuality means that by using strong duality the cuts are integrated in a
    bilinear way.
Default is BigM.

Note that for KKT and StrongDuality the subproblems become nonlinear, and thus
    an MINLP solver (e.g. Gurobi) has to be used to solve them.
    Moreover, in these cases, only zeros should be used as initialization
    method and only Lagrangian cuts should be determined, since the LP
    relaxation is no LP anymore and may yield useless results.
"""

mutable struct BigM <: AbstractCutProjectionRegime end
mutable struct SOS1 <: AbstractCutProjectionRegime end
mutable struct KKT <: AbstractCutProjectionRegime end
mutable struct StrongDuality <: AbstractCutProjectionRegime end

################################################################################
# BINARY APPROXIMATION
################################################################################
abstract type AbstractStateApproximationRegime end

"""
BinaryApproximation means that a dynamically refined binary approximation
    is used in the backward pass. It also implies that cuts have to be projected
    back to the original space.
NoStateApproximation means that all cuts are
    generated in the original space, and thus may not be tight.
Default is BinaryApproximation.

Note that so far, the binary_precision is the same for all stages and only
differs in the state variables.
"""

mutable struct BinaryApproximation <: AbstractStateApproximationRegime
    binary_precision::Dict{Symbol, Float64}
    cut_projection_regime::AbstractCutProjectionRegime
    function BinaryApproximation(;
        binary_precision = Dict{Symbol, Float64}(),
        cut_projection_regime = BigM(),
    )
        return new(binary_precision, cut_projection_regime)
    end
end

mutable struct NoStateApproximation <: AbstractStateApproximationRegime end

################################################################################
# DUAL NORMALIZATION
################################################################################
abstract type AbstractNormalizationRegime end

mutable struct L₁_Deep <: AbstractNormalizationRegime end
mutable struct L₂_Deep <: AbstractNormalizationRegime end
mutable struct L∞_Deep <: AbstractNormalizationRegime end
mutable struct L₁∞_Deep <: AbstractNormalizationRegime end
mutable struct Fischetti <: AbstractNormalizationRegime end #TODO
mutable struct ChenLuedtke <: AbstractNormalizationRegime end

mutable struct Core_Midpoint <: AbstractNormalizationRegime
    integer_relax::Bool
    function Core_Midpoint(;
        integer_relax = false,
        )
        return new(integer_relax)
    end
end

mutable struct Core_In_Out <: AbstractNormalizationRegime
    integer_relax::Bool
    function Core_In_Out(;
        integer_relax = false,
        )
        return new(integer_relax)
    end
end

mutable struct Core_Optimal <: AbstractNormalizationRegime
    integer_relax::Bool
    function Core_Optimal(;
        integer_relax = false,
        )
            return new(integer_relax)
    end
end

mutable struct Core_Relint <: AbstractNormalizationRegime end

mutable struct Core_Epsilon <: AbstractNormalizationRegime
    perturb::Float64
    integer_relax::Bool
    function Core_Epsilon(;
        perturb = 1e-6,
        integer_relax = false,
    )
            return new(perturb, integer_relax)
    end
end

"""
This AbstractType allows to use different normalization in the Lagrangian dual
problem if the unified Lagrangian framework is used.

In the first group of normalization approaches, a norm of the dual multipliers
is bounded.

L_1_Deep means that a normalization is used such that deepest cuts w.r.t. to
    the 1-norm are generated. This is equivalent to one normalization used
    by Chen & Luedtke and related to the normalization by Fischetti et al.
    if copy constraints are used.
L_2_Deep means that a normalization is used such that deepest cuts w.r.t.
    the 2-norm are generated.
L_Inf_Deep means that a normalization is used such that deepest cuts w.r.t.
    the supremum norm are generated.
L_1_Inf_Deep means that a normalization is used such that deepest cuts w.r.t.
    a linear combination of the 1-norm and sup-norm are generated. This can be
    interpreted as a linear approximation of the 2-norm.
ChenLuedtke means that the second normalization approach by Chen & Luedtke
    is used. Here, the dual multipliers have to be restricted to the span
    of Benders multipliers.

In the second group of normalization approaches, a linear function of the dual
multipliers (a linear pseudonorm) is bounded. In all but the first approach
(Fischetti), the coefficients of the linear pseudonorm are determined as the
direction between a core point in the epigraph and the current incumbent
(see Brandenberg & Stursberg for some theory behind this approach).
These approaches are also similar to the traditional strategy by Magnanti
and Wong to compute Pareto-optimal cuts.

Fischetti means that the classical normalization approach by Fischetti et al.
    (2010) is used for the Lagrangian dual.
Core_Midpoint means that the normalization is based on a core point which is
    the midpoint of the state space. This requires that all state variables
    are bounded from below and above.
Core_In_Out means that the normalization is based on a  core point which is
    a convex combination of the previous core point and the current incumbent
    (similar to the in-out-strategy for cutting-plane methods and the
    core point updates by Papadakos).
Core_Epsilon means that the normalization is based on a core point which is a
    a slight perturbation of the current incumbent. This ideas is similar to the
    approach by Sherali & Lunday to compute Pareto-optimal cuts.
Core_Optimal means that the normalization is based on a core point which is the
    optimal point of the current subproblem (with updated cuts and a non-fixed
    state variable).
Core_Relint means that the normalization is based on a core point which is
    determined by solving a special feasibility problem, which guarantees
    that the point lies in the relative interior of the epigraph of the
    closed convex envelope of the value function. This approach is based on
    a similar strategy by Conforti & Wolsey.

The parameter integer_relax is required because if we fix the state variables
    to the coordinates of the core point and if the problem contains integer
    original state or local variables, then these integer requirements may
    not be satisfiable, resulting in an infeasible problem to compute the
    y-coordinate of the core point. For instance, for Generation_Expansion,
    the sum of the binary state variables are fixed to non-integer values for
    the core point computation, but the original state variables still have to
    satisfy integer constraints. To avoid infeasibility, we can consider the LP
    relaxation of the subproblem. If the original state variables are continuous,
    this is not required (see SLDP_Example_1).

Default is L_1_Deep.
"""

################################################################################
# DUAL SPACE RESTRICTION
################################################################################
abstract type AbstractDualSpaceRegime end

mutable struct NoDualSpaceRestriction <: AbstractDualSpaceRegime end
mutable struct BendersSpanSpaceRestriction <: AbstractDualSpaceRegime
    K::Int64
    cut_type::Symbol #:Benders or :Lagrange
end

"""
NoDualSpaceRestriction means that no additional restriction is imposed
    on the dual space, so exact separation is possible.
BendersSpanSpaceRestriction means that only dual multipliers are considered
    which are in the span of the last K Benders multipliers (see Chen & Luedtke).
Default is NoDualSpaceRestriction.
"""

################################################################################
# COPY RESTRICTION
################################################################################
abstract type AbstractCopyRegime end

mutable struct StateSpaceCopy <: AbstractCopyRegime end
mutable struct ConvexHullCopy <: AbstractCopyRegime end
mutable struct NoBoundsCopy <: AbstractCopyRegime end

"""
StateSpaceCopy means that the copy variable has to satisfy the state space
    bounds and integer requirements. This is the classical approach in SDDP.jl.
ConvexHullCopy means that the copy variables has to be in the convex hull of
    the state space (e.g. [0,1] for binary state variables as in SDDiP).
    Since we cannot compute the convex hull for complicated state spaces,
    we assume that the state space is just box-constrained and maybe requires
    integer or binary states. Hence, using this regime the bounds are kept
    as they are, but the integer requirements are relaxed.
NoBoundsCopy means that the copy variables do not have to satisfy any state
    space constraints at all. To avoid unboundedness and infeasibility of the
    dual problem, we set the lower bound of the states to zero and
    the upper bound of the states to 1e9.
    This is equivalent to the approach of relaxing the linking constraint
    without introducing a copy constraint at all.
"""

################################################################################
# CUT FAMILY TO BE USED
################################################################################
abstract type AbstractDualityRegime end

mutable struct LagrangianDuality <: AbstractDualityRegime
    atol::Float64
    rtol::Float64
    iteration_limit::Int
    dual_initialization_regime::AbstractDualInitializationRegime
    dual_bound_regime::AbstractDualBoundRegime
    dual_solution_regime::AbstractDualSolutionRegime
    dual_choice_regime::AbstractDualChoiceRegime
    dual_status_regime::AbstractDualStatusRegime
    copy_regime::AbstractCopyRegime
    augmented::Bool
    function LagrangianDuality(;
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        dual_initialization_regime = ZeroDuals(),
        dual_bound_regime = BothBounds(),
        dual_solution_regime = Kelley(),
        dual_choice_regime = MinimalNormChoice(),
        dual_status_regime = Rigorous(),
        copy_regime = ConvexHullCopy(),
        augmented = false,
    )
        return new(atol, rtol, iteration_limit,
            dual_initialization_regime, dual_bound_regime,
            dual_solution_regime, dual_choice_regime, dual_status_regime,
            copy_regime, augmented)
    end
end

mutable struct LinearDuality <: AbstractDualityRegime end
mutable struct StrengthenedDuality <: AbstractDualityRegime end

mutable struct UnifiedLagrangianDuality <: AbstractDualityRegime
    atol::Float64
    rtol::Float64
    iteration_limit::Int
    dual_initialization_regime::AbstractDualInitializationRegime
    dual_bound_regime::AbstractDualBoundRegime
    dual_solution_regime::AbstractDualSolutionRegime
    dual_choice_regime::AbstractDualChoiceRegime
    dual_status_regime::AbstractDualStatusRegime
    normalization_regime::AbstractNormalizationRegime
    dual_space_regime::AbstractDualSpaceRegime
    copy_regime::AbstractCopyRegime
    user_dual_objective_bound::Union{Nothing,Float64}
    user_dual_multiplier_bound::Union{Nothing,Float64}

    function UnifiedLagrangianDuality(;
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        dual_initialization_regime = ZeroDuals(),
        dual_bound_regime = BothBounds(),
        dual_solution_regime = Kelley(),
        dual_choice_regime = MinimalNormChoice(),
        dual_status_regime = Rigorous(),
        normalization_regime = L₁_Deep(),
        dual_space_regime = NoDualSpaceRestriction(),
        copy_regime = ConvexHullCopy(),
        user_dual_objective_bound = nothing,
        user_dual_multiplier_bound = nothing,
    )
        return new(atol, rtol, iteration_limit,
            dual_initialization_regime, dual_bound_regime,
            dual_solution_regime, dual_choice_regime, dual_status_regime,
            normalization_regime, dual_space_regime, copy_regime,
            user_dual_objective_bound, user_dual_multiplier_bound)
    end
end


"""
LagrangianDuality means that the Lagrangian dual is (approximately) solved to obtain
    a cut. In this case, a lot of parameters have to be defined to configure
    the solution of the Lagrangian dual.
LinearDuality means that the LP relaxation is solved to obtain a cut.
StrengthenedDuality means that the Lagrangian relaxation is solved using the optimal
    dual multiplier of the LP relaxation to obtain a strengthened Benders cuts.
Default is LagrangianDuality.

The new addition UnifiedLagrangianDuality allows to determine cuts using the
unified framework originally proposed by Fischetti et al. (2010) and also used
by Chen & Luedtke (2021) for Lagrangian cuts in 2-stage stochastic programming.
This regime supports different choices for cut selection, e.g. deepest cuts,
facet-defining/Pareto-optimal cuts or the standard normalization by Fischetti
et al. It can be used in a multi-cut or a single-cut setting.
"""

################################################################################
# CUT AGGREGATION REGIME
################################################################################
abstract type AbstractCutAggregationRegime end

mutable struct SingleCutRegime <: AbstractCutAggregationRegime end
mutable struct MultiCutRegime <: AbstractCutAggregationRegime end

################################################################################
# CUT SELECTION
################################################################################
abstract type AbstractCutSelectionRegime end

mutable struct CutSelection <: AbstractCutSelectionRegime
    cut_deletion_minimum::Int
    function CutSelection(;
        cut_deletion_minimum = 1,
    )
        return new(cut_deletion_minimum)
    end
end

mutable struct NoCutSelection <: AbstractCutSelectionRegime end

"""
CutSelection means that a simple cut selection procedure is used to identify
    dominated cuts (as in SDDP).
NoCutSelection means that no such procedure is used, so all cuts are used
    in each iteration.
"""

################################################################################
# DEFINING STRUCT FOR DEFINITION OF CUT GENERATION REGIMES
################################################################################
"""
Instead of directly defining a duality_regime and a state_approximation_regime
in AlgoParams, we now allow to define so-called cut_generation_regimes first.
Each cut_generation_regime from type CutGenerationRegime contains a specific
duality_regime and a specific state_approximation_regime.
AlgoParams then contains a Vector of CutGenerationRegime.
This allows for generating different types of cuts in each iteration
(e.g. Lagrangian cuts and (strengthened) Benders cuts) instead of just one
of them. The additional parameters of cut_generation_regime allow to restrict
the corresponding regime to a subset of iterations.
"""

mutable struct CutGenerationRegime
    state_approximation_regime::AbstractStateApproximationRegime
    duality_regime::AbstractDualityRegime
    iteration_to_start::Int64
    iteration_to_stop::Union{Int64,Float64}
    gap_to_start::Float64       # not used so far
    gap_to_stop::Float64        # not used so far
    cut_away_approach::Bool

    function CutGenerationRegime(;
        state_approximation_regime = BinaryApproximation(),
        duality_regime = LagrangianDuality(),
        iteration_to_start = 1,
        iteration_to_stop = Inf,
        gap_to_start = Inf,
        gap_to_stop = 0.0,
        cut_away_approach = false,
    )
        return new(
            state_approximation_regime,
            duality_regime,
            iteration_to_start,
            iteration_to_stop,
            gap_to_start,
            gap_to_stop,
            cut_away_approach,
        )
    end
end


"""
iteration_to_start:     first iteration at which this regime is applied
iteration_to_stop:      last iteration at which this regime is applied
gap_to_start:           relative optimality gap at which this regime is first
                        applied (tricky for stochastic case)
gap_to_stop:            relative optimality gap at which this regime is last
                        applied (tricky for stochastic case)
cut_away_approach:      if true, a hierarchy of cuts is used, so that this
                        regime is only used if the incumbent is not cut away
                        by the previous cut already. This parameter should
                        not be true for the first CutGenerationRegime in
                        AlgoParams.

Note that using gap_to_start and gap_to_stop for stochastic problems may be
misleading as the upper bounds, and thus the gaps, are stochastic.
"""

################################################################################
# REGULARIZATION
################################################################################
abstract type AbstractNorm end

mutable struct L₁ <: AbstractNorm end
mutable struct L∞ <: AbstractNorm end

abstract type AbstractRegularizationRegime end

mutable struct Regularization <: AbstractRegularizationRegime
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    norm :: AbstractNorm
    norm_lifted :: AbstractNorm
    copy_regime :: AbstractCopyRegime

    function Regularization(;
        sigma = Float64[],
        sigma_factor = 5.0,
        norm = L₁(),
        norm_lifted = L₁(),
        copy_regime = ConvexHullCopy(),
    )
        return new(sigma, sigma_factor, norm, norm_lifted, copy_regime)
    end
end

mutable struct NoRegularization <: AbstractRegularizationRegime end

"""
Regularization means that in the forward pass some regularized value functions
    are considered. It also implies that a sigma test is conducted once the
    stopping criterion is satisfied. Furthermore, it can be exploited in combination
    with bounding the dual variables.
NoRegularization means that no regularization is used. This may be detrimental
    w.r.t. convergence. Note that it also makes it difficult to bound the
    dual multipliers and the bigM parameters appropriately.
    For bigM so far 1e4 is used. For dual multipliers no bound is applied
    even if BothBounds is chosen.
Default is Regularization.

Note that if a regularization is used, it is used in the forward and backward
    pass (in the latter to compute primal_obj and bound the dual variables).

sigma is the regularization parameter.
sigma_factor is a factor with which the regularization parameter is increased
    if required.
norm is the norm that is used as regularization function.
norm_lifted is the norm that is used as regularization function (in a weighted
    form) in the backward pass to bound the dual multipliers if a
    BinaryApproximation of the states is used in the cut generation
    process.
copy_regime defines the constraints that the copy variable z of the state x has
    to satisfy in the forward pass (or backward pass if no state approximation
    is used). For the backward pass with BinaryApproximation, this is separately
    defined by the duality_regime.
"""
#TODO: Maybe define the copy_regime one time in a completeley separate way.

################################################################################
# SIMULATION
################################################################################
# Sampling schemes (similar to the ones in SDDP.jl)
abstract type AbstractSamplingScheme end

mutable struct InSampleMonteCarlo <: AbstractSamplingScheme end

mutable struct OutOfSampleMonteCarlo <: AbstractSamplingScheme
    number_of_realizations :: Int
    simulation_seed :: Int

    function OutOfSampleMonteCarlo(;
        number_of_realizations = 10,
        simulation_seed = 121212,
    )
        return new(simulation_seed)
    end
end

# Simulation regimes
abstract type AbstractSimulationRegime end

mutable struct Simulation <: AbstractSimulationRegime
    sampling_scheme :: DynamicSDDiP.AbstractSamplingScheme
    number_of_replications :: Int

    function Simulation(;
        sampling_scheme = DynamicSDDiP.InSampleMonteCarlo,
        number_of_replications = 1000,
    )
        return new(sampling_scheme, number_of_replications)
    end
end

mutable struct NoSimulation <: AbstractSimulationRegime end

"""
Simulation means that after training the model we perform a simulation with
    number_of_replications. This can be either an in-sample simulation
    (SDDP.InSampleMonteCarlo) or an out-of-sample simulation
    (SDDP.OutOfSampleMonteCarlo). In the latter case, we have to provide a
    method to generate new scenario trees. Different sampling schemes from
    SDDP.jl such as HistoricalSampling or PSRSampling are not supported yet.
NoSimulation means that we do not perform a simulation after training the model,
    either because we do not want to or because we solve a determinist model.
Default is NoSimulation.
"""

################################################################################
# RESAMPLING
################################################################################
abstract type AbstractResamplingRegime end

mutable struct Resampling <: AbstractResamplingRegime
    resampling_limit :: Int

    function Resampling(;
        resampling_limit = 5,
    )
        return new(resampling_limit)
    end
end

mutable struct NoResampling <: AbstractResamplingRegime

"""
In the forward pass and in the simulation sometimes numerical issues occur so that
single subproblems become (or at least are identified) as infeasible or unbounded.
This means that the algorithm stops in such case or no statistical upper bound
is provided.

Resampling means that in both cases a resampling for the current stage is used
    to escape from the error. This procedure is repeated up to resampling_limit
    times before the algorithm actually stops with an error.
    In case of only a few infeasible/unbounded subproblems, this may help
    in still getting results from the algorithm.

NoResampling means that no resampling is applied, but in case of an error,
    this error is thrown immediately.

Note that catching every exception in the forward pass subproblems and to
    use a resampling is pretty risky, as it may also prevent us from identifying
    real infeasibilities in our problem.
    Therefore, it makes most sense to try to solve the problem first without
    resampling and only use this option if we really encounter numerical issues
    for some instances.
"""

################################################################################
# DEFINING STRUCT FOR SOLVERS TO BE USED
################################################################################
"""
For the Lagrangian subproblems a separate solver can be defined if for
    the LP/MILP solver numerical issues occur.
"""

struct AppliedSolvers
    LP :: Any
    MILP :: Any
    MIQCP :: Any
    MINLP :: Any
    NLP :: Any
    Lagrange :: Any

    function AppliedSolvers(;
        LP = "Gurobi",
        MILP = "Gurobi",
        MIQCP = "Gurobi",
        MINLP = "SCIP",
        NLP = "SCIP",
        Lagrange = "Gurobi",
        )
        return new(
            LP,
            MILP,
            MIQCP,
            MINLP,
            NLP,
            Lagrange
            )
    end
end

abstract type AbstractSolverApproach end

mutable struct GAMS_Solver <: AbstractSolverApproach end
mutable struct Direct_Solver <: AbstractSolverApproach end

################################################################################
# DEFINING STRUCT FOR CONFIGURATION OF ALGORITHM PARAMETERS
################################################################################
"""
Note that the parameters from risk_measure to refine_at_similar_nodes are
basic SDDP parameters which are required as we are using some functionality
from the package SDDP.jl. They should not be changed, though, as for different
choices the DynamicSDDiP algorithm will not work.
"""

mutable struct AlgoParams
    stopping_rules::Vector{SDDP.AbstractStoppingRule}
    regularization_regime::AbstractRegularizationRegime
    cut_aggregation_regime::AbstractCutAggregationRegime
    cut_selection_regime::AbstractCutSelectionRegime
    cut_generation_regimes::Vector{CutGenerationRegime}
    simulation_regime::AbstractSimulationRegime
    resampling_regime::AbstractResamplingRegime
    ############################################################################
    risk_measure::SDDP.AbstractRiskMeasure
    forward_pass::SDDP.AbstractForwardPass
    sampling_scheme::SDDP.AbstractSamplingScheme
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme
    parallel_scheme::SDDP.AbstractParallelScheme
    cut_type::SDDP.CutType
    refine_at_similar_nodes::Bool
    cycle_discretization_delta::Float64
    ############################################################################
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

    function AlgoParams(;
        stopping_rules = [DeterministicStopping()],
        regularization_regime = Regularization(),
        cut_aggregation_regime = SingleCutRegime(),
        cut_selection_regime = CutSelection(),
        cut_generation_regimes = [CutGenerationRegime()],
        simulation_regime = NoSimulation(),
        resampling_regime = NoResampling(),
        risk_measure = SDDP.Expectation(),
        forward_pass = SDDP.DefaultForwardPass(),
        sampling_scheme = SDDP.InSampleMonteCarlo(),
        backward_sampling_scheme = SDDP.CompleteSampler(),
        parallel_scheme = SDDP.Serial(),
        cut_type = SDDP.SINGLE_CUT,
        refine_at_similar_nodes = true,
        cycle_discretization_delta = 0.0,
        print_level = 2,
        log_frequency = 1,
        log_file = "DynamicSDDiP.log",
        run_numerical_stability_report = true,
        numerical_focus = false,
        silent = true,
        infiltrate_state = :none,
        seed = nothing,
        run_description = "",
        solver_approach = DynamicSDDiP.GAMS_Solver(),
    )
        return new(
            stopping_rules,
            regularization_regime,
            cut_aggregation_regime,
            cut_selection_regime,
            cut_generation_regimes,
            simulation_regime,
            resampling_regime,
            risk_measure,
            forward_pass,
            sampling_scheme,
            backward_sampling_scheme,
            parallel_scheme,
            cut_type,
            refine_at_similar_nodes,
            cycle_discretization_delta,
            print_level,
            log_frequency,
            log_file,
            run_numerical_stability_report,
            numerical_focus,
            silent,
            infiltrate_state,
            seed,
            run_description,
            solver_approach
        )
    end
end

################################################################################
# DEFINING STRUCT FOR CONFIGURATION OF TEST PROBLEM PARAMETERS
################################################################################
"""
Stores some parameters of the test problem that is solved.
This is mainly used for logging purposes.
"""

struct ProblemParams
    number_of_stages::Int
    number_of_realizations::Int
    tree_seed::Union{Nothing,Int}

    function ProblemParams(
        number_of_stages,
        number_of_realizations;
        tree_seed = nothing,
    )
        return new(
            number_of_stages,
            number_of_realizations,
            tree_seed
        )
    end
end

################################################################################
# DEFINING CUTS
################################################################################

abstract type Cut end

mutable struct NonlinearCut <: Cut
    intercept::Float64
    coefficients::Dict{Symbol,Float64}
    scaling_coeff::Float64
    ############################################################################
    trial_state::Dict{Symbol,Float64}
    anchor_state::Dict{Symbol,Float64}
    binary_state::Dict{Symbol,DynamicSDDiP.BinaryState}
    binary_precision::Dict{Symbol,Float64}
    ############################################################################
    sigma::Union{Nothing,Float64}
    ############################################################################
    cut_variables::Vector{JuMP.VariableRef}
    cut_constraints::Vector{JuMP.ConstraintRef}
    ############################################################################
    # obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    # belief_y::Union{Nothing,Dict{T,Float64} where {T}}
    ############################################################################
    non_dominated_count::Int
    ############################################################################
    iteration::Int64
    ############################################################################
    aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime
    duality_regime::DynamicSDDiP.AbstractDualityRegime
end

"""
The first two arguments define the cut coefficients (in the binary space).

The next four arguments store the trial_state from the forward pass, the
anchor_state at which the cut is constructed (and which may deviate from the
trial_state), the exact binary representation of the anchor_state and the
the binary precision at the moment the cut is constructed.
Actually, not all this information has to be stored. For example, the anchor_point
could as well be derived from binary_state and binary_precision.

The next argument stores the value of the regularization parameter sigma at
the moment the cut is created.

The next two arguments store references to the cut_variables and cut_constraints
(note that using the cut projection one cut refers to several such variables
and constraints).

The next two arguments are SDDP-specific and not required in my case.

Finally, the number of dominated cuts, which is required for the cut selection
process, and the current iteration number (to enumerate the cuts) is stored.

Note that if no binary approximation is used, trial_state and anchor_state
will always be the same and only one cut_constraint has to be stored.
"""

mutable struct LinearCut <: Cut
    intercept::Float64
    coefficients::Dict{Symbol,Float64}
    scaling_coeff::Float64
    ############################################################################
    trial_state::Dict{Symbol,Float64} # same as anchor state
    ############################################################################
    sigma::Union{Nothing,Float64}
    ############################################################################
    cut_constraint::Union{Nothing,JuMP.ConstraintRef}
    ############################################################################
    # obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    # belief_y::Union{Nothing,Dict{T,Float64} where {T}}
    ############################################################################
    non_dominated_count::Int
    ############################################################################
    iteration::Int64
    ############################################################################
    aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime
    duality_regime::DynamicSDDiP.AbstractDualityRegime
end

################################################################################
# DEFINING STORAGE FOR ITERATION RESULTS
################################################################################
"""
This is based on a similar struct in the SDDP package. It stores the results
corresponding to the current iteration.

Note that the upper_bound is statistical only if we consider stochastic problems.
Also note that current_sol only contains the values of the state variables
in the current solution.

Status refers to the number of iterations if I remember correctly.
Nonlinear_cuts is only required for logging, if I remember correctly.
"""

mutable struct IterationResult{T,S}
    lower_bound :: Float64
    upper_bound :: Float64
    current_sol :: Array{Dict{Symbol,Float64},1}
    scenario_path :: Vector{Tuple{T,S}}
    has_converged :: Bool
    status :: Symbol #NOTE
end

################################################################################
# DEFINING STORAGE FOR BACKWARD PASS RESULTS
################################################################################
"""
This is based on a similar struct in the SDDP package. It stores items
corresponding to the current backward pass.
"""

struct BackwardPassItems{T,U}
    "Given a (node, noise) tuple, index the element in the array."
    cached_solutions::Dict{Tuple{T,Any},Int}
    duals::Vector{Dict{Symbol,Float64}}
    dual_0_var::Vector{Float64}
    supports::Vector{U}
    nodes::Vector{T}
    probability::Vector{Float64}
    objectives::Vector{Float64}
    belief::Vector{Float64}
    bin_state::Vector{Dict{Symbol,BinaryState}}
    lag_iterations::Vector{Int}
    add_cut_flags::Vector{Bool}

    function BackwardPassItems(T, U)
        return new{T,U}(
            Dict{Tuple{T,Any},Int}(),
            Dict{Symbol,Float64}[],
            Float64[],
            U[],
            T[],
            Float64[],
            Float64[],
            Float64[],
            Dict{Symbol,Float64}[],
            Int[],
            Bool[]
        )
    end
end
