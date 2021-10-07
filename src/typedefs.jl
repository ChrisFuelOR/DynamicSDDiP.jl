# The structs
# > "NonlinearCut",
# > "InnerLoopIterationResult"
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

# Deterministic stopping rule
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

# Different regimes for state approximation
################################################################################
abstract type AbstractStateApproximationRegime end

mutable struct BinaryApproximation <: AbstractStateApproximationRegime
    binary_precision::Dict{Symbol, Float64} # so far, binary precision is the same for all stages
    function BinaryApproximation(;
        binary_precision = ...,
    )
        return new(binary_precision)
    end
end

mutable struct NoStateApproximation <: AbstractStateApproximationRegime end

"""
BinaryApproximation means that a dynamically refined binary approximation
    is used in the backward pass. It also implies that cuts have to be projected
    back to the original space.
NoStateApproximation means that all cuts are
    generated in the original space, and thus may not be tight.
Default is BinaryApproximation.
"""

# Different regimes for regularization
################################################################################
abstract type AbstractRegularizationRegime end

mutable struct Regularization <: AbstractRegularizationRegime
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    function BinaryApproximation(;
        sigma = 1,
        sigma_factor = 5
    )
        return new(sigma, sigma_factor)
    end
end

mutable struct NoRegularization <: AbstractRegularizationRegime end

"""
Regularization means that in the forward pass some regularized value functions
    are considered. It also implies that a sigma test is conducted once the
    stopping criterion is satisfied. Furthermore, it can be exploited in combination
    with bounding the dual variables.
NoRegularization means that no regularization is used. This may be detrimental
    w.r.t. convergence.
Default is Regularization.
"""

# Different regimes for initialization of dual multipliers
################################################################################
abstract type AbstractDualInitializationRegime end

mutable struct ZeroDuals <: AbstractDualInitializationRegime end
mutable struct LPDuals <: AbstractDualInitializationRegime end
mutable struct CPLEXFixed <: AbstractDualInitializationRegime end

"""
ZeroDuals means that the dual multipliers are initialized as zero.
LPDuals means that the LP relaxation is solved and the corresponding optimal
    dual multipliers are used as an initial solution.
CPLEXFixed can only be chosen if CPLEX is used to solve the subproblems.
    In contrast to Gurobi, CPLEX provides marginal values even when the primal
    MILP is solved. These dual values are determined by solving an LP after fixing
    the integer variables to their optimal values. They are hard to interpret,
    but still may be used for initialization.
Default is ZeroDuals.
"""

# Different regimes for solving the Lagrangian dual
################################################################################
abstract type AbstractDualSolutionRegime end

mutable struct Kelley <: AbstractDualSolutionRegime
    atol::Float64
    rtol::Float64
    iteration_limit::Int
    function LevelBundle(;
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        )
        return new(atol, rtol, iteration_limit)
    end
end

mutable struct LevelBundle <: AbstractDualSolutionRegime
    atol::Float64
    rtol::Float64
    iteration_limit::Int
    bundle_alpha::Float64
    bundle_factor::Float64
    level_factor::Float64
    function LevelBundle(;
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        bundle_alpha = 1.0,
        bundle_factor = 1.0,
        level_factor = 1.0,
        )
        return new(atol, rtol, iteration_limit, bundle_alpha, bundle_factor, level_factor)
    end
end

"""
Kelley means that a classical cutting-plane method is used to solve the dual
    problem.
LevelBundle means that a level bundle method with specified parameters is
    used to solve the dual problem.
Default is Kelley.
"""

# Different regimes for bounding the dual variables in the Lagrangian dual
################################################################################
abstract type AbstractDualBoundRegime end

mutable struct ValueBound <: AbstractDualBoundRegime end
mutable struct NormBound <: AbstractDualBoundRegime end
mutable struct BothBounds <: AbstractDualBoundRegime end

"""
ValueBound means that the optimal value of the Lagrangian dual is bounded from
    the beginning using the optimal value of the primal problem (as in the
    SDDP package). However, the dual multipliers are not bounded.
    This may result in infinitely steep cuts.
NormBound means that the dual multipliers in the Lagrangian dual are bounded
    in their norm. This makes most sense when using a regularization and choosing
    a dual bound related to the regularization parameter sigma.
BothBounds means that ValueBound and NormBound are both used.
Default is BothBounds.
"""

# Different regimes for handling the precision of the Lagrangian dual solution
################################################################################
abstract type AbstractDualStatusRegime end

mutable struct Rigorous <: AbstractDualStatusRegime end
mutable struct Lax <: AbstractDualStatusRegime end

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
This only makes a difference if the cut type is LagrangianCut.
Default is Rigorous.
"""

# Different regimes for choosing the best dual solution if degenerated
################################################################################
abstract type AbstractDualChoiceRegime end

mutable struct StandardChoice <: AbstractDualChoiceRegime end
mutable struct MagnantiWongChoice <: AbstractDualChoiceRegime end

"""
Standard choice means that the Lagrangian multipliers are used as determined
    by the cutting-plane or level bundle method.
MagnantiWongChoice means that it is attempted to determine pareto-optimal
    dual multipliers (as in the newer SDDP package version).
Default is MagnantiWongChoice.
"""

# Different regimes for the cut type
################################################################################
abstract type AbstractCutTypeRegime end

mutable struct LagrangianCut <: AbstractCutTypeRegime end
mutable struct BendersCut <: AbstractCutTypeRegime end
mutable struct StrengthenedCut <: AbstractCutTypeRegime end

"""
LagrangianCut means that the Lagrangian dual is (approximately) solved to obtain
    a cut.
BendersCut means that the LP relaxation is solved to obtain a cut.
StrengthenedCut means that the Lagrangian relaxation is solved using the optimal
    dual multiplier of the LP relaxation to obtain a strengthened Benders cuts.
Default is LagrangianCut.
"""

# Different regimes for the cut projection
################################################################################
abstract type AbstractCutProjectionRegime end

mutable struct BigM <: AbstractCutProjectionRegime end
mutable struct SOS1 <: AbstractCutProjectionRegime end
mutable struct Bilinear <: AbstractCutProjectionRegime end

"""
BigM means that the complementarity constraints of the KKT conditions of
    the cut projection closure are reformulated using a big-M approach.
SOS1 means that the complemnentarity constraints of the KKT conditions of
    the cut projection closure are reformulated using SOS1 constraints.
Bilinear means that by using strong duality the cuts are integrated in a
    bilinear way. This implies that MINLP solvers have to be used to solve
    the subproblems.
Default is BigM.
Note that this only makes a differences if a binary approximation is used,
    otherwise no projection is required at all.
"""

# Different regimes for the cut selection
################################################################################
abstract type AbstractCutSelectionRegime end

mutable struct CutSelection <: AbstractCutSelectionRegime end
mutable struct NoCutSelection <: AbstractCutSelectionRegime end

"""
CutSelection means that a simple cut selection procedure is used to identify
    dominated cuts (as in SDDP).
NoCutSelection means that no such procedure is used, so all cuts are used
    in each iteration.
"""

################################################################################
################################################################################

# Mutable struct for algorithmic parameters
mutable struct AlgoParams
    stopping_rules::Vector{SDDP.AbstractStoppingRule}
    state_approximation_regime::AbstractStateApproximationRegime
    regularization_regime::AbstractRegularizationRegime
    dual_initialization_regime::AbstractDualInitializationRegime
    dual_solution_regime::AbstractDualSolutionRegime
    dual_bound_regime::AbstractDualBoundRegime
    dual_status_regime::AbstractDualStatusRegime
    dual_choice_regime::AbstractDualChoiceRegime
    cut_type_regime::AbstractCutTypeRegime
    cut_projection_regime::AbstractCutProjectionRegime
    cut_selection_regime::AbstractCutSelectionRegime
    infiltrate_state::Symbol

    function AlgoParams(;
        stopping_rules = [DeterministicStopping()],
        state_approximation_regime = BinaryApproximation(),
        regularization_regime = Regularization(),
        dual_initialization_regime = ZeroDuals(),
        dual_solution_regime = Kelley(),
        dual_bound_regime = BothBounds(),
        dual_status_regime = Rigorous(),
        dual_choice_regime = MagnantiWongChoice(),
        cut_type_regime = LagrangianCut(),
        cut_projection_regime = BigM(),
        cut_selection_regime = CutSelection(),
        infiltrate_state = :None
        )
        return new(
            stopping_rules,
            state_approximation_regime,
            regularization_regime,
            dual_initialization_regime,
            dual_solution_regime,
            dual_bound_regime,
            dual_status_regime,
            dual_choice_regime,
            cut_type_regime,
            cut_projection_regime,
            cut_selection_regime,
            infiltrate_state,
            )
    end


end

# struct for solvers to be used
struct AppliedSolvers
    LP :: Any
    MILP :: Any
    MINLP :: Any
    NLP :: Any
    Lagrange :: Any # can be specified separately if numerical issues occur
end

# Struct to store information on a nonlinear cut
mutable struct NonlinearCut
    # cut coefficients (in binary space)
    ############################################################################
    intercept::Float64 # intercept of the cut (dual function value)
    coefficients::Dict{Symbol,Float64} # optimal dual variables in binary space
    # cut construction point
    ############################################################################
    # NOTE: not sure if all of these have to be stored (anchor_state can be
    # determined from other information for example)
    trial_state::Dict{Symbol,Float64} # trial point at which cut should have been created
    anchor_state::Dict{Symbol,Float64} # anchor point at which this cut was created
    binary_state::Dict{Symbol,BinaryState} # binary representation of anchor point
    binary_precision::Dict{Symbol,Float64} # binary precision at moment of creation
    # sigma at moment of creation
    ############################################################################
    sigma::Float64
    # references to variables and constraints of the cut projection closure
    ############################################################################
    cut_variables::Vector{JuMP.VariableRef}
    cut_constraints::Vector{JuMP.ConstraintRef}
    # SDDP-specific stuff
    ############################################################################
    obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    belief_y::Union{Nothing,Dict{T,Float64} where {T}}
    # number of non-dominated cuts (required for cut selection)
    ############################################################################
    non_dominated_count::Int # SDDP
    # iteration in which cut was created
    ############################################################################
    iteration::Int64
end

# mutable struct for iteration results
mutable struct IterationResult{T,S}
    lower_bound :: Float64
    upper_bound :: Float64 # statistical if we solve stochastic problems
    current_sol :: Array{Dict{Symbol,Float64},1} # current state (also required for binary refinement)
    scenario_path :: Vector{Tuple{T,S}}
    has_converged :: Bool
    status :: Symbol # solution status (i.e. number of iterations)
    nonlinear_cuts :: Dict{T, Vector{Any}}
    # NOTE: only required for logging, binary expansion
    # however, then also binary precision / K should be stored for interpretability
end

struct BackwardPassItems{T,U}
    "Given a (node, noise) tuple, index the element in the array."
    cached_solutions::Dict{Tuple{T,Any},Int}
    duals::Vector{Dict{Symbol,Float64}}
    supports::Vector{U}
    nodes::Vector{T}
    probability::Vector{Float64}
    objectives::Vector{Float64}
    belief::Vector{Float64}
    bin_state::Vector{Dict{Symbol,BinaryState}}
    lag_iterations::Vector{Int}
    lag_status::Vector{Symbol}

    function BackwardPassItems(T, U)
        return new{T,U}(
            Dict{Tuple{T,Any},Int}(),
            Dict{Symbol,Float64}[],
            U[],
            T[],
            Float64[],
            Float64[],
            Float64[],
            Dict{Symbol,Float64}[],
            Int[],
            Symbol[]
        )
    end
end
