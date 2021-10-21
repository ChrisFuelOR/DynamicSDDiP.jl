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
CPLEXFixed can only be chosen if CPLEX is used to solve the subproblems.
    In contrast to Gurobi, CPLEX provides marginal values even when the primal
    MILP is solved. These dual values are determined by solving an LP after fixing
    the integer variables to their optimal values. They are hard to interpret,
    but still may be used for initialization.
Default is ZeroDuals.
"""

mutable struct ZeroDuals <: AbstractDualInitializationRegime end
mutable struct LPDuals <: AbstractDualInitializationRegime end
#mutable struct CPLEXFixed <: AbstractDualInitializationRegime end

################################################################################
# SOLUTION METHOD FOR LAGRANGIAN DUAL
################################################################################
abstract type AbstractDualSolutionRegime end

"""
Kelley means that a classical cutting-plane method is used to solve the dual
    problem.
LevelBundle means that a level bundle method with specified parameters is
    used to solve the dual problem.
Default is Kelley.
"""

mutable struct Kelley <: AbstractDualSolutionRegime
    # atol::Float64
    # rtol::Float64
    # iteration_limit::Int
    # function Kelley(;
    #     atol = 1e-8,
    #     rtol = 1e-8,
    #     iteration_limit = 1000,
    #     )
    #     return new(atol, rtol, iteration_limit)
    # end
end

mutable struct LevelBundle <: AbstractDualSolutionRegime
    # atol::Float64
    # rtol::Float64
    # iteration_limit::Int
    level_factor::Float64
    function LevelBundle(;
        # atol = 1e-8,
        # rtol = 1e-8,
        # iteration_limit = 1000,
        level_factor = 1.0,
        )
        #return new(atol, rtol, iteration_limit, level_factor)
        return new(level_factor)
    end
end

################################################################################
# BOUNDS IN LAGRANGIAL DUAL
################################################################################
abstract type AbstractDualBoundRegime end

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
Standard choice means that the Lagrangian multipliers are used as determined
    by the cutting-plane or level bundle method.
MagnantiWongChoice means that it is attempted to determine pareto-optimal
    dual multipliers (as in the newer SDDP package version).
Default is MagnantiWongChoice.
"""

mutable struct StandardChoice <: AbstractDualChoiceRegime end
mutable struct MagnantiWongChoice <: AbstractDualChoiceRegime end

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
# REGULARIZATION
################################################################################
abstract type AbstractRegularizationRegime end

mutable struct Regularization <: AbstractRegularizationRegime
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    function Regularization(;
        sigma = Float64[],
        sigma_factor = 5.0
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
    function LagrangianDuality(;
        atol = 1e-8,
        rtol = 1e-8,
        iteration_limit = 1000,
        dual_initialization_regime = ZeroDuals(),
        dual_bound_regime = BothBounds(),
        dual_solution_regime = Kelley(),
        dual_choice_regime = MagnantiWongChoice(),
        dual_status_regime = Rigorous(),
    )
        return new(atol, rtol, iteration_limit,
            dual_initialization_regime, dual_bound_regime,
            dual_solution_regime, dual_choice_regime, dual_status_regime)
    end
end

mutable struct LinearDuality <: AbstractDualityRegime end
mutable struct StrengthenedDuality <: AbstractDualityRegime end

"""
LagrangianDuality means that the Lagrangian dual is (approximately) solved to obtain
    a cut. In this case, a lot of parameters have to be defined to configure
    the solution of the Lagrangian dual.
LinearDuality means that the LP relaxation is solved to obtain a cut.
StrengthenedDuality means that the Lagrangian relaxation is solved using the optimal
    dual multiplier of the LP relaxation to obtain a strengthened Benders cuts.
Default is LagrangianDuality.
"""

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
    state_approximation_regime::AbstractStateApproximationRegime
    regularization_regime::AbstractRegularizationRegime
    duality_regime::AbstractDualityRegime
    cut_selection_regime::AbstractCutSelectionRegime
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
    infiltrate_state::Symbol

    function AlgoParams(;
        stopping_rules = [DeterministicStopping()],
        state_approximation_regime = BinaryApproximation(),
        regularization_regime = Regularization(),
        duality_regime = LagrangianDuality(),
        cut_selection_regime = CutSelection(),
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
        infiltrate_state = :none,
    )
        return new(
            stopping_rules,
            state_approximation_regime,
            regularization_regime,
            duality_regime,
            cut_selection_regime,
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
            infiltrate_state,
        )
    end
end

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
    MINLP :: Any
    NLP :: Any
    Lagrange :: Any

    function AppliedSolvers(;
        LP = "Gurobi",
        MILP = "Gurobi",
        MINLP = "Gurobi",
        NLP = "SCIP",
        Lagrange = "Gurobi",
        )
        return new(
            LP,
            MILP,
            MINLP,
            NLP,
            Lagrange
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
    ############################################################################
    trial_state::Dict{Symbol,Float64}
    anchor_state::Dict{Symbol,Float64}
    binary_state::Dict{Symbol,DynamicSDDiP.BinaryState}
    binary_precision::Dict{Symbol,Float64}
    ############################################################################
    sigma::Float64
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
    ############################################################################
    trial_state::Dict{Symbol,Float64} # same as anchor state
    ############################################################################
    sigma::Float64
    ############################################################################
    cut_constraint::JuMP.ConstraintRef
    ############################################################################
    # obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    # belief_y::Union{Nothing,Dict{T,Float64} where {T}}
    ############################################################################
    non_dominated_count::Int
    ############################################################################
    iteration::Int64
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
    nonlinear_cuts :: Dict{T, Vector{Any}} #NOTE
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
