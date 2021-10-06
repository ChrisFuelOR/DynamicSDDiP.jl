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

# Mutable struct for bundle parameters if level bundle method is used
mutable struct BundleParams
    bundle_alpha :: Float64
    bundle_factor :: Float64
    level_factor :: Float64

# Mutable struct for algorithmic parameters
mutable struct AlgoParams
    # optimality tolerances for the MILP
    ############################################################################
    opt_rtol :: Float64
    opt_atol :: Float64
    # binary approximation parameters
    ############################################################################
    binary_approx :: Bool
    #NOTE: so far, binary precision is the same for all stages
    binary_precision :: Dict{Symbol, Float64}
    # regularization parameters
    ############################################################################
    regularization :: Bool
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    # lagrangian dual specific parameters
    ############################################################################
    lagrangian_atol :: Float64
    lagrangian_rtol :: Float64
    lagrangian_iteration_limit :: Int
    dual_initialization_method :: Symbol
    dual_solution_method :: Symbol
    dual_status_regime :: Symbol
    magnanti_wong :: Bool
    bundle_params :: Union{Nothing,DynamicSDDiP.BundleParams}
    # nonlinear cut parameters
    ############################################################################
    cut_projection_method :: Symbol
    cut_selection :: Bool
    # debugging parameter
    ############################################################################
    infiltrate_state :: Symbol
end

# struct for solvers to be used
struct AppliedSolvers
    LP :: Any
    MILP :: Any
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
