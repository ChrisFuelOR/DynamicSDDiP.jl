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
#import SDDP

# Struct for algorithmic parameters
# ------------------------------------------------------------------------------------------------------------------
# epsilon_innerLoop: optimality tolerance for inner loop
# binaryPrecision: epsilons for latest/current binary expansion (better use vector instead of dict?)
# sigma: parameters used to obtain the regularized problem (better vector?)
# ------------------------------------------------------------------------------------------------------------------
# possible extensions:
# maxcuts ::  Int64 # maximum number of cuts to be stored (for storage efficiency)
# dropcuts ::  Int64 # number of cuts dropped so far
# what about number of constraints and variables in the model?
# what about differences in the binary expansion of different components?

# Mutable struct for algorithmic parameters that may change during the iterations
# Vector{Float64} or Dict{Int64, Float64}?
mutable struct AlgoParams
    epsilon_innerLoop :: Float64 # optimality tolerance for inner loop (relative!)
    binaryPrecision :: Dict{Symbol, Float64}
    sigma :: Vector{Float64} # parameters used to obtain the regularized problem (better vector?)
    sigma_factor :: Float64
    infiltrate_state :: Symbol
    lagrangian_atol :: Float64
    lagrangian_rtol :: Float64
    lagrangian_iteration_limit :: Int
    dual_initialization_regime :: Symbol
    lagrangian_method :: Symbol
    bundle_alpha :: Float64
    bundle_factor :: Float64
    level_factor :: Float64
    cut_selection :: Bool
    lag_status_regime :: Symbol
end

# Struct for initial algorithmic parameters that remain fixed and characterize a model run
struct InitialAlgoParams
    epsilon_innerLoop :: Float64
    binaryPrecision :: Dict{Symbol, Float64}
    plaPrecision :: Array{Vector{Float64},1}
    sigma :: Vector{Float64}
    sigma_factor :: Float64
    lagrangian_atol :: Float64
    lagrangian_rtol :: Float64
    lagrangian_iteration_limit :: Int
    dual_initialization_regime :: Symbol
    lagrangian_method :: Symbol
    bundle_alpha :: Float64
    bundle_factor :: Float64
    level_factor :: Float64
    cut_selection :: Bool
    lag_status_regime :: Symbol
end

# struct for solvers to be used (maybe mutable)
struct AppliedSolvers
    LP :: Any
    MILP :: Any
    NLP :: Any
    Lagrange :: Any
end

# Struct to store information on a nonlinear cut
mutable struct NonlinearCut
    intercept::Float64 # intercept of the cut (Lagrangian function value)
    coefficients::Dict{Symbol,Float64} # optimal dual variables in binary space
    trial_state::Dict{Symbol,Float64} # point at which cut should have been created
    anchor_state::Dict{Symbol,Float64} # point at which this cut was created
    binary_state::Dict{Symbol,BinaryState} # point in binary space where cut was created
    binary_precision::Dict{Symbol,Float64} # binary precision at moment of creation
    sigma::Float64
    cutVariables::Vector{JuMP.VariableRef}
    cutConstraints::Vector{JuMP.ConstraintRef}
    obj_y::Union{Nothing,NTuple{N,Float64} where {N}} # SDDP
    belief_y::Union{Nothing,Dict{T,Float64} where {T}} # SDDP
    non_dominated_count::Int # SDDP
    iteration::Int64
end
    # TODO:
    # 1) So far, binary precision is the same for all components.
    # 2) Should we also store the expression of this cut in binary space?

# struct for inner loop iteration results
mutable struct InnerLoopIterationResult{T,S}
    # pid
    lower_bound :: Float64
    upper_bound :: Float64 # should be renamed as cumulative_value as in SDDP if we solve stochastic problems
    current_sol :: Array{Dict{Symbol,Float64},1} #Vector{Dict{Symbol, Float64}} # current solution of state variables (also required for binary refinement)
    scenario_path :: Vector{Tuple{T,S}}
    has_converged :: Bool
    status :: Symbol # solution status (i.e. number of iterations)
    nonlinearCuts :: Dict{T, Vector{Any}} # only required for logging, binary explanation
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
    #TODO: We could also store sigma and binary precision here possibly
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
