```@meta
EditURL = "cuts.jl"
```

# How cuts are represented in the code

In the code, the required information to handle cuts is stored in structs of type `Cut`, as known from SDDP.jl. However, the struct is slightly adjusted to satisfy the special form of the cuts above.

````@example cuts
using JuMP
using DynamicSDDiP

mutable struct LinearCut <: DynamicSDDiP.Cut
    intercept::Float64
    coefficients::Dict{Symbol,Float64}
    scaling_coeff::Float64
    trial_state::Dict{Symbol,Float64}
    sigma::Union{Nothing,Float64}
    cut_constraint::Union{Nothing,JuMP.ConstraintRef}
    non_dominated_count::Int
    iteration::Int64
    aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime
    duality_regime::DynamicSDDiP.AbstractDualityRegime
end
````

Its fields are defined as follows:

 * `intercept`: The constant intercept of the cut.
 * `coefficients`: This is the cut gradient vector. It can be computed using the values of certain dual multipliers.
 * `scaling_coeff`: This is the coefficient for the special dual multiplier $\pi_{n0}$ that is introduced in our paper. For standard Benders and Lagrangian cuts it is always 1.
 * `trial_state`: The incumbent where the cut is constructed.
 * `sigma`: The regularization parameter $\sigma_n$ at cut generation. This is only relevant if a Lipschitz regularization is applied, see [Setting algorithmic parameters](params.md).
 * `cut_constraint`: Reference to the JuMP constraint of the cut.
 * `non_dominated_count`: Number that is iteratively updated if `CutSelection` is used.
 * `iteration`: Stores the iteration in which the cut was generated. This is just for logging and debugging purposes.
 * `aggregation_regime`: Stores if this cut was generated using a `SINGLE_CUT` or `MULTI_CUT` approach. This is just for logging and debugging purposes.
 * `duality_regime`: Stores which type of cut this is, i.e. Benders cut, strengthened Benders cut, Lagrangian cut. For instance, this is relevant when `BendersSpanSpaceRestriction` is used.

## The special case of non-convex cuts

As discussed in [another preprint](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/), DynamicSDDiP.jl also allows for the generation of non-convex cuts. For these cuts, a separate type of struct is used because more information needs to be stored.

````@example cuts
using JuMP
using DynamicSDDiP

mutable struct NonlinearCut <: DynamicSDDiP.Cut
    intercept::Float64
    coefficients::Dict{Symbol,Float64}
    scaling_coeff::Float64
    trial_state::Dict{Symbol,Float64}
    anchor_state::Dict{Symbol,Float64}
    binary_state::Dict{Symbol,DynamicSDDiP.BinaryState}
    binary_precision::Dict{Symbol,Float64}
    sigma::Union{Nothing,Float64}
    cut_variables::Vector{JuMP.VariableRef}
    cut_constraints::Vector{JuMP.ConstraintRef}
    non_dominated_count::Int
    iteration::Int64
    aggregation_regime::DynamicSDDiP.AbstractCutAggregationRegime
    duality_regime::DynamicSDDiP.AbstractDualityRegime
end
````

Compared to `LinearCut` the additional or differing fields are:

 * `anchor_state`: This is the state at which the non-convex cut is constructed. By design, this may differ from the initial incumbent `trial_point` due to the binary approximation of the state space.
 * `binary_state`: This is the equivalent to the `anchor_state` but in the lifted binary space. In other words, this is the exact binary representation of `anchor_state`.
 * `binary_precision`: A dict of the precision that is used for each state component for the binary approximation.
 * `cut_variables`: The mixed-integer linear representation of the non-convex cut requires to introduce several auxiliary variables. The references are stored in this vector in order to keep track of which cut they belong to.
 * `cut_constraints`: The mixed-integer linear representation of the non-convex cut requires to introduce several (auxiliary) constraints. The references are stored in this vector in order to keep track of which cut they belong to.

!!! note "Remark"
    The dict `coefficients` will have a different dimension than for a `LinearCut`, as the coefficients refer to dual multipliers in the lifted binary space.

!!! note "Remark"
    Note that being able to deal with linear cuts (`LinearCut`) and non-convex cuts (`NonlinearCut`) required a lot of adjustments and re-definitions in the file `bellman.jl` compared to SDDP.jl.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

