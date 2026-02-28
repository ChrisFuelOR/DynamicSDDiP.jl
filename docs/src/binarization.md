```@meta
EditURL = "binarization.jl"
```

# Binarization and non-convex cuts

When solving a multistage stochastic mixed-integer linear program with our dynamic version of SDDiP, we have three options considering the state space:

 1. We consider the original state space. If it is not purely binary, then none of the cuts that are generated (see [Configuring the cut generation](cut_generation.md)) are guaranteed to be tight, so we lose any convergence guarantees.
 2. We manually approximate the state space with binary variables before running SDDiP (**static binarization**). This increases the state space, but allows for the generation of tight cuts. This is what is proposed in the original [SDDiP paper](https://link.springer.com/article/10.1007/s10107-018-1249-5).
 3. We apply a **temporary binary approximation** in the cut generation process. The cuts are then generated according to [Configuring the cut generation](cut_generation.md) but in a lifted space. The cuts are then projected back to the original state space, so that in the forward pass of SDDiP still the original state space can be considered. The binary approximation can be dynaimcally improved if it is required (this also explains the name of our implementation as **dynamic** SDDiP).

In our experiments for this paper, we only used linear cuts, and thus applied variants 1.) and 2.). If variant 2.) is used, then the state binarization has to be included in the model definition (see TODO). In both cases, the `state_approximation_regime` is set to `NoStateApproximation`.

Variant 3.) is implemented, as this repository was used for the experiments for [paper](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/)  as well. If variant 3.) is used, then `state_approximation_regime` is set to `BinaryApproximation` with some additional parameters (see below). Note that in this case, due to the projection to the original state space, the linear cuts created in the lifted space lead to non-convex approximations. In the linked paper, we call these non-convex cuts **cut projection closure (CPC)**.

## Binary approximation

The temporary and dynamic binary approximation mentioned in 3.) above can be controled using a `state_approximation_regime` in `AlgoParams`, see [Setting algorithmic parameters](params.md).

We have two options:

 * `BinaryApproximation`: This means that a temporary binary approximation is applied in the cut generation process. The parameter `binary_precision` can be used to initialize the precision of the approximation, which is later refined if needed. The parameter `cut_projection_regime` controls how the cuts are projected to the original state space.
 * `NoStateApproximation`: This means that cuts are generated in the original state space and either approach 1.) or 2.) from above is applied.

````@example binarization
using DynamicSDDiP
mutable struct BinaryApproximation <: DynamicSDDiP.AbstractStateApproximationRegime
    binary_precision::Dict{Symbol, Float64}
    cut_projection_regime::DynamicSDDiP.AbstractCutProjectionRegime
end

mutable struct NoStateApproximation <: DynamicSDDiP.AbstractStateApproximationRegime end
````

Note that so far, the implementation is done in such a way that `binary_precision` is the same for all stages.

Importantly, the binarization is applied to each state variable (i.e. each component of $x_{a(n)}^i$) separately.

The temporary binarization of the state space requires to relate the new binar variables to the original state variables. To this end, the `BinaryState` struct is introduced in the code.

````@example binarization
struct BinaryState
    value::Float64
    x_name::Symbol
    k::Int64
end
````

The field `value` allows to store the value of the corresponding state variable in the original space, i.e. a component of $x_{a(n)}^i$. This is required to restore the value later. The field `x_name` stores the name of the original state variable which the binary variable corresponds to. The field `k` stores the dimension of the corresponding binary variable. It depends on the `binary_precision`.

## Cut projection

The cut projection approach can be controled using parameter `cut_projection_regime`.

We have two options. Both are based on representing the CPC by KKT conditions.

 * `BigM`: This means that the complementarity constraints in the KKT conditions are reformulated using a Big-M approach. If a Lipschitz regularization is used, a natural bound for the Big-M parameters can be derived, see [paper](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/).
 * `SOS1`: This means that the complementarity constraints in the KKT conditions are reformulated using SOS-1 constraints.

````@example binarization
using DynamicSDDiP
mutable struct BigM <: DynamicSDDiP.AbstractCutProjectionRegime end
mutable struct SOS1 <: DynamicSDDiP.AbstractCutProjectionRegime end
````

!!! note "Remark"
    Note that for the generation of non-convex cuts with convergence guarantees, a Lipschitz regularization should be applied. For more details see [Setting algorithmic parameters](params.md).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

