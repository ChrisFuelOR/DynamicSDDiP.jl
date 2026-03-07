```@meta
EditURL = "index.jl"
```

## DynamicSDDiP.jl

This project implements a special version of **stochastic dual dynamic integer programming (SDDiP)** and is based on (an earlier version of) [SDDP.jl](https://github.com/odow/SDDP.jl) from Oscar Dowson.

As standard SDDiP, our code is tailored to solve multistage stochastic mixed-integer linear problems (MS-MILPs). However, in addition to standard cuts from the literature like Benders cuts, strengthened Benders cuts and Lagrangian cuts, our code allows to generate **deep Lagrangian cuts** (obtained by a norm-normalized dual problem) and **LN Lagrangian cuts** (obtained by a linearly normalized (LN) dual problem).
The theory behind these cuts as well as extensive computational experiments are presented in a [preprint](https://optimization-online.org/2024/08/a-new-framework-to-generate-lagrangian-cuts-in-multistage-stochastic-mixed-integer-programming/).

Note that this project also includes features like binary approximation of the state space, Lipschitz regularization and the generation of special non-convex cuts, obtained by projecting Lagrangian cuts from a lifted binary state space to the original state space (so-called **cut projection closure**). These topics are discussed in another [preprint](https://optimization-online.org/2024/08/on-lipschitz-regularization-and-lagrangian-cuts-in-multistage-stochastic-mixed-integer-linear-programming/).

````@example index
# License
````

`DynamicSDDiP.jl` is licensed under the [MPL 2.0 license](https://github.com/ChrisFuelOR/DynamicSDDiP.jl/blob/main/LICENSE).

````@example index
# Documentation
````

You can find the documentation at TODO.

````@example index
# Help
````

If you need help, please [open a GitHub issue](https://github.com/ChrisFuelOR/DynamicSDDiP.jl/issues).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

