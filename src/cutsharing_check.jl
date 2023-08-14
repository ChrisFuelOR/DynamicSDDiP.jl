using JuMP
using Infiltrator

"""
Function which constructs a symbolic expression for a part of the cut intercept
given some arguments.

T:  number of stages
p:  number of lags for the autoregressive process
t:  current stage for which a cut should be generated
τ:  part of the intercept which should be generated (later on we can iterate
    over τ to generate all parts)
"""

function cut_formula_tester(
    T::Int,
    p::Int,
    t::Int,
    τ::Int,
)

    if p == 0
        cut_formula_tester(T,t,τ)
    else

        model = JuMP.Model()
        JuMP.@variable(model, ξ[1:T])
        JuMP.@variable(model, ϕ[1:T])
        expr = JuMP.@NLexpression(model, 1)

        @assert(τ >= t)
        @assert(t > 1)

        # Iterate over number of factors
        for lag in 1:p
            # Current index (also refers to the first factor in the exponent)
            k = t - lag

            if k >= 1
                println("t: ", t ," , τ: ", τ, ", p: ", p, ", k: ", k)

                # Compute the basis exponent
                base = JuMP.@NLexpression(model, ϕ[k+1] * prod(ϕ[ℓ] for ℓ=t+1:τ))

                # Compute the subtrahend
                sub = JuMP.@NLexpression(model, prod(ϕ[ℓ] for ℓ=t+1:k+p))

                # Compute the exponent
                #@assert(k > t-p)

                if k > (τ-1)-p || τ == t
                    term = JuMP.@NLexpression(model, base)
                elseif k == t-p
                    term = JuMP.@NLexpression(model, base - 1)
                else # k in (t-p, (τ-1)-p]
                    term = JuMP.@NLexpression(model, base - sub)
                end

                # Update expression
                expr = JuMP.@NLexpression(model, expr * ξ[k]^term)

                Infiltrator.@infiltrate

            end
        end

        JuMP._nl_subexpression_string(JuMP.REPLMode, model)[expr.index]

    end

end

function cut_formula_tester_one(
    T::Int,
    t::Int,
    τ::Int,
)

    model = JuMP.Model()
    JuMP.@variable(model, ξ[1:T])
    JuMP.@variable(model, ϕ[1:T])

    @assert(τ >= t)
    @assert(t > 1)

    if τ == t
        expr = JuMP.@NLexpression(model, ξ[t-1]^(ϕ[t]))
    else
        expr = JuMP.@NLexpression(model, ξ[t-1]^(prod(ϕ[k] for k in t:τ)))
    end

    JuMP._nl_subexpression_string(JuMP.REPLMode, model)[expr.index]

end

# cut_formula_tester_one(4,4,4)
# cut_formula_tester(4,3,3)
# cut_formula_tester(4,3,4) #ERROR
# cut_formula_tester(4,2,2)
# cut_formula_tester(4,2,3)
# cut_formula_tester(4,2,4)

cut_formula_tester(4,1,2,3)
