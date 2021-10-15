using SDDP

abstract type Test end

print(Test)
print(typeof(Test))

mutable struct Subtest <: Test end

print(Subtest)
print(typeof(Subtest))
print(typeof(Subtest) == Subtest)
print(Subtest == Subtest)
lim = SDDP.IterationLimit
print(lim)
print(typeof(lim))
print(lim == SDDP.IterationLimit)


function return_test()

    obj = 1000.0
    dual_bound = Inf

    return (
        obj,
        dual_bound
    )

end

return_results = return_test()

typeof(return_results)






using JuMP
using GLPK
using GAMS
approx_model = JuMP.Model(GAMS.Optimizer)
set_optimizer(approx_model, optimizer_with_attributes(
    GAMS.Optimizer,
    "Solver"=>"Gurobi",
    "optcr"=>0.0,
    )
)
ν = JuMP.@variable(approx_model, [1:10], lower_bound=0)
typeof(ν)

@variable(approx_model, t <= 100)
@objective(approx_model, Max, t)

JuMP.unset_silent(approx_model)
JuMP.optimize!(approx_model)
JuMP.set_silent(approx_model)
JuMP.optimize!(approx_model)


# Create the dual variables
# Note that the real dual multipliers are split up into two non-negative
# variables here, which is required for the Magnanti Wong part later
@variable(approx_model, π⁺[1:5] >= 0)
approx_model[:π⁺]




function get_dual_solution(node::Node, lagrange::LagrangianDuality)
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal
    # duals.
    @assert JuMP.termination_status(node.subproblem) == MOI.OPTIMAL
    # Query the current MIP solution  here. This is used as a bound for the
    # cutting plane method.
    primal_obj = JuMP.objective_value(node.subproblem)
    # A sign bit that is used to avoid if-statements in the models.
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1
    # Storage for the cutting plane method.
    num_states = length(node.states)
    x_in_value = zeros(num_states)               # The original value of x.
    λ_k = zeros(num_states)                      # The current estimate for λ
    λ_star = zeros(num_states)                   # The best estimate for λ
    # The best estimate for the dual objective value, ignoring optimization
    # sense (bigger is better).
    L_star = -Inf
    h_expr = Vector{AffExpr}(undef, num_states)  # The expression for x̄ - x
    h_k = zeros(num_states)                      # The value of x̄_k - x
    # Start by relaxing the fishing constraint.
    for (i, (_, state)) in enumerate(node.states)
        # We're going to need this value later when we reset things.
        x_in_value[i] = JuMP.fix_value(state.in)
        h_expr[i] = @expression(node.subproblem, state.in - x_in_value[i])
        # Relax the constraint from the problem.
        JuMP.unfix(state.in)
        # We need bounds to ensure that the dual problem is feasible. However,
        # they can't be too tight. Let's use 1e9 as a default...
        lb = has_lower_bound(state.out) ? lower_bound(state.out) : -1e9
        ub = has_upper_bound(state.out) ? upper_bound(state.out) : 1e9
        JuMP.set_lower_bound(state.in, lb)
        JuMP.set_upper_bound(state.in, ub)
    end
    # Create the model for the cutting plane algorithm
    model = JuMP.Model(something(lagrange.optimizer, node.optimizer))
    @variable(model, λ⁺[1:num_states] >= 0)
    @variable(model, λ⁻[1:num_states] >= 0)
    @variable(model, t <= s * primal_obj)
    @expression(model, λ, λ⁺ .- λ⁻)
    @objective(model, Max, t)
    # Step 1: find an optimal dual solution and corresponding objective value.
    iter, t_k = 0, s * primal_obj
    while !isapprox(L_star, t_k, atol = lagrange.atol, rtol = lagrange.rtol)
        iter += 1
        if iter > lagrange.iteration_limit
            error("Iteration limit exceeded in Lagrangian subproblem.")
        end
        JuMP.optimize!(model)
        @assert JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(model)
        λ_k .= value.(λ)
        L_k = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
        JuMP.@constraint(model, t <= s * (L_k + h_k' * (λ .- λ_k)))
        if s * L_k >= L_star
            L_star = s * L_k
            λ_star .= λ_k
        end
    end
    # Step 2: given the optimal dual objective value, try to find the optimal
    # dual solution with the smallest L1-norm ‖λ‖.
    @objective(model, Min, sum(λ⁺) + sum(λ⁻))
    set_lower_bound(t, t_k)
    # The worst-case scenario in this for-loop is that we run through the
    # iterations without finding a new dual solution. However if that happens
    # we can just keep our current λ_star.
    for _ in (iter+1):lagrange.iteration_limit
        JuMP.optimize!(model)
        @assert JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        λ_k .= value.(λ)
        L_k = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
        if isapprox(L_star, L_k, atol = lagrange.atol, rtol = lagrange.rtol)
            # At this point we tried the smallest ‖λ‖ from the cutting plane
            # problem, and it returned the optimal dual objective value. No
            # other optimal dual vector can have a smaller norm.
            λ_star = λ_k
            break
        end
        JuMP.@constraint(model, t <= s * (L_k + h_k' * (λ .- λ_k)))
    end
    # Restore the fishing constraint x.in == x_in_value
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end
    λ_solution = Dict{Symbol,Float64}(
        name => λ_star[i] for (i, name) in enumerate(keys(node.states))
    )
    # Remember to correct the sign of the optimal dual objective value.
    return s * L_star, λ_solution
end


function _solve_primal_problem(
    model::JuMP.Model,
    λ::Vector{Float64},
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    h_k::Vector{Float64},
)
    primal_obj = JuMP.objective_function(model)
    JuMP.set_objective_function(
        model,
        @expression(model, primal_obj - λ' * h_expr),
    )
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == MOI.OPTIMAL
    h_k .= -JuMP.value.(h_expr)
    L_λ = JuMP.objective_value(model)
    JuMP.set_objective_function(model, primal_obj)
    return L_λ
end


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

v = [DeterministicStopping(), 4, 3.5]

DeterministicStopping in v

v[1] == DeterministicStopping()

isequal(v[1], DeterministicStopping)

typeof(DeterministicStopping)

a = v[1]

a

a == DeterministicStopping

isa(a, DeterministicStopping)

isa(, DeterministicStopping)



contains
