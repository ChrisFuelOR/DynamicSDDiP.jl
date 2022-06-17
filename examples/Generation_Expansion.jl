using SDDP
import GLPK
import Test
import Gurobi

const GRB_ENV_2 = Gurobi.Env()

function generation_expansion(duality_handler)
    build_cost = 1e4
    use_cost = 4
    num_units = 5
    capacities = ones(num_units)
    demand_vals =
        0.5 * [
            5 5 5 5 5 5 5 5
            4 3 1 3 0 9 8 17
            0 9 4 2 19 19 13 7
            25 11 4 14 4 6 15 12
            6 7 5 3 8 4 17 13
        ]
    # Cost of unmet demand
    penalty = 5e5
    # Discounting rate
    rho = 0.99
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        lower_bound = 0.0,
        optimizer = () -> Gurobi.Optimizer(GRB_ENV_2),
    ) do sp, stage
        @variable(
            sp,
            0 <= invested[1:num_units] <= 1,
            SDDP.State,
            Int,
            initial_value = 0
        )
        @variables(sp, begin
            generation >= 0
            unmet >= 0
            demand
        end)

        @constraints(
            sp,
            begin
                # Can't un-invest
                investment[i in 1:num_units], invested[i].out >= invested[i].in
                # Generation capacity
                sum(capacities[i] * invested[i].out for i in 1:num_units) >=
                generation
                # Meet demand or pay a penalty
                unmet >= demand - sum(generation)
                # For fewer iterations order the units to break symmetry, units are identical (tougher numerically)
                [j in 1:(num_units-1)], invested[j].out <= invested[j+1].out
            end
        )
        # Demand is uncertain
        SDDP.parameterize(ω -> JuMP.fix(demand, ω), sp, demand_vals[stage, :])

        @expression(
            sp,
            investment_cost,
            build_cost *
            sum(invested[i].out - invested[i].in for i in 1:num_units)
        )
        @stageobjective(
            sp,
            (investment_cost + generation * use_cost) * rho^(stage - 1) +
            penalty * unmet
        )

        JuMP.set_silent(sp)
    end
    if get(ARGS, 1, "") == "--write"
        # Run `$ julia generation_expansion.jl --write` to update the benchmark
        # model directory
        model_dir = joinpath(@__DIR__, "..", "..", "..", "benchmarks", "models")
        SDDP.write_to_file(
            model,
            joinpath(model_dir, "generation_expansion.sof.json.gz");
            test_scenarios = 100,
        )
        exit(0)
    end

    # det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer, time_limit = 600.0)
    # JuMP.set_silent(det_equiv)
    # JuMP.optimize!(det_equiv)
    # print(JuMP.objective_value(det_equiv))

    SDDP.train(
        model,
        iteration_limit = 20,
        log_frequency = 1,
        duality_handler = duality_handler,
        print_level = 2,
    )
    Test.@test SDDP.calculate_bound(model) ≈ 2.078860e6 atol = 1e3
    return
end

#generation_expansion(SDDP.ContinuousConicDuality())
generation_expansion(SDDP.LagrangianDuality(atol=1e-4,rtol=1e-4))
#generation_expansion(SDDP.StrengthenedConicDuality())
