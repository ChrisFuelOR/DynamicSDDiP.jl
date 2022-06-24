#  Copyright 2017-21, Oscar Dowson.                                     #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # SLDP: example 1

# This example is derived from Section 4.2 of the paper:
# Ahmed, S., Cabral, F. G., & da Costa, B. F. P. (2019). Stochastic Lipschitz
# Dynamic Programming. Optimization Online. [PDF](http://www.optimization-online.org/DB_FILE/2019/05/7193.pdf)

using SDDP, GLPK, Test, Gurobi
import JuMP

const GRB_ENV_2 = Gurobi.Env()

function sldp_example_one()
    model = SDDP.LinearPolicyGraph(
        stages = 8,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        JuMP.@variable(sp, x, SDDP.State, initial_value = 2.0, lower_bound = -5.0, upper_bound = 5.0)
        JuMP.@variables(sp, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)
        SDDP.@stageobjective(sp, 0.9^(t - 1) * (x⁺ + x⁻))
        JuMP.@constraints(sp, begin
            x.out == x.in + 2 * u - 1 + ω
            x⁺ >= x.out
            x⁻ >= -x.out
        end)
        points = [
            -0.3089653673606697,
            -0.2718277412744214,
            -0.09611178608243474,
            0.24645863921577763,
            0.5204224537256875,
        ]
        JuMP.set_silent(sp)
        return SDDP.parameterize(φ -> JuMP.fix(ω, φ), sp, [points; -points])
    end

    # det_equiv = SDDP.deterministic_equivalent(model, Gurobi.Optimizer, time_limit = 600.0)
    # JuMP.set_silent(det_equiv)
    # JuMP.optimize!(det_equiv)
    # print(JuMP.objective_value(det_equiv))

    duality_handler = SDDP.ContinuousConicDuality()
    #duality_handler = SDDP.StrengthenedConicDuality()
    #duality_handler = SDDP.LagrangianDuality(atol=1e-4,rtol=1e-4)

    SDDP.train(
        model,
        iteration_limit = 100,
        log_frequency = 1,
        duality_handler = duality_handler,
        print_level = 2,
    )

    Test.@test SDDP.calculate_bound(model) <= 1.1675
    return
end

Random.seed!(11111)
sldp_example_one()

Random.seed!(11112)
sldp_example_one()

Random.seed!(11113)
sldp_example_one()

Random.seed!(11114)
sldp_example_one()

Random.seed!(11115)
sldp_example_one()
