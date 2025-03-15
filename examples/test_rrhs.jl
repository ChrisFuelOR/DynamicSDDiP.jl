using GAMS
using JuMP

function test()

    model = JuMP.Model()
    eps = 0
    a = 12/(6+eps)^2

    JuMP.@variable(model, x)
    JuMP.@constraint(model, -a*x^2 + 12 <= 0)
    JuMP.@objective(model, Min, x)

    boxes_list = [[0,8]]#, [0,4], [4,8], [4,6], [5,6], [5.5,6], [5.75,6], [5.875,6], [5.9375,6]]

    for solver in ["minos", "ipopt", "knitro", "conopt", "snopt"]

        JuMP.set_optimizer(model, JuMP.optimizer_with_attributes(() -> GAMS.Optimizer(),"Solver"=>solver))
        JuMP.set_silent(model)

        for box in boxes_list
            JuMP.set_lower_bound(x, box[1])
            JuMP.set_upper_bound(x, box[2])

            JuMP.optimize!(model)

            println(solver)
            println(box)
            println(JuMP.termination_status(model))

            if JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
                println(JuMP.value(x))
                println(JuMP.objective_value(model))
            end
            println("----------------------------------------------------------")

        end
    end

end

test()
