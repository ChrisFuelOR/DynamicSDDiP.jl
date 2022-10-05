using GAMS
using JuMP

function test()

    model = JuMP.Model()
    eps = 0
    t = 2/3+1e-12

    JuMP.@variable(model, x)
    JuMP.@NLconstraint(model, x^3/3 - t*x^2 + (2*t-1)*x <= 0)
    JuMP.@objective(model, Min, -x)

    boxes_list = [[-2,0], [0,2], [-2,-1], [-1,0], [0,1], [1,2]]

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
