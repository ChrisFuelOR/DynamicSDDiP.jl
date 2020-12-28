using GAMS
using JuMP
using Gurobi
using MathOptInterface

m1 = JuMP.Model(GAMS.Optimizer)
JuMP.set_optimizer_attribute(m1, "Solver", "Gurobi")
#JuMP.set_silent(m1)
JuMP.@variable(m1, x>=0)
JuMP.@objective(m1, MOI.MIN_SENSE, x)
JuMP.@constraint(m1, x >= 2257.812325)
JuMP.optimize!(m1)

m2 = JuMP.Model(Gurobi.Optimizer)
#JuMP.set_silent(m2)
JuMP.@variable(m2, y>=0)
JuMP.@objective(m2, MOI.MIN_SENSE, y)
JuMP.@constraint(m2, y >= 2257.812325)

JuMP.optimize!(m2)

println()
println("Optimal point calling Gurobi via GAMS:")
println(JuMP.value(x))
println("Optimal point calling Gurobi directly:")
println(JuMP.value(y))
