# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Function which identifies the correct solver for a given solution step.
"""
function identify_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
)

    solver = nothing

    if algorithmic_step in [:main_subproblem, :forward_pass, :backward_pass]
        solver = applied_solvers.solver_subproblem
    elseif algorithmic_step in [:lagrange_relax]
        if !isnothing(applied_solvers.solver_Lagrange_relax)
            solver = applied_solvers.solver_Lagrange_relax
        else
            solver = applied_solvers.solver_subproblem
        end
    elseif algorithmic_step in [:LP_relax]
        if !isnothing(applied_solvers.solver_LP_relax)
            solver = applied_solvers.solver_LP_relax
        else
            solver = applied_solvers.solver_subproblem
        end
    elseif algorithmic_step in [:reg]
        if !isnothing(applied_solvers.solver_reg)
            solver = applied_solvers.solver_reg
        else
            solver = applied_solvers.solver_subproblem
        end
    elseif algorithmic_step in [:norm]
        if !isnothing(applied_solvers.solver_norm)
            solver = applied_solvers.solver_norm
        else
            solver = applied_solvers.solver_subproblem
        end
    elseif algorithmic_step in [:kelley]
        solver = applied_solvers.solver_Lagrange_approx
    elseif algorithmic_step in [:level_bundle]
        if !isnothing(applied_solvers.solver_Lagrange_Bundle)
            solver = applied_solvers.solver_Lagrange_Bundle
        else
            solver = applied_solvers.solver_Lagrange_approx
        end
    elseif algorithmic_step in [:subgradient]
        solver = applied_solvers.solver_Lagrange_subgradient
    elseif algorithmic_step in [:cut_selection]
        solver = applied_solvers.solver_cut_selection
    end

    # TODO: If non-convex cuts have been created for which Strong Duality
    # or KKT cut_projection_regime is used, then we have to make sure that
    # a nonlinear solver is used for the subproblems from there on.

    return solver

end

"""
Function which assigns the correct solvers to be used given that we
use GAMS.jl as a framework for solver communication.
"""
function set_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    solver::String,
    solver_approach::DynamicSDDiP.GAMS_Solver,
)
    numerical_focus = algo_params.numerical_focus ? 1 : 0
    tolerance = applied_solvers.solver_tol

    # Set solver with tolerance
    ############################################################################
    if solver == "CPLEX"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            () -> GAMS.Optimizer(),
            "Solver"=>solver,
            "optcr"=>tolerance,
            "numericalemphasis"=>numerical_focus)
        )
    elseif solver == "Gurobi"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            () -> GAMS.Optimizer(),
            "Solver"=>solver,
            "optcr"=>tolerance,
            "NumericFocus"=>numerical_focus)
        )
    else
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            () -> GAMS.Optimizer(),
            "Solver"=>solver,
            "optcr"=>tolerance)
        )
    end

    # Numerical focus warning
    ############################################################################
    if numerical_focus == 1 && !(solver in ["Gurobi", "CPLEX"])
        @warn("Numerical focus only works with Gurobi or CPLEX.")
    end

    # Silence solver
    ############################################################################
    if algo_params.silent
        JuMP.set_silent(subproblem)
    else
        JuMP.unset_silent(subproblem)
    end

    return
end

"""
Function which assigns the correct solvers to be used given that we
use the solvers directly in JuMP.
"""
function set_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    solver::String,
    solver_approach::DynamicSDDiP.Direct_Solver,
)
    numerical_focus = algo_params.numerical_focus ? 1 : 0
    tolerance = applied_solvers.solver_tol
    time_limit = applied_solvers.solver_time

    # Set solver with tolerance
    ############################################################################
    if solver in ["CPLEX", "BARON"]
        error("Solver can only be used with our GAMS license")
    elseif solver == "Gurobi"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GURB_ENV[]),
            "MIPGap"=>tolerance,
            "TimeLimit"=>time_limit,
            "NumericFocus"=>numerical_focus)
            )
    elseif solver == "SCIP"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            SCIP.Optimizer(),
            "display_verblevel"=>0,
            "limits_gap"=>tolerance)
            )
    else
        error("Solver has to set-up first.")
    end

    # Numerical focus warning
    ############################################################################
    if numerical_focus == 1 && !(solver in ["Gurobi", "CPLEX"])
        @warn("Numerical focus only works with Gurobi or CPLEX.")
    end

    # Silence solver
    ############################################################################
    if algo_params.silent
        JuMP.set_silent(subproblem)
    else
        JuMP.unset_silent(subproblem)
    end

    return
end


"""
Function which sets the solver for the first time when a specific problem
is constructed.
"""
function set_solver_initially!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
    solver_approach::Union{DynamicSDDiP.GAMS_Solver,DynamicSDDiP.Direct_Solver}
)

    solver = identify_solver!(subproblem, algo_params, applied_solvers, algorithmic_step)
    set_solver!(subproblem, algo_params, applied_solvers, solver, solver_approach)

    return
end

"""
Function which re-sets the solver for a specific phase of the algorithm
if applied_solvers.general_solver is set to false.
"""
function reset_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
    solver_approach::Union{DynamicSDDiP.GAMS_Solver,DynamicSDDiP.Direct_Solver}
)

    if applied_solvers.general_solver
        return
    else
        solver = identify_solver!(subproblem, algo_params, applied_solvers, algorithmic_step)
        set_solver!(subproblem, algo_params, applied_solvers, solver, solver_approach)
    end

    return
end



"""
Sometimes using a lot of cutting-planes some numerical issues occur, and the
subproblem is no longer bounded or feasible. However, from my experience CPLEX
handles this way better than Gurobi, so that CPLEX may still yield a feasible
solution for the subproblem. On the contrary, using CPLEX by means of GAMS.jl
is way slower than directly using Gurobi. Therefore, we revert to CPLEX
in case that a subproblem is infeasible/unbounded only.
"""
function elude_numerical_issues!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
)

    numerical_focus = algo_params.numerical_focus ? 1 : 0
    JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
        () -> GAMS.Optimizer(),
        "Solver"=>"CPLEX",
        "optcr"=>1e-4,
        "numericalemphasis"=>numerical_focus)
        )
    JuMP.optimize!(subproblem)

    @assert JuMP.termination_status(subproblem) == MOI.OPTIMAL

    return
end
