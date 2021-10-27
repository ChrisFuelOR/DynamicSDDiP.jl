# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Function which assigns the correct solvers to be used in the specific part
of DynamicSDDiP
"""
function set_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
)

    cut_projection_regime = algo_params.state_approximation_regime.cut_projection_regime

    # CHOOSE THE CORRECT TYPE OF SOLVER AND SOLVER
    ############################################################################
    if algorithmic_step in [:forward_pass, :backward_pass]
        if isa(cut_projection_regime, DynamicSDDiP.KKT)
            solver = applied_solvers.MINLP
        elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
            solver = applied_solvers.MIQCP
        else
            solver = applied_solvers.MILP
        end
    elseif algorithmic_step in [:lagrange_relax]
        if isa(cut_projection_regime, DynamicSDDiP.KKT)
            solver = applied_solvers.MINLP
        elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
            solver = applied_solvers.MIQCP
        else
            solver = applied_solvers.Lagrange
        end
    elseif algorithmic_step in [:level_bundle]
        solver = applied_solvers.NLP
    elseif algorithmic_step in [:LP_relax, :kelley, :cut_selection]
        # TODO: What about nonlinear cut projections here?
        solver = applied_solvers.LP
    end

    # CHECK NUMERICAL FOCUS
    ############################################################################
    if algo_params.numerical_focus
        numerical_focus = 1
    else
        numerical_focus = 0
    end

    # SET THE CORRECT SOLVER WITH THE REQUIRED PROPERTIES
    ############################################################################
    if solver == "CPLEX"
        set_optimizer(subproblem, optimizer_with_attributes(
            GAMS.Optimizer,
            "Solver"=>solver,
            "optcr"=>0.0,
            "numericalemphasis"=>numerical_focus)
            )
    elseif solver == "Gurobi" && isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
        set_optimizer(subproblem, optimizer_with_attributes(
            GAMS.Optimizer,
            "Solver"=>solver,
            "optcr"=>0.0,
            "NumericFocus"=>numerical_focus,
            "nonConvex"=>2)
            )
    elseif solver == "Gurobi"
        set_optimizer(subproblem, optimizer_with_attributes(
            GAMS.Optimizer,
            "Solver"=>solver,
            "optcr"=>0.0,
            "NumericFocus"=>numerical_focus)
            )
    else
    set_optimizer(subproblem, optimizer_with_attributes(
        GAMS.Optimizer,
        "Solver"=>solver,
        "optcr"=>0.0)
        )

        if numerical_focus == 1
            @warn("Numerical focus only works with Gurobi or CPLEX.")
        end

    end

    if algo_params.silent
        JuMP.set_silent(subproblem)
    else
        JuMP.unset_silent(subproblem)
    end

    return
end
