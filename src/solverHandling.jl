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

    model = subproblem.ext[:sddp_policy_graph]
    iteration = model.ext[:iteration]

    # CHOOSE THE CORRECT TYPE OF SOLVER
    ########################################################################
    if algorithmic_step in [:level_bundle]
        solver = applied_solvers.NLP
    elseif algorithmic_step in [:LP_relax, :kelley, :cut_selection]
        # TODO: What about nonlinear cut projections here?
        solver = applied_solvers.LP
    elseif algorithmic_step in [:l₂]
        solver = "scip"
    else
        # For all other algorithmic steps, we have to consider the projection
        # of the non-convex cuts if BinaryApproximation is used.
        # We loop over all cut_generation_regimes to see if such non-convex
        # cuts have been created so far.

        for cut_generation_regime in algo_params.cut_generation_regimes
            # NON-CONVEX CUTS HAVE BEEN CREATED SO FAR
            if isa(cut_generation_regime.state_approximation_regime, DynamicSDDiP.BinaryApproximation) && cut_generation_regime.iteration_to_start <= iteration

                cut_projection_regime = cut_generation_regime.state_approximation_regime.cut_projection_regime

                if algorithmic_step in [:forward_pass, :backward_pass]
                    if isa(cut_projection_regime, DynamicSDDiP.KKT)
                        solver = applied_solvers.MINLP
                    elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
                        solver = applied_solvers.MINLP
                        #solver = applied_solvers.MIQCP
                    else
                        solver = applied_solvers.MILP
                    end
                elseif algorithmic_step in [:lagrange_relax]
                    if isa(cut_projection_regime, DynamicSDDiP.KKT)
                        solver = applied_solvers.MINLP
                    elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
                        solver = applied_solvers.MINLP
                        #solver = applied_solvers.MIQCP
                    else
                        solver = applied_solvers.Lagrange
                    end
                end

            # NO NON-CONVEX CUTS HAVE BEEN CREATED SO FAR
            else
                if algorithmic_step in [:forward_pass, :backward_pass]
                    solver = applied_solvers.MILP
                elseif algorithmic_step in [:lagrange_relax]
                    solver = applied_solvers.Lagrange
                end
            end
        end
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
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            GAMS.Optimizer,
            "Solver"=>solver,
            "optcr"=>0.0,
            "numericalemphasis"=>numerical_focus)
            )
    elseif solver == "Gurobi"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            GAMS.Optimizer,
            "Solver"=>solver,
            "optcr"=>0.0,
            "NumericFocus"=>numerical_focus)
            )
    else
    JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
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
