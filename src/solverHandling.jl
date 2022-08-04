# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Function which assigns the correct solvers to be used in the specific part
of DynamicSDDiP, given that we use GAMS.jl as a framework for solver communication.
"""
function set_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
    solver_approach::DynamicSDDiP.GAMS_Solver,
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
        solver = applied_solvers.MINLP
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
            () -> GAMS.Optimizer(),
            "Solver"=>solver,
            "optcr"=>1e-4,
            "numericalemphasis"=>numerical_focus)
            )
    elseif solver == "Gurobi"
        JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
            () -> GAMS.Optimizer(),
            "Solver"=>solver,
            "optcr"=>1e-4,
            "NumericFocus"=>numerical_focus)
            )
    else
    JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
        () -> GAMS.Optimizer(),
        "Solver"=>solver,
        "optcr"=>1e-4)
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


"""
Function which assigns the correct solvers to be used in the specific part
of DynamicSDDiP, given that we use solvers directly in JuMP.
"""
function set_solver!(
    subproblem::JuMP.Model,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    algorithmic_step::Symbol,
    solver_approach::DynamicSDDiP.Direct_Solver,
)

    model = subproblem.ext[:sddp_policy_graph]
    iteration = model.ext[:iteration]

    # CHOOSE THE CORRECT TYPE OF SOLVER
    ########################################################################
    # if algorithmic_step in [:level_bundle]
    #     solver = applied_solvers.NLP
    # elseif algorithmic_step in [:LP_relax, :kelley, :cut_selection]
    #     # TODO: What about nonlinear cut projections here?
    #     solver = applied_solvers.LP
    # elseif algorithmic_step in [:l₂]
    #     solver = applied_solvers.MINLP
    # else
    #     # For all other algorithmic steps, we have to consider the projection
    #     # of the non-convex cuts if BinaryApproximation is used.
    #     # We loop over all cut_generation_regimes to see if such non-convex
    #     # cuts have been created so far.
    #
    #     for cut_generation_regime in algo_params.cut_generation_regimes
    #         # NON-CONVEX CUTS HAVE BEEN CREATED SO FAR
    #         if isa(cut_generation_regime.state_approximation_regime, DynamicSDDiP.BinaryApproximation) && cut_generation_regime.iteration_to_start <= iteration
    #
    #             cut_projection_regime = cut_generation_regime.state_approximation_regime.cut_projection_regime
    #
    #             if algorithmic_step in [:forward_pass, :backward_pass]
    #                 if isa(cut_projection_regime, DynamicSDDiP.KKT)
    #                     solver = applied_solvers.MINLP
    #                 elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
    #                     solver = applied_solvers.MINLP
    #                     #solver = applied_solvers.MIQCP
    #                 else
    #                     solver = applied_solvers.MILP
    #                 end
    #             elseif algorithmic_step in [:lagrange_relax]
    #                 if isa(cut_projection_regime, DynamicSDDiP.KKT)
    #                     solver = applied_solvers.MINLP
    #                 elseif isa(cut_projection_regime, DynamicSDDiP.StrongDuality)
    #                     solver = applied_solvers.MINLP
    #                     #solver = applied_solvers.MIQCP
    #                 else
    #                     solver = applied_solvers.Lagrange
    #                 end
    #             end
    #
    #         # NO NON-CONVEX CUTS HAVE BEEN CREATED SO FAR
    #         else
    #             if algorithmic_step in [:forward_pass, :backward_pass]
    #                 solver = applied_solvers.MILP
    #             elseif algorithmic_step in [:lagrange_relax]
    #                 solver = applied_solvers.Lagrange
    #             end
    #         end
    #     end
    # end
    #
    # # CHECK NUMERICAL FOCUS
    # ############################################################################
    # if algo_params.numerical_focus
    #     numerical_focus = 1
    # else
    #     numerical_focus = 0
    # end
    #
    # # SET THE CORRECT SOLVER WITH THE REQUIRED PROPERTIES
    # ############################################################################
    # if solver == "CPLEX" || solver == "BARON"
    #     error("Solver can only be used with our GAMS license")
    # elseif solver == "Gurobi"
    #     JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
    #         () -> Gurobi.Optimizer(GRB_ENV[]),
    #         "MIPGap"=>1e-4,
    #         "NumericFocus"=>numerical_focus)
    #         )
    # elseif solver == "SCIP"
    #     JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
    #         SCIP.Optimizer(),
    #         "display_verblevel"=>0,
    #         "limits_gap"=>0.05)
    #         )
    # elseif solver == "GLPK"
    #     JuMP.set_optimizer(subproblem, JuMP.optimizer_with_attributes(
    #         GLPK.Optimizer())
    #         )
    # else
    #     error("Solver has to set-up first.")
    # end
    #
    # if numerical_focus == 1 && !(solver == "Gurobi" || solver == "CPLEX")
    #     @warn("Numerical focus only works with Gurobi or CPLEX.")
    # end
    # 
    # if algo_params.silent
    #     JuMP.set_silent(subproblem)
    # else
    #     JuMP.unset_silent(subproblem)
    # end

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


""" A possible alternative is to attempt a solution using the numerical_focus
parameter. However, from my experience this does not help in most cases.
See previously used code below.
"""

# TRY RECOVERING FROM NUMERICAL ISSUES USING numerical_focus
########################################################################
# if (JuMP.termination_status(approx_model) != MOI.OPTIMAL) && !algo_params.numerical_focus
#     algo_params.numerical_focus = true
#     set_solver!(node.subproblem, algo_params, applied_solvers, :lagrange_relax, algo_params.solver_approach)
#     JuMP.optimize!(approx_model)
#     @assert JuMP.termination_status(approx_model) == MOI.OPTIMAL
#     algo_params.numerical_focus = false
# else (JuMP.termination_status(approx_model) != MOI.OPTIMAL) && algo_params.numerical_focus
#     @assert JuMP.termination_status(approx_model) == MOI.OPTIMAL
# end

""" We could also use the SDDP numerical_recovery approach, see below.
However, from my experience this does not help in most of our cases. """

#MOI.Utilities.reset_optimizer(subproblem)
#JuMP.optimize!(subproblem)
