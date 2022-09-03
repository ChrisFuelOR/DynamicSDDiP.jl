"""
Solving aggregated dual in the backward pass.
"""
function solve_aggregated_dual(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    i::Int64,
    items::BackwardPassItems,
    belief::Float64,
    # belief_state,
    # objective_state,
    outgoing_state::Dict{Symbol,Float64},
    epi_state::Float64,
    primal_obj::Float64,
    add_cut_flag::Bool,
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    scenario_path,
    algo_params::DynamicSDDiP.AlgoParams,
    cut_generation_regime::DynamicSDDiP.CutGenerationRegime,
    applied_solvers::DynamicSDDiP.AppliedSolvers,
    dual_solution_regime::DynamicSDDiP.Kelley,
) where {T}

duality_regime = cut_generation_regime.duality_regime

############################################################################
# DO THIS FOR EACH FOLLOWING MARKOV STATE / STAGE
############################################################################
length_scenario_path = length(scenario_path)
for child in node.children
    if isapprox(child.probability, 0.0, atol = 1e-6)
        continue
    end
    child_node = model[child.term]

    ########################################################################
    # SOME INITIALIZATIONS
    ########################################################################
    # storages for return of dual values and binary state values (trial point)
    # note that with NoStateApproximation bin_state will just remain empty
    dual_values = Dict{Symbol,Float64}()
    bin_state = Dict{Symbol, BinaryState}()
    number_of_noise = length(SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node))
    number_of_states_per_noise = get_number_of_states(child_node, cut_generation_regime.state_approximation_regime)

    # storages for information on Lagrangian dual
    lag_obj = 0
    lag_iterations = 0
    lag_status = :none

    ########################################################################
    # INITIALIZE DUALS
    ########################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "dual_initialization" begin
        dual_vars_agg = zeros(number_of_states_per_noise)
        dual_0_var = 1.0 # (0 is not suitable for relatively complete recourse)
    end

    ########################################################################
    # GET BOUNDS FOR LAGRANGIAN DUAL
    ########################################################################
    bound_results = get_dual_bounds(child_node, node_index+1, algo_params, primal_obj, duality_regime.dual_bound_regime)
    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    ############################################################################
    # GET CORE POINT RESP. NORMALIZATION COEFFICIENTS
    ############################################################################
    """ If we do not use a normalization that requires a core point, we simply
    return nothing """
    normalization_coeff = get_normalization_coefficients_agg(
        node,
        number_of_states_per_noise,
        number_of_noise,
        epi_state,
        algo_params,
        backward_sampling_scheme,
        cut_generation_regime.state_approximation_regime,
        duality_regime.normalization_regime,
        duality_regime.copy_regime
    )

    try
        ####################################################################
        # SOLVE THE AGGREGATED LAGRANGIAN DUAL PROBLEM
        ####################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
            results = solve_aggregated_lagrangian_dual(
                child_node,
                node_index+1,
                outgoing_state,
                scenario_path,
                i,
                epi_state,
                normalization_coeff,
                primal_obj,
                dual_vars_agg,
                dual_0_var,
                bound_results,
                backward_sampling_scheme,
                algo_params,
                cut_generation_regime,
                applied_solvers,
                duality_regime.dual_solution_regime
            )
        end

        lag_obj = results.lag_obj
        lag_iterations = results.iterations
        lag_status = results.lag_status
        dual_0_var = results.dual_0_var

        model.ext[:agg_lag_iterations] += results.iterations

        # Counter to compare only number of iterations for converged cases
        if lag_status in (:opt, :conv, :sub, :iter, :mn_opt, :mn_iter)
            model.ext[:corr_lag_iterations] += results.iterations
            model.ext[:corr_realizations] += 1
        end

        ####################################################################
        # CHECK STATUS FOR ABNORMAL BEHAVIOR
        ####################################################################
        # if status is not as intended, the algorithm terminates with an error
        lagrangian_status_check(model, lag_status, duality_regime.dual_status_regime)

        Infiltrator.@infiltrate algo_params.infiltrate_state in [:all, :lagrange]

    catch e
        #SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end

    ############################################################################
    # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
    ############################################################################
    if isapprox(lag_obj, 0.0, atol=1e-8)

        add_cut_flag = false

    elseif dual_0_var == 0
        """ There should be no feasibility cuts in our setting."""

        add_cut_flag = false

    elseif !isnothing(normalization_coeff) && all(normalization_coeff.ω .== 0.0) && isapprox(normalization_coeff.ω₀, 0.0, atol=1e-8)
         """ If the linear pseudonorm is used, but all coefficients are zero,
         then the Lagrangian dual is not normalized, but unbounded. Analogously,
         a zero function is optimized over the reverse polar set, which can yield
         any point in this unbounded set. Therefore, we are not guaranteed to
         obtain any meaningful cut.

         Note that we artificially bound the dual objective, so that
         we obtain valid, but possibly very large cut coefficients.
         However, experiments show that the cuts still tend to be invalid due
         to numerical issues for these large coefficients. Therefore, in this
         case we do not construct a new cut at all (or at least restrict
         to redundant cut).
         """

         add_cut_flag = false

     else
         # We have to correct the intercept. We do this at this point, as (in
         # contrast to some other corrections) it is not required for Benders cuts.
         lag_obj = lag_obj + epi_state * dual_0_var
     end

    store_dual_values!(child_node, dual_values, dual_vars_agg, bin_state, cut_generation_regime.state_approximation_regime)

    #Infiltrator.@infiltrate

    # Store output in items
    push!(items.duals, dual_values)
    push!(items.dual_0_var, dual_0_var)
    push!(items.supports, SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node)[1]) # not required in aggregated case
    push!(items.nodes, child_node.index)
    push!(items.probability, child.probability * 1.0) # not required in aggregated case
    push!(items.objectives, lag_obj)
    push!(items.belief, 0.0) # not required in my case
    push!(items.bin_state, bin_state)
    push!(items.lag_iterations, lag_iterations)
    push!(items.add_cut_flags, add_cut_flag)

    ########################################################################
    # RECHANGE STATE SPACE
    ########################################################################
    TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
        rechangeStateSpace!(child_node, child_node.subproblem, outgoing_state, cut_generation_regime.state_approximation_regime)
    end

    Infiltrator.@infiltrate algo_params.infiltrate_state in [:all]

end

############################################################################
# DROP SCENARIO PATH'S LAST ENTRY
############################################################################
if length(scenario_path) == length_scenario_path
    # No-op. There weren't any children to solve.
else
    # Drop the last element (i.e., the one we added).
    pop!(scenario_path)
end

return

end
