"""
Solving aggregated dual in the backward pass.
"""
function solve_aggregated_dual(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    node_index::Int64,
    items::BackwardPassItems,
    belief::Float64,
    # belief_state,
    # objective_state,
    outgoing_state::Dict{Symbol,Float64},
    epi_state::Float64,
    primal_obj::Float64,
    backward_sampling_scheme::SDDP.AbstractBackwardSamplingScheme,
    scenario_path,
    algo_params::DynamicSDDiP.AlgoParams,
    applied_solvers::DynamicSDDiP.AppliedSolvers
) where {T}

    duality_regime = algo_params.duality_regime

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
        number_of_states_per_noise = get_number_of_states(child_node, algo_params.state_approximation_regime)
        number_of_states = number_of_states_per_noise * number_of_noise

        # storages for information on Lagrangian dual
        lag_obj = 0
        lag_iterations = 0
        lag_status = :none

        ########################################################################
        # INITIALIZE DUALS
        ########################################################################
        dual_vars = zeros(number_of_states)
        dual_0_var = 1.0 # shouldn't be zero if we do not consider feasibility cuts

        ########################################################################
        # GET BOUNDS FOR LAGRANGIAN DUAL
        ########################################################################
        bound_results = get_dual_bounds(child_node, node_index+1, algo_params, primal_obj, duality_regime.dual_bound_regime)
        @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        try
            ####################################################################
            # SOLVE THE AGGREGATED LAGRANGIAN DUAL PROBLEM
            ####################################################################
            TimerOutputs.@timeit DynamicSDDiP_TIMER "solve_lagrange" begin
                results = solve_aggregated_lagrangian(
                    child_node,
                    node_index+1,
                    outgoing_state,
                    scenario_path,
                    epi_state,
                    primal_obj,
                    dual_vars,
                    dual_0_var,
                    bound_results,
                    backward_sampling_scheme,
                    algo_params,
                    applied_solvers,
                    duality_regime.dual_solution_regime
                )
            end

            lag_obj = results.lag_obj
            lag_iterations = results.iterations
            lag_status = results.lag_status
            dual_0_var = results.dual_0_var

            ####################################################################
            # CHECK STATUS FOR ABNORMAL BEHAVIOR
            ####################################################################
            # if status is not as intended, the algorithm terminates with an error
            lagrangian_status_check(lag_status, duality_regime.dual_status_regime)

            @infiltrate algo_params.infiltrate_state in [:all, :lagrange]

        catch e
            #SDDP.write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
            rethrow(e)
        end

        ########################################################################
        # SET DUAL VARIABLES AND STATES CORRECTLY FOR RETURN
        ########################################################################
        # This includes a division by dual_0_var so that we can use the
        # original formulation of the Bellman function
        store_dual_values!(child_node, dual_values, dual_vars, dual_0_var, bin_state, algo_params.state_approximation_regime)

        @infiltrate

        # Store output in items
        push!(items.duals, dual_values)
        push!(items.supports, SDDP.sample_backward_noise_terms(backward_sampling_scheme, child_node)[1]) # not required in aggregated case
        push!(items.nodes, child_node.index)
        push!(items.probability, child.probability * 1.0) # not required in aggregated case
        push!(items.objectives, 0.0) # not required in my case
        push!(items.belief, 0.0) # not required in my case
        push!(items.bin_state, bin_state)
        push!(items.lag_iterations, lag_iterations)
        push!(items.lag_status, lag_status)

        ########################################################################
        # RECHANGE STATE SPACE
        ########################################################################
        TimerOutputs.@timeit DynamicSDDiP_TIMER "space_change" begin
            rechangeStateSpace!(child_node, child_node.subproblem, outgoing_state, algo_params.state_approximation_regime)
        end

        @infiltrate algo_params.infiltrate_state in [:all]

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


function store_dual_values!(node::SDDP.Node, dual_values::Dict{Symbol, Float64},
    dual_vars::Vector{Float64}, dual_0_var::Float64, bin_state::Dict{Symbol, BinaryState},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation)

    old_rhs = node.ext[:backward_data][:old_rhs]

    for (i, name) in enumerate(keys(node.ext[:backward_data][:bin_states]))
        dual_values[name] = dual_vars[i] / dual_0_var

        value = old_rhs[i]
        x_name = node.ext[:backward_data][:bin_x_names][name]
        k = node.ext[:backward_data][:bin_k][name]
        bin_state[name] = BinaryState(value, x_name, k)
    end

    return
end

function store_dual_values!(node::SDDP.Node, dual_values::Dict{Symbol, Float64},
    dual_vars::Vector{Float64}, dual_0_var::Float64, bin_state::Dict{Symbol, BinaryState},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation)

    for (i, name) in enumerate(keys(node.states))
        dual_values[name] = dual_vars[i] / dual_0_var
    end

    return
end
