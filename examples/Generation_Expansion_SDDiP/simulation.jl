import SDDP
import DynamicSDDiP

include("scenario_tree.jl")

function simulate(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    problem_params::DynamicSDDiP.ProblemParams,
    simulation_regime::DynamicSDDiP.Simulation
    )

    simulate(model, algo_params, problem_params, simulation_regime.sampling_scheme)

    return
end


function simulate(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    problem_params::DynamicSDDiP.ProblemParams,
    simulation_regime::DynamicSDDiP.NoSimulation
    )

    return
end


function simulate(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    problem_params::DynamicSDDiP.ProblemParams,
    sampling_scheme::DynamicSDDiP.InSampleMonteCarlo
    )

    # SIMULATE THE MODEL
    ############################################################################
    simulations = SDDP.simulate(model, algo_params.simulation_regime.number_of_replications)

    # OBTAINING BOUNDS AND CONFIDENCE INTERVAL
    ############################################################################
    objectives = map(simulations) do simulation
        return sum(stage[:stage_objective] for stage in simulation)
    end

    μ, ci = SDDP.confidence_interval(objectives)

    # get last lower bound again
    lower_bound, _ = DynamicSDDiP.calculate_bound(model)

    # LOGGING OF SIMULATION RESULTS
    ############################################################################
    log_simulation_results(algo_params, μ, ci, lower_bound)

    return
end

function simulate(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    problem_params::DynamicSDDiP.ProblemParams,
    sampling_scheme::DynamicSDDiP.OutOfSampleMonteCarlo
    )

    ############################################################################
    # GET A SAMPLE PATH USING THE EXISTING DISTRIBUTIONS
    ############################################################################
    sampling_scheme = SDDP.OutOfSampleMonteCarlo(model, use_insample_transition = true) do node
        return get_scenario_path(algo_params, node)
    end

    # SIMULATE THE MODEL
    ############################################################################
    simulations = SDDP.simulate(model, algo_params.simulation_regime.number_of_replications)

    # OBTAINING BOUNDS AND CONFIDENCE INTERVAL
    ############################################################################
    objectives = map(simulations) do simulation
        return sum(stage[:stage_objective] for stage in simulation)
    end

    μ, ci = SDDP.confidence_interval(objectives)

    # get last lower bound again
    lower_bound, _ = DynamicSDDiP.calculate_bound(model)

    # LOGGING OF SIMULATION RESULTS
    ############################################################################
    log_simulation_results(algo_params, μ, ci, lower_bound)

    return
end


function log_simulation_results(
    algo_params::DynamicSDDiP.AlgoParams,
    μ::Float64,
    ci::Float64,
    lower_bound::Float64
)

    log_file_handle = open(algo_params.log_file, "a")
    DynamicSDDiP.print_helper(DynamicSDDiP.print_simulation, log_file_handle, algo_params, μ, ci, lower_bound)
    close(log_file_handle)

    return
end