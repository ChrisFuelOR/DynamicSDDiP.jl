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

function simulate(
    model::SDDP.PolicyGraph,
    algo_params::DynamicSDDiP.AlgoParams,
    problem_params::DynamicSDDiP.ProblemParams,
    number_of_servers::Int,
    number_of_clients::Int,
    sampling_scheme::DynamicSDDiP.HistoricalSample,
    )

    number_of_realizations = problem_params.number_of_realizations
    number_of_stages = problem_params.number_of_stages
    number_of_total_scenarios = number_of_realizations^(number_of_stages - 1)

    ############################################################################
    # GET A SAMPLE PATH USING THE EXISTING DISTRIBUTIONS
    ############################################################################
    historical_sample = get_historical_sample(algo_params, problem_params, number_of_servers, number_of_clients)

    # SIMULATE THE MODEL
    ############################################################################
    simulations = SDDP.simulate(model, number_of_total_scenarios, sampling_scheme = SDDP.Historical(historical_sample))

    # OBTAINING BOUNDS AND CONFIDENCE INTERVAL
    ############################################################################
    objectives = map(simulations) do simulation
        return sum(stage[:stage_objective] for stage in simulation)
    end

    upper_bound = sum(objectives)/number_of_total_scenarios

    # get last lower bound again
    lower_bound, _ = DynamicSDDiP.calculate_bound(model)

    # LOGGING OF SIMULATION RESULTS
    ############################################################################
    log_simulation_results_historical(algo_params, upper_bound, lower_bound)

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

function log_simulation_results_historical(
    algo_params::DynamicSDDiP.AlgoParams,
    upper_bound::Float64,
    lower_bound::Float64
)

    log_file_handle = open(algo_params.log_file, "a")
    DynamicSDDiP.print_helper(DynamicSDDiP.print_simulation_historical, log_file_handle, algo_params, upper_bound, lower_bound)
    close(log_file_handle)

    return
end
