import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
import Random
import CSV
import DataFrames

function get_recombining_scenario_tree(algo_params::DynamicSDDiP.AlgoParams, problem_params::DynamicSDDiP.ProblemParams)

    # Set the seed from algo_params
    Random.seed!(problem_params.tree_seed)

    # Define arrays for support and probabilities
    support_array = Vector{Array{NamedTuple{(:xi1, :xi2, :xi3),Tuple{Float64,Float64,Float64}},1}}(undef,problem_params.number_of_stages)
    probabilities_array = Vector{Array{Float64,1}}(undef,problem_params.number_of_stages)

    """ We use the scenario tree data from clsp_scenarios_synthetic.csv,
    which was provided by Filipe Cabral. This includes deterministic data
    for stage 1 and 20 realizations for all the other stages (up to 118).
    """

    # Read noise data from .csv
    noise_df = CSV.read("clsp_scenarios_synthetic.csv", DataFrames.DataFrame)

    # Group noise data by stage
    noise_list = DataFrames.groupby(noise_df, :stage)

    # STORE STAGEWISE DATA IN VECTORS
    ############################################################################
    # Stage 1
    # Get noise values for all three products
    support_1 = Vector(noise_list[1].xi1)
    support_2 = Vector(noise_list[1].xi2)
    support_3 = Vector(noise_list[1].xi3)

    # Combination
    support = [(xi1 = support_1[1], xi2 = support_2[1], xi3 = support_3[1])]
    probability = [1.0]

    # Add both to a list of support and probablities for all stages
    support_array[1] = support
    probabilities_array[1] = probability

    # Further stages
    for t in 2:problem_params.number_of_stages

        # Get noise values for all three products
        support_1 = Vector(noise_list[t].xi1)
        support_2 = Vector(noise_list[t].xi2)
        support_3 = Vector(noise_list[t].xi3)

        # Combination
        support = [(xi1 = support_1[i], xi2 = support_2[i], xi3 = support_3[i]) for i in 1:problem_params.number_of_realizations]
        # We take an SAA approach
        probability = 1/problem_params.number_of_realizations * ones(problem_params.number_of_realizations)

        ########################################################################
        # Add both to a list of support and probablities for all stages
        support_array[t] = support
        probabilities_array[t] = probability

    end

    return (support_array = support_array, probabilities_array = probabilities_array)

end


function get_scenario_path(algo_params::DynamicSDDiP.AlgoParams, stage::Int)

    # Set the seed from algo_params
    Random.seed!(algo_params.simulation_regime.sampling_scheme.simulation_seed)

    Error("Out of sample simulation not implemented for this problem yet.")

    return
end


function get_historical_sample(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)

    # Preparation
    number_of_realizations = problem_params.number_of_realizations
    number_of_stages = problem_params.number_of_stages
    number_of_total_scenarios = number_of_realizations^(number_of_stages - 1)
    historical_samples_all = Vector{Array{Tuple{Int64,NamedTuple{(:xi1, :xi2, :xi3),Tuple{Float64,Float64,Float64}}},1}}(undef, number_of_total_scenarios)

    # Iterate over number_of_total_scenarios
    for s in 1:number_of_total_scenarios
        historical_sample = Vector{Tuple{Int64,NamedTuple{(:xi1, :xi2, :xi3),Tuple{Float64,Float64,Float64}}}}(undef, number_of_stages)

        # Iterate over stages
        for t in 1:number_of_stages
            k = 1

            # Number of equivalent stage-t scenarios for each (stage-T) scenario
            equiv_number = number_of_realizations^(number_of_stages - 1) / number_of_realizations^(t - 1)

            # Get correct index for stage-t scenario for overall scenario s
            # Index in traditional scenario tree
            while k * equiv_number < s
                k += 1
            end

            # Index in recombining tree
            pos = k - div(k, number_of_realizations) * number_of_realizations
            if pos == 0
                pos = number_of_realizations
            end

            # Add stage realization to historical_sample
            historical_sample[t] = (t, scenario_tree.support_array[t][pos])
        end

        # Add this scenario path to historical_samples_all
        historical_samples_all[s] = historical_sample
    end

    return historical_samples_all
end
