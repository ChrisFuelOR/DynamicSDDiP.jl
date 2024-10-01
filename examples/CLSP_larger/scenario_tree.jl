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
    support_array = Vector{Array{NamedTuple{(:xi1, :xi2, :xi3, :xi4, :xi5, :xi6, :xi7, :xi8, :xi9, :xi10),Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}},1}}(undef,problem_params.number_of_stages)
    probabilities_array = Vector{Array{Float64,1}}(undef,problem_params.number_of_stages)

    """ We use the scenario tree data from clsp_scenarios_synthetic.csv,
    which was provided by Filipe Cabral. This includes deterministic data
    for stage 1 and 20 realizations for all the other stages (up to 118).
    """

    # Read noise data from .csv
    noise_df = CSV.read("clsp_scenarios_larger.csv", DataFrames.DataFrame)

    # Group noise data by stage
    noise_list = DataFrames.groupby(noise_df, :stage)

    # STORE STAGEWISE DATA IN VECTORS
    ############################################################################
    # Stage 1
    # Get noise values for all three products
    support_1 = Vector(noise_list[1][:,"1"])
    support_2 = Vector(noise_list[1][:,"2"])
    support_3 = Vector(noise_list[1][:,"3"])
    support_4 = Vector(noise_list[1][:,"4"])
    support_5 = Vector(noise_list[1][:,"5"])
    support_6 = Vector(noise_list[1][:,:"6"])
    support_7 = Vector(noise_list[1][:,:"7"])
    support_8 = Vector(noise_list[1][:,:"8"])
    support_9 = Vector(noise_list[1][:,:"9"])
    support_10 = Vector(noise_list[1][:,:"10"])

    # Combination
    support = [(xi1 = support_1[1], xi2 = support_2[1], xi3 = support_3[1], xi4 = support_4[1], xi5 = support_5[1], xi6 = support_6[1], xi7 = support_7[1], xi8 = support_8[1], xi9 = support_9[1], xi10 = support_10[1])]
    probability = [1.0]

    # Add both to a list of support and probablities for all stages
    support_array[1] = support
    probabilities_array[1] = probability

    # Further stages
    for t in 2:problem_params.number_of_stages

        # Get noise values for all three products
        support_1 = Vector(noise_list[t][:,"1"])
        support_2 = Vector(noise_list[t][:,"2"])
        support_3 = Vector(noise_list[t][:,"3"])
        support_4 = Vector(noise_list[t][:,"4"])
        support_5 = Vector(noise_list[t][:,"5"])
        support_6 = Vector(noise_list[t][:,:"6"])
        support_7 = Vector(noise_list[t][:,:"7"])
        support_8 = Vector(noise_list[t][:,:"8"])
        support_9 = Vector(noise_list[t][:,:"9"])
        support_10 = Vector(noise_list[t][:,:"10"])

        # Combination
        support = [(xi1 = support_1[i], xi2 = support_2[i], xi3 = support_3[i], xi4 = support_4[i], xi5 = support_5[i], xi6 = support_6[i], xi7 = support_7[i], xi8 = support_8[i], xi9 = support_9[i], xi10 = support_10[i]) for i in 1:problem_params.number_of_realizations]
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
    historical_samples_all = Vector{Array{Tuple{Int64,NamedTuple{(:xi1, :xi2, :xi3, :xi4, :xi5, :xi6, :xi7, :xi8, :xi9, :xi10),Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}},1}}(undef, number_of_total_scenarios)

    # Iterate over number_of_total_scenarios
    for s in 1:number_of_total_scenarios
        historical_sample = Vector{Tuple{Int64,NamedTuple{(:xi1, :xi2, :xi3, :xi4, :xi5, :xi6, :xi7, :xi8, :xi9, :xi10),Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}}}(undef, number_of_stages)
        
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
