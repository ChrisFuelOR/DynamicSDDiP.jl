import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
import Random
import CSV
import DataFrames

function get_scenarios(algo_params::DynamicSDDiP.AlgoParams, problem_params::DynamicSDDiP.ProblemParams, number_of_clients::Int, data::Vector{String})

    # Get probability data
    probability = parse.(Float64,split(chop(data[8],head=1),","))
    @assert problem_params.number_of_realizations == length(probability)

    # Get realization data
    h_data = DelimitedFiles.readdlm(IOBuffer(data[7]), ',', ';')

    @assert number_of_clients == length(h_data)
    h = Array{Float64}(undef, number_of_clients, problem_params.number_of_realizations)

    for i in 1:number_of_clients
        if i == 1
            parsed_line = parse.(Float64,split(chop(h_data[i],head=1,tail=0)))
        elseif i == number_of_clients
            parsed_line = parse.(Float64,split(chop(h_data[i],head=0,tail=1)))
        else
            parsed_line = parse.(Float64,split(h_data[i]))
        end

        h[i,:] = parsed_line
        i += 1
    end

    # Get realization data into right format (we need a vector of vectors)
    support = Vector{Vector{Float64}}(undef, problem_params.number_of_realizations)

    for i in 1:length(support)
        support[i] = h[:,i]
    end

    return (support_array = support, probabilities_array = probability)

end


function get_scenario_path(algo_params::DynamicSDDiP.AlgoParams, stage::Int)

    # Set the seed from algo_params
    Random.seed!(algo_params.simulation_regime.sampling_scheme.simulation_seed)

    Error("Out of sample simulation not implemented for this problem yet.")

    return
end


function get_historical_sample(problem_params::DynamicSDDiP.ProblemParams, number_of_servers::Int, number_of_clients::Int)

    # In the 2-stage case we can just use the original scenarios here

    # Read data from file
    file_path = string("sslp1_", number_of_servers, "_", number_of_clients, "_", problem_params.number_of_realizations, ".txt")
    data_file = open(file_path)
    data = readlines(data_file)
    close(data_file)

    scenarios = get_scenarios(algo_params, problem_params, number_of_clients, data)

    return scenarios
end
