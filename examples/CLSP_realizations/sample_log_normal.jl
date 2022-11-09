import CSV
import Distributions
import Infiltrator
import DataFrames
import Random


function generate_demand_distribution(tree_seed::Int, number_of_stages::Int, number_of_realizations::Int)

    distribution_list = [
        Distributions.LogNormal(0.13, 0.32),
        Distributions.LogNormal(0.11, 0.31),
        Distributions.LogNormal(0.09, 0.3),
    ]

    # Set the seed
    Random.seed!(tree_seed)

    # Create a dataframe for the data
    df = DataFrames.DataFrame(stage = Int[], scenario = Int[], probablity = Float64[])
    for state in 1:length(distribution_list)
        df[:, Symbol(state)] = Float64[]
    end

    # Iterate over stages
    for stage in 1:number_of_stages

        # generate a deterministic demand for first stage by taking one sample only
        if stage == 1
            data_list = Any[]
            push!(data_list, stage)
            push!(data_list, 1)
            push!(data_list, 1.0)

            # Iterate over state variables
            for state in 1:length(distribution_list)
                distribution = distribution_list[state]
                push!(data_list, rand(distribution, 1)[1])
            end

            push!(df, data_list)

        # otherwise generate several samples
        else
            # Iterate over realizations
            for outcome in 1:number_of_realizations
                data_list = Any[]
                push!(data_list, stage)
                push!(data_list, outcome)
                push!(data_list, 1/number_of_realizations)

                # Iterate over state variables
                for state in 1:length(distribution_list)
                    distribution = distribution_list[state]
                    push!(data_list, rand(distribution, 1)[1])
                end

                push!(df, data_list)
            end

        end

    end

    # Write the data frame to csv
    CSV.write("clsp_scenarios_realizations.csv", df)

end

generate_demand_distribution(45678, 20, 50)
