import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
import Random
import Distributions
import CSV
import DataFrames
import StatsBase


function get_recombining_scenario_tree(algo_params::DynamicSDDiP.AlgoParams, problem_params::DynamicSDDiP.ProblemParams, itineraries::Any, classes::Any, days::Vector{Int64})

    # Set the seed from algo_params
    Random.seed!(problem_params.tree_seed)

    """
    We use the data from data.csv which is given in the paper by Möller et al.
    """
    all_data_df = CSV.read("data.csv", DataFrames.DataFrame)

    # CONSTRUCT THE SCENARIO TREE DATA
    ############################################################################
    support_df = DataFrames.DataFrame(i = Any[], j = Any[], t = Int64[], support = Vector{Vector{Float64}}())

    # Iterate over each itinerary-class combination
    for i in itineraries
        for j in classes
            # Get row for given itinerary-class combination
            current_row = filter(row -> row.i == String(i.sym) && row.j == String(j.sym), all_data_df)

            # Get parameters for the corresponding distributions
            k_gamma = current_row.k[1]
            θ_gamma = current_row.theta[1]
            a_beta = current_row.a[1]
            b_beta = current_row.b[1]

            # Sample from the gamma distribution 1000 times
            # This gives u 1000 different numbers of cumulative bookings over the whole horizon
            gamma_distribution = Distributions.Gamma(k_gamma, θ_gamma)
            G_list = rand(gamma_distribution, 1000) / (1 - i.cancellation_rate)

            # Stage 1 demand. We assume a deterministic demand of 0.
            push!(support_df, (i, j, 1, [0.0]))

            # Iterate over the different time slots
            for t in 2:problem_params.number_of_stages

                stage_request_realizations = Vector{Int64}()

                # Iterate over the realizations of G and determine the cumulative
                # bookings up to this time for each case
                for G in G_list
                    # Use a beta distribution to do this. Note that a beta
                    # distribution is defined on [0,1], so this has to be
                    # adapted accordingly
                    beta_distribution = Distributions.Beta(a_beta, b_beta)
                    stage_requests = round(G * (Distributions.cdf(beta_distribution, 1-days[t]/days[1]) - Distributions.cdf(beta_distribution, 1-days[t-1]/days[1])))
                    push!(stage_request_realizations, stage_requests)

                end

                # Sample from the stage_request_realization number_of_realizations
                # times to get the realizations for the stagewise independent
                # process
                stage_support = StatsBase.sample(stage_request_realizations, problem_params.number_of_realizations)

                #println((i.sym, j.sym, t, mean(stage_support), mean(stage_request_realizations)))

                push!(support_df, (i, j, t, stage_support))
                #push!(prob_df, (i, j, t, stage_prob))
            end
        end
    end

    # RESHAPE THE DATA INTO THE REQUIRED FORM
    ############################################################################
    support_vector_all = Vector{Any}()
    prob_vector_all = Vector{Any}()

    # Iterate over the stages
    for t in 1:problem_params.number_of_stages
        # Get stage data
        support_stage_df = filter(row -> row.t == t, support_df)

        support_vector_stage = Vector{Any}()
        prob_vector_stage = Vector{Float64}()

        # Iterate over realizations

        if t == 1
            support_vector_stage_realization = Vector{Any}()

            # Iterate over itinerary-class combinations, i.e. data frame rows
            for row in eachrow(support_stage_df)
                # Create a tuple from this row data using only the l-th realization
                support_data = (row.i, row.j, row.support[1])
                push!(support_vector_stage_realization, support_data)
            end

            push!(support_vector_stage, support_vector_stage_realization)
            push!(prob_vector_stage, 1.0)

        else
            for l in 1:problem_params.number_of_realizations

                support_vector_stage_realization = Vector{Any}()

                # Iterate over itinerary-class combinations, i.e. data frame rows
                for row in eachrow(support_stage_df)
                    # Create a tuple from this row data using only the l-th realization
                    support_data = (row.i, row.j, row.support[l])
                    push!(support_vector_stage_realization, support_data)
                end

                push!(support_vector_stage, support_vector_stage_realization)
                push!(prob_vector_stage, 1/problem_params.number_of_realizations)

            end
        end

        push!(support_vector_all, support_vector_stage)
        push!(prob_vector_all, prob_vector_stage)
    end

    return (support_vector = support_vector_all, prob_vector = prob_vector_all)

end


function get_scenario_path(algo_params::DynamicSDDiP.AlgoParams, stage::Int)

    # Set the seed from algo_params
    Random.seed!(algo_params.simulation_regime.sampling_scheme.simulation_seed)

    Error("Out of sample simulation not implemented for this problem yet.")

    return
end
