import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise
import Random
import Distributions

function get_recombining_scenario_tree(algo_params::DynamicSDDiP.AlgoParams, problem_params::DynamicSDDiP.ProblemParams)

    # Set the seed from algo_params
    Random.seed!(problem_params.tree_seed)

    # Define arrays for support and probabilities
    support_array = Vector{Array{Float64,1}}(undef,problem_params.number_of_stages)
    probabilities_array = Vector{Array{Float64,1}}(undef,problem_params.number_of_stages)

    # DETERMINISTIC VALUES FOR STAGE 1
    ############################################################################
    # NOTE: Note that these values are not used in the model anyway
    support = [0]
    probability = [1.0]

    # Add both to a list of support and probablities for all stages
    support_array[1] = support
    probabilities_array[1] = probability

    # SAMPLE FOR EACH STAGE SEPARATELY (STAGEWISE INDEPENDENCE)
    ############################################################################
    for t in 2:problem_params.number_of_stages

        # Use the same fixed values for each stage
        # TODO: This should be changed later on.
        points = [
            -0.3089653673606697,
            -0.2718277412744214,
            -0.09611178608243474,
            0.24645863921577763,
            0.5204224537256875,
        ]

        support = [points; -points]

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
