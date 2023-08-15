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

    initial_demand = 0.57e6
    initial_gas_price = 9.37

    # Define arrays for support and probabilities
    support_array = Vector{Array{NamedTuple{(:dem, :gas),Tuple{Float64,Float64}},1}}(undef,problem_params.number_of_stages)
    probabilities_array = Vector{Array{Float64,1}}(undef,problem_params.number_of_stages)

    # DETERMINISTIC VALUES FOR STAGE 1
    ############################################################################
    support_demand = 0
    support_gas = initial_gas_price

    support = [(dem = support_demand, gas = support_gas)]
    probability = [1.0]

    # Add both to a list of support and probablities for all stages
    support_array[1] = support
    probabilities_array[1] = probability

    # SAMPLE FOR EACH STAGE SEPARATELY (STAGEWISE INDEPENDENCE)
    ############################################################################
    for t in 2:problem_params.number_of_stages
        # Mean value
        mean = initial_demand * 1.00727^(t-1)

        # Demand follows a uniform distribution
        # (we allow for a deviation of t percent)
        demand_distribution = Distributions.Uniform(mean*(1-t/100), mean*(1+t/100)) #TODO

        # Sample from this distribution
        support_demand = (rand(demand_distribution, problem_params.number_of_realizations) .- initial_demand) / 8760.0

        # Remove negative demand values
        support_demand[support_demand .< 0] .= 0

        ########################################################################
        # Mean value
        mean = initial_gas_price * 1.041188^(t-1)
        std_dev = 0.8*1.05^(t-1) #own assumption

        # Gas prices follow a truncated normal distribution
        gas_price_distribution = Distributions.truncated(Distributions.Normal(mean,std_dev), 5, 30) #TODO

        # Sample from this distribution
        support_gas = rand(gas_price_distribution, problem_params.number_of_realizations)

        ########################################################################
        # Combination
        support = [(dem = support_demand[i], gas = support_gas[i]) for i in 1:problem_params.number_of_realizations]
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

    initial_demand = 0.57e6
    initial_gas_price = 9.37
    number_of_realizations = algo_params.simulation_regime.sampling_scheme.number_of_realizations

    if stage == 0
        # DETERMINISTIC VALUES FOR STAGE 1
        ########################################################################
        support_demand = 0
        support_gas = initial_gas_price

        support = [(dem = support_demand, gas = support_gas)]
        probability = [1.0]

    else
        # SAMPLE FROM THE UNDERLYING DISTRIBUTION
        ########################################################################
        t = stage

        # Mean value
        mean = initial_demand * 1.00727^(t-1)

        # Demand follows a uniform distribution (we allow for a deviation of t percent)
        demand_distribution = Distributions.Uniform(mean*(1-t/100), mean*(1+t/100))

        # Sample from this distribution
        support_demand = (rand(demand_distribution, number_of_realizations) .- initial_demand) / 8760.0

        # Remove negative demand values
        support_demand[support_demand .< 0] .= 0

        ########################################################################
        # Mean value
        mean = initial_gas_price * 1.041188^(t-1)
        std_dev = 0.8*1.05^(t-1) #own assumption

        # Gas prices follow a truncated normal distribution
        gas_price_distribution = Distributions.truncated(Distributions.Normal(mean,std_dev), 5, 30)

        # Sample from this distribution
        support_gas = rand(gas_price_distribution, number_of_realizations)

        ########################################################################
        # Combination
        support = [(dem = support_demand[i], gas = support_gas[i]) for i in 1:number_of_realizations]
        # We take an SAA approach
        probability = 1/number_of_realizations * ones(number_of_realizations)

    end

    return [SDDP.Noise(ω, p) for (ω, p) in zip(support, probability)]
end
