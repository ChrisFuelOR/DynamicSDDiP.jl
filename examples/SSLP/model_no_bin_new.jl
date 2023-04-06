import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
import DataFrames
import DelimitedFiles
import CSV
using Revise

include("scenario_tree.jl")


function model_no_bin_definition(problem_params::DynamicSDDiP.ProblemParams, number_of_clients::Int, number_of_servers::Int, data::Vector{String}, scenarios)

    """
    This model is the same as analyzed in the paper by Chen and Luedtke (2022)
    and first proposed by Ntaimo and Sen (2005).

    This is the model version without binary approximation of the state variables.
    """

    # cost vector for servers at given locations (length of number_of_servers)
    c = parse.(Float64,split(chop(data[2],head=1),","))

    # revenues from client i being served by server at location j (matrix of size number_of_clients * number_of_servers)
    q_data = DelimitedFiles.readdlm(IOBuffer(data[3]), ',', ';')
    @assert number_of_clients == length(q_data)
    q = Array{Float64}(undef, number_of_clients, number_of_servers)

    for i in 1:number_of_clients
        if i == 1
            parsed_line = parse.(Float64,split(chop(q_data[i],head=1,tail=0)))
        elseif i == number_of_clients
            parsed_line = parse.(Float64,split(chop(q_data[i],head=0,tail=1)))
        else
            parsed_line = parse.(Float64,split(q_data[i]))
        end

        q[i,:] = parsed_line
        i += 1
    end

    # shortage cost (deterministic, length of number_of_servers)
    q0 = parse.(Float64,split(chop(data[4],head=1),","))

    # demand data (matrix of size number_of_clients * number_of_servers)
    d_data = DelimitedFiles.readdlm(IOBuffer(data[5]), ',', ';')
    @assert number_of_clients == length(d_data)
    d = Array{Float64}(undef, number_of_clients, number_of_servers)

    for i in 1:number_of_clients
        if i == 1
            parsed_line = parse.(Float64,split(chop(d_data[i],head=1,tail=0)))
        elseif i == number_of_clients
            parsed_line = parse.(Float64,split(chop(d_data[i],head=0,tail=1)))
        else
            parsed_line = parse.(Float64,split(d_data[i]))
        end

        d[i,:] = parsed_line
    end

    # server capacity (same for all servers)
    u = parse(Float64, data[6])

    ############################################################################

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # Add state variable
        JuMP.@variable(subproblem, x[1:number_of_servers], SDDP.State, Bin, initial_value = 0)

        if t == 1
            # Stage objective
            SDDP.@stageobjective(subproblem, sum(c[j] * x[j].out for j in 1:number_of_servers))

            # Bound constraints
            #JuMP.@constraint(subproblem, sum(x[j].out for j in 1:number_of_servers) <= ν)
            # JuMP.@constraint(subproblem, sum(x[j] for j in 1:number_of_servers) >= 0)

        elseif t == 2
            # Add local variables
            JuMP.@variable(subproblem, y[1:number_of_clients,1:number_of_servers], Bin)
            JuMP.@variable(subproblem, overflow[1:number_of_servers] >= 0)

            # Random variable
            JuMP.@variable(subproblem, h[1:number_of_clients])

            # Stage objective
            SDDP.@stageobjective(subproblem, -sum(q[i,j] * y[i,j] for i in 1:number_of_clients, j in 1:number_of_servers) + sum(q0[j] * overflow[j] for j in 1:number_of_servers))

            # Add constraints
            JuMP.@constraint(subproblem, cons1[j in 1:number_of_servers], sum(d[i,j] * y[i,j] for i in 1:number_of_clients) - overflow[j] <= u * x[j].in)
            JuMP.@constraint(subproblem, cons2[i in 1:number_of_clients], sum(y[i,j] for j in 1:number_of_servers) == h[i])

            # PARAMETERIZE THE RANDOM VARIABLES
            ########################################################################
            # Get the support and probability for the current stage
            support = scenarios.support_array
            probability = scenarios.probabilities_array

            # Parameterize the demand
            SDDP.parameterize(subproblem, support, probability) do ω
                for i in 1:length(h)
                    JuMP.fix(h[i], ω[i])
                end
            end
        end

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_no_bin_set_up(
    number_of_servers::Int,
    number_of_clients::Int,
    number_of_realizations::Int;
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
    tree_seed::Int = 12345
)

    # Read data from file
    file_path = string("sslp1_", number_of_servers, "_", number_of_clients, "_", number_of_realizations, ".txt")
    data_file = open(file_path)
    data = readlines(data_file)
    close(data_file)

    # Get parameters number_of_servers, number_of_clients, number_of_realizations from file
    # Make sure that they are equal to the ones chosen so that the correct data is used
    param_vector = parse.(Int,split(chop(data[1],head=1),","))
    @assert param_vector[1] == number_of_servers
    @assert param_vector[2] == number_of_clients
    @assert param_vector[3] == number_of_realizations

    number_of_stages = 2

    ############################################################################
    # DEFINE PROBLEM PARAMS
    ############################################################################
    problem_params = DynamicSDDiP.ProblemParams(number_of_stages, number_of_realizations, tree_seed = tree_seed)

    ############################################################################
    # GET FINITE SCENARIO TREE FOR MODEL
    ############################################################################
    scenarios = get_scenarios(algo_params, problem_params, number_of_clients, data)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_no_bin_definition(problem_params, number_of_clients, number_of_servers, data, scenarios)

    return (model = model, problem_params = problem_params)
end
