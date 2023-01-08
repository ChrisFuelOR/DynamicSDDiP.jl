import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using Revise

include("scenario_tree.jl")


struct Leg
    sym::Symbol
end

struct Itinerary
    sym::Symbol
    cancellation_rate::Float64
    legs::Vector{Symbol}
end

struct Compartment
    sym::Symbol
    capacity::Int
end

struct FareClass
    sym::Symbol
    compartment::Symbol
end


function model_no_bin_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree, itineraries::Vector{Itinerary}, classes::Vector{FareClass})

    """
    This model is based on a paper by Möller, Römisch and Weber (2008).
    It is also used in the original SDDiP paper for experiments.
    Some parts of the implementation and data are based on the code available
    in msppy.py.
    """

    legs = [
        Leg(:AH),
        Leg(:HA),
        Leg(:BH),
        Leg(:HB),
        Leg(:CH),
        Leg(:HC),
    ]

    compartments = [
        Compartment(:B, 24),
        Compartment(:E, 216),
    ]

    """ We use the data from data.csv which is given in the paper by Möller et al."""
    all_data_df = CSV.read("data.csv", DataFrames.DataFrame)

    fare = JuMP.Containers.DenseAxisArray{Int64}(undef, itineraries, classes)

    for i in itineraries
        for j in classes
            fare[i,j] = filter(row -> row.i == String(i.sym) && row.j == String(j.sym), all_data_df).fare[1]
        end
    end

    model = SDDP.LinearPolicyGraph(
        stages = problem_params.number_of_stages,
        lower_bound = -500000,
        #lower_bound = 0.0,
        optimizer = GAMS.Optimizer,
        sense = :Min,
    ) do subproblem, t

        num_itineraries = length(itineraries)
        num_classes = length(classes)
        num_legs = length(legs)
        num_compartments = length(compartments)

        # ORIGINAL STATE VARIABLES
        ########################################################################
        # Cumulative number of fulfilled bookings
        JuMP.@variable(
            subproblem,
            0 <= cum_book[i in itineraries, j in classes] <= sum(k.capacity * (1 + problem_params.number_of_stages * i.cancellation_rate) for k in compartments if k.sym == j.compartment),
            SDDP.State,
            Int,
            initial_value = 0,
        )

        # Cumulative number of fulfilled cancellation
        JuMP.@variable(
            subproblem,
            0 <= cum_canc[i in itineraries, j in classes] <= 1, #TODO: Upper Bound
            SDDP.State,
            Int,
            initial_value = 0
        )

        # LOCAL VARIABLES
        ########################################################################
        # Number of fulfilled bookings on this stage
        JuMP.@variable(
            subproblem,
            0 <= book[i in itineraries, j in classes],
            Int
        )

        # Number of fulfilled cancellations on this stage
        JuMP.@variable(
            subproblem,
            0 <= canc[i in itineraries, j in classes],
            Int
        )

        # RANDOM VARIABLES
        ########################################################################
        # RHS uncertainty of demand
        JuMP.@variable(
            subproblem,
            demand[i in itineraries, j in classes]
        )

        # CONSTRAINTS
        ########################################################################
        # State equation bookings
        JuMP.@constraint(
            subproblem,
            state_bookings[i in itineraries, j in classes],
            cum_book[i,j].out == cum_book[i,j].in + book[i,j]
        )

        # State equation cancellations
        JuMP.@constraint(
            subproblem,
            state_cancellations[i in itineraries, j in classes],
            cum_canc[i,j].out == cum_canc[i,j].in + canc[i,j]
        )

        # Cancellation constraints
        JuMP.@constraint(
            subproblem, [i in itineraries, j in classes],
            cum_canc[i,j].out <= i.cancellation_rate * cum_book[i,j].out + 0.5
        )

        JuMP.@constraint(
            subproblem, [i in itineraries, j in classes],
            cum_canc[i,j].out >= i.cancellation_rate * cum_book[i,j].out - 0.5
        )

        # Capacity constraint
        JuMP.@constraint(
            subproblem, [l in legs, k in compartments],
            sum(cum_book[i,j].out - cum_canc[i,j].out for i in itineraries for j in classes if (l.sym in i.legs && j.compartment == k.sym))
            <= k.capacity
        )

        # Demand constraint
        JuMP.@constraint(
            subproblem,
            demand_eq[i in itineraries, j in classes],
            book[i,j] <= demand[i,j]
        )

        # PARAMETERIZE THE RANDOM VARIABLES
        ########################################################################
        # Get the support and probability for the current stage
        support_vector_stage = scenario_tree.support_vector[t]
        prob_vector_stage = scenario_tree.prob_vector[t]

        # Parameterize the model using the uncertain demand
        SDDP.parameterize(subproblem, support_vector_stage, prob_vector_stage) do ω
            for l = 1:length(ω)
                i = ω[l][1]
                j = ω[l][2]
                JuMP.fix(demand[i,j], ω[l][3])
            end
        end

        # STAGE OBJECTIVE
        ########################################################################
        SDDP.@stageobjective(
            subproblem,
            -sum(fare[i,j] * book[i,j] - fare[i,j] * canc[i,j] for i in itineraries for j in classes)
        )

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_no_bin_set_up(
    number_of_stages::Int,
    number_of_realizations::Int;
    algo_params::DynamicSDDiP.AlgoParams = DynamicSDDiP.AlgoParams(),
    applied_solvers::DynamicSDDiP.AppliedSolvers = DynamicSDDiP.AppliedSolvers(),
    tree_seed::Int = 12345
)

    ############################################################################
    # DEFINE PROBLEM PARAMS
    ############################################################################
    problem_params = DynamicSDDiP.ProblemParams(number_of_stages, number_of_realizations, tree_seed = tree_seed)

    ############################################################################
    # DEFINE MODEL DATA THAT IS ALSO REQUIRED FOR THE SCENARIO TREE
    ############################################################################
    classes = [
        FareClass(:B1, :B),
        FareClass(:B2, :B),
        FareClass(:E1, :E),
        FareClass(:E2, :E),
        FareClass(:E3, :E),
        FareClass(:E4, :E),
    ]

    itineraries = [
        Itinerary(:AH, 0.1, [:AH]),
        Itinerary(:HA, 0.1, [:HA]),
        Itinerary(:BH, 0.1, [:BH]),
        Itinerary(:HB, 0.1, [:HB]),
        Itinerary(:CH, 0.05, [:CH]),
        Itinerary(:HC, 0.05, [:HC]),
        Itinerary(:AHB, 0.05, [:AH, :HB]),
        Itinerary(:BHA, 0.05, [:BH, :HA]),
        Itinerary(:AHC, 0.0, [:AH, :HC]),
        Itinerary(:CHA, 0.0, [:CH, :HA]),
        Itinerary(:BHC, 0.0, [:BH, :HC]),
        Itinerary(:CHB, 0.0, [:CH, :HB]),
    ]

    days = [182, 126, 84, 56, 35, 21, 14, 10, 7, 5, 3, 2, 1, 0]

    ############################################################################
    # GET FINITE SCENARIO TREE FOR MODEL
    ############################################################################
    scenario_tree = get_recombining_scenario_tree(algo_params, problem_params, itineraries, classes, days)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_no_bin_definition(problem_params, scenario_tree, itineraries, classes)

    return (model = model, problem_params = problem_params)
end


# function model_no_bin_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)
#
#     """
#     This model is based on a paper by Möller, Römisch and Weber (2008).
#     It is also used in the original SDDiP paper for experiments.
#     Some parts of the implementation and data are based on the code available
#     in msppy.py.
#     """
#
#     legs = [
#         Leg(:AH),
#         Leg(:HA),
#         Leg(:BH),
#         Leg(:HB),
#         Leg(:CH),
#         Leg(:HC),
#     ]
#
#     itineraries = [
#         Itinerary(:AH, 0.1, [:AH]),
#         Itinerary(:HA, 0.1, [:HA]),
#         Itinerary(:BH, 0.1, [:BH]),
#         Itinerary(:HB, 0.1, [:HB]),
#         Itinerary(:CH, 0.05, [:CH]),
#         Itinerary(:HC, 0.05, [:HC]),
#         Itinerary(:AHB, 0.05, [:AH, :HB]),
#         Itinerary(:BHA, 0.05, [:BH, :HA]),
#         Itinerary(:AHC, 0.0, [:AH, :HC]),
#         Itinerary(:CHA, 0.0, [:CH, :HA]),
#         Itinerary(:BHC, 0.0, [:BH, :HC]),
#         Itinerary(:CHB, 0.0, [:CH, :HB]),
#     ]
#
#     compartments = [
#         Compartment(:B, 24),
#         Compartment(:E, 216),
#     ]
#
#     classes = [
#         FareClass(:B1, :B),
#         FareClass(:B2, :B),
#         FareClass(:E1, :E),
#         FareClass(:E2, :E),
#         FareClass(:E3, :E),
#         FareClass(:E4, :E),
#     ]
#
#     #TODO: Read Fare Prices from .csv
#
#     model = SDDP.LinearPolicyGraph(
#         stages = problem_params.number_of_stages,
#         lower_bound = 0.0,
#         optimizer = GAMS.Optimizer,
#         sense = :Min,
#     ) do subproblem, t
#
#         num_itineraries = length(itineraries)
#         num_classes = length(classes)
#         num_legs = length(legs)
#         num_compartments = length(compartments)
#
#         # ORIGINAL STATE VARIABLES
#         ########################################################################
#         # Cumulative number of fulfilled bookings
#         JuMP.@variable(
#             subproblem,
#             0 <= cum_book[i in 1:num_itineraries, j in 1:num_classes] <= 1, #TODO: Upper Bound
#             SDDP.State,
#             Int,
#             initial_value = 0
#         )
#
#         # Cumulative number of fulfilled cancellation
#         JuMP.@variable(
#             subproblem,
#             0 <= cum_canc[i in 1:num_itineraries, j in 1:num_classes] <= 1, #TODO: Upper Bound
#             SDDP.State,
#             Int,
#             initial_value = 0
#         )
#
#         # LOCAL VARIABLES
#         ########################################################################
#         # Number of fulfilled bookings on this stage
#         JuMP.@variable(
#             subproblem,
#             0 <= book[i in 1:num_itineraries, j in 1:num_classes],
#             Int
#         )
#
#         # Number of fulfilled cancellations on this stage
#         JuMP.@variable(
#             subproblem,
#             0 <= canc[i in 1:num_itineraries, j in 1:num_classes],
#             Int
#         )
#
#         # RANDOM VARIABLES
#         ########################################################################
#         # RHS uncertainty of demand
#         JuMP.@variable(
#             subproblem,
#             demand
#         )
#
#         # CONSTRAINTS
#         ########################################################################
#         # State equation bookings
#         JuMP.@constraint(
#             subproblem,
#             state_bookings[i in 1:num_itineraries, j in 1:num_classes],
#             cum_book[i,j].out == cum_book[i,j].in + book[i,j]
#         )
#
#         # State equation cancellations
#         JuMP.@constraint(
#             subproblem,
#             state_cancellations[i in 1:num_itineraries, j in 1:num_classes],
#             cum_canc[i,j].out == cum_canc[i,j].in + canc[i,j]
#         )
#
#         # Cancellation constraints
#         JuMP.@constraint(
#             subproblem, [i in 1:num_itineraries, j in 1:num_classes],
#             cum_canc[i,j].out <= itineraries[i].cancellation_rate * cum_book[i,j].out + 0.5
#         )
#
#         JuMP.@constraint(
#             subproblem, [i in 1:num_itineraries, j in 1:num_classes],
#             cum_canc[i,j].out >= itineraries[i].cancellation_rate * cum_book[i,j].out - 0.5
#         )
#
#         # Capacity constraint
#         JuMP.@constraint(
#             subproblem, [l in legs, k in compartments],
#             sum(cum_book[i,j].out - cum_canc[i,j].out for i in 1:num_itineraries for j in 1:num_classes if (l.sym in itineraries[i].legs && classes[j].compartment == k.sym))
#             <= k.capacity
#         )
#
#         # Demand constraint
#         JuMP.@constraint(
#             subproblem,
#             demand_eq,
#             sum(book[i] for i in 1:num_itineraries) <= demand
#         )
#
#         # Switch the model to silent mode
#         JuMP.set_silent(subproblem)
#
#         return
#     end
#
#     return model
# end
