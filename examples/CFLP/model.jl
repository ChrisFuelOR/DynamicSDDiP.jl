import SDDP
import DynamicSDDiP
import JuMP
import Infiltrator
using DataFramesMeta
using Revise


function model_definition(problem_params::DynamicSDDiP.ProblemParams, scenario_tree)

    """
    This model is a multistage stochastic capacitated facility location problem 
    under facility disruption uncertainty.
    The problem and the data were provided for us by Bonn Kleiford Seranilla and Nils Löhndorf.
    The problem is based on a non-stochastic single-stage problem from Beasly's OR library.
    """

    # Number of stages
    T = problem_params.number_of_stages

    # Number of facilities and customers
    number_of_facilities = 16
    number_of_customers = 50
    
    # Capacity (constant for each facility and stage)
    C = 50

    # Demand per customer (constant for each stage)
    D = [1, 0, 6, 13, 0, 5, 23, 10, 0, 0, 54, 9, 14, 1, 6, 5, 2, 30, 2, 1, 0, 8, 5, 3, 8, 3, 43, 5, 4, 4, 2, 3, 6, 129, 3, 3, 36, 22, 7, 3, 16, 11, 2, 5, 22, 7, 2, 0, 14, 2]

    # Fixed cost of opening a facility (constant for each stage)
    f = [75000, 75000, 75000, 75000, 75000, 75000, 75000, 75000, 75000, 75000, 0, 75000, 75000, 75000, 75000, 75000]

    # Transportation cost from facility to customer (constant for each stage)
    d = [[67.0, 103.0, 76.0, 52.0, 57.0, 66.0, 43.0, 38.0, 64.0, 53.0, 52.0, 41.0, 73.0, 50.0, 103.0, 60.0],
        [32.0, 54.0, 38.0, 23.0, 26.0, 32.0, 18.0, 22.0, 31.0, 25.0, 22.0, 17.0, 51.0, 21.0, 53.0, 28.0],
        [49.0, 264.0, 196.0, 138.0, 91.0, 149.0, 218.0, 353.0, 151.0, 236.0, 98.0, 193.0, 574.0, 111.0, 229.0, 154.0],
        [323.0, 299.0, 210.0, 296.0, 212.0, 200.0, 642.0, 801.0, 259.0, 692.0, 230.0, 487.0, 1351.0, 405.0, 605.0, 529.0],
        [17.0, 21.0, 15.0, 10.0, 12.0, 13.0, 15.0, 9.0, 13.0, 17.0, 11.0, 10.0, 20.0, 13.0, 25.0, 18.0],
        [64.0, 237.0, 161.0, 103.0, 74.0, 123.0, 158.0, 272.0, 124.0, 177.0, 70.0, 139.0, 454.0, 69.0, 203.0, 109.0],
        [819.0, 284.0, 431.0, 657.0, 588.0, 485.0, 1386.0, 1552.0, 531.0, 1473.0, 569.0, 1023.0, 2595.0, 964.0, 1319.0, 1183.0],
        [333.0, 265.0, 63.0, 167.0, 135.0, 88.0, 515.0, 579.0, 109.0, 573.0, 127.0, 335.0, 1057.0, 322.0, 600.0, 465.0],
        [20.0, 24.0, 18.0, 13.0, 15.0, 16.0, 18.0, 12.0, 15.0, 20.0, 13.0, 12.0, 19.0, 16.0, 28.0, 21.0],
        [14.0, 19.0, 14.0, 8.0, 10.0, 11.0, 11.0, 5.0, 11.0, 14.0, 9.0, 7.0, 19.0, 11.0, 23.0, 15.0],
        [1410.0, 2059.0, 1041.0, 126.0, 460.0, 661.0, 1983.0, 2202.0, 581.0, 2415.0, 252.0, 975.0, 4619.0, 1061.0, 2884.0, 1854.0],
        [176.0, 320.0, 153.0, 84.0, 12.0, 90.0, 327.0, 413.0, 92.0, 375.0, 24.0, 197.0, 792.0, 167.0, 419.0, 285.0],
        [382.0, 424.0, 153.0, 158.0, 115.0, 51.0, 626.0, 712.0, 115.0, 704.0, 103.0, 384.0, 1352.0, 366.0, 775.0, 558.0],
        [19.0, 50.0, 40.0, 34.0, 22.0, 35.0, 55.0, 83.0, 35.0, 59.0, 24.0, 48.0, 131.0, 32.0, 46.0, 42.0],
        [171.0, 360.0, 253.0, 162.0, 158.0, 211.0, 73.0, 217.0, 207.0, 125.0, 146.0, 119.0, 399.0, 111.0, 324.0, 143.0],
        [256.0, 356.0, 251.0, 157.0, 184.0, 212.0, 184.0, 81.0, 204.0, 233.0, 162.0, 126.0, 350.0, 192.0, 398.0, 249.0],
        [70.0, 104.0, 63.0, 25.0, 39.0, 47.0, 68.0, 71.0, 44.0, 87.0, 30.0, 27.0, 174.0, 50.0, 130.0, 82.0],
        [784.0, 925.0, 366.0, 274.0, 235.0, 230.0, 1263.0, 1413.0, 114.0, 1448.0, 228.0, 740.0, 2740.0, 752.0, 1594.0, 1148.0],
        [94.0, 124.0, 77.0, 35.0, 50.0, 60.0, 102.0, 102.0, 56.0, 124.0, 41.0, 56.0, 217.0, 80.0, 162.0, 117.0],
        [85.0, 141.0, 105.0, 72.0, 78.0, 91.0, 54.0, 63.0, 88.0, 67.0, 73.0, 58.0, 95.0, 66.0, 135.0, 77.0],
        [15.0, 20.0, 13.0, 6.0, 9.0, 10.0, 16.0, 16.0, 10.0, 19.0, 7.0, 9.0, 30.0, 12.0, 26.0, 18.0],
        [231.0, 487.0, 360.0, 261.0, 234.0, 304.0, 149.0, 338.0, 306.0, 111.0, 223.0, 209.0, 560.0, 173.0, 401.0, 159.0],
        [120.0, 198.0, 96.0, 38.0, 22.0, 58.0, 196.0, 238.0, 62.0, 228.0, 0.0, 107.0, 469.0, 101.0, 268.0, 177.0],
        [48.0, 128.0, 108.0, 89.0, 67.0, 94.0, 119.0, 181.0, 94.0, 116.0, 71.0, 108.0, 282.0, 71.0, 77.0, 79.0],
        [240.0, 396.0, 246.0, 110.0, 136.0, 189.0, 201.0, 246.0, 177.0, 271.0, 105.0, 52.0, 592.0, 148.0, 445.0, 264.0],
        [41.0, 121.0, 81.0, 56.0, 29.0, 58.0, 116.0, 171.0, 59.0, 128.0, 32.0, 90.0, 294.0, 56.0, 131.0, 87.0],
        [2814.0, 4067.0, 3258.0, 2532.0, 2647.0, 2944.0, 2113.0, 2107.0, 2894.0, 2396.0, 2481.0, 2222.0, 1240.0, 2388.0, 3925.0, 2615.0],
        [110.0, 221.0, 114.0, 55.0, 24.0, 74.0, 192.0, 254.0, 75.0, 223.0, 16.0, 116.0, 496.0, 90.0, 265.0, 169.0],
        [85.0, 224.0, 141.0, 74.0, 66.0, 107.0, 115.0, 204.0, 108.0, 140.0, 59.0, 89.0, 370.0, 29.0, 205.0, 99.0],
        [124.0, 254.0, 221.0, 190.0, 155.0, 198.0, 239.0, 340.0, 199.0, 234.0, 160.0, 222.0, 504.0, 161.0, 80.0, 174.0],
        [37.0, 111.0, 82.0, 58.0, 46.0, 66.0, 64.0, 118.0, 66.0, 61.0, 44.0, 68.0, 187.0, 36.0, 85.0, 32.0],
        [46.0, 133.0, 78.0, 43.0, 28.0, 56.0, 96.0, 146.0, 57.0, 113.0, 23.0, 69.0, 267.0, 39.0, 133.0, 79.0],
        [134.0, 351.0, 259.0, 173.0, 152.0, 211.0, 168.0, 328.0, 213.0, 156.0, 146.0, 193.0, 531.0, 109.0, 278.0, 73.0],
        [3726.0, 2291.0, 2033.0, 3228.0, 2613.0, 2299.0, 6812.0, 8105.0, 2580.0, 7283.0, 2617.0, 5126.0, 13615.0, 4514.0, 6447.0, 5758.0],
        [97.0, 180.0, 120.0, 71.0, 75.0, 98.0, 57.0, 106.0, 96.0, 85.0, 63.0, 48.0, 213.0, 59.0, 178.0, 95.0],
        [120.0, 181.0, 114.0, 53.0, 73.0, 88.0, 108.0, 102.0, 83.0, 139.0, 59.0, 41.0, 269.0, 84.0, 214.0, 136.0],
        [976.0, 736.0, 595.0, 833.0, 659.0, 569.0, 1852.0, 2220.0, 730.0, 1987.0, 659.0, 1372.0, 3787.0, 1199.0, 1749.0, 1553.0],
        [607.0, 635.0, 273.0, 309.0, 204.0, 150.0, 977.0, 1145.0, 247.0, 1095.0, 205.0, 635.0, 2090.0, 583.0, 1144.0, 874.0],
        [544.0, 651.0, 521.0, 403.0, 444.0, 472.0, 476.0, 343.0, 468.0, 536.0, 427.0, 392.0, 409.0, 483.0, 733.0, 557.0],
        [71.0, 86.0, 64.0, 78.0, 52.0, 56.0, 152.0, 197.0, 70.0, 163.0, 56.0, 119.0, 326.0, 94.0, 135.0, 123.0],
        [380.0, 707.0, 395.0, 158.0, 164.0, 279.0, 498.0, 626.0, 269.0, 619.0, 100.0, 225.0, 1331.0, 231.0, 812.0, 474.0],
        [397.0, 529.0, 322.0, 136.0, 204.0, 245.0, 411.0, 408.0, 228.0, 506.0, 161.0, 206.0, 919.0, 314.0, 696.0, 475.0],
        [161.0, 207.0, 156.0, 110.0, 125.0, 137.0, 125.0, 78.0, 133.0, 148.0, 115.0, 100.0, 131.0, 130.0, 231.0, 156.0],
        [169.0, 325.0, 233.0, 162.0, 162.0, 198.0, 93.0, 207.0, 195.0, 95.0, 145.0, 130.0, 320.0, 124.0, 275.0, 125.0],
        [1686.0, 2107.0, 1692.0, 1319.0, 1446.0, 1537.0, 1347.0, 1093.0, 1505.0, 1492.0, 1360.0, 1248.0, 899.0, 1424.0, 2244.0, 1605.0],
        [571.0, 667.0, 531.0, 409.0, 453.0, 480.0, 525.0, 390.0, 469.0, 588.0, 426.0, 405.0, 459.0, 491.0, 759.0, 610.0],
        [155.0, 184.0, 143.0, 106.0, 120.0, 128.0, 142.0, 101.0, 125.0, 161.0, 111.0, 105.0, 149.0, 131.0, 212.0, 167.0],
        [25.0, 39.0, 30.0, 22.0, 23.0, 26.0, 17.0, 19.0, 26.0, 20.0, 21.0, 18.0, 20.0, 20.0, 37.0, 23.0],
        [340.0, 342.0, 244.0, 313.0, 219.0, 208.0, 690.0, 866.0, 272.0, 743.0, 238.0, 519.0, 1479.0, 429.0, 648.0, 565.0],
        [70.0, 119.0, 78.0, 41.0, 48.0, 63.0, 49.0, 64.0, 60.0, 68.0, 40.0, 26.0, 149.0, 45.0, 126.0, 74.0]]

    model = SDDP.LinearPolicyGraph(
        stages = T,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
        sense = :Min,
    ) do subproblem, t

        # Add state variable
        JuMP.@variable(subproblem, u[1:number_of_facilities], SDDP.State, Bin, initial_value = 1)

        # Add local decision variables
        JuMP.@variable(subproblem, x[1:number_of_facilities, 1:number_of_customers] >= 0, Int) #TODO: Should this be bounded from below?
        JuMP.@variable(subproblem, y[1:number_of_facilities], Bin)

        # Add random variable (disruption)
        JuMP.@variable(subproblem, disrupt[1:number_of_facilities]) #TODO: Bin

        # Add coupling constraint 
        # JuMP.@constraint(subproblem, coupling[i=1:number_of_facilities], u[i].out == u[i].in * (1 - disrupt[i]) + y[i])
        """ 
        Even though both variables are fixed before we solve the subproblem, 
        the original coupling constraints are interpreted as nonlinear if we consider ξ and u.in variables. 
        As both are binary, we can exactly reformulate their product, though, to draw on LP solvers. 
        """
        JuMP.@variable(subproblem, auxiliary[1:number_of_facilities], Bin)
        JuMP.@constraint(subproblem, coupling[i=1:number_of_facilities], u[i].out == auxiliary[i] + y[i])
        JuMP.@constraint(subproblem, [i=1:number_of_facilities], auxiliary[i] <= u[i].in)
        JuMP.@constraint(subproblem, [i=1:number_of_facilities], auxiliary[i] <= 1 - disrupt[i])  
        JuMP.@constraint(subproblem, [i=1:number_of_facilities], auxiliary[i] >= u[i].in + (1-disrupt[i]) - 1)  

        # Add further constraints
        JuMP.@constraint(subproblem, cap[i=1:number_of_facilities], sum(x[i,j] for j in 1:number_of_customers) <= C * u[i].out)
        JuMP.@constraint(subproblem, dem[j=1:number_of_customers], sum(x[i,j] for i in 1:number_of_facilities) == D[j])

        # Add stage objective
        scaling_factor = 1/1000
        # scaling_factor = 1
        SDDP.@stageobjective(subproblem, scaling_factor * sum(f[i] * y[i] for i in 1:number_of_facilities) + scaling_factor * sum(d[j][i] * x[i,j] for j in 1:number_of_customers for i in 1:number_of_facilities))

        # PARAMETERIZE THE RANDOM VARIABLES
        ########################################################################
        # Get the support and probability for the current stage
        support = Vector{Vector{Float64}}()

        if t == 1
            considered_realizations = 1
        else
            considered_realizations = problem_params.number_of_realizations
        end 

        # Iterate over scenarios to get the right values
        for s in 1:considered_realizations
            support_data = DataFramesMeta.@rsubset(scenario_tree, :t == t, :s == s)[!, Not([:row_number, :s, :t])]
            support_data = [values(support_data[1,:])...]
            push!(support, support_data)
        end

        # Parameterize the disruptions
        SDDP.parameterize(subproblem, support) do ω
            for i in 1:number_of_facilities
               JuMP.fix(disrupt[i], ω[i], force=true)
            end
        end

        # Switch the model to silent mode
        JuMP.set_silent(subproblem)

        return
    end

    return model
end


function model_set_up(
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
    # GET FINITE SCENARIO TREE FOR MODEL
    ############################################################################
    scenario_tree = get_scenario_data(100, 16, 20, :medium)

    ############################################################################
    # DEFINE MODEL
    ############################################################################
    model = model_definition(problem_params, scenario_tree)

    return (model = model, problem_params = problem_params)
end
