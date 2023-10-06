import Infiltrator
import CSV
import DataFrames
import DataFramesMeta

function read_scenario_data(file_name::String)
    return CSV.read(file_name, DataFrames.DataFrame; delim=",", header=false)
end

function get_scenario_data(T::Int, F::Int, S::Int, case::Symbol)

    if case == :high
        file_name = "RANDOM_T100_F16_HIGH.csv" 
    elseif case == :medium
        file_name = "RANDOM_T100_F16_MED.csv" 
    elseif case == :low
        file_name = "RANDOM_T100_F16_LOW.csv" 
    end
    scenario_df = read_scenario_data(file_name)

    scenario_df = deleteat!(scenario_df, 1:2) # first two rows is for some reason not used
    scenario_array = Matrix(scenario_df)
    scenario_array = reshape(transpose(scenario_array), (S*T*F,1))
    scenario_array = transpose(reshape(scenario_array, (F,S*T)))

    column_headers = ["F_$i" for i in 1:F]
    scenario_new_df = DataFrames.DataFrame(scenario_array, column_headers)
    scenario_new_df[!, :row_number] = 1:size(scenario_new_df,1)

    # Add scenario index row
    scenario_new_df[!, :s] = scenario_new_df[!, :row_number]
    scenario_new_df = DataFramesMeta.@transform(scenario_new_df, :s = :row_number/T)
    scenario_new_df = DataFramesMeta.@eachrow scenario_new_df begin :s = ceil(Int, :s) end
    scenario_new_df[!, :s] = convert.(Int, scenario_new_df[!,:s])

    # Add time index row
    scenario_new_df[!, :t] = scenario_new_df[!, :row_number]
    scenario_new_df = DataFramesMeta.@eachrow scenario_new_df begin :t = mod(:row_number, T) end
    scenario_new_df = DataFramesMeta.@eachrow scenario_new_df begin 
        if :t == 0    
            :t = T
        end
    end

    return scenario_new_df
end

