using Infiltrator
using CSV
using DataFrames


function get_key_numbers(file_name::String)

    open(file_name) do file
        counter_separator = 0
        line_number = 0
        first_line_number = 0
        last_line_number = 0
        number_of_stages = 0
        number_of_realizations = 0

        while !eof(file)
            line = readline(file)
            line_number += 1

            if line == "────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────"
                counter_separator += 1
            end

            if counter_separator == 2
                first_line_number = line_number + 1
                counter_separator = 0
            end

            if startswith(line, "Number of stages: ")
                number_of_stages = parse(Int, replace(line, "Number of stages: " => ""))
            end

            if startswith(line, "Number of realizations per stage: ")
                number_of_realizations = parse(Int, replace(line, "Number of realizations per stage: " => ""))
            end

            if startswith(line, "Terminating DynamicSDDiP with status")
                last_line_number = line_number - 2
            end

            if eof(file) && last_line_number == 0
                last_line_number = line_number
            end

        end

        lines_to_skip_at_end = line_number - last_line_number

        return (first_line_number = first_line_number, lines_to_skip_at_end = lines_to_skip_at_end, number_of_stages = number_of_stages, number_of_realizations = number_of_realizations)
    end
end

function read_data(file_name::String, first_line_number::Int, lines_to_skip_at_end::Int)

    df = CSV.read(file_name, header=false, skipto=first_line_number, footerskip=lines_to_skip_at_end, delim=" ", ignorerepeated=true, DataFrame)
    return df
end


function read_txt_file(file_name::String)

    # Get line number of first and last row of data
    key_numbers = get_key_numbers(file_name)

    total_number_of_nodes = (key_numbers.number_of_stages - 1) * key_numbers.number_of_realizations

    # Read main data from table into dataframe
    df = read_data(file_name, key_numbers.first_line_number, key_numbers.lines_to_skip_at_end)

    # Select required columns
    df = select!(df, Not([:Column2, :Column3, :Column5, :Column8, :Column9, :Column10, :Column11, :Column12, :Column13, :Column14, :Column15]))

    # Renaming columns
    rename!(df, :Column1 => :Iteration)
    rename!(df, :Column4 => :LB)
    rename!(df, :Column6 => :Total_Time)
    rename!(df, :Column7 => :It_Time)
    rename!(df, :Column16 => :Lag_It)

    # Add column for number of Lagrangian iterations per node
    df[:, :Lag_It_per_node] = df[:, :Lag_It] / total_number_of_nodes

    # Remove previous column
    df = select!(df, Not(:Lag_It))

    # Handling of correction for Lagrangian dual iterations
    # This is only logged for newer runs, so not included in older log files
    if "Column17" in names(df) || "Column18" in names(df)

        # Rename the columns
        rename!(df, :Column17 => :Lag_It_Corr)
        rename!(df, :Column18 => :Realizations_Corr)

        # Add column for this corrected value
        df[:, :Lag_It_per_node_corr] = df[:, :Lag_It_Corr] ./ df[:, :Realizations_Corr]
        df[:, :Lag_It_per_node_corr] = replace!(df[:, :Lag_It_per_node_corr], NaN => 0.0)

        # Remove columns which are no longer required
        df = select!(df, Not([:Lag_It_Corr, :Realizations_Corr]))

    else
        # Otherwise just use the same information as for the previous column
        df[:, :Lag_It_per_node_corr] = df[:, :Lag_It_per_node]
        @warn("No correction of number of Lagrangian dual iterations possible. Not logged for this run.")
    end

    return df
end


read_txt_file("C:/Users/cg4102/Documents/julia_logs/to_plot/clsp_binary_L1_16.log")
