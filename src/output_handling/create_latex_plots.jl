struct PlotConfig
    file_path::String
    plottype::Symbol
    color::String
	opacity::Float64
	mark::String
    property::String
    legend::String
	description::String
    restrict_it::Union{Nothing,Int}
end

include("read_txt_file.jl")


function create_latex_plots()

	############################################################################
	############################################################################
	############################################################################
	# TO BE ADAPTED FOR EACH NEW PLOT
	############################################################################

    file_path_latex = "C:/Users/cg4102/Documents/julia_plots/clsp_large_16_no_bin_lb_time.tex"
	plottype = :lb_time

    # Create an array of struct instances which define the plots to be created
    plots_to_create = (
        PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Single.log", plottype, "red", 1.0, "none", "solid", "Single-cut", "SDDiP Single-cut", nothing),
        PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", "Multi-cut", "SDDiP Multi-cut", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"$ \ell_1 $", "L1", nothing),
        PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1sup.log", plottype, "green", 1.0, "none", "solid", raw"$ \ell_{1 \infty} $", "L1sup", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Lsup_MNC.log", plottype, "green!70!black", 1.0, "none", "solid", raw"$ \ell_\infty$", "Lsup", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", "Mid", "Mid", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", "Eps", "Eps", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", "Relint", "Relint", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Opt.log", plottype, "violet", 1.0, "none", "solid", "Opt", "Relint", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_InOut.log", plottype, "orange", 1.0, "none", "solid", "In-Out", "Relint", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_SB_Multi.log", plottype, "gray", 1.0, "none", "solid", "SB multi", "SB multi", nothing),
    )

    # Create header of LaTeX file
	create_latex_header(file_path_latex, [0.0, 70.0, 0.0, 30000.0])

	############################################################################
	############################################################################
	############################################################################

    # For each plot in plots_to_create add the plot data in LaTeX
	for plot in plots_to_create
		create_plot(file_path_latex, plot)
	end

    # Create footer of LaTeX file
	create_latex_footer(file_path_latex)

end

function create_latex_header(file_path_latex::String, axis_limits::Vector{Float64})

	io = open(file_path_latex, "w")

	println(io, "\\documentclass[tikz]{standalone}")
	println(io, "\\usetikzlibrary{intersections}")
	println(io, "\\usepackage{pgfplots}")
	println(io, "\\usepackage{etoolbox}")
	println(io)

	println(io, "\\begin{document}")
	println(io)

	println(io, "\\begin{tikzpicture}[scale=1,line cap=round,every mark/.append style={mark size=1pt}]")
	println(io, "%Styles")
	println(io, "\\pgfplotsset{")
	println(io, "axis line style={gray!30!black},")
	println(io, "every axis label/.append style ={gray!30!black},")
	println(io, "every tick label/.append style={gray!30!black}")
	println(io, "}")
	println(io)

	println(io, "\\begin{axis}")
	println(io, "[")
	println(io, "axis on top = true,")
	println(io, "axis lines = left,")
	println(io, "xmin=", axis_limits[1], ", xmax=", axis_limits[2], ", ymin=", axis_limits[3], ", ymax=", axis_limits[4], ",")
	println(io, "ticklabel style={fill=white},")
	println(io, "xlabel={Time [min]}, ylabel={LB},")
	println(io, "legend columns=3,") #TODO
	println(io, "legend style={at={(0.5,-0.25)},anchor=north, /tikz/every even column/.append style={column sep=0.5cm}, font=\\scriptsize, draw=none},")
	println(io, "clip=false] ")
	println(io)
	println(io)

	close(io)

end

function create_plot(file_path_latex::String, plot_config::PlotConfig)

	io = open(file_path_latex, "a")

	# PLOT HEADER
	############################################################################
	println(io, "%", plot_config.description)
	println(io, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

	print(io, "\\addplot[color=")
	print(io, plot_config.color)
	print(io, ",mark=")
	print(io, plot_config.mark)
	print(io, ",opacity=")
	print(io, plot_config.opacity)
	if !isnothing(plot_config.property)
		print(io, ", ", plot_config.property)
	end
	print(io, "] coordinates {")
	println(io)

	# PLOT DATA
	############################################################################
	# Get data from file
	df = read_txt_file(plot_config.file_path)

	# Check which columns are required
	if plot_config.plottype == :lb_time
		df = select!(df, [:Total_Time, :LB])
		# Convert to minuts
		df[:, :Total_Time] = df[:, :Total_Time] / 60.0
	elseif plot_config.plottype == :lb_it
		df = select!(df, [:Iteration, :LB])
	elseif plot_config.plottype == :it_time
		df = select!(df, [:Iteration, :It_Time])
	elseif plot_config.plottype == :lag_it
		df = select!(df, [:Iteration, :Lag_It_per_node])
	elseif plot_config.plottype == :lag_it_corr
		df = select!(df, [:Iteration, :Lag_It_per_node_corr])
	else
		Error("Plot type not defined.")
	end

	# Plot data points for each row
	for row in eachrow(df)

		if !isnothing(plot_config.restrict_it) && rownumber(row) <= plot_config.restrict_it
			print(io, "(   ")
			print(io, row[1])
			print(io, "   ,   ")
			print(io, row[2])
			print(io, ")   ")
			println(io)
		elseif isnothing(plot_config.restrict_it)
			print(io, "(   ")
			print(io, row[1])
			print(io, "   ,   ")
			print(io, row[2])
			print(io, ")   ")
			println(io)
		end
	end

	# PLOT FOOTER
	############################################################################
	println(io, "};")
	print(io, "\\addlegendentry{")
	print(io, plot_config.legend)
	print(io, "}")
	println(io)
	println(io)
	println(io)

	close(io)

end

function create_latex_footer(file_path_latex::String)

	io = open(file_path_latex, "a")

	println(io)
	println(io, "\\end{axis}")
	println(io, "\\end{tikzpicture}")
	println(io, "\\end{document}")

	close(io)

end

create_latex_plots()
