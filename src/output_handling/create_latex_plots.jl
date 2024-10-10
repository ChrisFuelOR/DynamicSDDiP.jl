struct PlotConfig
    file_path::String
    plottype::Symbol
    color::String
	opacity::Float64
	mark::String
	property::String
	mark_options::String
    mark_repeat::Int
	mark_size::String
	legend::String
    restrict_it::Union{Nothing,Int}
end

include("read_txt_file.jl")


function create_latex_plots()

	############################################################################
	############################################################################
	############################################################################
	# TO BE ADAPTED FOR EACH NEW PLOT
	############################################################################

	input_folder_path = "C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CFLP_scaled/"
    #input_folder_path = "C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP"
    #input_folder_path = "C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP_Large"
    file_path_latex = "C:/Users/cg4102/Documents/julia_plots/Paper 2024/CFLP/CFLP_S_SB_lb_time2.tex"
	plottype = :lb_time

    # Create header of LaTeX file
	create_latex_header(file_path_latex, [0.0, 300.0, 50000.0, 80000.0], plottype)

    # Create an array of struct instances which define the plots to be created
	plots_to_create = (
		###################################################################### CLSP - 4 - 20 - 3 #############################################################################################################################
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_Single_SB.log", plottype, "black", 1.0, "none", "dashed", raw"Single", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"Multi", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_Mid_SB.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_16_20_InOut_SB.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", nothing),

		# PlotConfig(input_folder_path * "CLSP_RR2N_4_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_6_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_10_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_10_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_10_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_10_20_CL_Span.log", plottype, "yellow", 1.0, "none", "solid", raw"Span", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_CL_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_10_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		
			
		# ###################################################################### CLSP - 16 - 20 - 10 #############################################################################################################################
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_Single.log", plottype, "black", 1.0, "none", "dashed", raw"Single", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2_Large_O_16_20_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2_Large_O_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR2N_Large_O_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_Conv_90.log", plottype, "teal", 1.0, "none", "solid", raw"Conv", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_SB.log", plottype, "red", 1.0, "none", "dashed", raw"SB", nothing),
		# PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "dashed", raw"BMulti", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_SB_Multi.log", plottype, "red", 1.0, "none", "solid", raw"SBMulti", nothing),
		
		#PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_Large_O_16_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_CL_Span.log", plottype, "yellow", 1.0, "none", "solid", raw"Span", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_Large_O_16_20_CL_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2_Large_O_16_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR2N_Large_O_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),
		#PlotConfig(input_folder_path * "CLSP_RR_Large_O_16_20_Conv_90_CL.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", nothing),

		# ###################################################################### CFLP #############################################################################################################################
		PlotConfig(input_folder_path * "CFLP_S_Single_SB.log", plottype, "black", 1.0, "o", "dotted", "solid", 2, "1pt", raw"Single", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Multi_SB.log", plottype, "black", 1.0, "o", "solid", "solid", 2, "1pt", raw"Multi", nothing),
		PlotConfig(input_folder_path * "CFLP_S_L1_SB.log", plottype, "blue", 1.0, "asterisk", "solid", "solid", 2, "1pt", raw"L1", nothing),
		PlotConfig(input_folder_path * "CFLP_S_L1sup_SB.log", plottype, "green!70!black", 1.0, "mercedes star", "solid", "solid", 2, "1pt", raw"L1sup", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Lsup_SB.log", plottype, "green", 1.0, "10-pointed star", "solid", "solid", 2, "1pt", raw"Lsup", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Mid_N_sbopt_SB.log", plottype, "yellow!50!brown", 1.0, "square", "solid", "solid", 2, "1pt", raw"Mid", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Eps_N_sbopt_SB.log", plottype, "magenta", 1.0, "triangle", "solid", "solid", 2, "1pt", raw"Eps", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Relint_N_sbopt_SB.log", plottype, "cyan", 1.0, "diamond", "solid", "solid", 2, "1pt", raw"Relint", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Conv_50_N_SB.log", plottype, "violet", 1.0, "pentagon", "solid", "solid", 2, "1pt", raw"Conv(50)", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Conv_75_N_SB.log", plottype, "violet", 1.0, "pentagon", "dashed", "solid", 2, "1pt", raw"Conv(75)", nothing),
		PlotConfig(input_folder_path * "CFLP_S_Conv_90_N_SB.log", plottype, "violet", 1.0, "pentagon", "dotted", "solid", 2, "1pt", raw"Conv(90)", nothing),
		PlotConfig(input_folder_path * "CFLP_S_B.log", plottype, "olive", 1.0, "oplus", "dotted", "solid", 20, "1pt", raw"B", nothing),
		PlotConfig(input_folder_path * "CFLP_S_SB.log", plottype, "red", 1.0, "otimes", "dotted", "solid", 20, "1pt", raw"SB", nothing),
		PlotConfig(input_folder_path * "CFLP_S_B_Multi.log", plottype, "olive", 1.0, "oplus", "solid", "solid", 20, "1pt", raw"BMulti", nothing),
		PlotConfig(input_folder_path * "CFLP_S_SB_Multi.log", plottype, "red", 1.0, "otimes", "solid", "solid", 20, "1pt", raw"SBMulti", nothing),


	)

	############################################################################
	############################################################################
	############################################################################

    # For each plot in plots_to_create add the plot data in LaTeX
	for plot in plots_to_create
		create_plot(file_path_latex, plot)
		println("done")
	end

    # Create footer of LaTeX file
	create_latex_footer(file_path_latex)

end

function create_latex_header(file_path_latex::String, axis_limits::Vector{Float64}, plottype::Symbol)

	io = open(file_path_latex, "w")

	println(io, "\\documentclass[tikz]{standalone}")
	println(io, "\\usetikzlibrary{intersections}")
	println(io, "\\usepackage{pgfplots}")
	println(io, "\\usepackage{etoolbox}")
	println(io, "\\usetikzlibrary{matrix}")
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

	if plottype == :lb_time
		println(io, "xlabel={Time [min]}, ylabel={LB},")
	elseif plottype == :lb_it
		println(io, "xlabel={It}, ylabel={LB},")
	elseif plottype == :it_time
		println(io, "xlabel={It}, ylabel={Time [sec]},")
	elseif plottype == :lag_it
		println(io, "xlabel={It}, ylabel={Lag It},")
	elseif plottype == :lag_it_corr
		println(io, "xlabel={It}, ylabel={Lag It},")
	else
		error("Plot type not defined.")
	end

	println(io, "clip=false] ")
	println(io)
	println(io)

	close(io)

end

function create_plot(file_path_latex::String, plot_config::PlotConfig)

	io = open(file_path_latex, "a")

	# PLOT HEADER
	############################################################################
	println(io, "%", plot_config.legend)
	println(io, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

	print(io, "\\addplot[color=")
	print(io, plot_config.color)
	print(io, ",opacity=")
	print(io, plot_config.opacity)
	print(io, ",mark=")
	print(io, plot_config.mark)
	if !isnothing(plot_config.property)
		print(io, ", ", plot_config.property)
	end
	if !isnothing(plot_config.mark_options)
		print(io, ",mark options={")
		print(io, plot_config.mark_options)
		print(io, "}, ")
	end
	print(io, ",mark repeat=")
	print(io, plot_config.mark_repeat)
	if !isnothing(plot_config.mark_size)
		print(io, ",mark size=")
		print(io, plot_config.mark_size)
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
		error("Plot type not defined.")
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
	print(io, "\\label{")
	print(io, plot_config.legend)
	print(io, "}")
	println(io)
	println(io)
	println(io)

	close(io)

end

function create_latex_footer(file_path_latex::String)

	io = open(file_path_latex, "a")

	println(io, "\\end{axis}")
	println(io)
	println(io, "\\matrix[")
	println(io, "matrix of nodes,")
	println(io, "anchor=north,")
	println(io, "inner sep=0.2em,")
	println(io, "column 1/.style={nodes={anchor=center}},")
	println(io, "column 2/.style={nodes={anchor=west},font=\\scriptsize},")
	println(io, "column 3/.style={nodes={anchor=center}},")
	println(io, "column 4/.style={nodes={anchor=west},font=\\scriptsize},")
	println(io, "column 5/.style={nodes={anchor=center}},")
	println(io, "column 6/.style={nodes={anchor=west},font=\\scriptsize},")
	println(io, "]")
	println(io, "at(current bounding box.south){")
	
	println(io, raw"\ref{Single} & \texttt{SB + L} (single) & \ref{L1} & \texttt{SB +} $\ell_1$-deep & \ref{Mid} & \texttt{SB + LN-Mid} \\ ")
	println(io, raw"\ref{Multi} & \texttt{SB + L} (multi) & \ref{Lsup} & \texttt{SB +} $\ell_\infty$-deep & \ref{Eps} & \texttt{SB + LN-Eps} \\ ")
	println(io, raw"\ref{B} & \texttt{B} (single) & \ref{L1sup} & \texttt{SB +} $\ell_{1 \infty}$-deep & \ref{Relint} & \texttt{SB + LN-Relint} \\ ")
	println(io, raw"\ref{BMulti} & \texttt{B} (multi) & & & \ref{Conv(50)} & \texttt{SB + LN-Conv(50)} \\ ")
	println(io, raw"\ref{SB} & \texttt{SB} (single) & & & \ref{Conv(75)} & \texttt{SB + LN-Conv(75)} \\ ")
	println(io, raw"\ref{SBMulti} & \texttt{SB} (multi) & & & \ref{Conv(90)} & \texttt{SB + LN-Conv(90)} \\ ")
	
	# println(io, raw"\ref{Single} & \texttt{L} (single) & \ref{L1} & $\ell_1$-deep & \ref{Mid} & \texttt{LN-Mid} \\ ")
	# println(io, raw"\ref{Multi} & \texttt{L} (multi) & \ref{Lsup} & $\ell_\infty$-deep & \ref{Eps} & \texttt{LN-Eps} \\ ")
	# println(io, raw"\ref{B} & \texttt{B} (single) & \ref{L1sup} & $\ell_{1 \infty}$-deep & \ref{Relint} & \texttt{LN-Relint} \\ ")
	# println(io, raw"\ref{BMulti} & \texttt{B} (multi) & & & \ref{Conv(50)} & \texttt{LN-Conv(50)} \\ ")
	# println(io, raw"\ref{SB} & \texttt{SB} (single) & & & \ref{Conv(75)} & \texttt{LN-Conv(75)} \\ ")
	# println(io, raw"\ref{SBMulti} & \texttt{SB} (multi) & & & \ref{Conv(90)} & \texttt{LN-Conv(90)} \\ ")
	
	# println(io, raw"\ref{Single} & \texttt{L} (single) & \ref{L1} & $\ell_1$-deep & \ref{Mid} & \texttt{LN-Mid} \\ ")
	# println(io, raw"\ref{Multi} & \texttt{L} (multi) & \ref{LsupMNC} & $\ell_\infty$-deep & \ref{Eps} & \texttt{LN-Eps} \\ ")
	# println(io, raw"\ref{B} & \texttt{B} (single) & \ref{L1sup} & $\ell_{1 \infty}$-deep & \ref{Relint} & \texttt{LN-Relint} \\ ")
	# println(io, raw"\ref{BMulti} & \texttt{B} (multi) & & & \ref{InOut} & \texttt{LN-In-Out} \\ ")
	# println(io, raw"\ref{SB} & \texttt{SB} (single) \\ ")
	# println(io, raw"\ref{SBMulti} & \texttt{SB} (multi) \\ ")
	
	println(io, "};")
	println(io)
	println(io)

	println(io, "\\end{tikzpicture}")
	println(io, "\\end{document}")

	close(io)

end

create_latex_plots()
