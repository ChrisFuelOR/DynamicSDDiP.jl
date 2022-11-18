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

    file_path_latex = "C:/Users/cg4102/Documents/julia_plots/CLSP_Large_O_16_20_lb_time.tex"
	plottype = :lb_time

    # Create an array of struct instances which define the plots to be created
    plots_to_create = (
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_ON_16_20_Single.log", plottype, "red", 1.0, "none", "solid", "Single-cut", "SDDiP Single-cut", nothing),
        #PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_ON_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", "Multi-cut", "SDDiP Multi-cut", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Single.log", plottype, "red", 1.0, "none", "solid", "Single-cut", "SDDiP Single-cut", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", "Multi-cut", "SDDiP Multi-cut", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"$ \ell_1$", "L1", nothing),
        PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1sup.log", plottype, "green", 1.0, "none", "solid", raw"$ \ell_{1 \infty}$", "L1sup", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Lsup_MNC.log", plottype, "green!70!black", 1.0, "none", "solid", raw"$ \ell_\infty MNC$", "Lsup", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", "Mid", "Mid", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", "Eps", "Eps", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", "Relint", "Relint", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt.txt", plottype, "violet", 1.0, "none", "solid", "Opt", "Opt", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut.txt", plottype, "orange", 1.0, "none", "solid", "In-Out", "In-Out", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_B.log", plottype, "purple", 1.0, "none", "solid", "B", "B", 5000),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_B_Multi.log", plottype, "purple", 1.0, "none", "dashed", "B Multi", "B Multi", 5000),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_SB.log", plottype, "olive", 1.0, "none", "solid", "SB", "SB", 5000),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_SB_Multi.log", plottype, "olive", 1.0, "none", "dashed", "SB Multi", "SB Multi", 5000),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Mid_Bd.txt", plottype, "yellow", 1.0, "none", "dash dot", "Mid Bd", "Mid Bd", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Eps_Bd.txt", plottype, "magenta", 1.0, "none", "dash dot", "Eps Bd", "Eps Bd", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Relint_Bd.txt", plottype, "cyan", 1.0, "none", "dash dot", "Relint Bd", "Relint Bd", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt_Bd.txt", plottype, "violet", 1.0, "none", "dash dot", "Opt Bd", "Opt Bd", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut_Bd.txt", plottype, "orange", 1.0, "none", "dash dot", "In-Out Bd", "In-Out Bd", 50),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Single_CS.log", plottype, "red", 1.0, "none", "dashed", "Single-cut CS", "SDDiP Single-cut CS", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Multi_CS.log", plottype, "black", 1.0, "none", "dashed", "Multi-cut CS", "SDDiP Multi-cut CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_L1_NumFocus.log", plottype, "blue", 0.5, "none", "solid", raw"$ \ell_1 $ NF", "L1 NF", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_L1sup_CS.log", plottype, "green", 1.0, "none", "dashed", raw"$ \ell_{1 \infty} $ CS", "L1sup CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Lsup_CS.log", plottype, "green!70!black", 1.0, "dashed", "solid", raw"$ \ell_\infty$ CS", "Lsup CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Mid_CS.log", plottype, "yellow", 1.0, "none", "dashed", "Mid CS", "Mid CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Eps_CS.log", plottype, "magenta", 1.0, "none", "dashed", "Eps CS", "Eps CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_16_20_Relint_CS.log", plottype, "cyan", 1.0, "none", "dashed", "Relint CS", "Relint CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt_CS.txt", plottype, "violet", 1.0, "none", "dashed", "Opt CS", "Opt _CS", 50),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut_CS.txt", plottype, "orange", 1.0, "none", "dashed", "In-Out CS", "In-Out CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Single_MNC.txt", plottype, "red", 1.0, "none", "dotted", "Single-cut MNC", "SDDiP Single-cut MNC", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Multi_MNC.txt", plottype, "black", 1.0, "none", "dotted", "Multi-cut MNC", "SDDiP Multi-cut MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1_MNC.txt", plottype, "blue", 1.0, "none", "dotted", raw"$ \ell_1 $ MNC", "L1 MNC", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1sup_MNC.txt", plottype, "green", 1.0, "none", "dotted", raw"$ \ell_{1 \infty} $ MNC", "L1sup MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Lsup_MNC.txt", plottype, "green!70!black", 1.0, "dotted", "solid", raw"$ \ell_\infty$ MNC", "Lsup MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Mid_MNC.txt", plottype, "yellow", 1.0, "none", "dotted", "Mid MNC", "Mid MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Eps_MNC.txt", plottype, "magenta", 1.0, "none", "dotted", "Eps MNC", "Eps MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Relint_MNC.txt", plottype, "cyan", 1.0, "none", "dotted", "Relint MNC", "Relint MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt_MNC.txt", plottype, "violet", 1.0, "none", "dotted", "Opt MNC", "Opt MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut_MNC.txt", plottype, "orange", 1.0, "none", "dotted", "In-Out MNC", "In-Out MNC", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Single_Bundle_CS.txt", plottype, "red", 1.0, "none", "dotted", "Single-cut Bun CS", "SDDiP Single-cut Bun", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Multi_Bundle_CS.txt", plottype, "black", 1.0, "none", "dotted", "Multi-cut Bun CS", "SDDiP Multi-cut Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_L1_bd_Bundle_CS.txt", plottype, "blue", 1.0, "none", "dotted", raw"$ \ell_1 $ Bun CS", "L1 Bun", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_L1sup_bd_Bundle_CS.txt", plottype, "green", 1.0, "none", "dotted", raw"$ \ell_{1 \infty} $ Bun CS", "L1sup Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Lsup_bd_Bundle_CS_MNC.txt", plottype, "green!70!black", 1.0, "none", "solid", raw"$ \ell_\infty$ MNC Bun CS", "Lsup MNC Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Mid_bd_Bundle_CS.txt", plottype, "yellow", 1.0, "none", "dotted", "Mid Bd Bun CS", "Mid Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Eps_bd_Bundle_CS.txt", plottype, "magenta", 1.0, "none", "dotted", "Eps Bd Bun CS", "Eps Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_16_20_Relint_bd_Bundle_CS.txt", plottype, "cyan", 1.0, "none", "dotted", "Relint Bd Bun CS", "Relint Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt_Bundle2.txt", plottype, "violet", 1.0, "none", "dotted", "Opt Bun", "Opt Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut_Bundle2.txt", plottype, "orange", 1.0, "none", "dotted", "In-Out Bun", "In-Out Bun", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Single_Bundle_B.txt", plottype, "red", 1.0, "none", "dashed", "Single-cut Bun + SB", "SDDiP Single-cut Bun + SB", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Multi_Bundle_B2.txt", plottype, "black", 1.0, "none", "dashed", "Multi-cut Bun + SB", "SDDiP Multi-cut Bun + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1_Bundle_B2.txt", plottype, "blue", 1.0, "none", "dashed", raw"$ \ell_1 $ Bun + SB", "L1 Bun + SB", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1sup_Bundle_B2.txt", plottype, "green", 1.0, "none", "dashed", raw"$ \ell_{1 \infty} $ Bun + SB", "L1sup Bun + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Lsup_MNC_Bundle_B2.txt", plottype, "green!70!black", 1.0, "none", "dashed", raw"$ \ell_\infty$ MNC Bun + SB", "Lsup MNC Bun + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Mid_Bundle2_B2.txt", plottype, "yellow", 1.0, "none", "dashed", "Mid Bun + SB", "Mid Bun + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Eps_bd_Bundle2_B2.txt", plottype, "magenta", 1.0, "none", "dashed", "Eps Bun + SB", "Eps Bun + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Relint_bd_Bundle2_B2.txt", plottype, "cyan", 1.0, "none", "dashed", "Relint Bun + SB", "Relint Bun + SB", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_L1_CS_SB.log", plottype, "blue", 1.0, "none", "solid", raw"CL + $ \ell_1 $", "CL L1", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_Lsup_CS_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"$CL + \ell_\infty$", "CL + Lsup MNC", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_Mid_CS_SB.log", plottype, "yellow", 1.0, "none", "solid", "CL + Mid", "CL + Mid", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_Span_CS_SB.log", plottype, "teal", 1.0, "none", "solid", "CL + Span", "CL + Span", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_Eps_CS_SB.log", plottype, "magenta", 1.0, "none", "solid", "CL + Eps", "CL + Eps", nothing),
		##PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_CL_Relint_CS_SB.log", plottype, "cyan", 1.0, "none", "solid", "CL + Relint", "CL + Relint", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Single_Bundle2_CS_B (ws02).txt", plottype, "red", 1.0, "none", "solid", "Single-cut Bun CS + SB", "SDDiP Single-cut CS", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Multi_Bundle2_CS_B (ws02).txt", plottype, "black", 1.0, "none", "solid", "Multi-cut Bun CS + SB", "SDDiP Multi-cut CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Lsup_Bundle2_MNC_CS_B (ws02).txt", plottype, "green!70!black", 1.0, "none", "solid", raw"$ \ell_\infty$ Bun CS", "Lsup CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1_Bundle2_CS_B (ws02).txt", plottype, "blue", 1.0, "none", "solid", raw"$ \ell_1 $ Bun CS + SB", "L1 CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_L1sup_Bundle2_CS_B (ws02).txt", plottype, "green", 1.0, "none", "solid", raw"$ \ell_{1 \infty} $ Bun CS + SB", "L1sup CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Mid_Bundle2_CS_B (ws02).txt", plottype, "yellow", 1.0, "none", "solid", "Mid Bun CS + SB", "Mid CS + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Eps_Bundle2_CS_B (ws02).txt", plottype, "magenta", 1.0, "none", "solid", "Eps Bun CS + SB", "Eps CS + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Relint_Bundle2_CS_B (ws02).txt", plottype, "cyan", 1.0, "none", "solid", "Relint Bun CS + SB", "Relint CS + SB", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_InOut_Bundle2_CS_B (ws02).txt", plottype, "orange", 1.0, "none", "solid", "In-Out Bun CS", "In-Out CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_10_20_Opt_Bundle2_CS_B (ws02).txt", plottype, "violet", 1.0, "none", "solid", "Opt Bun CS", "Opt CS", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Single_MNC.log", plottype, "red", 1.0, "none", "solid", "Single-cut", "SDDiP Single-cut", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Multi_MNC.log", plottype, "black", 1.0, "none", "solid", "Multi-cut", "SDDiP Multi-cut", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1_MNC.log", plottype, "blue", 1.0, "none", "solid", raw"$ \ell_1 $", "L1", nothing),
        # PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_L1sup_MNC.log", plottype, "green", 1.0, "none", "solid", raw"$ \ell_{1 \infty} $", "L1sup", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Lsup_MNC.log", plottype, "green!70!black", 1.0, "none", "solid", raw"$ \ell_\infty$", "Lsup", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Mid_MNC.log", plottype, "yellow", 1.0, "none", "solid", "Mid", "Mid", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Eps_MNC.log", plottype, "magenta", 1.0, "none", "solid", "Eps", "Eps", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_Relint_MNC.log", plottype, "cyan", 1.0, "none", "solid", "Relint", "Relint", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_B.log", plottype, "purple", 1.0, "none", "solid", "B", "B", 100),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_B_Multi.log", plottype, "purple", 1.0, "none", "dashed", "B Multi", "B Multi", 100),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_SB.log", plottype, "olive", 1.0, "none", "solid", "SB", "SB", 100),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/CLSP_Large_O_16_20_SB_Multi.log", plottype, "olive", 1.0, "none", "dashed", "SB Multi", "SB Multi", 100),
    )

    # Create header of LaTeX file
	create_latex_header(file_path_latex, [0.0, 70.0, 0.0, 30000.0], plottype)

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

function create_latex_header(file_path_latex::String, axis_limits::Vector{Float64}, plottype::Symbol)

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
