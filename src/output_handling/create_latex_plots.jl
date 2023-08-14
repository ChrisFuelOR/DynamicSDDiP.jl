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

    file_path_latex = "C:/Users/cg4102/Documents/julia_plots/CLSP_RR_16_20_Large_lag_it.tex"
	plottype = :lag_it

    # Create header of LaTeX file
	create_latex_header(file_path_latex, [0.0, 200.0, 0.0, 25000.0], plottype)

    # Create an array of struct instances which define the plots to be created
	plots_to_create = (
		###################################################################### CLSP - 4 - 20 - 3 #############################################################################################################################
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Single.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_InOut.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_SB.log", plottype, "red", 1.0, "none", "dashed", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_SB_Multi.log", plottype, "red", 1.0, "none", "solid", raw"SB Multi", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Single_SB.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_InOut_SB.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Updated_with_Bin/CLSP_16_20_prim_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_CL_Span.log", plottype, "green!70!black", 1.0, "none", "solid", raw"Span", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR_16_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_CL_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_RR2_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),

		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Single.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_InOut.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Single_SB.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_InOut_SB.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_SB.log", plottype, "red", 1.0, "none", "dashed", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_SB_Multi.log", plottype, "red", 1.0, "none", "solid", raw"SB Multi", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_CL2_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR_Large_O_16_20_CL2_Span.log", plottype, "green!70!black", 1.0, "none", "solid", raw"Span", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_CL2_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_CL2_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_CL2_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP Large/CLSP_RR2_Large_O_16_20_CL2_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_4_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint (21)", "", nothing),
		#
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint with LB (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB with LB (21) + Single (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB with LB (21) + Multi (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB with LB (21) + L1 (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB with LB (21) + L1sup (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB with LB (21) + Lsup (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB with LB (21) + Mid (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB with LB (21) + Eps (41)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB with LB (21) + Relint (41)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps (21) with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB2_4_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint (21) with LB (31)", "", nothing),
		# #
		# # #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_B.log", plottype, "olive", 1.0, "none", "solid", raw"B with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_SB.log", plottype, "purple", 1.0, "none", "solid", raw"SB with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_4_20_SB_Multi.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi with LB (21)", "", nothing),
		# # #
		# ###################################################################### CLSP - 16 - 20 - 3 #############################################################################################################################
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Single_CS.log", plottype, "red", 1.0, "none", "solid", raw"Single CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Multi_CS.log", plottype, "black", 1.0, "none", "solid", raw"Multi CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_L1_CS2.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_L1sup_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Lsup_CS2.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Mid_CS2.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_InOut_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"InOut bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Eps_CS2.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Opt_CS.log", plottype, "black", 1.0, "none", "solid", raw"Opt bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_Relint_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_B_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_SB_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_B_Multi_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_SB_Multi_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd Kel" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid bd Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd Kel", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_SB.log", plottype, "purple", 1.0, "none", "dashed", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "dotted", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_SB_Multi.log", plottype, "purple", 1.0, "none", "dotted", raw"SB Multi", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O2_16_20_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"B + L1 bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Span.log", plottype, "black", 1.0, "none", "solid", raw"B + Span bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"B + Lsup MNC bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"B + Mid bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"B + Eps bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"B + Relint bd (21) with CL-20", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_L1_CS.log", plottype, "blue", 1.0, "none", "solid", raw"B + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Span_CS.log", plottype, "black", 1.0, "none", "solid", raw"B + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Lsup_CS.log", plottype, "green", 1.0, "none", "solid", raw"B + Lsup MNC bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Mid_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"B + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Eps_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"B + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_ON_16_20_CL_Relint_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"B + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"B + L1 (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_Span.log", plottype, "black", 1.0, "none", "solid", raw"B + Span (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"B + Lsup MNC (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"B + Mid (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"B + Eps (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_O_16_20_prim_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"B + Relint (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps with LB (31)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_LB_16_20_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint with LB (31)", "", nothing),
		# #
		# ###################################################################### CLSP - 16 - 20 - 3 MORE ANALYSES ################################################################################################################
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_SB1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 -- SB(1),L(1)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_SB21.log", plottype, "blue!50!white", 1.0, "none", "solid", raw"L1 -- SB(1),L(21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_B1.log", plottype, "teal", 1.0, "none", "solid", raw"L1 -- B(1),L(1)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_B21.log", plottype, "cyan", 1.0, "none", "solid", raw"L1 -- B(1),L(21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_B_10.log", plottype, "blue", 1.0, "none", "dashed", raw"L1 -- B(1),L(11) -- CL(10)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_SB_10.txt", plottype, "blue!50!white", 1.0, "none", "dashed", raw"L1 -- SB(1),L(11) -- CL(10)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_B_20.log", plottype, "teal", 1.0, "none", "dashed", raw"L1 -- B(1),L(21) -- CL(20)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_SB_20.log", plottype, "cyan", 1.0, "none", "dashed", raw"L1 -- SB(1),L(21) -- CL(20)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_B_30.log", plottype, "green!50!black", 1.0, "none", "dashed", raw"L1 -- B(1),L(31) -- CL(30)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_CL_SB_30.log", plottype, "green", 1.0, "none", "dashed", raw"L1 -- SB(1),L(31) -- CL(30)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21.log", plottype, "green!50!black", 1.0, "none", "solid", raw"L1 -- LB(10,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_11.log", plottype, "green", 1.0, "none", "solid", raw"L1 -- LB(10,11)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31.log", plottype, "lime", 1.0, "none", "solid", raw"L1 -- LB(10,31)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_5_21.log", plottype, "yellow", 1.0, "none", "solid", raw"L1 -- LB(5,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_5_11.log", plottype, "pink", 1.0, "none", "solid", raw"L1 -- LB(5,11)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_5_31.log", plottype, "lightgray", 1.0, "none", "solid", raw"L1 -- LB(5,31)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_subopt.log", plottype, "gray", 1.0, "none", "solid", raw"L1 -- LB(10,21) subopt", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_Lag50.log", plottype, "darkgray", 1.0, "none", "solid", raw"L1 -- LB(10,21) Lag50", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_Tol.log", plottype, "black", 1.0, "none", "solid", raw"L1 -- LB(10,21) tol", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31_Lag50.log", plottype, "olive", 1.0, "none", "solid", raw"L1 -- LB(10,31) Lag50", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31_Lag50_onlyLB.log", plottype, "brown", 1.0, "none", "solid", raw"L1 -- LB(10,31) Lag50 after LB", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31_Tol.log", plottype, "orange", 1.0, "none", "solid", raw"L1 -- LB(10,31) tol", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_SB_31.log", plottype, "red", 1.0, "none", "solid", raw"L1 -- SB(1),L(31) -- LB(10,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31_SB_21.log", plottype, "purple", 1.0, "none", "solid", raw"L1 -- SB(1),L(21) -- LB(10,31)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_SB_1.log", plottype, "magenta", 1.0, "none", "solid", raw"L1 -- SB(1),L(1) -- LB(10,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_31_SB_1.log", plottype, "violet", 1.0, "none", "solid", raw"L1 -- SB(1),L(1) -- LB(10,31)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_SB_31_CL.log", plottype, "black", 1.0, "none", "dashed", raw"L1 -- SB(1),L(31) -- CL(10) -- LB(10,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_11_SB_31_CL_v2.log", plottype, "brown", 1.0, "none", "dashed", raw"L1 -- SB(1),L(31) -- CL(10) -- LB(10,11)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_21_B_31_CL.log", plottype, "red", 1.0, "none", "dashed", raw"L1 -- B(1),L(31) -- CL(10) -- LB(10,21)", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_OT_16_20_L1_LB_10_11_B_31_CL_v2.log", plottype, "purple", 1.0, "none", "dashed", raw"L1 -- B(1),L(31) -- CL(20) -- LB(10,11)", "", nothing),

		# ###################################################################### CLSP - 16 - 20 - 10 #############################################################################################################################
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O2_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_B.log", plottype, "olive", 1.0, "none", "solid", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_SB.log", plottype, "purple", 1.0, "none", "solid", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_SB_Multi.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21)" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Single_SB_CS.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Multi_SB_CS.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1_SB_CS.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1sup_SB_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Lsup_SB_CS.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Mid_SB_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Eps_SB_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Relint_SB_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_B_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_SB_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_B_Multi_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_SB_Multi_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Single_MNC_CS.log", plottype, "red", 1.0, "none", "solid", raw"Single MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Multi_MNC_CS.log", plottype, "black", 1.0, "none", "solid", raw"Multi MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1_MNC_CS.log", plottype, "blue", 1.0, "none", "solid", raw"L1 MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_L1sup_MNC_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Lsup_MNC_CS.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Mid_MNC_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Eps_MNC_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_Relint_MNC_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint MNC + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_L1_CS.log", plottype, "blue", 1.0, "none", "solid", raw"B + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Span_CS.log", plottype, "black", 1.0, "none", "solid", raw"B + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Lsup_CS.log", plottype, "green", 1.0, "none", "solid", raw"B + Lsup MNC bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Mid_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"B + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Eps_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"B + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Relint_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"B + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_L1_CS_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Span_CS_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Lsup_CS_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Mid_CS_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Eps_CS_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Relint_CS_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Span_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Lsup_SB_MNC.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_O_16_20_CL_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Large_LB_16_20_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint with LB (21)", "", nothing),
		# #
		# ###################################################################### CLSP - 16 - 50 - 3 #############################################################################################################################
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Single_2h.log", plottype, "red", 1.0, "none", "solid", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Multi_2h.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1_2h.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1sup_2h.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Lsup_MNC_2h.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Mid_2h.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Eps_2h.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Relint_2h.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_SB.log", plottype, "purple", 1.0, "none", "dashed", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_B_Multi.log", plottype, "olive", 1.0, "none", "dotted", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_SB_Multi.log", plottype, "purple", 1.0, "none", "dotted", raw"SB Multi", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Lsup_MNC_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (1)" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (1)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (1)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Lsup_MNC_SB2.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21)" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1_SB2.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1sup_SB2.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Mid_SB2.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Eps_SB2.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Relint_SB2.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Single_CS.log", plottype, "red", 1.0, "none", "solid", raw"Single + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Multi_CS.log", plottype, "black", 1.0, "none", "solid", raw"Multi + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1_CS.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1sup_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Lsup_MNC_CS.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Mid_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid bd + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Eps_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Relint_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_B_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_SB_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_B_Multi_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_SB_Multi_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Single_SB_CS.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_MultiSB__CS.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1_SB_CS.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_L1sup_SB_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Lsup_MNC_SB_CS.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Mid_SB_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Eps_SB_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_Relint_SB_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Span_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_L1_SB_CS.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Span_SB_CS.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Lsup_SB_CS.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20 + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Mid_SB_CS.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Eps_SB_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_O_16_50_CL_Relint_SB_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_16_50_prim_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Mid.log", plottype, "yellow", 1.0, "none", "solid", raw"Mid (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint (21) with LB (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Mid_SB.log", plottype, "yellow", 1.0, "none", "solid", raw"SB + Mid (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps (21) with LB (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/New after handling numerical issues (RR Re-run)/CLSP/CLSP_Real_LB_16_50_prim_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint (21) with LB (21)", "", nothing),
	)

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

	println(io, "legend columns=2,") #TODO
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
	println(io, "%", plot_config.legend)
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
