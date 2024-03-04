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

    file_path_latex = "C:/Users/cg4102/Documents/julia_plots/Paper 2024/CLSP_RR_Relint_lb_time.tex"
	plottype = :lb_time

    # Create header of LaTeX file
	create_latex_header(file_path_latex, [0.0, 500.0, 0.0, 7000.0], plottype)

    # Create an array of struct instances which define the plots to be created
	plots_to_create = (
		###################################################################### CLSP - 4 - 20 - 3 #############################################################################################################################
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_Single_SB.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_Mid_SB.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_InOut_SB.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_4_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_6_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_10_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_10_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_10_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2N_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"LsupMNC", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_CL_Span.log", plottype, "yellow", 1.0, "none", "solid", raw"Span", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_CL_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_10_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Single.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_InOut.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_B.log", plottype, "olive", 1.0, "none", "dashed", raw"B", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_SB.log", plottype, "red", 1.0, "none", "dashed", raw"SB", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"BMulti", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_10_20_SB_Multi.log", plottype, "red", 1.0, "none", "solid", raw"SBMulti", "", nothing),

		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Single_SB.log", plottype, "black", 1.0, "none", "dashed", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Mid_SB.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_InOut_SB.log", plottype, "orange", 1.0, "none", "solid", raw"InOut", "", nothing),

		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		#PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Updated_with_Bin/CLSP_16_20_prim_CL_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_CL_Span.log", plottype, "green!70!black", 1.0, "none", "solid", raw"Span", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR_16_20_CL_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_CL_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_CL_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_RR2_16_20_CL_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),

		
		# ###################################################################### CLSP - 16 - 20 - 10 #############################################################################################################################
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Single.log", plottype, "red", 1.0, "none", "solid", raw"Single", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Multi.log", plottype, "black", 1.0, "none", "solid", raw"Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1 bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC bd" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O2_16_20_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint bd", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_B.log", plottype, "olive", 1.0, "none", "solid", raw"B", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_SB.log", plottype, "purple", 1.0, "none", "solid", raw"SB", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_B_Multi.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_SB_Multi.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Single_SB.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Multi_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1sup_SB.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Lsup_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21)" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Mid_SB.log", plottype, "violet", 1.0, "none", "solid", raw"SB + Mid bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21)", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21)", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Single_SB_CS.log", plottype, "red", 1.0, "none", "solid", raw"SB + Single (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Multi_SB_CS.log", plottype, "black", 1.0, "none", "solid", raw"SB + Multi (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1_SB_CS.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1sup_SB_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"SB + L1sup bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Lsup_SB_CS.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Mid_SB_CS.log", plottype, "violet", 1.0, "none", "solid", raw"SB + Mid bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Eps_SB_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Relint_SB_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_B_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_SB_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_B_Multi_CS.log", plottype, "olive", 1.0, "none", "solid", raw"B Multi + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_SB_Multi_CS.log", plottype, "purple", 1.0, "none", "solid", raw"SB Multi + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Single_MNC_CS.log", plottype, "red", 1.0, "none", "solid", raw"Single MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Multi_MNC_CS.log", plottype, "black", 1.0, "none", "solid", raw"Multi MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1_MNC_CS.log", plottype, "blue", 1.0, "none", "solid", raw"L1 MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_L1sup_MNC_CS.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Lsup_MNC_CS.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC + CS" "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Mid_MNC_CS.log", plottype, "violet", 1.0, "none", "solid", raw"Mid MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Eps_MNC_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps MNC + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_Relint_MNC_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint MNC + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_L1.log", plottype, "blue", 1.0, "none", "solid", raw"L1", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_L1sup.log", plottype, "green!70!black", 1.0, "none", "solid", raw"L1sup", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_Lsup.log", plottype, "green", 1.0, "none", "solid", raw"Lsup MNC", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_Mid.log", plottype, "violet", 1.0, "none", "solid", raw"Mid", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_Eps.log", plottype, "magenta", 1.0, "none", "solid", raw"Eps", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_prim_Relint.log", plottype, "cyan", 1.0, "none", "solid", raw"Relint", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_L1_CS.log", plottype, "blue", 1.0, "none", "solid", raw"B + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Span_CS.log", plottype, "black", 1.0, "none", "solid", raw"B + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Lsup_CS.log", plottype, "green", 1.0, "none", "solid", raw"B + Lsup MNC bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Mid_CS.log", plottype, "violet", 1.0, "none", "solid", raw"B + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Eps_CS.log", plottype, "magenta", 1.0, "none", "solid", raw"B + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Relint_CS.log", plottype, "cyan", 1.0, "none", "solid", raw"B + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_L1_CS_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Span_CS_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Lsup_CS_SB.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Mid_CS_SB.log", plottype, "violet", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Eps_CS_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20 + CS", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Relint_CS_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20 + CS", "", nothing),
		# #
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_L1_SB.log", plottype, "blue", 1.0, "none", "solid", raw"SB + L1 bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Span_SB.log", plottype, "black", 1.0, "none", "solid", raw"SB + Span bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Lsup_SB_MNC.log", plottype, "green", 1.0, "none", "solid", raw"SB + Lsup MNC bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Mid_SB.log", plottype, "violet", 1.0, "none", "solid", raw"SB + Mid bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Eps_SB.log", plottype, "magenta", 1.0, "none", "solid", raw"SB + Eps bd (21) with CL-20", "", nothing),
		# PlotConfig("C:/Users/cg4102/Documents/julia_logs/Restructuring/Final paper/CLSP/CLSP_Large_O_16_20_CL_Relint_SB.log", plottype, "cyan", 1.0, "none", "solid", raw"SB + Relint bd (21) with CL-20", "", nothing),

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
	println(io, raw"\ref{Single} & \texttt{L} (single) & \ref{L1} & $\ell_1$-deep & \ref{Mid} & \texttt{LN-Mid} \\ ")
	println(io, raw"\ref{Multi} & \texttt{L} (multi) & \ref{LsupMNC} & $\ell_\infty$-deep & \ref{Eps} & \texttt{LN-Eps} \\ ")
	println(io, raw"\ref{B} & \texttt{B} (single) & \ref{L1sup} & $\ell_{1 \infty}$-deep & \ref{Relint} & \texttt{LN-Relint} \\ ")
	println(io, raw"\ref{BMulti} & \texttt{B} (multi) & & & \ref{InOut} & \texttt{LN-In-Out} \\ ")
	println(io, raw"\ref{SB} & \texttt{SB} (single) \\ ")
	println(io, raw"\ref{SBMulti} & \texttt{SB} (multi) \\ ")
	println(io, "};")
	println(io)
	println(io)

	println(io, "\\end{tikzpicture}")
	println(io, "\\end{document}")

	close(io)

end

create_latex_plots()
