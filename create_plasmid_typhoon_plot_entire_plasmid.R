# Clear the workspace
rm(list = ls())
graphics.off()

library(nucleR)

# Load the functions
source("/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/get_peak_signal.function.R")
source("/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/create_gene_schematic_on_typhoon_plot.function.R")
source("/data/illumina_pipeline/scripts/create_typhoon_plots/create_typhoon_plot_functions.R")

# Set the dataset
dataset_id = 18:24
dataset = c("2677", "2678", "2679", "2680", "2681", "2682", "2683")

# Iterate through each dataset
for(i in 1:length(dataset_id)){

	# Set the plot
	#	If want transparent, add 'bg = "transparent" '
	png(file = paste("/data/collaborator_data/bell/azmi_ishara/2015_07_16_azmi_ishara/", 
					 "output_plots/entire_plasmid/", dataset[i], "_plasmid.png", sep = ""),
		width = 15, height = 7.5, units = "in", res = 400) # bg = "transparent")

	# Set the axes
	feature_chr = "pUC19_ARS1_3807bp"
	feature_chr = "pARS_Belsky"
	x_start = 0
	x_end = 4000
	bam_file_name = paste("/data/collaborator_data/aligned_experiments/CD", dataset_id[i], 
						  "/CD", dataset_id[i], ".bam", sep = ""
						 )

	# Setup the screen
	scr.m = matrix(c(0, 1, 0.95, 1,
					 0, 1, 0.85, 0.95,
					 0, 1, 0.5, 0.85,
					 0, 1, 0, 0.5),
				   ncol = 4, byrow = T
				  )

	# Close the screens
	close.screen(all.screens = T)

	# Split the screens
	s_num.v = split.screen(scr.m)

	# Set the x_left parameter
	x_left = 4.5
	x_right = 4.5

	# Set up the plot
	screen(s_num.v[1])
	par(mar = c(0, x_left, 0, x_right))

	# Make the gene schematic
	#	If want transparent, add 'bg = "transparent" '
	plot(0, 0, type = "n", bty = "n", # bg = "transparent",
		 xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
		 ylim = c(0, 1), yaxs = "i", yaxt = "n",
		 ann = F
		)

	text(x = x_start + 10, y = 0.5, adj = 0, labels = dataset[i], cex = 1.5)

	text(x = 916, y = 0.5, labels = "ACS", srt = 90, col = "darkgreen", cex = 1.125)
	# text(x = 945, y = 0.5, labels = "B1", srt = 90, col = "darkgreen", cex = 1.125)
	# text(x = 980, y = 0.5, labels = "B2", srt = 90, col = "darkgreen", cex = 1.125)
	# text(x = 1023, y = 0.5, labels = "B3", srt = 90, col = "darkgreen", cex = 1.125)

	# Set up the chromatin schematic
	screen(s_num.v[2])
	par(mar = c(0, x_left, 0, x_right))

	set_chromatin_schematic(x_start, x_end)

	# Get the nucleosome feature density
	nuc.v = get_mnase_feature_density(bam_file_name, feature_chr, x_start, x_end, 125, 175, 20)

	# Get the subnuc peaks
	nuc_peaks.df = get_mod_peaks(nuc.v, (x_start + x_end)/2, min_thresh = 0.25, peak_width = 75)
	 
	# Plot the nucleosome
	plot_nucleosome(nuc_peaks.df, 1.5, 0.5, 0.2)

	# Plot the nucleosome track
	screen(s_num.v[3])
	par(mar = c(2, x_left, 2, x_right), cex = 1.25)

	plot(x_start:x_end, nuc.v, type = "l", col = "red",
		 xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
		 ylim = c(0, max(nuc.v)),
		 ylab = "Nuc density", xlab = NA
		)
	make_typhoon_plot_axes(x_start, x_end, x_step = 250, x_digits = 2, addXAxis = T, addYAxis = F)

	# Get the nuc_pos relative to the start site
	abline(v = nuc_peaks.df$pos, lty = 2)
	axis(3, at = nuc_peaks.df$pos, cex.axis = 0.75, las = 3)

	# Get the typhoon plot
	screen(s_num.v[4])
	par(mar = c(5, x_left, 2, x_right), cex = 1.25)

	# Get the mat
	mat.m = get_typhoon_plot_mat_java(bam_file_name, feature_chr, x_start, x_end)
	colnames(mat.m) = x_start:x_end

	# Normalize to 10E6 reads
	mat.m = mat.m * (10E6) / get_total_read_number(bam_file_name, feature_chr)

	# Make the typhoon plot
	dens_dot_plot(mat.m, z_max = 5000, 
				  x_label = "Plasmid coordinate (kb)", y_label = "Fragment length (bp)", plot_title = "",
				  x_axt = "n", y_axt = "n" 
				 )
	make_typhoon_plot_axes(x_start, x_end, x_step = 250, x_digits = 2, addXAxis = T, addYAxis = T)

	# Close the screens
	close.screen(all.screens = T)

	# Close the device
	dev.off()

}
