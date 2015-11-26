# Clear the workspace
rm(list = ls())
graphics.off()

library(nucleR)

# Load the functions
source("/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/get_peak_signal.function.R")
source("/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/create_gene_schematic_on_typhoon_plot.function.R")
source("/data/illumina_pipeline/scripts/create_typhoon_plots/create_typhoon_plot_functions.R")
source("/data/illumina_pipeline/scripts/process_bam_scripts/process_bam_functions.R")

# Set the dataset
dataset = paste("CD", 3:10, sep = "")
dataset_name = c("7567_ISW1a", "7568_ISW1b", "7569_ISW2", "7570_INO80", 
				 "7571_SWI-SNF", "7572_RSC", "7573_CHD1", "7574_Nap1"
				)

# Set the plot
#	If want transparent, add 'bg = "transparent" '
png(file = paste("/data/collaborator_data/bell/azmi_ishara/2015_07_16_azmi_ishara/output_plots/",
				 "plasmid_center_acs/compare_nucleosome_structure.png", sep = ""),	
    width = 10, height = 7.5, units = "in", res = 400) # bg = "transparent")

# Set the axes
feature_chr = "pUC19_ARS1_3807bp"
x_start = 916 - 350
x_end = 916 + 350

# bam_file_name = paste("/data/collaborator_data/bell_azmi_ishara/2015_07_16_azmi_ishara/output_plots/"
#					  "plasmid_center_acs/", dataset, ".bam", sep = ""
#					 )

# Setup the screen
scr.m = matrix(c(0, 1, 0.95, 1,
				 0, 1, 0, 0.95
				),
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

# Label the origin
text(x = 916, y = 0.5, labels = "ACS", srt = 90, col = "darkgreen", cex = 1.125)
text(x = 945, y = 0.5, labels = "B1", srt = 90, col = "darkgreen", cex = 1.125)
text(x = 980, y = 0.5, labels = "B2", srt = 90, col = "darkgreen", cex = 1.125)
text(x = 1023, y = 0.5, labels = "B3", srt = 90, col = "darkgreen", cex = 1.125)

# Plot the nucleosome track
screen(s_num.v[2])
par(mar = c(4, x_left, 2, x_right), cex = 1.25)

# Set up the plot
plot(0, 0, type = "n",
	 xlim = c(-350, 350), xaxs = "i", 
	 ylim = c(0, 2.5), yaxs = "i",
	 xlab = "Relative distance from ACS of ARS1 (bp)",
	 ylab = "Nuc density"
	)

# Set the color vector
color_vector.v = c("black", "red", "darkgreen", "blue", "maroon", "darkgoldenrod", "darkcyan", "darkorange")

# Iterate through each bam_file
for(b in 1:length(dataset)){

	# Get the bam_file_name
	bam_file_name = paste("/data/collaborator_data/aligned_experiments/", dataset[b], "/", dataset[b], ".bam", sep = "")

	# Get the nucleosome feature density
	nuc.v = get_mnase_feature_density(bam_file_name, feature_chr, x_start, x_end, 125, 175, 20)

	# Get the total number of reads
	read_number = get_total_read_number(bam_file_name, feature_chr)

	# Normalize the nuc.v to 200,000 reads
	# nuc.v = nuc.v * (2.0E5)/read_number

	# Enter into the plot
	lines(-350:350, nuc.v, col = color_vector.v[b])

}

# Enter the legend
legend("top", legend = dataset_name, bty = "n", lwd = 2, col = color_vector.v, ncol = 2)

# Close the screen
close.screen(all.screens = T)

# Close the device
dev.off()
