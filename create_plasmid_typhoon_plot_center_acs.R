# Clear the workspace
rm(list = ls())
graphics.off()

# Load the libraries
# library("TyphoonPlot")
# library("nucleR")

# Set the dataset
dataset_id.v = paste("CD", 3:17, sep = "")
names(dataset_id.v) = c("ISW1a", "ISW1b", "ISW2", "INO80", "SWI/SNF", "RSC", "CHD1", "Nap1",
                        "ORC + ISW1a", "ORC + SWI/SNF", "ORC + RSC", "ORC + ISW1a",
                        "ORC + MCM + ISW1a", "ORC + MCM + SWI/SNF", "ORC + MCM + RSC"
                        )


# Set the plot_dir
plot_dir = paste("/home/jab112/alchemy_data2_jab112/2015_bell_lab_ishara/output_plots/",
				 "plasmid_center_acs/plots_500bp_win", sep = ""
				)

# Set the BAM File dir
bam_file_dir = "/home/jab112/collaborator_data/aligned_experiments"

# Set the axes
feature_chr = ""
feature_pos = 916
win = 500

# Set the x-coordinate positions
x_left = feature_pos - win
x_right = feature_pos + win
x_step = win / 2






# Iterate through each dataset
for(i in 1:length(dataset_id.v)){
# for(i in 1){

  cat("Processing ", dataset_id.v[i], " (", names(dataset_id.v)[i], ")...\n", sep = "")

	# Get the BAM file
	bam_file_name = paste(bam_file_dir, "/", dataset_id.v[i], "/", dataset_id.v[i], ".bam", sep = "")

	# Set the feature_chr
	feature_chr = names(scanBamHeader(bam_file_name)[[1]]$targets)

	# Get the nucleosome feature density
	nuc.v = GetMNaseFeatureDensity(bam_file_name, feature_chr, x_left, x_right, 125, 175, 20)

	# Get the nucleosome peaks
	nuc_peaks.df = GetDensityPeaks(cov.v = nuc.v, peak_width = 75, isPeakMax = FALSE)

	# Subset on the 2 nuc_peaks that are closest to the feature position
  upstream_nuc_peaks.df = nuc_peaks.df[which(nuc_peaks.df$pos < feature_pos),]
  downstream_nuc_peaks.df = nuc_peaks.df[which(nuc_peaks.df$pos > feature_pos),]

  # Update the nuc_peaks.df
  nuc_peaks.df = rbind(upstream_nuc_peaks.df[nrow(upstream_nuc_peaks.df),],
                       downstream_nuc_peaks.df[1,]
                      )




	# Set the plot
	#	If want transparent, add 'bg = "transparent" '
	png(file = paste(plot_dir, "/", dataset_id.v[i], ".png", sep = ""),
		width = 10, height = 7.5, units = "in", res = 300) # bg = "transparent")

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
	mar_left = 4.5
	mar_right = 4.5

	# Set up the schematic plot
	screen(s_num.v[1])

		# Set the plot margin parameters
		par(mar = c(0, mar_left, 0, mar_right))

		# Make the gene schematic
		#	If want transparent, add 'bg = "transparent" '
		plot(0, 0, type = "n", bty = "n", # bg = "transparent",
			 xlim = c(x_left, x_right), xaxs = "i", xaxt = "n",
			 ylim = c(0, 1), yaxs = "i", yaxt = "n",
			 ann = F
			)

		# Enter the plot name
		text(x = x_left + 10, y = 0.5, adj = 0, labels = names(dataset_id.v)[i], cex = 1.5)

		# Label the ACS sequence elements
		text(x = 916, y = 0.5, labels = "ACS", srt = 90, col = "darkgreen", cex = 1.125)
		text(x = 945, y = 0.5, labels = "B1", srt = 90, col = "darkgreen", cex = 1.125)
		text(x = 980, y = 0.5, labels = "B2", srt = 90, col = "darkgreen", cex = 1.125)
		text(x = 1023, y = 0.5, labels = "B3", srt = 90, col = "darkgreen", cex = 1.125)







	# Set up the chromatin schematic
	screen(s_num.v[2])

		# Set the plot margin parameters
		par(mar = c(0, mar_left, 0, mar_right))

		# Set up a blank plot schematic
		SetChromatinSchematic(x_left, x_right)

		# Plot the nucleosome
		PlotNucleosome(nuc_peaks.df, y_max = max(nuc_peaks.df$sig), 0.5, 0.2)




	# Plot the nucleosome track
	screen(s_num.v[3])

		# Set the plot margin parameters
		par(mar = c(2, mar_left, 2, mar_right), cex = 1.25)

		# Set up the density nucleosome plot
		plot(x_left:x_right, nuc.v, type = "l", col = "red",
			 xlim = c(x_left, x_right), xaxs = "i", xaxt = "n",
			 ylim = range(nuc.v),
			 ylab = "Nuc density", xlab = NA
			)

		# Add in the x-axis label
		axis(1, at = seq(x_left + x_step/2, x_right - x_step/2, x_step), labels = F, tcl = 0.625 * par()$tcl)
		axis(1, at = seq(x_left, x_right, x_step), labels = seq(-win, win, x_step))

		# Indicate the nucleosome positions on the plot
		abline(v = nuc_peaks.df$pos, lty = 2)
		axis(3, at = nuc_peaks.df$pos, labels = nuc_peaks.df$pos - feature_pos)




if(0){
	# Get the typhoon plot
	screen(s_num.v[4])

		# Set the plot margin parameters
		par(mar = c(5, mar_left, 2, mar_right), cex = 1.25)

		# Get the mat
		mat.m = GetTyphoonPlotMat(bam_file_name, feature_chr, x_left, x_right)

		# Make the typhoon plot
		DensDotPlot(mat.m, z_max = quantile(as.vector(mat.m), probs = 0.99),
					x_label = "Relative distance from ACS of ARS1 (bp)",
					y_label = "Fragment length (bp)", plot_title = "",
					x_axt = "n", y_axt = "n"
				   )

		# Set the Axes
		MakeTyphoonPlotAxes(-win, win, addXAxis = F, addYAxis = T)
		axis(1, at = seq(x_left, x_right, x_step), labels = seq(-win, win, x_step))
		axis(1, at = seq(x_left + x_step/2, x_right - x_step/2, x_step), labels = F, tcl = 0.625 * par()$tcl)

}

  # Get the fragment length distribution
	screen(s_num.v[4])

		# Set the plot margin parameters
		par(mar = c(5, mar_left, 2, mar_right), cex = 1.25)

		# Obtain the paired-end reads over the region
		reads.gr = ConvertPairedReadBAMToGR(bam_file_name, feature_chr, x_left, x_right)

		# Create a GenomicRanges object of the nucleosome positions
		nucleosome_pos.gr = GRanges(seqnames = feature_chr,
		                            ranges = IRanges(start = nuc_peaks.df$pos, width = 1)
		                           )

		# Find the overlaps
		overlaps.l = as.list(findOverlaps(nucleosome_pos.gr, reads.gr))

		# Get the read lenght distributions for each position
		frag_length_dist.l = lapply(overlaps.l, function(x){return(width(reads.gr[x]))})

		# Get the inter-quartile range for each dataset
		BoxplotIQR = function(x.v){

		  # Get the 25th and 75th percentiles
		  boxplot_quantiles.v = quantile(x.v, probs = c(0.25, 0.75))

		  # Get the interquartile range
		  iqr = diff(boxplot_quantiles.v)

		  # Get the IQR distance from each quantile
		  min_y = boxplot_quantiles.v[1] - 1.5 * iqr
		  max_y = boxplot_quantiles.v[2] + 1.5 * iqr

		  # Return the extent of the boxplot range
		  return(c(min_y, max_y))

		}

		# Get the IQR for the boxplot
		boxplot_range.l = lapply(frag_length_dist.l, BoxplotIQR)

		# Set up the plot
		plot(0, 0, type = "n",
		     xlim = c(x_left, x_right), xaxs = "i", xaxt = "n",
		     ylim = range(unlist(boxplot_range.l)), yaxt = "n",
		     xlab = "Relative distance from ACS of ARS1 (bp)",
		     ylab = "Fragment length (bp)"
		)

		# Add in the boxplot
		boxplot(frag_length_dist.l, outline = F, at = nuc_peaks.df$pos,
		        pars = list(boxwex = 75), names = NA,
		        xaxt = "n", yaxt = "n",
		        add = T
		)

		# pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5)

		# Set the axes
		axis(1, at = seq(x_left, x_right, x_step), labels = seq(-win, win, x_step))
		axis(1, at = seq(x_left + x_step/2, x_right - x_step/2, x_step), labels = F, tcl = 0.625 * par()$tcl)

		axis(2, at = seq(0, 250, 50), labels = T)
		axis(2, at = seq(25, 225, 50), labels = F, tcl = 0.625 * par()$tcl)

	# Close the screens
	close.screen(all.screens = T)

	# Close the device
	dev.off()

}
