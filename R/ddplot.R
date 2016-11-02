#' @title Create a heatmap showing the correlations in two conditions.
#' @description This function orders the differences in correlations between conditions by the median strength of correlation differences for each gene and plots a heatmap of the correlations in each condition (lower = condition A, upper = condition B) using the heatmap.2 function from the gplots package.
#' @param dcObject A differential correlation object from which correlation and differential correlation matrices will be extracted. Optional; can also input the correlation matrices and differential correlation matrix individually.
#' @param corMatA Optional, correlation matrix from condition A. Will be plotted in the lower left triangle.
#' @param corMatB Optional, correlation matrix from condition B. Will be plotted in the upper right triangle.
#' @param zDiff Optional, difference measure of correlations between conditions A and B.
#' @param color_palette Colors for plotting the heatmap. If not specified, defaults to a color-blind palette where blue corresponds to a negative correlation and orange/red corresponds to a positive one.
#' @param customize_heatmap Option to remove some default options in the heatmap plot, to allow users to add custom options.
#' @param heatmapClassic Option to make the heatmap more granular (e.g., not showing the individual gene symbols) and more of a "classic" type of heatmap. Overrides most other heatmap options.
#' @param corPower The power to raise the correlations to before plotting the classic heatmap. Larger correlation powers emphasize larger correlation values relatively more compared to smaller correlation values.
#' @param flip Switch the ordering of z-differences to be inverse. Default = TRUE
#' @param ... Additional plotting arguments to the heatmap.2 function.
#' @return The sorted difference in z-score matrix in both conditions, which you can use to create your own plot if you'd prefer.
#' @export
ddplot <- function(dcObject = NULL, corMatA = NULL, corMatB = NULL,
	zDiff = NULL, flip = TRUE, color_palette = NULL,
	customize_heatmap = FALSE, heatmapClassic = FALSE, corPower = 2, ...){

	if (!requireNamespace("gplots", quietly = TRUE)) {
		stop("gplots is needed for the heatmap.2 function to work. Please install it.",
			call. = FALSE)
	}

	if(!is.null(dcObject) & !is.null(corMatA)){
		stop("Must input either a dcObject *or* the individual correlation matrices and the correlation differences z-score matrix.")
	}

	if(!is.null(dcObject)){
		corMatA = slot(dcObject, "corA")
		corMatB = slot(dcObject, "corB")
		zDiff = slot(dcObject, "ZDiff")
	}

  zDiff[lower.tri(zDiff)] = t(zDiff)[lower.tri(t(zDiff))]

	#calculate the median differential correlation for each gene
	zDiff_row_medians = matrixStats::rowMedians(zDiff, na.rm = TRUE)
	names(zDiff_row_medians) = rownames(zDiff)

	#order based on the median of z differences
	if(!flip){
		z_diff_order = names(sort(zDiff_row_medians))
	} else {
		z_diff_order = names(sort(zDiff_row_medians, decreasing = TRUE))
	}

	corMatA_order = corMatA[z_diff_order, z_diff_order]
	corMatB_order = corMatB[z_diff_order, z_diff_order]

	corMat_total = corMatA_order
	corMat_total[upper.tri(corMat_total)] = corMatB_order[upper.tri(corMatB_order)]
	diag(corMat_total) = NA

	if(is.null(color_palette)){
		color_palette = colorRampPalette(c("#0072B2", "#56B4E9", "white", "#E69F00", "#E34234"))(n = 200)
	}

	breaks_hm = seq(-1, 1, 0.01)

	if(!heatmapClassic){
		#also see http://stackoverflow.com/a/34608069/560791
		message("Plotting heatmap using heatmap.2...")
		if(!customize_heatmap){
			gplots::heatmap.2(corMat_total, dendrogram = "none", Rowv = FALSE,
				Colv = FALSE, trace = "none", col = color_palette, na.color = "black",
		    scale = "none", margins = c(6, 2), symkey = FALSE, symbreaks = TRUE,
		    density.info = 'histogram', denscol = "black", breaks = breaks_hm,
		    keysize = 1, srtCol = 45, key.par = list(mar = c(3.5, 0, 3, 0)),
		    lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(2.5, 5), lwid = c(1, 10, 1),
				sepcolor = "black",
				colsep = 0:ncol(corMat_total), rowsep = 0:nrow(corMat_total),
				...)
		} else {
			gplots::heatmap.2(corMat_total, dendrogram = "none", Rowv = FALSE,
				Colv = FALSE, trace = "none", col = color_palette, na.color = "black",
				scale = "none", margins = c(6, 2), symkey = FALSE, symbreaks = TRUE,
				density.info = "histogram", denscol = "black", breaks = breaks_hm,
				srtCol = 45, key.par = list(mar = c(3.5, 0, 3, 0)),
				lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(2.5, 5), lwid = c(1, 10, 1),
				...)
		}
	}

	if(heatmapClassic){

		#used to flip matrix about the x rather than just the diagonal axis (as in transpose)
		rotate <- function(x) t(apply(x, 2, rev))

		diag(corMat_total) = NA
		corMat_total = abs(corMat_total)
		corMat_total = corMat_total^corPower
		stats::heatmap(rotate(corMat_total), col = rev(heat.colors(100)),
			Rowv = NA, Colv = NA, labRow = FALSE, labCol = FALSE)

	}

	return(corMat_total)

}
