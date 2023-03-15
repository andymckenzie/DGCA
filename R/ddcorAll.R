#' @title Calls the DGCA pairwise pipeline.
#' @description Runs the full discovery of differential correlation (ddcor) section for comparing pairwise correlations across conditions in the Differential Gene Correlation Analysis (DGCA) package.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param inputMatB Optional, secondary input matrix that allows you to calculate correlation and differential correlation for the rows between inputMat and imputMatB. Default = NULL.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param splitSet Optional character vector that splits the first matrix into two matrices and calculates differential correlation across these matrices. Common use case is when you want the differential correlation of a small set of identifiers (e.g., one), compared with all of the other identifiers in the matrix in each condition. Cannot be used when a second matrix is inputted -- setting both of arguments to non-NULL values will result in an error.
#' @param impute A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman". Default = "pearson".
#' @param nPairs Either a number, specifying the number of top differentially correlated identifier pairs to display in the resulting table, or a the string "all" specifying that all of the pairs should be returned. If splitSet is specified, this is reset to the number of non-splitSet identifiers in the input matrix, and therefore will not be evaluated.
#' @param sortBy Character string specifying the way by which you'd like to sort the resulting table.
#' @param adjust Allows for resulting p-values to be corrected for multiple hypothesis tests, optional. Some non-default choices require the "fdrtool" package or the "qvalue". Default = "perm", which use permutation samples. Other options include "perm" which means that no p-value adjustment is performed, methods in ?p.adjust (i.e., "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), and methods in ?fdrtool (i.e., "fndr", "pct0", "locfdr").
#' @param nPerms Number of permutations to generate. If NULL, permutation testing will not be performed. Default = "10".
#' @param plotFdr Allows for plotting of fdrtool p-value adjustment result OR empirical FDR q-value adjustment technique, if either of these are chosen. Requires fdrtool package OR qvalue package. Default = FALSE.
#' @param classify Binary value specifying whether the correlation values in each condition and differential correlation scores should be used to classifying the resulting identifiers into groups. Default = TRUE
#' @param sigThresh If classify = TRUE, this numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.
#' @param corSigThresh If classify = TRUE, this numeric value specifies the p-value threshold at which a correlation p-value is deemed significant. Default = 0.05.
#' @param heatmapPlot Option indicating whether a heatmap of the differential correlations between the two conditions should be plotted. Default = TRUE.
#' @param customize_heatmap Option to remove some default options in the heatmap plot, to allow users to add custom options.
#' @param color_palette Color palette for plotting the heatmap. If not specified, the heatmap defaults to a red-green color-blind palette with bluish green indicating negative correlations and vermillion indicating positive correlations. Default = NULL
#' @param heatmapClassic Option to make the heatmap more granular (e.g., not showing the individual gene symbols) and more of a "classic" type of heatmap. Overrides most other heatmap options.
#' @param corPower The power to raise the correlations to before plotting the classic heatmap. Larger correlation powers emphasize larger correlation values relatively more compared to smaller correlation values.
#' @param verbose Option indicating whether the program should give more frequent updates about its operations. Default = FALSE.
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Note that this does NOT affect the underlying correlation values, but does affect the z-score difference of correlation calculation in the dcTopPairs table. Default = 0.99
#' @param getDCorAvg Logical, specifying whether the average difference in correlation between groups should be calculated. Default = FALSE
#' @param dCorAvgType Character vector specifying the type of average differential correlation calculation that should be performed. Only evaluated if dCorAge is TRUE. Types = c("gene_average", "total_average", "both"). gene_average calculates whether each genes' differential correlation with all others is more than expected via permutation samples (and empirical FDR adjustment, in the case of > 1 gene), while total_average calculates whether the total average differential correlation is higher than expected via permutation samples. "both" performs both of these. If splitSet is specified, then only genes in the splitSet have their average gene differential correlation calculated if gene_average is chosen.
#' @param dCorAvgMethod Character vector specifying the method for calculating the "average" differential correlation calculation that should be used. Options = "median", "mean".
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none") prior to calculating differential correlation. This could be used if, e.g., you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @param oneSidedPVal If the dCorAvgType test is total_average, this option specifies whether a one-sided p-value should be reported, as opposed to a two-sided p-value. That is, if the average difference of z-scores is greater than zero, test whether the permutation average difference of z-scores are less than that average to get the p-value, and vice versa for the case that the average difference of z-scores is less than 0. Otherwise, test whether the absolute value of the average difference in z-scores is greater than the absolute values of the permutation average difference in z-scores. Default = FALSE.
#' @param ... Additional plotting arguments if heatmapPlot = TRUE.
#' @return Typically, the returned object is a data frame of the table of differential correlations between conditions. In the case that dCorAvg is calculated, the returned object is instead a list containing that table as well as the object summarizing the difference in average correlation for the specified portion of the data set.
#' @examples
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' ddcor_res = ddcorAll(inputMat = darmanis_subset, design = design_mat,
#' 	compare = c("oligodendrocyte", "neuron"))
#' @export
ddcorAll <- function(inputMat, design, compare, inputMatB = NULL, splitSet = NULL,
	impute = FALSE, corrType = "pearson", nPairs = "all", sortBy = "zScoreDiff",
	adjust = "perm", nPerms = 10, classify = TRUE, sigThresh = 1,
	corSigThresh = 0.05, heatmapPlot = FALSE, color_palette = NULL, verbose = FALSE, plotFdr = FALSE,
	corr_cutoff = 0.99, signType = "none", getDCorAvg = FALSE, dCorAvgType = "gene_average",
	dCorAvgMethod = "median", oneSidedPVal = FALSE, customize_heatmap = FALSE,
	heatmapClassic = FALSE, corPower = 2, ...){

	################################
	# check inputs

	if(!is.null(splitSet)){
		if(!mode(splitSet) == "character") stop("splitSet must be character type.\n")
	}

	if(!is.null(splitSet) & !is.null(inputMatB)){
		stop("You cannot input both a character vector to split the first input matrix into two as well as a second matrix of identifiers.")
	}

	if(is.null(rownames(inputMat))){
		stop("inputMat must specify identifiers as rownames.")
	}

	if(adjust == "perm" & nPerms == 0){
		stop("If you choose permutation for p-value adjustment, then you need to generate at least one permutation sample by setting nPerms > 0.")
	}

	if(adjust != "perm" & nPerms > 0 & getDCorAvg == FALSE){
		warning("If you are not choosing permutation for p-value adjustment or calculating the differential correlation average, then you may be wasting time by generating permutation samples. Consider setting nPerms to 0.")
	}

	if(!(is.numeric(nPairs) | nPairs == "all")){
		stop("nPairs must either be numeric or be a character vector \"all\" to specify that all pairs should be returned.")
	}

	if(nPairs == "all" & is.null(inputMatB)){
		nPairs = (nrow(inputMat)^2)/2 - nrow(inputMat)/2
	}

	if(nPairs == "all" & !is.null(inputMatB)){
		nPairs = nrow(inputMat) * nrow(inputMatB)
	}

	if(!dCorAvgMethod %in% c("median", "mean")){
		stop("The differential correlation average method chosen must be one of median or mean.")
	}

	##############################
	#set SAF to FALSE while restoring to default when the function is finished
	SAF = getOption("stringsAsFactors", FALSE)
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

	#############################
	# subset the input matrix as necessary

	if(!is.null(splitSet)){
		splitSetFound = splitSet %in% rownames(inputMat)
		#make sure that at least some of the splitSet entries are found in rownames
		nPairs = (nrow(inputMat) - sum(splitSetFound)) * sum(splitSetFound)
		if(sum(splitSetFound) == 0){
			stop("None of the splitSet identifiers were found in the rownames of the input matrix.")
		} else if(sum(splitSetFound) > 0){
			splitSetRows = rownames(inputMat) %in% splitSet
			inputMatB = inputMat[splitSetRows, , drop = FALSE]
			#switch the input matrix to not have the splitSet elements
			inputMat = inputMat[!splitSetRows, , drop = FALSE]
			if(verbose){
			message(sum(splitSetRows), " row(s) corresponding to the ", length(splitSet),
				" identifier(s) found in splitSet were found in the input matrix.")
			}
		}
	}

	secondMat = FALSE
	if(!is.null(splitSet) | !is.null(inputMatB)){
		secondMat = TRUE
	}

	##############################
	# calculate correlation and differential correlation

	ddcor_res = getDCors(inputMat = inputMat, design = design,
		compare = compare, inputMatB = inputMatB,
		impute = impute, corrType = corrType, corr_cutoff = corr_cutoff,
		signType = signType)

	if(nPerms > 0){
		ddcor_perm = getDCorPerm(inputMat = inputMat, design = design,
			compare = compare, inputMatB = inputMatB,
			impute = impute, corrType = corrType, nPerms = nPerms,
			corr_cutoff = corr_cutoff, signType = signType)
	}

	##############################
	# extract the differential correlation table as requested

	if(adjust != "perm"){
		ddcor_table = dcTopPairs(dcObject = ddcor_res, nPairs = nPairs,
			adjust = adjust, plotFdr = plotFdr, classify = classify, compare = compare,
			sigThresh = sigThresh, corSigThresh = corSigThresh, verbose = verbose,
			secondMat = secondMat)
	}

	if(adjust == "perm"){
		ddcor_table = dcTopPairs(dcObject = ddcor_res, nPairs = nPairs,
			adjust = adjust, plotFdr = plotFdr, classify = classify, compare = compare,
			sigThresh = sigThresh, corSigThresh = corSigThresh,
			zScorePerm = ddcor_perm, verbose = verbose, secondMat = secondMat)
	}

	#############################
	# return the table

	ddcor_table = ddcor_table[order(-abs(ddcor_table[ , sortBy])), ]

	###############################
	# plot the differential correlation object using a heatmap

	if(heatmapPlot){
		ddplot_plot = ddplot(dcObject = ddcor_res, color_palette = color_palette,
			customize_heatmap = customize_heatmap, heatmapClassic = heatmapClassic,
			corPower = corPower, ...)
	}

	############################
	# if requested, find the average of correlation differences between conditions

	if(getDCorAvg){
		ZDiffs = slot(ddcor_res, "ZDiff")
		if(!dCorAvgType == "both"){
			avg_dcor = dCorAvg(zDiff = ZDiffs, zDiffPerm = ddcor_perm,
				dCorAvgType = dCorAvgType, secondMat = secondMat, dCorAvgMethod = dCorAvgMethod)
			ddcor_res = list(dcPair = ddcor_table, avg_dcor = avg_dcor)
		} else {
			gene_avg_dcor = dCorAvg(zDiff = ZDiffs, zDiffPerm = ddcor_perm,
				dCorAvgType = "gene_average", secondMat = secondMat, dCorAvgMethod = dCorAvgMethod)
			total_avg_dcor = dCorAvg(zDiff = ZDiffs, zDiffPerm = ddcor_perm,
				dCorAvgType = "total_average", secondMat = secondMat, dCorAvgMethod = dCorAvgMethod)
			ddcor_res = list(dcPair = ddcor_table, gene_avg_dcor = gene_avg_dcor,
				total_avg_dcor = total_avg_dcor)
		}
	}

	if(!getDCorAvg){
		ddcor_res = ddcor_table
	}

	#############################
	# return results, in different formats depending on what options were chosen

	return(ddcor_res)

}
