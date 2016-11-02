#' @title Calculate pairwise differential correlations
#' @description Find the differential correlation between two conditions.
#' @param corMatsObj A class object containing a named list of lists of matrices, one list per condition, with each list containing a correlation matrix, a correlation significance p-values matrix, and a "number of samples used to calculate the correlation" matrix.
#' @param compare Vector of two character strings, each corresponding to one name in the list of correlation matrices that should be compared.
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none"). This could be used if you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @return A dcPair class object, containing the difference in z scores for each comparison, the p-values of that differences, and the original correlation matrices and significances for subsequent classification steps.
#' @export
pairwiseDCor <- function(corMatsObj, compare, corr_cutoff = 0.99, corrType = "pearson",
	secondMat = FALSE, signType = "none"){

	SAF = getOption("stringsAsFactors")
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

	inputLists = slot(corMatsObj, "corMatList")

	if(!mode(compare) == "character") stop("Compare argument is not a character vector.")

	if(!all(compare %in% names(inputLists))) stop("Comparison names are not in the input lists.\n")

	matA = inputLists[[compare[1]]]
	matB = inputLists[[compare[2]]]

	AB_res = dCorMats(matA[["corrs"]], matA[["nsamps"]],
		matB[["corrs"]], matB[["nsamps"]], corr_cutoff = corr_cutoff,
		corrType = corrType, secondMat = secondMat, signType = signType)

	#may need these dimnames in the Zdiff matrix for dCorAvg
	colnames(AB_res$diffs) = colnames(matA$corrs)
	rownames(AB_res$diffs) = rownames(matA$corrs)

	dcPair = list(matA$corrs, matA$pvals, matB$corrs, matB$pvals,
		AB_res$diffs, AB_res$pvals)
	names(dcPair) = c("corA", "corPvalA", "corB", "corPvalB",
		"ZDiff", "PValDiff")

	dcPair_res = new("dcPair",
		corA = matA$corrs,
		corPvalA = matA$pvals,
		corB = matB$corrs,
		corPvalB = matB$pvals,
		ZDiff = AB_res$diffs,
		PValDiff = AB_res$pvals)

	return(dcPair_res)

}
