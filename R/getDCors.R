#' @title Get groupwise correlations and pairwise differential correlations.
#' @description Takes input and methods to perform getCor as well as group-specific pairwiseDCor.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see ?model.matrix for more information.
#' @param inputMatB Optional, secondary input matrix that allows you to calculate correlation and differential correlation for the rows between inputMat and imputMatB.
#' @param compare Vector of two character strings, each corresponding to one name in the list of correlation matrices that should be compared.
#' @param impute A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes.
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none"). This could be used if you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @return A dcPair class object, containing the difference in z-scores for each comparison, the p-values of that differences, and the original correlation matrices and significances for subsequent classification steps.
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' dcors_res = getDCors(inputMat = darmanis_subset, design = design_mat, compare = c("oligodendrocyte", "neuron"))
#'
#' @export
getDCors <- function(inputMat, design, compare, inputMatB = NULL, impute = FALSE, corrType = "pearson",
	corr_cutoff = 0.99, signType = "none"){

	corMats_res = getCors(inputMat = inputMat, design = design, inputMatB = inputMatB, corrType = corrType,
		impute = impute)
	if(is.null(inputMatB)){
		dcPairs_res = pairwiseDCor(corMats_res, compare = compare, corr_cutoff = corr_cutoff,
			corrType = corrType, signType = signType)
	} else {
		dcPairs_res = pairwiseDCor(corMats_res, compare = compare, corr_cutoff = corr_cutoff,
			corrType = corrType, secondMat = TRUE, signType = signType)
	}
	return(dcPairs_res)

}
