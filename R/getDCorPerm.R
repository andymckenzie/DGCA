#' @title Get permuted groupwise correlations and pairwise differential correlations
#' @description Takes input and methods and randomly permutes the data to do getCor as well as group-specific pairwiseDCor.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param inputMatB Optional, secondary input matrix that allows you to calculate correlation and differential correlation for the rows between inputMat and imputMatB.
#' @param compare Vector of two character strings, each corresponding to one name in the list of correlation matrices that should be compared.
#' @param impute A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param nPerms Number of permutations to generate.
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Default = 0.99
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none"). This could be used if you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @return An array of permuted differences in z-scores calculated between conditions, with the third dimension corresponding to the number of permutations performed.
#' @export
getDCorPerm <- function(inputMat, design, compare, inputMatB = NULL, impute = FALSE,
	nPerms = 10, corrType = "pearson", corr_cutoff = 0.99, signType = "none"){

	secondMat = FALSE
	if(!is.null(inputMatB)){
		zPermMat = array(dim = c(nrow(inputMat), nrow(inputMatB), nPerms))
		secondMat = TRUE
	} else {
		zPermMat = array(dim = c(nrow(inputMat), nrow(inputMat), nPerms))
	}

	for(i in 1:nPerms){
		message("Calculating permutation number ", i, ".")
		inputMat_perm = inputMat[ , sample(ncol(inputMat)), drop = FALSE]
		if(secondMat){
			inputMatB_perm = inputMatB[ , sample(ncol(inputMatB)), drop = FALSE]
			corMats_res = getCors(inputMat_perm, design = design,
				inputMatB = inputMatB_perm, corrType = corrType, impute = impute)
		} else {
			corMats_res = getCors(inputMat_perm, design = design,
				corrType = corrType, impute = impute)
		}
		dcPairs_res = pairwiseDCor(corMats_res, compare, corr_cutoff = corr_cutoff,
			secondMat = secondMat, signType = signType)
		zscores = slot(dcPairs_res, "ZDiff")
		zPermMat[ , , i] = zscores
	}

	return(zPermMat)

}
