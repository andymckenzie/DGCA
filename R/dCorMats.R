#' @title Finds differential correlations between matrices. 
#' @description Takes two corresponding correlation and nsamp matrices and returns matrices for the scaled difference in correlation as well as the p-value of that difference.
#' @param matA Correlation matrix with numeric entries.
#' @param nmatA Number of samples (nsamp) matrix with numeric entries, corresponding to the number of samples used for each of the correlations calculated in matA.
#' @param matB Correlation matrix with numeric entries.
#' @param nmatB Number of samples (nsamp) matrix with numeric entries, corresponding to the number of samples used for each of the correlations calculated in matB.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not. If no, then computations will only be performed the upper triangle of the input matrices.
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none"). This could be used if you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @return A list of two differential correlation matrices: one for the difference in z-scores and one for the corresponding p-values.
#' @export
dCorMats <- function(matA, nmatA, matB, nmatB,
	corr_cutoff = 0.99, corrType = "pearson", secondMat = FALSE, signType = "none"){

	if(!identical(dim(matA), dim(nmatA))){
		stop("Matrices are not the same size.")
	}
	if(!identical(dim(matA), dim(matB))){
		stop("Matrices are not the same size.")
	}
	if(!identical(dim(matA), dim(nmatB))){
		stop("Matrices are not the same size.")
	}

	n_col = ncol(matA)
	n_row = nrow(matA)

	matA[matA > corr_cutoff] = corr_cutoff
	matB[matB > corr_cutoff] = corr_cutoff

	if(signType == "positive"){
		matA[matA < 0] = 0
		matB[matB < 0] = 0
	}

	if(signType == "negative"){
		matA[matA > 0] = 0
		matB[matB > 0] = 0
	}

	if(!secondMat){

		matA = matA[upper.tri(matA)]
		nmatA = nmatA[upper.tri(nmatA)]
		matB = matB[upper.tri(matB)]
		nmatB = nmatB[upper.tri(nmatB)]

		dcorr_res = dCorrs(as.vector(matA), as.vector(nmatA),
			as.vector(matB), as.vector(nmatB), corrType = corrType)

		diffs = matrix(NA, ncol = n_col, nrow = n_row)
		diffs[upper.tri(diffs)] = dcorr_res
		rownames(diffs) = rownames(matA)
		colnames(diffs) = colnames(matA)

		pvals_res = 2 * (pnorm(abs(dcorr_res), mean = 0.00, sd = 1.00,
			lower.tail = FALSE, log.p = FALSE))

		pvals = matrix(NA, ncol = n_col, nrow = n_row)
		pvals[upper.tri(pvals)] = pvals_res
		rownames(pvals) = rownames(matA)
		colnames(pvals) = colnames(matA)

	} else {

		dcorr_res = dCorrs(as.vector(matA), as.vector(nmatA),
			as.vector(matB), as.vector(nmatB), corrType = corrType)

		diffs = matrix(dcorr_res, ncol = n_col, nrow = n_row)
		rownames(diffs) = rownames(matA)
		colnames(diffs) = colnames(matA)

		pvals = 2 * (pnorm(abs(dcorr_res), mean = 0.00, sd = 1.00,
			lower.tail = FALSE, log.p = FALSE))

		pvals = matrix(pvals, ncol = n_col, nrow = n_row)
		rownames(pvals) = rownames(matA)
		colnames(pvals) = colnames(matA)

	}

	return(list("diffs" = diffs, "pvals" = pvals))

}
