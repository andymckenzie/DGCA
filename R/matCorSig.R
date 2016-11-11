#' @title Calculate correlation matrix p-values.
#' @description Calculate two-sided p-values from a pairwise correlations matrix and a corresponding "number of samples" matrix.
#' @param corrs Computed correlation matrix.
#' @param nsamp Computed number of samples used per call in the correlation matrix.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @return A matrix of p-values.
#' @references HMisc R package https://cran.r-project.org/web/packages/Hmisc/index.html
#' @export
matCorSig <- function(corrs, nsamp, secondMat = FALSE){

	if(!secondMat){

		pvals = matrix(NA, nrow = nrow(corrs), ncol = ncol(corrs))
		corrs = corrs[upper.tri(corrs)]
		nsamp = nsamp[upper.tri(nsamp)]
		df = nsamp - 2
		pvals_upper = 2 * (1 - pt(abs(corrs) * sqrt(df) / sqrt(1 - corrs^2), df))
		pvals[upper.tri(pvals)] = pvals_upper
		pvals[abs(corrs) == 1] = 0

	} else {

		df = nsamp - 2
		pvals = matrix(2 * (1 - pt(abs(corrs) * sqrt(df) /
			sqrt(1 - corrs^2), df)), nrow = nrow(corrs))
		pvals[abs(corrs) == 1] = 0

	}

	return(pvals)

}
