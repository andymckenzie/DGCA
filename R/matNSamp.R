#' @title Find the number of non-missing values.
#' @description This function calculates the pairwise number of non-missing values in a matrix.
#' @param matA Input data matrix with numeric entries.
#' @param impute Binary value; if true, indicates that imputation was performed previously, and so checking for NAs is not necessary.
#' @param matB Optional input data matrix with which the comparison with matA will be made.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @return A number of samples (nsamp) matrix.
#' data(darmanis); darmanis_subset = as.matrix(darmanis[1:30, ])
#' nsamp_res = matNSamp(darmanis_subset)
#' darmanis_subset[1, 1] = NA
#' nsamp_res_na = matNSamp(darmanis_subset)
#' @export
matNSamp <- function(matA, impute = FALSE, matB = NULL, secondMat = FALSE){

	if(!secondMat){
		#if there are no NA values in the matrix, then nsamp is directly calculated
		if(impute | sum(is.na(matA)) == 0){
			nsamp = matrix(nrow(matA), nrow = ncol(matA), ncol = ncol(matA))
		} else {
			#else count the number of NA values in each comparison
		  matA[!is.na(matA)] = 1
		  matA[is.na(matA)] = 0
		  nsamp = t(matA) %*% matA
		}
	}

	if(secondMat == TRUE){
		#if there are no NA values in the matrix, then nsamp is directly calculated
		if(impute | sum(is.na(matA)) == 0 & sum(is.na(matB)) == 0){
			nsamp = matrix(nrow(matA), nrow = ncol(matA), ncol = ncol(matB))
		} else {
			#else count the number of NA values in each comparison
		  matA[!is.na(matA)] = 1
		  matA[is.na(matA)] = 0
		  matB[!is.na(matB)] = 1
		  matB[is.na(matB)] = 0
		  nsamp = t(matA) %*% matB
		}
	}

  return(nsamp)

}
