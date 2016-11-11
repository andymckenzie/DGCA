#' @title Calculate a correlation matrix.
#' @description This function takes one or two input matrices and calculates a correlation matrix from it using the speed-optimized correlation function from WGCNA.
#' @param matA Input data matrix with numeric entries.
#' @param corrType The type of correlation to be performed. Either "pearson" or "spearman".
#' @param matB Optional input data matrix with which the comparison with matA will be made.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @param use The "use" method for performing the correlation calculation. See ?cor for more information. Default = "pairwise.complete.obs" (which is one of the speed-optimized versions; see ?WGCNA::cor for more).
#' @return A correlation matrix.
#' data(darmanis); darmanis_subset = darmanis[1:30, ]
#' matcor_res = matCorr(matA = darmanis_subset, corrType = "pearson")
#' @export
matCorr <- function(matA, corrType, use = "pairwise.complete.obs", matB = NULL, secondMat = FALSE){

	if(!secondMat){
		if(corrType %in% "pearson"){
			corrs = WGCNA::cor(matA, use = use)
		}
		if(corrType %in% "spearman"){
			corrs = WGCNA::cor(matA, use = use, method = "spearman")
		}
	}

	if(secondMat){
		if(corrType %in% "pearson"){
			corrs = WGCNA::cor(matA, matB, use = use)
		}
		if(corrType %in% "spearman"){
		  corrs = WGCNA::cor(matA, matB, use = use,
				method = "spearman")
		}
	}

	return(corrs)

}
