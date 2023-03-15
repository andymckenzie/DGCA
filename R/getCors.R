#' @title Compute matrices necessary for differential correlation calculation.
#' @description As a first step in the standard DGCA workflow, this function reads in a design matrix and an input matrix, splits the input matrix by the named groups defined in the design matrix, and outputs a list of matrices (correlation matrix, correlation significance matrix, and number of samples matrix) to be used by downstream parts of the analysis.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param inputMatB Optional, secondary input matrix that allows you to calculate correlation and differential correlation for the rows between inputMat and imputMatB.
#' @param impute A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman". Default = "pearson".
#' @return A corMats S4 class object, containing a list of matrices from each group, the design matrix, and a character vector of options.
#' @examples
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' cors_res = getCors(inputMat = darmanis_subset, design = design_mat)
#' @export
getCors <- function(inputMat, design, inputMatB = NULL, impute = FALSE, corrType = "pearson"){

	##############################
	#set SAF to FALSE while restoring to default when the function is finished
	SAF = getOption("stringsAsFactors", FALSE)
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

	##################################
	#check inputs
	if(!corrType %in% c("pearson", "spearman")) stop("corrType should be one of \"pearson\" or \"spearman\".\n")

	#check the design matrix
	if(!mode(design) == "numeric") stop("Design matrix must be numeric.\n")

	#design matrix should be made up of only 0's and 1's
	if(!all(design %in% c(0, 1))) stop("Design matrix must be made up of 0's and 1's.\n")

	#check the correspondence between the input and design matrices
	if(nrow(design) != ncol(inputMat)) stop("The number of rows in the design matrix must be equal to the number of columns in the input matrix.")

	if(!is.null(inputMatB)){
		#check the correspondence between the input and design matrices
		if(nrow(design) != ncol(inputMatB)) stop("The number of rows in the design matrix must be equal to the number of columns in the input matrix.")

	}

	###################################
	# if necessary, perform imputation on the whole matrix (prior to subsetting)
	if(impute == TRUE){

		if (!requireNamespace("impute", quietly = TRUE)) {
			stop("The R package impute is needed for the impute knn function to work. Please install it.",
				call. = FALSE)
		}

		if(any(is.na(inputMat))){

			imputed = impute::impute.knn(as.matrix(inputMat))
			inputMat = imputed$data

			} else {message("No NA values in input, so ignoring impute = TRUE.\n")}

			if(!is.null(inputMatB)){
				if(any(is.na(inputMatB))){

					imputedB = impute::impute.knn(as.matrix(inputMatB))
					inputMatB = imputedB$data

					} else {message("No NA values in secondary input, so ignoring impute = TRUE.\n")}

			}

	}

	####################################
	#if there is only one input matrix to be compared with itself, entrywise
	if(is.null(inputMatB)){

		####################################
		#use the design matrix to split the input matrix a list of sub-matrices
		designRes = getGroupsFromDesign(inputMat, design)
		groupList = designRes[[1]]

		#create list of lists to store the results
		groupMatLists = as.list(rep(NA, length(groupList)))
		names(groupMatLists) = designRes[[2]]

		###################################
		#calculate the correlation-related matrices for each of the sub-matrices

		for(i in 1:length(groupList)){

			#for extensibility, the first two fxns below can accept multiple matrices
			corr = matCorr(t(groupList[[i]]), corrType = corrType)
			nsamp = matNSamp(t(groupList[[i]]), impute = impute)
			pval = matCorSig(corr, nsamp)

			groupMatLists[[i]] = list(corrs = corr, pvals = pval, nsamps = nsamp)

		}

	}

	####################################
	#if there are two input matrices to be compared
	if(!is.null(inputMatB)){

		####################################
		#use the design matrix to split the input matrix a list of sub-matrices
		designRes = getGroupsFromDesign(inputMat = inputMat, design = design,
			inputMatB = inputMatB, secondMat = TRUE)
		groupListA = designRes[[1]]
		groupListB = designRes[[2]]

		#create list of lists to store the results
		groupMatLists = as.list(rep(NA, length(groupListA)))
		names(groupMatLists) = designRes[[3]]

		###################################
		#calculate the correlation-related matrices for each of the sub-matrices

		for(i in 1:length(groupListA)){

			#for extensibility, the first two fxns below can accept multiple matrices
			corr = matCorr(matA = t(groupListA[[i]]), corrType = corrType,
				matB = t(groupListB[[i]]), secondMat = TRUE)

			nsamp = matNSamp(matA = t(groupListA[[i]]), impute = impute,
				matB = t(groupListB[[i]]), secondMat = TRUE)

			pval = matCorSig(corr, nsamp, secondMat = TRUE)

			groupMatLists[[i]] = list(corrs = corr, pvals = pval, nsamps = nsamp)

		}

	}

	#################################
	#create the class object that will be outputted and used by subsequent analysis steps
	#list of lists of the matrices, named by condition
	#options, e.g. corrType, +/- impute, and design matrix (for differential expression wrapper...)

	options = c(corrType)
	names(options) = c("corrType")

	corMatsObj = new("corMats", corMatList = groupMatLists,
		design = design, options = options)

	return(corMatsObj)

}
