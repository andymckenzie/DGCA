#' @title Split input matrix(es) based on the design matrix.
#' @description This function splits the input matrix(es) based on a design matrix, into a named list of subsetted matrices. If the design matrix has no names, this function will create names for the resulting list of matrices.
#' @param inputMat Input (e.g., expression) matrix which will be subsetted.
#' @param inputMatB Optional input (e.g., expression) matrix which will be subsetted in the same way.
#' @param design Standard design matrix, must specify at least two conditions. For more info, see ?model.matrix
#' @param secondMat Logical value indicating whether there is a second input matrix to be subsetted.
#' @return A list whose first element is a list of subsetted matrices and whose second element is a list of group names.
#' @examples
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' groups_from_design = getGroupsFromDesign(inputMat = darmanis_subset, design = design_mat)
#' str(groups_from_design)
#' @export
getGroupsFromDesign <- function(inputMat, design, inputMatB = NULL, secondMat = FALSE){

	nGroups = qr(design)$rank

	if(nGroups <= 1) stop("Design matrix must specify at least 2 conditions.")

	listGroupMats = list()
	listGroupNames = list()

	#if the design matrix is not named, assign it basic names
	if(!is.character(colnames(design))){
		colnames(design) = toupper(letters[1:ncol(design)])
	}

	if(!secondMat){
		for(i in 1:nGroups){
			tmpCols = as.numeric(which(design[,i] == 1))
			listGroupMats[[i]] = inputMat[ , tmpCols, drop = FALSE]
			listGroupNames[[i]] = colnames(design)[i]
		}
		return(list(listGroupMats, listGroupNames))
	}

	if(secondMat){
		listGroupMatsB = list()
		for(i in 1:nGroups){
			tmpCols = as.numeric(which(design[,i] == 1))
			listGroupMats[[i]] = inputMat[ , tmpCols, drop = FALSE]
			listGroupMatsB[[i]] = inputMatB[, tmpCols, drop = FALSE]
			listGroupNames[[i]] = colnames(design)[i]
		}
		return(list(listGroupMats, listGroupMatsB, listGroupNames))
	}

}
