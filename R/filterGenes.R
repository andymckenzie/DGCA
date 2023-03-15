#' @title Filter rows out of a matrix.
#' @description Filter out rows in an input matrix that are not above a certain percentile with respect to a central tendency and/or dispersion measure. To be used, e.g, prior to differential correlation testing with the function ddcorall.
#' @param inputMat The matrix (or data.frame) of all numeric values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers that you are interested in protecting from the filter, if any.
#' @param filterTypes Vector containing up to two character strings, specifying the methods that should be used for filtering genes. Options include "central" and "dispersion" for filtering based on the measures of central tendency and dispersion, respectively. To use both, set this to c("central", "dispersion").
#' @param keepRows Optional character vector, specifying rownames (i.e., symbols) that should not be filtered out of the matrix even if they are found to be below the quantile specified for either the central tendency or dispersion, as applicable.
#' @param filterCentralType Method to be used for filtering for the central tendency of the input matrix. Options = "mean" (for arithmetic mean) and "median".
#' @param filterDispersionType Method to be used for filtering for the dispersion of the input matrix. Options = "dispersion_index", "cv" (for coefficient of variation), and "variance".
#' @param filterCentralPercentile If central tendency filtering is used, the quantile of the central tendency below which rows will be filtered out.
#' @param filterDispersionPercentile If dispersion filtering is used, the quantile of the dispersion measure below which rows will be filtered out.
#' @param sequential If both central tendency and dispersion measures and used for filtering the input matrix, then sequential is a logical flag indicating whether the central tendency filtering steps should be performed prior to the dispersion filtering step (and quantile cutoff specification; if sequential = TRUE), or independently (if sequential = FALSE).
#' @param allGroups Logical for whether genes need to pass the filter in all of the groups specified in the design matrix.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @return A filtered matrix.
#' @examples
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' filtered_mat = filterGenes(inputMat = darmanis_subset, filterTypes = "central")
#' filtered_mat_both = filterGenes(inputMat = darmanis_subset,
#'  filterTypes = c("central", "dispersion"), filterCentralType = "mean",
#'  filterDispersionPercentile = 0.1)
#' filtered_mat_all_groups = filterGenes(inputMat = darmanis_subset,
#'  design = design_mat, filterTypes = "dispersion", allGroups = TRUE)
#' @export
filterGenes <- function(inputMat, filterTypes = "central",
  keepRows = NULL, filterCentralType = "median",
  filterDispersionType = "dispersion_index", filterCentralPercentile = 0.25,
  filterDispersionPercentile = 0.25, sequential = FALSE, allGroups = FALSE,
  design = NULL){

  ##############################
	#set SAF to FALSE while restoring to default when the function is finished
	SAF = getOption("stringsAsFactors", FALSE)
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("The library matrixStats is needed for the filtering function to work. Please install it.",
      call. = FALSE)
  }

  if(!allGroups){

    rows_above_central = rep(TRUE, nrow(inputMat))
    if("central" %in% filterTypes){
      if(filterCentralType == "mean"){
        row_central = rowMeans(data.matrix(inputMat))
      }
      if(filterCentralType == "median"){
        row_central = matrixStats::rowMedians(data.matrix(inputMat))
      }
      row_central_cutoff = quantile(row_central, filterCentralPercentile)
      rows_above_central = row_central > row_central_cutoff
    }

    if(!is.null(keepRows)){
      rows_above_central = rows_above_central | rownames(inputMat) %in% keepRows
    }

    if(sequential){
      inputMat = inputMat[rows_above_central, ]
    }

    rows_above_dispersion = rep(TRUE, nrow(inputMat))
    if("dispersion" %in% filterTypes){
      row_vars = matrixStats::rowVars(data.matrix(inputMat))
      if(filterDispersionType == "dispersion_index"){
        row_means = rowMeans(data.matrix(inputMat))
        row_dispersion = row_vars/abs(row_means)
      }
      if(filterDispersionType == "cv"){
        row_means = rowMeans(data.matrix(inputMat))
        row_dispersion = sqrt(row_vars)/abs(row_means)
      }
      if(filterDispersionType == "variance"){
        row_dispersion = row_vars
      }
      row_dispersion_cutoff = quantile(row_dispersion, filterDispersionPercentile)
      rows_above_dispersion = row_dispersion > row_dispersion_cutoff
    }

    if(!is.null(keepRows)){
      rows_above_dispersion = rows_above_dispersion | rownames(inputMat) %in% keepRows
    }

    if(sequential){
      rows_to_keep = rows_above_dispersion
    } else {
      rows_to_keep = rows_above_central & rows_above_dispersion
    }

  } else {

    if(is.null(design)) stop("If allGroups == TRUE, then you must input a design matrix as well.")

    designRes = getGroupsFromDesign(inputMat, design)
    groupList = designRes[[1]]

    for(i in 1:length(groupList)){

      rows_above_central = rep(TRUE, nrow(groupList[[i]]))
      if("central" %in% filterTypes){
        if(filterCentralType == "mean"){
          row_central = rowMeans(data.matrix(groupList[[i]]))
        }
        if(filterCentralType == "median"){
          row_central = matrixStats::rowMedians(data.matrix(groupList[[i]]))
        }
        row_central_cutoff = quantile(row_central, filterCentralPercentile)
        rows_above_central = row_central > row_central_cutoff
      }

      if(!is.null(keepRows)){
        rows_above_central = rows_above_central | rownames(groupList[[i]]) %in% keepRows
      }

      if(sequential){
        groupList[[i]] = groupList[[i]][rows_above_central, ]
      }

      rows_above_dispersion = rep(TRUE, nrow(groupList[[i]]))
      if("dispersion" %in% filterTypes){
        row_vars = matrixStats::rowVars(data.matrix(groupList[[i]]))
        if(filterDispersionType == "dispersion_index"){
          row_means = rowMeans(data.matrix(groupList[[i]]))
          row_dispersion = row_vars/abs(row_means)
        }
        if(filterDispersionType == "cv"){
          row_means = rowMeans(data.matrix(groupList[[i]]))
          row_dispersion = sqrt(row_vars)/abs(row_means)
        }
        if(filterDispersionType == "variance"){
          row_dispersion = row_vars
        }
        row_dispersion_cutoff = quantile(row_dispersion, filterDispersionPercentile)
        rows_above_dispersion = row_dispersion > row_dispersion_cutoff
      }

      if(!is.null(keepRows)){
        rows_above_dispersion = rows_above_dispersion | rownames(groupList[[i]]) %in% keepRows
      }

      if(i == 1){
        if(sequential){
          rows_to_keep = rows_above_dispersion
        } else {
          rows_to_keep = rows_above_central & rows_above_dispersion
        }
      } else {
        if(sequential){
          rows_to_keep = rows_above_dispersion & rows_to_keep
        } else {
          rows_to_keep = rows_above_central & rows_above_dispersion & rows_to_keep
        }
      }

    }

  }

  filteredMat = inputMat[rows_to_keep, ]

  return(filteredMat)

}
