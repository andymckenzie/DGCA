#' @title Plot gene pair correlations in multiple conditions.
#' @description Takes the original input matrix, a design matrix, and two gene symbols to plot the corelation in the conditions specified.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman". Default = "pearson".
#' @param geneA The first gene symbol.
#' @param geneB The second gene symbol.
#' @param oneRow Coerce all of the conditions to be plotted on the same row (as opposed to wrapping to multiple rows; relevant if there are >3 conditions).
#' @param log Logical, indicating whether the data should be log2-transformed prior to plotting (after adding a small constant of 0.5 to avoid problems with the log transform).
#' @param smooth Whether to perform lm-based smoothing of the trend in each condition and add this to the plot.
#' @param ylab Override the y-axis label to one of your choice.
#' @param xlab Override the x-axis label to one of your choice.
#' @return A ggplot2 object that can be plotted, further modified, and/or saved.
#' @export
plotCors <- function(inputMat, design, compare, corrType = "pearson",
  geneA, geneB, oneRow = FALSE, smooth = TRUE, log = FALSE, ylab = NULL, xlab = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The R package ggplot2 is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors", FALSE)
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  if(!all(compare %in% colnames(design))) stop("Comparison group names must be column names in the design matrix.")

  if(!geneA %in% rownames(inputMat)) stop("geneA was not found in the rownames of the input matrix.")

  if(!geneB %in% rownames(inputMat)) stop("geneB was not found in the rownames of the input matrix.")

  if(!corrType %in% c("pearson", "spearman")) stop("corrType should be one of \"pearson\" or \"spearman\".\n")

  #check the design matrix
  if(!mode(design) == "numeric") stop("Design matrix must be numeric.\n")

  #design matrix should be made up of only 0's and 1's
  if(!all(design %in% c(0, 1))) stop("Design matrix must be made up of 0's and 1's.\n")

  #check the correspondence between the input and design matrices
  if(nrow(design) != ncol(inputMat)) stop("The number of rows in the design matrix must be equal to the number of columns in the input matrix.")

  ####################################
  #use the design matrix to split the input matrix a list of sub-matrices
  designRes = getGroupsFromDesign(inputMat, design)
  groupList = designRes[[1]]

  ###################################
  #find the paired values for the two genes in each of the conditions

  data_mat = matrix(ncol = 3)

  for(i in 1:length(compare)){

    #get the group number to extract
    group = which(designRes[[2]] %in% compare[i])

    tmp = groupList[[group]]
    tmpGeneA = tmp[rownames(tmp) %in% geneA, ]
    tmpGeneB = tmp[rownames(tmp) %in% geneB, ]
    cond = rep(compare[i], ncol(tmp))

    tmp_mat = data.frame(as.numeric(tmpGeneA), as.numeric(tmpGeneB), as.character(cond))
    colnames(tmp_mat) = c("GeneA", "GeneB", "Condition")

    if(i == 1){
      data_mat = tmp_mat
    } else {
      data_mat = rbind(data_mat, tmp_mat)
    }

  }

  data_mat$Condition = factor(data_mat$Condition, levels = unique(data_mat$Condition))

  if(log){
    data_mat$GeneA = log(data_mat$GeneA + 0.5, 2)
    data_mat$GeneB = log(data_mat$GeneB + 0.5, 2)
  }

  dc_plot = ggplot2::ggplot(data_mat, ggplot2::aes(x = GeneA, y = GeneB, color = Condition)) +
    ggplot2::geom_point() + ggplot2::theme_bw()

  if(oneRow){
    dc_plot = dc_plot + ggplot2::facet_wrap(~ Condition, nrow = 1)
  } else {
    dc_plot = dc_plot + ggplot2::facet_wrap(~ Condition)
  }

  if(smooth){
    dc_plot = dc_plot + ggplot2::geom_smooth(method = 'lm', formula = y~x)
  }

  if(!is.null(ylab)){
    dc_plot = dc_plot + ggplot2::ylab(ylab)
  } else {
    dc_plot = dc_plot + ggplot2::ylab(paste0("Expression ", geneB))
  }

  if(!is.null(xlab)){
    dc_plot = dc_plot + ggplot2::xlab(xlab)
  } else {
    dc_plot = dc_plot + ggplot2::xlab(paste0("Expression ", geneA))
  }

  return(dc_plot)

}
