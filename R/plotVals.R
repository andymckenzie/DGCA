#' @title Creates a dotplot of the overall values for an individual gene in multiple conditions.
#' @description Takes the original input matrix, a design matrix, and one gene symbols (row name of the original matrix) to plot its values in the conditions specified, using a dotplot, +/- a summary bar. Will remove NAs prior to plotting.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose values you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param gene The gene symbol (row identifier).
#' @param log Logical, indicating whether the data should be log2-transformed prior to plotting (after adding a small constant of 0.5 to avoid problems with the log transform and stabilize the variance with respect to the mean).
#' @param add_summary_bar Logical indicating whether to include a summary bar for each group.
#' @param summary_bar If summary bar included, type of bar to use to calculate. Options = "mean", "median"
#' @param summary_width Horizontal width of summary crossbar included.
#' @param dotplot_width The width of the dots in the dotplot. See ?geom_dotplot for more information.
#' @param dotplot_binwidth The binwidth size for the dots in the dotplot. If NULL, will be calculated automatically.
#' @param dotplot_size The diameter of the dots in the dotplot. Affects the chart relative to the dotplot binwidth. If NULL, will be calculated automatically.
#' @param ylab Override the y-axis label to one of your choice.
#' @param xlab Override the x-axis label to one of your choice.
#' @return A ggplot2 object that can be plotted, further modified, and/or saved.
#' @export
plotVals <- function(inputMat, design, compare, gene,
    log = FALSE, ylab = NULL, xlab = NULL, add_summary_bar = TRUE,
    summary_bar = "mean", summary_width = 0.75,
    dotplot_width = 0.5, dotplot_size = NULL, dotplot_binwidth = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The R package ggplot2 is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The R package cowplot is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors", FALSE)
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  if(!all(compare %in% colnames(design))) stop("Comparison group names must be column names in the design matrix.")

  if(!gene %in% rownames(inputMat)) stop("gene was not found in the rownames of the input matrix.")

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
    tmpGene = tmp[rownames(tmp) %in% gene, ]
    cond = rep(compare[i], ncol(tmp))

    tmp_mat = data.frame(as.numeric(tmpGene), as.character(cond))
    colnames(tmp_mat) = c("Gene", "Condition")

    if(i == 1){
      data_mat = tmp_mat
    } else {
      data_mat = rbind(data_mat, tmp_mat)
    }

  }

  data_mat = data_mat[!is.na(data_mat$Gene), ]

  data_mat$Condition = factor(data_mat$Condition, levels = unique(data_mat$Condition))

  if(log){
    data_mat$Gene = log(data_mat$Gene + 0.5, 2)
  }

  if(is.null(dotplot_binwidth)){
    dotplot_binwidth = (max(data_mat$Gene) - min(data_mat$Gene)) / 60
  }

  if(is.null(dotplot_size)){
    cond_table = table(data_mat$Condition)
    largest_cond = names(sort(cond_table, decreasing = TRUE))[1]
    data_mat_largest_cond = data_mat[data_mat$Condition == largest_cond, ]
    dotplot_size = 10 * sqrt(1/nrow(data_mat_largest_cond))
  }

  message(paste0("The dotplot_binwidth is ",
    prettyNum(dotplot_binwidth, digits = 3, format = "g"), " and the dotplot_size is ",
    prettyNum(dotplot_size, digits = 3, format = "g"), "."))

  vals_plot = ggplot2::ggplot(data_mat, ggplot2::aes(x = Condition, y = Gene, fill = Condition)) +
    ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",
      width = dotplot_width, dotsize = dotplot_size, binwidth = dotplot_binwidth) +
      cowplot::theme_cowplot()

  if(add_summary_bar){
    if(summary_bar == "median"){
      vals_plot = vals_plot + ggplot2::stat_summary(fun.y = median,
        fun.ymin = median, fun.ymax = median, geom = "crossbar", width = summary_width)
    }
    if(summary_bar == "mean"){
      vals_plot = vals_plot + ggplot2::stat_summary(fun.y = mean,
        fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = summary_width)
    }
  }

  if(!is.null(ylab)){
    vals_plot = vals_plot + ggplot2::ylab(ylab)
  } else {
    vals_plot = vals_plot + ggplot2::ylab(paste0("Expression ", gene))
  }

  if(!is.null(xlab)){
    vals_plot = vals_plot + ggplot2::xlab(xlab)
  } else {
    vals_plot = vals_plot + ggplot2::xlab("")
  }

  return(vals_plot)

}
