#' DGCA: An R package for Differential Gene Correlation Analysis
#'
#' The DGCA package provides three major classes of functions:
#' 1) Functions to calculate correlations, correlation significances, and number of samples in each correlation calculation based on an input matrix and a design matrix.
#' 2) Functions to calculate differential correlations between regions, which in the current package can only be pairwise (i.e., one condition vs another).
#' 3) Functions to extract, sort information about the differential correlation calculations in a convenient format.
#' The first two functions comprise the discovery of differential correlation (ddcor) portion of the package, which is why the names of the functions and object names often begin with ddcor.
#' Note that DGCA makes use of the 	SAF = getOption("stringsAsFactors", FALSE); on.exit(options(stringsAsFactors = SAF)); options(stringsAsFactors = FALSE) design pattern many times in order to avoid errors related to stringsAsFactors in porting code to new environments. This should not affect the stringsAsFactors options in your environment; however, you may want to be aware of this.
#'
#' @docType package
#' @name DGCA
#' @import WGCNA matrixStats methods
#' @importFrom grDevices colorRampPalette heat.colors
#' @importFrom graphics par plot
#' @importFrom stats median model.matrix p.adjust pnorm pt quantile runif
#' @importFrom utils capture.output globalVariables head str
NULL
#> NULL

globalVariables(c("GeneA", "GeneB", "Gene", "Condition", "cols", "rows", "vals"))

#' @title Brain sample ages vector.
#' @description A vector specifying the ages of the brain cell types from the single-cell RNA-seq study. The total data set can be downloaded by following the links in the original paper.
#' @references Darmanis S, Sloan SA, Zhang Y, et al. A survey of human brain transcriptome diversity at the single cell level. Proc Natl Acad Sci USA. 2015;112(23):7285-90.
"ages_darmanis"

#' @title Design matrix of cell type specifications of the single-cell RNA-seq samples.
#' @description Cell type specifications were performed by the authors. The total data set can be downloaded by following the links in the original paper.
#' @references Darmanis S, Sloan SA, Zhang Y, et al. A survey of human brain transcriptome diversity at the single cell level. Proc Natl Acad Sci USA. 2015;112(23):7285-90.
"design_mat"

#' @title Single-cell gene expression data from different brain cell types.
#' @description This data set has been filtered to include the genes that are in the 95th percentile or above in both oligodendrocytes and neurons. The total data set can be downloaded by following the links in the original paper.
#' @references Darmanis S, Sloan SA, Zhang Y, et al. A survey of human brain transcriptome diversity at the single cell level. Proc Natl Acad Sci USA. 2015;112(23):7285-90.
"darmanis"
