#' @title Ranks genes by their total number of differentially correlated gene pairs.
#' @description Returns list of lists for the top differentially correlated gene pairs in each direction and/or class.
#' @param ddcor_res The table of differential correlations outputted from ddcor. Expected to have pValDiff or pValDiff_adj columns as well as zScoreDiff, Gene1, +/- Classes columns.
#' @param pval_gene_thresh p-value threshold to call a gene as having significant differential correlation or not.
#' @param adjusted Logical indicating whether adjusted p-values from the differential correlation table (i.e., column "pValDiff_adj", when adjusted = TRUE) or unadjusted p-values (i.e., column "pValDiff", when adjusted = FALSE) should be used to subset the table into significant and non-significant portions.
#' @param geneNameCol Character vector specifying the name of the columns that are used to extract the gene symbols. Note that the default is c("Gene1", "Gene2"), but this only makes sense in the context of a full DGCA experiment. In the case of a splitSet, you may want to use "Gene1" to avoid counting the splitSet names in all of the categories.
#' @param nGenes Number of genes to display in the resulting table. Default = "all", but also can be restricted to a particular number.
#' @param classes Gets the number of differentially correlated gene pairs associated with each of the differential correlation classes.
#' @return A data frame with corresponding lists of genes most associated with each of the directions and/or correlation classes.
#' @examples
#' data(darmanis); data(design_mat); darmanis_subset = darmanis[1:30, ]
#' ddcor_res = ddcorAll(inputMat = darmanis_subset, design = design_mat,
#'  compare = c("oligodendrocyte", "neuron"))
#' top_genes = topDCGenes(ddcor_res)
#' @export
topDCGenes <- function(ddcor_res, adjusted = FALSE, pval_gene_thresh = 0.05,
  geneNameCol = c("Gene1", "Gene2"), nGenes = "all", classes = TRUE){

  .getTop <- function(ddcor_res){
    number_links = c(ddcor_res$Gene1, ddcor_res$Gene2)
    ddcor_res_genes = table(number_links)
    ddcor_res_genes = sort(ddcor_res_genes, decreasing = TRUE)
    if(!nGenes == "all"){
      top_genes = head(ddcor_res_genes, nGenes)
    } else {
      top_genes = ddcor_res_genes
    }
    return(top_genes)
  }

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors", FALSE)
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  if(!adjusted) {
    ddcor_res_sig = ddcor_res[ddcor_res$pValDiff < pval_gene_thresh, ]
  } else {
    ddcor_res_sig = ddcor_res[ddcor_res$pValDiffAdj < pval_gene_thresh, ]
  }

  res_list = list()

  res_list[["all"]] = .getTop(ddcor_res_sig)
  res_list[["zdiff_up"]] = .getTop(ddcor_res_sig[ddcor_res_sig$zScoreDiff > 0, ])
  res_list[["zdiff_down"]] = .getTop(ddcor_res_sig[ddcor_res_sig$zScoreDiff < 0, ])

  if(classes){
    res_list[["class_+/+"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "+/+", ])
    res_list[["class_+/0"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "+/0", ])
    res_list[["class_+/-"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "+/-", ])
    res_list[["class_0/+"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "0/+", ])
    res_list[["class_0/0"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "0/0", ])
    res_list[["class_0/-"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "0/-", ])
    res_list[["class_-/+"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "-/+", ])
    res_list[["class_-/0"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "-/0", ])
    res_list[["class_-/-"]] = .getTop(ddcor_res_sig[ddcor_res_sig$Classes == "-/-", ])
  }

  return(res_list)

}
