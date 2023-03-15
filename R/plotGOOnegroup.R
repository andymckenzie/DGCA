#' @title Plot results from a hypergeometric enrichment test for one condition.
#' @description Uses ggplot2 to create a horizontal bar plot of the p-values (or odds-ratios) from enrichment tests of GO terms derived from differentially correlated gene sets (or any gene sets inputted into upstream functions). Note that the first column of each data frame is removed to allow for row binding, and otherwise the column names should match.
#' @param dfList A named list of data frames corresponding to different GO term enrichments.
#' @param nTerms The number of most-enriched terms to plot from each GO term type.
#' @param minSize The number of genes above which a gene set should be removed from analysis (e.g., because it is so small as to be overly specific and untrustworthy).
#' @param maxSize The number of genes above which a gene set should be removed from analysis (e.g., because it is so big as to be overly generic and relatively uninteresting).
#' @param dataCol Column of the input matrix to be plotted in the bar plot. If "Pvalue", it will be -log10 transformed prior to plotting. If not "Pvalue", the x-axis label should be changed manually following the function call.
#' @param namesCol The column specifying the name of the GO terms to be plotted.
#' @param labelsCol The column specifying the fill labels of the GO terms to be plotted.
#' @param legendTitle The title for the legend in the resulting plot.
#' @param adjustPVals Logical, indicating whether or not p-values should be adjusted by the Benjamini-Hochberg method.
#' @return An ggplot2 object that can be called in order to plot it, saved, or modified.
#' @export
plotGOOneGroup <- function(dfList, nTerms = 5, minSize = 50, maxSize = 500, dataCol = "Pvalue",
  namesCol = "Term", labelsCol = "Ontology", legendTitle = "GO Type", adjustPVals = FALSE){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The R package ggplot2 is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors", FALSE)
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  for(i in 1:length(dfList)){
    tmp = dfList[[i]]
    if(adjustPVals){
      tmp$Pvalue = p.adjust(tmp$Pvalue, method = "BH")
    }
    tmp = tmp[tmp$Size > minSize, ]
    tmp = tmp[tmp$Size < maxSize, ]
    #remove the GO ID column to enable row binding
    tmp = tmp[ , -1]
    tmp = tmp[order(tmp[ , dataCol], decreasing = FALSE), ]
    tmp = head(tmp, nTerms)
    tmp = tmp[order(tmp[ , dataCol], decreasing = TRUE), ]
    dfList[[i]] = tmp
  }

  data_df = rbind(dfList[[1]], dfList[[2]], dfList[[3]])

  if(dataCol == "Pvalue"){
    data_df[ , dataCol] = -log(data_df[ , dataCol], 10)
  }

  data_df[ , namesCol] = factor(data_df[ , namesCol], levels = data_df[ , namesCol])
  data_df[ , labelsCol] = factor(data_df[ , labelsCol], levels = data_df[ , labelsCol])

  cbPalette = c("#CC79A7", "#E69F00", "#56B4E9")
  go_plot = ggplot2::ggplot(data_df, ggplot2::aes(x = get(namesCol), y = get(dataCol), fill = get(labelsCol))) +
    ggplot2::geom_bar(stat = 'identity') + ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = cbPalette, guide = ggplot2::guide_legend(title = "", reverse = TRUE)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = legendTitle))

  if(dataCol == "Pvalue"){
    go_plot = go_plot + ggplot2::ylab(expression(paste("-Log"[10], " p-Value", sep = " ")))
  }

  return(go_plot)

}
