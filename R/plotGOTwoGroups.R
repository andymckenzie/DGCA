#' @title Plot results from a hypergeometric enrichment test to compare two conditions.
#' @description Uses plotrix to create a pyramid plot of the odds-ratios from enrichment tests of GO terms derived from differentially correlated gene sets (or any gene sets inputted into upstream functions) in two conditions. Note that the first column of each data frame is removed to allow for row binding, and otherwise the column names should match.
#' @param dfList1 A named list of data frames corresponding to different GO term enrichments from the first condition or differential correlation class.
#' @param dfList2 A named list of data frames corresponding to different GO term enrichments from the second condition or differential correlation class.
#' @param nTerms The number of most difference-in-enrichment (i.e., difference adjusted p-value for the difference in log OR) terms to plot from each GO term type.
#' @param minSize The number of genes above which a gene set should be removed from analysis (e.g., because it is so small as to be overly specific and untrustworthy).
#' @param maxSize The number of genes above which a gene set should be removed from analysis (e.g., because it is so big as to be overly generic and relatively uninteresting).
#' @param pValCutoff p-Value cutoff to specify how significant each term enrichment must be in at least one of the groups to be considered for the comparison between the conditions. If no cutoff is desired, set to 1.
#' @param GOTermTypes Character vector for the GO term types to be plotted.
#' @param labelsCol The column specifying the fill labels of the GO terms to be plotted.
#' @param filterSignificant Logical, indicating whether or not the p-value for the difference in ORs between groups should be required to be less than a particular threshold prior to downstream analysis.
#' @param filterSigThresh Number indicating the threshold of the p-value for the difference in ORs between groups required if filterSignificant is TRUE.
#' @param adjustPVals Logical, indicating whether or not p-values for the difference in log odds between conditions should be adjusted by the Benjamini-Hochberg method.
#' @param plotrix_gap Parameter specifying the size of the gap between the two sides of the ORs.
#' @param labels Character vector specifying labels for the pyramid plot; first entry = left label, second entry = middle label, third entry = right label.
#' @param fill_zero_cats Logical, indicating whether counts of zero in some groups should be included in the comparisons. Adds 0.1 to the ORs and in the denominator of the count for the calculation of the SE in the cases where there are zero counts identified, as a convservative measure to prevent finding infinite differences between groups.
#' @return A data frame with the values used to create the plot.
#' @export
plotGOTwoGroups <- function(dfList1, dfList2, nTerms = 5, minSize = 40, maxSize = 1000,
  labelsCol = "Ontology", adjustPVals = TRUE, plotrix_gap = 20,
  GOTermTypes = c("BP", "CC", "MF"), pValCutoff = 0.05, filterSignificant = FALSE,
  filterSigThresh = 0.05, labels = c("Corr Class 1", "GO Term Name", "Corr Class 2"), fill_zero_cats = FALSE){

  if (!requireNamespace("plotrix", quietly = TRUE)) {
    stop("The R package plotrix is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  data_df = data.frame()

  for(i in 1:length(GOTermTypes)){
    tmp_res_merge = merge(dfList1[[GOTermTypes[i]]],
      dfList2[[GOTermTypes[i]]], by = paste0("GO", GOTermTypes[i], "ID"), all = TRUE)

    tmp_res_merge = tmp_res_merge[ , !colnames(tmp_res_merge) %in% paste0("GO", GOTermTypes[i], "ID")]

    tmp_res_merge$OddsRatio.x_plot = tmp_res_merge$OddsRatio.x
    tmp_res_merge$OddsRatio.y_plot = tmp_res_merge$OddsRatio.y

    if(fill_zero_cats == TRUE){
      for(i in 1:length(tmp_res_merge)){
        if(is.na(tmp_res_merge$Count.x[i])){
          tmp_res_merge$Count.x[i] = 0
          tmp_res_merge$OddsRatio.x[i] = 0.1
          tmp_res_merge$OddsRatio.x_plot[i] = 0
          tmp_res_merge$Pvalue.x[i] = 1
          tmp_res_merge$Ontology.x[i] = tmp_res_merge$Ontology.y[i]
          tmp_res_merge$Term.x[i] = tmp_res_merge$Term.y[i]
          tmp_res_merge$Size.x[i] = tmp_res_merge$Size.y[i]
          tmp_res_merge$gene_set_size.x[i] = mean(tmp_res_merge$gene_set_size.x, na.rm = TRUE)
          tmp_res_merge$universe_size.x[i] = mean(tmp_res_merge$universe_size.x, na.rm = TRUE)
          n11 = tmp_res_merge$Count.x[i]
          n21 = tmp_res_merge$Size.x[i] - tmp_res_merge$Count.x[i]
          n12 = tmp_res_merge$gene_set_size.x[i] - tmp_res_merge$Count.x[i]
          n22 = tmp_res_merge$universe_size.x[i] - n11 - n12 - n21
          tmp_res_merge$LogOddsRatioSE.x[i] = sqrt((1/(n11 + 0.1)) + (1/n12) + (1/n21) + (1/n22))
        }
        if(is.na(tmp_res_merge$Count.y[i])){
          tmp_res_merge$Count.y[i] = 0
          tmp_res_merge$OddsRatio.y[i] = 0.1
          tmp_res_merge$OddsRatio.y_plot[i] = 0
          tmp_res_merge$Pvalue.y[i] = 1
          tmp_res_merge$Ontology.y[i] = tmp_res_merge$Ontology.x[i]
          tmp_res_merge$Term.y[i] = tmp_res_merge$Term.x[i]
          tmp_res_merge$Size.y[i] = tmp_res_merge$Size.x[i]
          tmp_res_merge$gene_set_size.y[i] = mean(tmp_res_merge$gene_set_size.y, na.rm = TRUE)
          tmp_res_merge$universe_size.y[i] = mean(tmp_res_merge$universe_size.y, na.rm = TRUE)
          n11 = tmp_res_merge$Count.y[i]
          n21 = tmp_res_merge$Size.y[i] - tmp_res_merge$Count.y[i]
          n12 = tmp_res_merge$gene_set_size.y[i] - tmp_res_merge$Count.y[i]
          n22 = tmp_res_merge$universe_size.y[i] - n11 - n12 - n21
          tmp_res_merge$LogOddsRatioSE.y[i] = sqrt((1/(n11 + 0.1)) + (1/n12) + (1/n21) + (1/n22))
        }
      }
    }

    #filter to only include groups with the appropriate characteristics
    tmp_res_merge_filt = tmp_res_merge[tmp_res_merge$Size.x > minSize & tmp_res_merge$Size.y > minSize, ]
    tmp_res_merge_filt = tmp_res_merge_filt[tmp_res_merge_filt$Size.x < maxSize & tmp_res_merge_filt$Size.y < maxSize, ]
    tmp_res_merge_filt = tmp_res_merge_filt[tmp_res_merge_filt$Pvalue.x < pValCutoff | tmp_res_merge_filt$Pvalue.y < pValCutoff, ]

    #get the difference in log odds
    tmp_res_merge_filt$diffLogOdds = log(tmp_res_merge_filt$OddsRatio.x) - log(tmp_res_merge_filt$OddsRatio.y)

    #get the significance of the difference in log odds
    tmp_res_merge_filt$diffLogOddsZScore = abs(tmp_res_merge_filt$diffLogOdds /
      sqrt(tmp_res_merge_filt$LogOddsRatioSE.x^2 + tmp_res_merge_filt$LogOddsRatioSE.y^2))
    tmp_res_merge_filt$diffOddsPValue = pnorm(tmp_res_merge_filt$diffLogOddsZScore, lower.tail = FALSE)

    #adjust the difference in ORs P-Values for MHT
    columnSel = "diffOddsPValue"
    if(adjustPVals){
      tmp_res_merge_filt$diffOddsPValue_adj = p.adjust(tmp_res_merge_filt$diffOddsPValue, method = "BH")
      columnSel = "diffOddsPValue_adj"
    }

    #optionally, filter for only those categories with a significant difference
    if(filterSignificant){
      tmp_res_merge_filt = tmp_res_merge_filt[tmp_res_merge_filt[ , columnSel] < filterSigThresh, ]
    }

    #find the top N terms from each of the groups
    tmp_res_merge_filt = tmp_res_merge_filt[order(tmp_res_merge_filt[ , columnSel], decreasing = FALSE), ]
    tmp_res_merge_filt_up1 = tmp_res_merge_filt[tmp_res_merge_filt$diffLogOdds > 0, ]
    tmp_df1 = head(tmp_res_merge_filt_up1, nTerms)
    tmp_res_merge_filt_up2 = tmp_res_merge_filt[tmp_res_merge_filt$diffLogOdds < 0, ]
    tmp_df2 = head(tmp_res_merge_filt_up2, nTerms)

    tmp_df = rbind(tmp_df1, tmp_df2)
    tmp_df = tmp_df[order(tmp_df[ , columnSel], decreasing = FALSE), ]

    tmp_df = tmp_df[!is.na(tmp_df$Pvalue.x), ]

    if(i == 1){
      data_df = tmp_df
    } else {
      data_df = rbind(data_df, tmp_df)
    }
  }

  if(nrow(data_df) == 0) stop("There are no GO categories to print that fit the specified criteria.")

  print(data_df)

  data_df = data_df[!duplicated(data_df$Term.x), ]

  labels_OR = 1:ceiling(max(c(data_df$OddsRatio.x, data_df$OddsRatio.y)))

  colors = c(rep("#E69F00", sum(data_df$Ontology.y == "MF")),
    rep("#56B4E9", sum(data_df$Ontology.y == "CC")),
    rep("#009E73", sum(data_df$Ontology.y == "BP")))
  par(mar = plotrix::pyramid.plot(rev(data_df$OddsRatio.x_plot), rev(data_df$OddsRatio.y_plot),
    labels = rev(data_df$Term.x), gap = plotrix_gap, show.values = TRUE, unit = "Odds Ratio",
    top.labels = labels, laxlab = labels_OR, raxlab = labels_OR,
    lxcol = colors, rxcol = colors))

  return(data_df)

}
