#' @title Get average empirical differential correlations.
#' @description Finds the average (median or mean) of differential correlations for either one gene versus all others or for all gene pairs in the input matrix of differences in z-scores between conditions, with significance calculated via a comparison with the permutation samples.
#' @param zDiff Matrix containing the actual difference of z-scores for each gene-gene pair.
#' @param zDiffPerm Matrix containing the differences of z-scores for each gene-gene pair in simulated data.
#' @param dCorAvgType Character vector specifying the type of average differential correlation calculation that should be performed. Types = c("gene_average", "total_average"). gene_average calculates whether each genes' differential correlation with all others is more than expected via permutation samples (and subsequent empirical FDR adjustment, in the case of > 1 gene), while total_average calculates whether the total average differential correlation is higher than expected via permutation samples. If splitSet is specified, then only genes in the splitSet have their average gene differential correlation calculated if gene_average is chosen.
#' @param dCorAvgMethod Character vector specifying the method for calculating the "average" differential correlation calculation that should be used. Options = "median", "mean".
#' @param oneSidedPVal If the dCorAvgType test is total_average, this option specifies whether a one-sided p-value should be reported, as opposed to a two-sided p-value. That is, if the average difference of z-scores is greater than zero, test whether the permutation average difference of z-scores are less than that average to get the p-value, and vice versa for the case that the average difference of z-scores is less than 0. Otherwise, test whether the absolute value of the average difference in z-scores is greater than the absolute values of the permutation average difference in z-scores. Default = FALSE.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @return A list containing the average difference(s) in z-score and the empirical p-value of that statistic calculated using the permutation samples.
#' @export
dCorAvg <- function(zDiff, zDiffPerm, dCorAvgType, oneSidedPVal = FALSE, secondMat = FALSE,
  dCorAvgMethod = "median"){

  message("Calculating the differential correlation average.")

  if(!dCorAvgMethod %in% c("median", "mean")){
    stop("The differential correlation average method chosen must be one of median or mean.")
  }

  if(dCorAvgType == "gene_average"){
    if(!secondMat){
      #since the lower.tri is full of NAs in the case of no secondMat, this is necessary to calculate the values for each gene
      for(i in 1:dim(zDiffPerm)[3]){
        tmp = zDiffPerm[,,i]
        tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        zDiffPerm[,,i] = tmp
      }
      zDiff[lower.tri(zDiff)] = t(zDiff)[lower.tri(t(zDiff))]
    }
    zdiff_medians = numeric(ncol(zDiff))
    empirical_pval = numeric(ncol(zDiff))
    nGenes = ncol(zDiff)
    for(i in 1:nGenes){
      zdiffs_gene_nonself = zDiff[-i, i]
      if(dCorAvgMethod == "median"){
        zdiff_medians[i] = median(zdiffs_gene_nonself)
      } else if(dCorAvgMethod == "mean"){
        zdiff_medians[i] = mean(zdiffs_gene_nonself)
      }
      zdiff_perm_gene_nonself = zDiffPerm[-i, i, ]
      if(dCorAvgMethod == "median") {
        zdiff_perm_gene_medians = matrixStats::colMedians(zdiff_perm_gene_nonself)
      } else if(dCorAvgMethod == "mean") {
        zdiff_perm_gene_medians = colMeans(zdiff_perm_gene_nonself)
      }
      empirical_pval[i] = 1 - sum(abs(zdiff_medians[i]) > abs(zdiff_perm_gene_medians))/
        length(zdiff_perm_gene_medians)
    }
    avg_dcor_df = data.frame(Gene = colnames(zDiff), avgZDiff = zdiff_medians,
      empirical_pVal = empirical_pval)
    avg_dcor_df$pVal_adj = p.adjust(avg_dcor_df$empirical_pVal, method = "BH")
    avg_dcor_df = avg_dcor_df[order(abs(avg_dcor_df$avgZDiff), decreasing = TRUE), ]
    return(avg_dcor_df)
  }

  if(dCorAvgType == "total_average"){
    if(dCorAvgMethod == "median"){
      zdiff_median = median(zDiff, na.rm = TRUE)
      zdiff_median_permwise = apply(zDiffPerm, 3, median, na.rm = TRUE)
    } else if(dCorAvgMethod == "mean"){
      zdiff_median = mean(zDiff, na.rm = TRUE)
      zdiff_median_permwise = apply(zDiffPerm, 3, mean, na.rm = TRUE)
    }
    nPerms = dim(zDiffPerm)[3]
    if(oneSidedPVal){
      if(zdiff_median > 0){
        pVal = 1 - (sum(zdiff_median < zdiff_median_permwise)/nPerms)
      }
      if(zdiff_median < 0){
        pVal = 1 - sum(zdiff_median > zdiff_median_permwise)/nPerms
      }
    } else if(!oneSidedPVal){
      pVal = 1 - sum(abs(zdiff_median) > abs(zdiff_median_permwise))/nPerms
    }
    return(list(total_zdiff = zdiff_median, pVal = pVal))
  }

}
