#' @title Calculate q-values from DGCA class objects based on permutation-based empirical null statistics.
#' @description First, estimate empirical p-values based on a comparison of the actual and permuted test statistics. Next, estimate the proportion of true null hypotheses using the qvalue package as well as qvalues from the empirical p-values, using this value. If the estimated pi0 <= 0, then sequentially recalculates using increasingly conservative set of lambda values, until lambda = 0.5.
#' @param dcObject The original S4 class object containing the test statistics to be extracted.
#' @param permObject The array of matrices containing the null test statistics.
#' @param secondMat Logical, indicating whether a second matrix was used in the construction of this dcObject and permObject. If FALSE, the upper.tri of both are extracted to avoid double counting test statistics.
#' @param testSlot The slot of the dcObject to be removed for use as the actual test statistic.
#' @param verbose Whether summaries of the q-value operations should be reported.
#' @param plotFdr Allows for plotting of fdrtool p-value adjustment result OR empirical FDR q-value adjustment technique, if either of these are chosen. Requires fdrtool package OR qvalue package. Default = FALSE.
#' @return A list containing a vectof of empirical p-values and a vector of q-values, both of the same length as the original actual test statistics.
#' @export
permQValue <- function(dcObject, permObject, secondMat, testSlot,
  verbose = FALSE, plotFdr = FALSE){

  test_stat_actual = slot(dcObject, testSlot)
  if(!secondMat) test_stat_actual = test_stat_actual[upper.tri(test_stat_actual)]
  test_stat_actual = as.numeric(test_stat_actual)

  if(!secondMat) permObject = permObject[apply(permObject, 3, upper.tri)]
  perm_stats = as.numeric(permObject)

  message("Calculating empirical p-values using the permutation sample statistics.")

  pvalues = bigEmpPVals(stat = abs(test_stat_actual), stat0 = abs(perm_stats))

  message("Calculating qvalues from the empirical p-values.")

  qobj = tryCatch(
    {
      WGCNA::qvalue(p = pvalues, lambda = seq(0.05, 0.95, 0.05))
    }, error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      cat("\n")
      message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
      return(NA)
    })
  #if the qvalue computation returned without error, then its format should be a list; if not, there was an error.
  if(!is.list(qobj)){
    qobj = tryCatch(
      {
        WGCNA::qvalue(p = pvalues, lambda = seq(0.1, 0.9, 0.05))
      }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        cat("\n")
        message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
        return(NA)
      })
  }
  if(!is.list(qobj)){
    qobj = tryCatch(
      {
        WGCNA::qvalue(p = pvalues, lambda = seq(0.2, 0.8, 0.05))
      }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        cat("\n")
        message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
        return(NA)
      })
  }
  if(!is.list(qobj)){
    qobj = tryCatch(
      {
        WGCNA::qvalue(p = pvalues, lambda = seq(0.3, 0.7, 0.05))
      }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        cat("\n")
        message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence... if this doesn't work, will report the empirical p-values and the adjusted q-values as NA values to indicate that q-value adjustment did not work.")
        qobj$qvalues = rep(NA, length(pvalues))
        return(qobj)
      })
  }

  if(verbose) summary(qobj)
  if(plotFdr) plot(qobj)

  return(list(empPVals = pvalues, pValDiffAdj = qobj$qvalues))

}
