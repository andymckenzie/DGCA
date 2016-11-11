#' @title Use speed-optimized sorting to calculate p-values observed and simulated null test statistic using a reference pool distribution.
#' @description A reimplementation of qvalue::empPvals designed to work faster and require less memory in the average case. Unlike qvalue::empPvals, *requires* the use of a reference pool distribution rather than having this as an option. Another difference of this function that the original is that it handles ties for test statistics equal to 0, to handle cases where test statistics are thresholded and may be zero more commonly than expected by chance (i.e., very rarely).
#' @param stat The vector of test statistics for which p-values should be returned.
#' @param stat0 The vector or matrix of simulated null test statistics (if matrix, will be coerced to a vector).
#' @param increasing Logical indicating whether the test statistics because more extreme as they increase in value from 0. Negative numbers are not allowed as inputs (i.e., the test statistic must be monotonic).
#' @references Please see ?qvalue::empPVals for more; from which this function was adapted.
#' @return A vector of p-values adjusted for the null statistics calculated, to be used as an input to the qvalue function.
#' @author John Storey, Andrew McKenzie
#' @examples
#' test_stat = rnorm(100, 1, 1)
#' test_stat0 = rnorm(1000, 0, 1)
#' emp_pvals = bigEmpPVals(test_stat, test_stat0)
#' @export
bigEmpPVals <- function(stat, stat0, increasing = TRUE){

  m = length(stat)
  m0 = length(stat0)

  if(is.matrix(stat0)) stat0 = as.vector(stat0)

  message("Sorting the combination of the actual and permuted test statistics.")
  #create and order a logical vector for the presence of actual vs empirical test stats
  v = c(rep(TRUE, m), rep(FALSE, m0))
  if(increasing){
    v = v[sort.list(c(stat, stat0), decreasing = TRUE, method = "quick", na.last = NA)]
  } else {
    v = v[sort.list(c(stat, stat0), decreasing = FALSE, method = "quick", na.last = NA)]
  }

  message("Finding the proportion of each actual test statistic greater than the permuted test statistics.")
  #at each position of the real test stats, find the proportion of real test stats that are greater than the null test stats
  u = 1:length(v)
  w = 1:m
  p = (u[v == TRUE] - w) / m0

  #handle cases where the test stat equals zero by making p-values uniform in this case
  if(min(stat, na.rm = TRUE) == 0){
    n_zeros = sum(stat == 0, na.rm = TRUE)
    first_zero = p[m-n_zeros+1]
    p[(m-n_zeros+1):m] = runif(n_zeros, min = first_zero, max = 1)
  }

  #match to the original test stats; randomly in the case of ties to avoid double-counting and introducing redundancy
  if(increasing){
    p = p[rank(-stat, ties.method = "random")]
  } else {
    p = p[rank(stat, ties.method = "random")]
  }

  #correct numbers so that the actual test stats that are greater than all of the nulls can only have a min p-value of at most 1/#null test stats
  p = pmax(p, 1/m0)

  return(p)

}
