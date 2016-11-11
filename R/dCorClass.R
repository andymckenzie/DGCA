#' @title Classify differential correlations. 
#' @description Classifies identifiers (e.g., genes) into one of the different categories pairwise-differential correlation classes. These categories are one of the Cartesian products of "Up Correlation", "No Correlation", and "Down Correlation" in each of the conditions, as well as a category for "no significant differential correlation".
#' @param corsA Numeric vector of correlations between gene pairs in condition A.
#' @param corsB Numeric vector of correlations between gene pairs in condition B.
#' @param pvalsA Numeric vector of the significance of correlation calls between gene pairs in condition A.
#' @param pvalsB Numeric vector of the significance of correlation calls between gene pairs in condition B.
#' @param dCorPVals Numeric vector of the differential correlation p-value calls.
#' @param sigThresh If classify = TRUE, this numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.
#' @param corSigThresh Threshold at which the correlation p-values must be below in order to be called "significant". Default = 0.05.
#' @param convertClasses Logical indicating whether the returned classes should be in numeric (factor) format or character format indicating the "actual" class.
#' @return A numeric vector of classes derived from each of the input vectors.
#' @examples
#' rho1 = runif(100, -1, 1); rho2 = runif(100, -1, 1)
#' pvalsA = runif(100, 0, 1); pvalsB = runif(100, 0, 1); dcor_pvals = runif(100, 0, 1)
#' cor_classes = dCorClass(rho1, pvalsA, rho2, pvalsB, dcor_pvals)
#' cor_classes = dCorClass(rho1, pvalsA, rho2, pvalsB, dcor_pvals, convertClasses = TRUE)
#' @export
dCorClass <- function(corsA, pvalsA, corsB, pvalsB, dCorPVals, sigThresh = 1,
	corSigThresh = 0.05, convertClasses = FALSE){

  if(!(all.equal(length(corsA), length(corsB), length(pvalsA),
		length(pvalsB), length(dCorrs)))) stop("All of the input vectors should be the same length.")

	classes = rep(0, length(corsA))

	sigs = (dCorPVals < sigThresh)
	pvA = (pvalsA < corSigThresh)
	pvB = (pvalsB < corSigThresh)
	cAup = (corsA > 0)
	cBup = (corsB > 0)

	#UpUp
	classes[which(sigs & pvA & pvB & cAup & cBup)] = 1
	#UpNon
	classes[which(sigs & pvA & !pvB & cAup)] = 2
	#UpDown
	classes[which(sigs & pvA & pvB & cAup & !cBup)] = 3
	#NonUp
	classes[which(sigs & !pvA & pvB & cBup)] = 4
	#NonNon
	classes[which(sigs & !pvA & !pvB)] = 5
	#NonDown
	classes[which(sigs & !pvA & pvB & !cBup)] = 6
	#DownUp
	classes[which(sigs & pvA & pvB & !cAup & cBup)] = 7
	#DownNon
	classes[which(sigs & pvA & !pvB & !cAup)] = 8
	#DownDown
	classes[which(sigs & pvA & pvB & !cAup & !cBup)] = 9

	if(convertClasses){
		classes = factor(classes, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
			labels = c("NonSig", "+/+", "+/0", "+/-",
			"0/+", "0/0", "0/-", "-/+", "-/0", "-/-"))
	}

	return(classes)

}
