#' @title Creates a data frame for the top differentially correlated gene pairs in your data set.
#' @description Reads in a dcPair object and outputs a table of all gene pairs (or just the top n pairs), sorted by their unadjusted differential correlation p-value.
#' @param dcObject The dcPair class object which you'd like to convert into a table.
#' @param nPairs The number of gene pairs to display in the resulting table.
#' @param adjust Allows for resulting p-values to be corrected for multiple hypothesis tests, optional. Some non-default choices require the "fdrtool" package or the "qvalue". Default = "none", which means that no p-value adjustment is performed. Other options include "perm" to use permutation samples, methods in ?p.adjust (i.e., "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), and methods in ?fdrtool (i.e., "fndr", "pct0", "locfdr").
#' @param plotFdr Allows for plotting of p-value adjustment result, if this is chosen. Requires fdrtool or qvalue package. Default = FALSE.
#' @param classify Binary value specifying whether the correlation values in each condition and differential correlation scores should be used to classifying the resulting identifiers into groups. Default = TRUE
#' @param sigThresh If classify = TRUE, this numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.
#' @param corSigThresh If classify = TRUE, this numeric value specifies the p-value threshold at which a correlation p-value is deemed significant. Default = 0.05.
#' @param zScorePerm A matrix of values with z-scores from permutation tests to be used to generate empirical p-values. Default = NULL.
#' @param verbose Whether summaries of the operations should be reported.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @return A table containing columns for each name in the considered gene pair (the order of which is arbitrary), correlation values in each condition, differences in z-score of the correlation, and p-values for that z-score difference.
#' @export
dcTopPairs <- function(dcObject, nPairs, adjust = "none", plotFdr = FALSE,
	classify = TRUE, sigThresh = 1, corSigThresh = 0.05, zScorePerm = NULL,
	verbose = FALSE, compare = NULL, secondMat = FALSE) {

	SAF = getOption("stringsAsFactors", FALSE)
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

	if(!is(dcObject, 'dcPair')) stop("Input to dcTopPairs must be a dcPair object.")

	###################################
	# sort unadjusted p-values
	pValDiff = slot(dcObject, "PValDiff")

	#sort ignores NAs so, since only the upper.tri is non-NA, do not need to extract the upper.tri
	sorted = sort(pValDiff, index.return = 1)
	sorted_pvals = sorted$x
	sorted_indx = sorted$ix[1:nPairs]

	###################################
	# extract and create the top pairs data frame
	if(!secondMat){
		indx = which(upper.tri(pValDiff, diag = FALSE), arr.ind = TRUE)
		corA = slot(dcObject, "corA")
		gene1 = rownames(corA)[indx[,1]][sorted_indx]
		gene2 = colnames(corA)[indx[,2]][sorted_indx]
	} else {
		#get indices of the minimum values
		dIndx = arrayInd(sorted_indx, dim(pValDiff))
		corA = slot(dcObject, "corA")
		gene1 = rownames(corA)[dIndx[,1]]
		gene2 = colnames(corA)[dIndx[,2]]
	}

	extractAndOrderSlots <- function(extractSlot){
		tmp = slot(dcObject, extractSlot)
		if(!secondMat) tmp = tmp[upper.tri(tmp)]
		tmp = tmp[sorted_indx]
		return(tmp)
	}

	corA = extractAndOrderSlots("corA")
	corApVal = extractAndOrderSlots("corPvalA")
	corB = extractAndOrderSlots("corB")
	corBpVal = extractAndOrderSlots("corPvalB")
	zScoreDiff = extractAndOrderSlots("ZDiff")
	pValDiff_unadj = extractAndOrderSlots("PValDiff")

	if(is.null(compare)){
		colnames_topPairs = c("Gene1", "Gene2", "corA", "corA_pVal", "corB",
			"corB_pVal", "zScoreDiff", "pValDiff")
	} else {
		colnames_topPairs = c("Gene1", "Gene2",
			paste0(compare[1], "_cor"), paste0(compare[1], "_pVal"),
			paste0(compare[2], "_cor"), paste0(compare[2], "_pVal"), "zScoreDiff", "pValDiff")
	}

	topPairs = data.frame(gene1, gene2, corA, corApVal, corB, corBpVal,
		zScoreDiff, pValDiff_unadj, stringsAsFactors = FALSE)

	#####################################
	# adjust p-values
	if(!adjust == "none"){
		if(!adjust == "perm"){
			pValDiffAdj_sorted = adjustPVals(pValDiff_unadj, adjust = adjust) #has already been ordered
			topPairs = data.frame(topPairs, pValDiffAdj_sorted, stringsAsFactors = FALSE)
			colnames_topPairs = c(colnames_topPairs, "pValDiff_adj")
		} else {
			pValDiffAdj = permQValue(dcObject, zScorePerm, secondMat = secondMat, testSlot = "ZDiff",
				verbose = verbose, plotFdr = plotFdr)
			empPVals_sorted = pValDiffAdj[["empPVals"]][sorted_indx]
			pValDiffAdj_sorted = pValDiffAdj[["pValDiffAdj"]][sorted_indx]
			topPairs = data.frame(topPairs, empPVals_sorted,
				pValDiffAdj_sorted, stringsAsFactors = FALSE)
			colnames_topPairs = c(colnames_topPairs, "empPVals", "pValDiff_adj")
		}
	}

	message("Classifying the differential correlation calls.")

	if(classify == TRUE){
		if(adjust == "none"){
			class_list = dCorClass(corA, corApVal, corB, corBpVal, pValDiff_unadj,
				sigThresh = sigThresh, corSigThresh = corSigThresh)
		} else {
			class_list = dCorClass(corA, corApVal, corB, corBpVal, pValDiffAdj_sorted,
				sigThresh = sigThresh, corSigThresh = corSigThresh)
		}
		class_factor = factor(class_list, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
		labels = c("NonSig", "+/+", "+/0", "+/-",
			"0/+", "0/0", "0/-", "-/+", "-/0", "-/-"))
		topPairs = data.frame(topPairs, class_factor, stringsAsFactors = FALSE)
		colnames_topPairs = c(colnames_topPairs, "Classes")
	}

	colnames(topPairs) = colnames_topPairs
	return(topPairs)

}
