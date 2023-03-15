#' @title Integration function to use MEGENA to perform network analyses of DGCA results.
#' @description Takes a table of results from a DGCA analysis and inputs it into the MEGENA package pipeline.
#' @param ddcor_res The table of differential correlations outputted from ddcor. Expected to have pValDiff or pValDiff_adj columns as well as zScoreDiff, Gene1, +/- Classes columns.
#' @param pval_gene_thresh p-value threshold to call a gene as having significant differential correlation or not.
#' @param adjusted Logical indicating whether adjusted p-values from the differential correlation table (i.e., column "pValDiff_adj", when adjusted = TRUE) or unadjusted p-values (i.e., column "pValDiff", when adjusted = FALSE) should be used to subset the table into significant and non-significant portions.
#' @param evalCompactness Logical indicating whether or not the resulting modules should be filtered for compactness. For inputs with relatively small numbers of significant gene pairs, this may not be desirable. Note that if this option is not chosen, all of the modules will be returned, but some of the module-specific results will not be available for all of these modules.
#' @param nPerm The number of permutations to use in evaluating module hubs and module compactness in do.MEGENA.
#' @param modulePVal The p-value threshold used to include or disclude modules following module compactness evaluation in do.MEGENA.
#' @param hubPVal The p-value threshold used to classify a gene as a hub within a module.
#' @param minModSize The minimum module size.
#' @param maxModSize The minimum module size.
#' @param saveOutput Whether the output of MEGENA should be saved in the current directory. Default = FALSE.
#' @param parallelize Logical indicating whether or not multiple cores should be utilized as a form of parallel processing. Requires the doMC package.
#' @param nCores If parallelize is TRUE, the number of cores to use in the processing. Ignored if parallelize is FALSE.
#' @param ... Additional arguments to do.MEGENA from the MEGENA R package.
#' @return A list containing a the planar filter network, the data frame of identified differentially correlated modules, as well as various other objects including module-specific hub genes, depending on the parameters chosen.
#' @export
ddMEGENA <- function(ddcor_res, adjusted = TRUE, pval_gene_thresh = 0.05,
  evalCompactness = TRUE, nPerm = 100, hubPVal = 0.05, modulePVal = 0.05,
  minModSize = 10, maxModSize = 1000, saveOutput = FALSE, parallelize = FALSE,
  nCores = 4, ...){

  ##############################
	#set SAF to FALSE while restoring to default when the function is finished
	SAF = getOption("stringsAsFactors", FALSE)
	on.exit(options(stringsAsFactors = SAF))
	options(stringsAsFactors = FALSE)

  if(!adjusted) {
    ddcor_res_sig = ddcor_res[ddcor_res$pValDiff < pval_gene_thresh, ]
  } else {
    ddcor_res_sig = ddcor_res[ddcor_res$pValDiff_adj < pval_gene_thresh, ]
  }

  if (!requireNamespace("MEGENA", quietly = TRUE)) {
    stop("The R package MEGENA is needed for the ddMEGENA function to work. Please install it.",
      call. = FALSE)
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The R package igraph is needed for the ddMEGENA function to work. Please install it.",
      call. = FALSE)
  }

  if(parallelize){
    if (!requireNamespace("doMC", quietly = TRUE)) {
      stop("The R package doMC is needed for the ddMEGENA function to be parallelized. If you are on Windows, it and therefore parallelization of MEGENA through this function is unfortunately not available. Otherwise, please install it.",
        call. = FALSE)
    }
  }

  ddcor_res_megena = ddcor_res_sig[ , colnames(ddcor_res_sig) %in% c("Gene1", "Gene2", "zScoreDiff")]
  ddcor_res_megena$zScoreDiff = abs(ddcor_res_megena$zScoreDiff)

  pfn_res = MEGENA::calculate.PFN(ddcor_res_megena, doPar = parallelize, num.cores = nCores)

  #normalize to follow the convention that weights are in [0, 1]
  pfn_res$weight = (pfn_res$weight/max(pfn_res$weight)) * 0.999999999

  g = igraph::graph.data.frame(pfn_res, directed = FALSE)

  MEGENA.output = MEGENA::do.MEGENA(g, mod.pval = modulePVal, hub.pval = hubPVal,
    remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
    doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)

  if(!evalCompactness){
    MEGENA_modules = MEGENA.output$module.outpu$modules
    MEGENA_modules_df = data.frame(Genes = unlist(MEGENA_modules),
      Modules = rep(names(MEGENA_modules), sapply(MEGENA_modules, length)))
    megena_output = list(modules = MEGENA_modules_df, full = MEGENA.output)
  }

  if(evalCompactness){
    output = MEGENA::MEGENA.ModuleSummary(MEGENA.output, mod.pval = modulePVal, hub.pval = hubPVal,
      min.size = minModSize, max.size = maxModSize, output.sig = TRUE)
    MEGENA_modules_df = utils::stack(output$modules)
    colnames(MEGENA_modules_df) = c("Genes", "Modules")
    megena_output = list(modules = MEGENA_modules_df, summary = output, full = MEGENA.output)
  }

  return(megena_output)

}
