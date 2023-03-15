#' @title Perform module GO-trait correlation
#' @description Takes input vectors of gene symbols, labels of corresponding modules, and a universe gene set and leverages the GOstats package to perform GO enrichment analysis.
#' @param genes A character vector specifying gene symbols, present as rows in the inputMat, corresponding to each module label in the labels argument.
#' @param labels A character vector specifying module label names, one for each gene symbol in the genes argument, with overlap allowed (i.e., each gene can be in more than one module).
#' @param universe Character vector of gene symbols which should be used as the background in the hypergeomtric test. If using this in the context of a DGCA experiment, this gene list most likely should be the gene set post-filtering, but prior to differential correlation analysis.
#' @param HGNC_clean Logical indicating whether the input gene symbols should be switched to clean HGNC symbols using the checkGeneSymbols function from the R package HGNChelper. Only applies if HGNC symbols are inputted.
#' @param pval_GO_cutoff Cutoff for the unadjusted p-values of gene ontology terms in the enrichment tests that should be displayed in the resulting table.
#' @param HGNC_switch Logical indicating whether or not the input gene symbols need to be switched from HGNC to Ensembl, the latter of which is required for GOstats enrichment test. Note that this is done by selecting the first Ensembl symbol that maps to a particular HGNC symbol, which is not always unique. If you need more precision on the conversion, you should do this outside of the function and insert the Ensembl list to the function.
#' @param gene_ontology A string specifying the branch of GO that should be used for enrichment analysis. One of "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all". If "all" is chosen, then this function finds the enrichment for all of the terms and combines them into one table. Default = "all"
#' @param annotation The library indicating the GO annotation database from which the Go terms should be mapped to gene symbols. Default = "org.Hs.eg.db", which is the table for Homo sapiens. Other common choices include "org.Mm.eg.db", "org.Rn.eg.db". The corresponding annotation library needs to be installed.
#' @param calculateVariance Optionally, find the variance of the odds ratio for each enrichment test. In particular, this finds the standard error of the log odds ratio, which converges to a normal distribution much more quickly than the non-log OR.
#' @param conditional Logical specifying whether the GO analysis should be done conditionally to take into account the hierarchical structure of the GO database in making sense of the gene set enrichments.
#' @return A list of lists of df's, one corresponding to each module, containing GO enrichment information for each module in each of the GO categories selected.
#' @export
moduleGO <- function(genes, labels, universe, HGNC_clean = TRUE, HGNC_switch = TRUE,
  gene_ontology = "all", pval_GO_cutoff = 1, annotation = "org.Hs.eg.db",
  conditional = FALSE, calculateVariance = FALSE){

  if (!requireNamespace("GOstats", quietly = TRUE)) {
    stop("The R package GOstats is needed for the moduleGO function to work. Please install it.",
      call. = FALSE)
  }
  if(HGNC_clean){ if (!requireNamespace("HGNChelper", quietly = TRUE)) {
      stop("The R package HGNChelper is needed for cleaning HGNC symbols with the HGNC_clean option you chose. Please install it.",
        call. = FALSE)
  }}
  if(HGNC_switch){ if (!requireNamespace(annotation, quietly = TRUE)) {
      stop(paste0("The R package ", annotation, " is needed for the moduleGO function to work with the annotation argument you chose. Please install it."),
        call. = FALSE)
  }}

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors", FALSE)
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  if(HGNC_clean){
    universe = switchGenesToHGCN(universe)
    genes = switchGenesToHGCN(genes)
  }
  if(HGNC_switch){
    univmap = AnnotationDbi::select(get(annotation), universe, "ENTREZID", "SYMBOL")
    universe = univmap[!duplicated(univmap[,1]), 2]
  }
  universe = universe[!is.na(universe)]
  universe = universe[!duplicated(universe)]


  labels_names = unique(labels)
  moduleGO_list = list()

  for(i in 1:length(labels_names)){
    message(paste0("Finding GO enrichments for module #", i, ", which is called ", labels_names[i]))
    genes_vector = genes[labels == labels_names[i]]

    if(length(genes_vector) == 0){
      name_df = paste0("enrichment_", labels_names[i])
      moduleGO_list[[name_df]] = list()
      next
    }

    if(HGNC_switch){
      univmap = AnnotationDbi::select(get(annotation), genes_vector, "ENTREZID", "SYMBOL")
      genes_vector = univmap[!duplicated(univmap[,1]), 2]
    }
    genes_vector = genes_vector[!is.na(genes_vector)]
    genes_vector = genes_vector[!duplicated(genes_vector)]

    moduleGO_df_tmp = findGOTermEnrichment(gene_vector = genes_vector,
      universe = universe,
      pval_GO_cutoff = pval_GO_cutoff, HGNC_clean = HGNC_clean,
      HGNC_switch = HGNC_switch, gene_ontology = gene_ontology,
      conditional = conditional, annotation = annotation)

    if(calculateVariance){
      if(length(gene_ontology) > 1 | gene_ontology == "all"){
        if(length(moduleGO_df_tmp) > 0){
          for(j in 1:length(moduleGO_df_tmp)){
            moduleGO_df_tmp_go = moduleGO_df_tmp[[j]]
            n11 = moduleGO_df_tmp_go$Count
            n21 = moduleGO_df_tmp_go$Size - moduleGO_df_tmp_go$Count
            n12 = moduleGO_df_tmp_go$gene_set_size - moduleGO_df_tmp_go$Count
            n22 = moduleGO_df_tmp_go$universe_size - n11 - n12 - n21
            moduleGO_df_tmp_go$LogOddsRatioSE = sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
            moduleGO_df_tmp[[j]] = moduleGO_df_tmp_go
          }
        }
      } else {
        n11 = moduleGO_df_tmp_go$Count
        n21 = moduleGO_df_tmp_go$Size - moduleGO_df_tmp_go$Count
        n12 = length(genes_vector) - moduleGO_df_tmp_go$Count
        n22 = length(universe) - n11 - n12 - n21
        moduleGO_df_tmp$LogOddsRatioSE = sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
      }
    }
    name_df = paste0("enrichment_", labels_names[i])
    moduleGO_list[[name_df]] = moduleGO_df_tmp
  }

  return(moduleGO_list)

}

#' @title Extract results from the module GO analysis
#' @description Turns the list of lists from the moduleGO function into a more comprehensible data frame. Note that if a GO term enrichment does not exist for that module, it is set as NA.
#' @param moduleGO_list The list of list of data frames for each module to be turned into a data frame.
#' @param labels Optional, a list of module names. Optional; if not inputted, these will be extracted from the list of module GO enrichments.
#' @return A data frame summarizing the GO term enrichments from each group, with columns on the ordered by the minimum p-value for OR term enrichment in any group.
#' @export
extractModuleGO <- function(moduleGO_list, labels = NULL){

  n = length(moduleGO_list)

  if(is.null(labels)){
    labels = names(moduleGO_list)
    labels = sapply(strsplit(labels, "_", fixed = TRUE), "[[", 2)
  } else {
    if(!length(labels) == n) stop("The length of the module names must be the same as the number of modules in the GO list result.")
  }

  for(i in 1:n){
    tmp_list = moduleGO_list[[i]]
    #rbind each data frame
    tmp_df = NULL
    for(j in 1:length(tmp_list)){
      #first column should be the GO ID
      if(!grepl("ID", colnames(tmp_list[[j]][1]))) stop("First column of each data frame should be the GO ID type.")
      colnames(tmp_list[[j]])[1] = "GOID"
      if(nrow(tmp_list[[j]]) > 0){
        tmp_df = rbind(tmp_df, tmp_list[[j]])
      }
    }
    colnames(tmp_df) = gsub("OddsRatio", paste0("OR_", labels[i]), colnames(tmp_df))
    colnames(tmp_df) = gsub("Pvalue", paste0("pVal_", labels[i]), colnames(tmp_df))
    module_colnames_keep = c("GOID", paste0("OR_", labels[i]), paste0("pVal_", labels[i]))
    if(i == 1){
      all_modules_GO = tmp_df[ , colnames(tmp_df) %in% c("Term", "Ontology",
        "Size", module_colnames_keep)]
      all_modules_GO = all_modules_GO[ , c("Term", "GOID", "Ontology", "Size",
        paste0("OR_", labels[i]), paste0("pVal_", labels[i]))]
    } else {
      tmp_df = tmp_df[ , colnames(tmp_df) %in% c("Term", "GOID", "Ontology", "Size",
        paste0("OR_", labels[i]), paste0("pVal_", labels[i]))]
      all_modules_GO = merge(all_modules_GO, tmp_df, by = "GOID", all = TRUE)
      #replace the NAs
      all_modules_GO[is.na(all_modules_GO$Term.x), c("Term.x", "Ontology.x", "Size.x")] =
        all_modules_GO[is.na(all_modules_GO$Term.x), c("Term.y", "Ontology.y", "Size.y")]
      all_modules_GO = all_modules_GO[ , !colnames(all_modules_GO) %in% c("Term.y", "Ontology.y", "Size.y")]
      col_names = colnames(all_modules_GO); col_names = gsub("Ontology.x", "Ontology", col_names)
      col_names = gsub("Term.x", "Term", col_names); col_names = gsub("Size.x", "Size", col_names)
      colnames(all_modules_GO) = col_names
    }
  }

  all_modules_GO$min_pval = apply(as.matrix(all_modules_GO[ , grepl("pVal_", colnames(all_modules_GO))]), 1, min)

  all_modules_GO = all_modules_GO[order(all_modules_GO$min_pval, decreasing = FALSE), ]
  return(all_modules_GO)

}

#' @title Plot extracted results from module-based GO enrichment analysis using ggplot2.
#' @description Takes a data frame of enrichment results in multiple modules and plots the results. Note that if a GO term enrichment does not exist for that module, it is set as 0 for an OR or 1 for a p-value.
#' @param df The data frame of term enrichments to be plotted.
#' @param nTerms The number of terms for each module whose GO terms with the minimum enrichment p-values for that group should be plotted.
#' @param termVector Optional character vector of GO term strings to plot, overriding other options.
#' @param plotOR Logical, indicating whether odds ratios should be plotted on the heatmap, instead of -log10 p-values (the default).
#' @param modules Optional, a list of module names to plot. Optional; if not inputted, all of the module names in the data frame will be used.
#' @param heatmapColor Optional specification of the heatmap colors. If not specified, ?heat.colors will be used.
#' @param text_size Text size of axes and legend in plot.
#' @param axis_text_col Color of axis text.
#' @param axis_x_text_angle Angle of x-axis text.
#' @param guide_title Optionally, specify the title of legend.
#' @param coord_flip Whether the coordinates should be flipped.
#' @param adjust If p-values are plotted, whether or not the enrichment p-values from each module should be adjusted by the Benjamini-Hochberg method.
#' @return A data frame summarizing the GO term enrichments from each group, with columns on the ordered by the minimum p-value for OR term enrichment in any group.
#' @export
plotModuleGO <- function(df, nTerms = 5, termVector = NULL, modules = NULL,
  heatmapColor = NULL, plotOR = FALSE, axis_text_col = "black",
  axis_x_text_angle = 45, text_size = 10,
  guide_title = NULL, coord_flip = FALSE, adjust = TRUE){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The R package ggplot2 is needed for this function to work. Please install it.",
      call. = FALSE)
  }

  if(is.null(guide_title)){
    if(!plotOR) guide_title = expression(paste("-Log"[10], " ", italic(p), sep = ""))
    if(plotOR) guide_title = "Odds Ratio"
  }

  if(!plotOR){
    matrix_vals = df[ , grepl("pVal_", colnames(df))]
  } else {
    matrix_vals = df[ , grepl("OR_", colnames(df))]
  }

  matrix_vals = as.matrix(matrix_vals)

  if(plotOR) matrix_vals[is.na(matrix_vals)] = 0
  if(!plotOR) matrix_vals[is.na(matrix_vals)] = 1

  if(!plotOR & adjust){
    matrix_vals = t(apply(matrix_vals, 1, p.adjust, method = "BH"))
  }

  colnames_modules = colnames(df)[grepl("pVal_", colnames(df))]
  module_names = sapply(strsplit(colnames_modules, "_", fixed = TRUE), "[[", 2)
  rownames(matrix_vals) = make.unique(df$Term)
  colnames(matrix_vals) = make.unique(module_names)

  if(!is.null(modules)){
    matrix_vals = matrix_vals[ , colnames(matrix_vals) %in% modules]
  }

  if(!is.null(termVector)){
    matrix_vals = matrix_vals[rownames(matrix_vals) %in% termVector, ]
  } else {
    termsKeep = vector()
    for(i in 1:length(module_names)){
      tmp = matrix_vals[, i]
      if(!plotOR) tmp = sort(tmp, decreasing = FALSE)
      if(plotOR) tmp = sort(tmp, decreasing = TRUE)
      terms = head(names(tmp), nTerms)
      termsKeep = c(termsKeep, terms)
    }
    termsKeep = unique(termsKeep)
    matrix_vals = matrix_vals[rev(termsKeep), ]
  }

  dat = data.frame(rows = rep(rownames(matrix_vals), ncol(matrix_vals)),
    cols = rep(colnames(matrix_vals), each = nrow(matrix_vals)),
    vals = as.vector(matrix_vals), stringsAsFactors = FALSE)

  if(!plotOR) dat$vals = -log((dat$vals), 10)

  if(is.null(heatmapColor)) heatmapColor = rev(heat.colors(50))

  #keep ordering
  dat$rows = factor(dat$rows, levels = unique(dat$rows), ordered = TRUE)
  dat$cols = factor(dat$cols, levels = unique(dat$cols), ordered = TRUE)

  heatmap_plot = ggplot2::ggplot(dat, ggplot2::aes(factor(cols), factor(rows))) +
    ggplot2::geom_tile(ggplot2::aes(fill = vals), colour = 'black') +
    ggplot2::scale_fill_gradientn(colours = heatmapColor,
      guide = ggplot2::guide_colourbar(title = guide_title,
      label.theme = ggplot2::element_text(size = text_size, angle = 0))) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, title.theme = ggplot2::element_text(size= text_size))) +
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = axis_x_text_angle,
    size = text_size, hjust = 1, vjust = 1, colour = axis_text_col),
    axis.text.y = ggplot2::element_text(size = text_size, colour = axis_text_col))

  if(coord_flip) heatmap_plot = heatmap_plot + ggplot2::coord_flip()

  return(heatmap_plot)

}
