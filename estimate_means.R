
library(edgeR)
library(boot)

#parameters
bootstrap_runs = 1
level = "exon"

setwd("/Users/amckenz/Dropbox/zhang/brain_gnxp/")

#functions
source("/Users/amckenz/Dropbox/zhang/brain_gnxp/cpm_uq.R")

std_error <- function(x) sd(x)/sqrt(length(x))

#for bootstrap
samplemean <- function(x, d) {
  return(mean(x[d]))
}

source("/Users/amckenz/Dropbox/zhang/brain_gnxp/cell_mean_estimation_pipeline_faster.R")

darmanis_gnxp = readRDS(paste0("processed_data/darmanis_gnxp_matched_", level, ".rds")) #darmanis_gnxp_no_outliers_
darmanis_meta = readRDS(paste0("processed_data/darmanis_meta_matched_", level, ".rds"))
darmanis_features = readRDS(paste0("processed_data/darmanis_features_", level, ".rds"))
# darmanis_features$geneSymbol = make.unique(darmanis_features$geneSymbol)
darmanis_df_full = cell_mean_estimation_pipeline_faster(darmanis_gnxp, darmanis_meta, darmanis_features, p = 0.99, bootstrap = FALSE)

#save one file per gene symbol
for(i in 1:length(unique(darmanis_features$Geneid))){
  ensembl = unique(darmanis_features$Geneid)[i]
  exons_tmp = darmanis_df_full[darmanis_features$Geneid == ensembl, ]
  symbol_tmp = darmanis_features[darmanis_features$Geneid == ensembl, "geneSymbol"][1]
  # saveRDS(exons_tmp, paste0("cell_means/darmanis/", symbol_tmp, "_", ensembl, "_", level, ".rds"))
  saveRDS(exons_tmp, paste0("cell_means/darmanis/", level, "/", symbol_tmp, "_", level, ".rds"))
}

#zeisel
zeisel_gnxp = readRDS(paste0("processed_data/zeisel _gnxp_matched_", level, ".rds")) #zeisel _gnxp_no_outliers_
zeisel_meta = readRDS(paste0("processed_data/zeisel _meta_matched_", level, ".rds"))
zeisel_features = readRDS(paste0("processed_data/zeisel _features_", level, ".rds"))
# zeisel _features$geneSymbol = make.unique(zeisel _features$geneSymbol)
zeisel_df_full = cell_mean_estimation_pipeline_faster(zeisel_gnxp, zeisel_meta, zeisel_features, p = 0.99, bootstrap = FALSE)

for(i in 1:length(unique(darmanis_features$Geneid))){
  ensembl = unique(darmanis_features$Geneid)[i]
  exons_tmp = darmanis_df_full[darmanis_features$Geneid == ensembl, ]
  symbol_tmp = darmanis_features[darmanis_features$Geneid == ensembl, "geneSymbol"][1]
  # saveRDS(exons_tmp, paste0("cell_means/darmanis/", symbol_tmp, "_", ensembl, "_", level, ".rds"))
  saveRDS(exons_tmp, paste0("cell_means/darmanis/", level, "/main_cells/", symbol_tmp, "_", level, ".rds"))
}

zeisel_meta$celltypes = zeisel_meta$cell_types
zeisel_meta$celltypes = zeisel_meta$cell_types
zeisel_meta$celltypes = gsub("ClauPyramidal", "(none)", zeisel_meta$celltypes)
zeisel_df_full_subtypes = cell_mean_estimation_pipeline_faster(zeisel_gnxp, zeisel_meta, zeisel_features, p = 0.99, bootstrap = FALSE)

for(i in 1:length(unique(darmanis_features$Geneid))){
  ensembl = unique(darmanis_features$Geneid)[i]
  exons_tmp = darmanis_df_full[darmanis_features$Geneid == ensembl, ]
  symbol_tmp = darmanis_features[darmanis_features$Geneid == ensembl, "geneSymbol"][1]
  # saveRDS(exons_tmp, paste0("cell_means/darmanis/", symbol_tmp, "_", ensembl, "_", level, ".rds"))
  saveRDS(exons_tmp, paste0("cell_means/darmanis/", level, "/sub_cells/", symbol_tmp, "_", level, ".rds"))
}


# #zeisel data
#for zeisel, need to write the
# zeisel_gnxp = readRDS(paste0("processed_data/zeisel_gnxp_no_outliers_", level, ".rds"))
# zeisel_meta = readRDS(paste0("processed_data/zeisel_meta_no_outliers_", level, ".rds"))
# zeisel_features = readRDS(paste0("processed_data/zeisel_features_", level, ".rds"))
# zeisel_features$geneSymbol = zeisel_features$Symbol
# # zeisel_gnxp_sub = zeisel_gnxp[1:2000, ]
# # zeisel_features_sub = zeisel_features[1:2000, ]
# # zeisel_df = cell_mean_estimation_pipeline(gnxp = zeisel_gnxp_sub, meta = zeisel_meta, features = zeisel_features_sub, p = 0.99, bootstrap = TRUE)
# zeisel_df_full = cell_mean_estimation_pipeline(zeisel_gnxp, zeisel_meta, zeisel_features, p = 0.99, bootstrap = FALSE)
# saveRDS(zeisel_df_full, paste0("cell_means/zeisel_", level, ".rds"))
#
# #tasic data
# tasic_gnxp = readRDS(paste0("processed_data/tasic_gnxp_no_outliers_", level, ".rds"))
# tasic_meta = readRDS(paste0("processed_data/tasic_meta_no_outliers_", level, ".rds"))
# tasic_meta$celltypes = tasic_meta$sub_class
# tasic_meta$celltypes[grepl("L", tasic_meta$celltypes)] = "Neuron"
# tasic_meta$celltypes[grepl("Sst", tasic_meta$celltypes)] = "Neuron"
# tasic_meta$celltypes = gsub("Vip", "Neuron", tasic_meta$celltypes)
# tasic_meta$celltypes = gsub("Th", "Neuron", tasic_meta$celltypes)
# tasic_meta$celltypes = gsub("Pvalb", "Neuron", tasic_meta$celltypes)
# tasic_meta$celltypes = gsub("Ndnf", "Neuron", tasic_meta$celltypes)
# tasic_meta$celltypes = make.names(tasic_meta$celltypes)
# tasic_features = readRDS(paste0("processed_data/tasic_features_", level, ".rds"))
# tasic_df_full = cell_mean_estimation_pipeline(tasic_gnxp, tasic_meta, tasic_features,
#   p = 0.99, bootstrap = FALSE)
# saveRDS(tasic_df_full, paste0("cell_means/tasic_", level, ".rds"))

#z14 data
# z14_gnxp = readRDS("processed_data/z14_gnxp_matched_gene.rds")
# z14_meta = readRDS("processed_data/z14_meta_matched_gene.rds")
# z14_features = readRDS("processed_data/z14_features_gene.rds")
# z14_features$geneSymbol = z14_features$Symbol
# # z14_gnxp_sub = z14_gnxp[1:2000, ]
# # z14_df = cell_mean_estimation_pipeline(z14_gnxp, z14_meta, z14_features, p = 0.99, bootstrap = TRUE)
# z14_df_full_se = cell_mean_estimation_pipeline(z14_gnxp, z14_meta, z14_features, p = 0.75, bootstrap = FALSE)
# saveRDS(z14_df_full_se, "cell_means/z14_gene.rds")
#
# # #z15 data
# z15_gnxp = readRDS("processed_data/z15_gnxp_matched_gene.rds")
# z15_meta = readRDS("processed_data/z15_meta_matched_gene.rds")
# z15_features = readRDS("processed_data/z15_features_gene.rds")
# # z15_gnxp_sub = z15_gnxp[1:2000, ]
# # z15_df = cell_mean_estimation_pipeline(z15_gnxp, z15_meta, z15_features, p = 0.99, bootstrap = TRUE)
# z15_df_full_se = cell_mean_estimation_pipeline(z15_gnxp, z15_meta, z15_features, p = 0.75, bootstrap = FALSE)
# saveRDS(z15_df_full_se, "cell_means/z15.rds")
