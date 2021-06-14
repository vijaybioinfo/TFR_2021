#!/usr/bin/R

library(Seurat)
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/clustering_utilities.R")
source("/home/ciro/scripts/handy_functions/devel/utilities.R")
source("/home/ciro/scripts/handy_functions/devel/filters.R")
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
load("/home/ciro/scripts/handy_functions/data/signatures_vijaylab.rdata")
source("/home/ciro/scripts/functions/gene_name_convertion.R")

dir.create("/home/ciro/large/simon/results/signatures/SiEs09_module_score", recursive = TRUE)
mycells = readRDS("/home/ciro/large/simon/results/clustering/SiEs09_tfr_xdoublt_xclust4/.object_stem_seurat_mean0.01_pct15.rds")

signatures_list <- lapply(signatures_vijaylab, function(x) human2mouse(x) ) # takes a long time
str(signatures_vijaylab)
str(signatures_list)
signatures_subset <- clean_feature_list(
 mat = mycells@assays$RNA@data, features = signatures_list, filterby = "p~2", verbose = TRUE
)

mycells <- signature_scoring(
 object = mycells,
 prefix = "/home/ciro/large/simon/results/signatures/SiEs09_module_score/",
 lsignatures = signatures_subset,
 confounders = c("orig.disease", "RNA_snn_res.0.4"),
 verbose = TRUE
)

gsea_results <- gsea_matrix(
  mat = expm1(mycells@assays$RNA@data),
  groups = "RNA_snn_res.0.4",
  metadata = mycells@meta.data,
  metric = "Signal2Noise",
  gsea_list = signatures_subset,
  method = "fgsea",
  path = "/home/ciro/large/simon/results/signatures/SiEs09_gsea/",
  plot_it = TRUE,
  classical_plot = !TRUE,
  verbose = TRUE
)
x <- gsea_plot_summary(
  tests_list = gsea_results,
  # pathways = grep(pattern = "cd8", x = names(signatures_subset), value = TRUE, invert = TRUE),
  path = "/home/ciro/large/simon/results/signatures/SiEs09_gsea/",
  padjthr = 0.05,
  nesthr = 0
)
saveRDS(signatures_subset, file = "/home/ciro/large/simon/results/signatures/gene_lists.rds")
