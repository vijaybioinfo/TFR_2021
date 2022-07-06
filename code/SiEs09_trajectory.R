#!/usr/bin/R

#############
# Monocle 3 #
#############

# This script performs trajectory analysis of single-cell data using Monocle 3
# If it's a Seurat object

# # DelayedMatrixStats mess... they recommend version >=1.8, I have 1.4
# # I can't find vesion 1.8
# devtools::install_github("PeteHaitch/DelayedMatrixStats")
# # This depends on DelayedArray >= 0.15.3, then S4Vectors >= 0.27,
# # and this in turn depends on BiocGenerics >= 0.36.0
# # but these are not availabel on R < 4... Then I realised I needed Bioconductor 3.10!!!!
# cd /home/ciro/R/x86_64-pc-linux-gnu-library/3.6
# wget https://bioconductor.org/packages/3.10/bioc/src/contrib/DelayedMatrixStats_1.8.0.tar.gz
# wget https://bioconductor.org/packages/3.10/bioc/src/contrib/BiocGenerics_0.32.0.tar.gz
# wget https://bioconductor.org/packages/3.10/bioc/src/contrib/S4Vectors_0.24.4.tar.gz
# R CMD INSTALL BiocGenerics_0.32.0.tar.gz
# R CMD INSTALL S4Vectors_0.24.4.tar.gz
# R CMD INSTALL DelayedMatrixStats_1.8.0.tar.gz
# # Trying to install monocle3 on R 4
# # configure: GDAL: 1.11.4
# # checking GDAL version >= 2.0.1... no
# # configure: error: sf is not compatible with GDAL versions below 2.0.1
# conda install -c conda-forge gdal
# # Dind't work:
# # configure: Install failure: compilation and/or linkage problems.
# # configure: error: GDALAllRegister not found in libgdal.
# conda remove gdal
# # from source
# cd ~/bin
# wget https://download.osgeo.org/gdal/2.4.2/gdal-2.4.2.tar.gz
# tar -zvxf gdal-2.4.2.tar.gz
# cd gdal-2.4.2
# sh autogen.sh
# # following https://tldp.org/HOWTO/Software-Building-HOWTO-3.html
# ./configure
# make
# make install

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
rm(.Last)
library(monocle3)
library(Seurat)
is_pristine <- function(x) { !is(x@seed, "DelayedOp") }
sessionInfo()

### Functions
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(
  cds,
  time_bin = 1
){
  time_bin_i <- names(which.max(table(cds[[time_bin]])))
  cell_ids <- which(colData(cds)[, time_bin] == time_bin_i)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  tvar <- as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[tvar]
}

### Input
redu = c("UMAP", "tSNE", "PCA", "LSI", "Aligned")[1]
# CD4T24
edataf = "/home/ciro/large/simon/results/clustering/SiEs09_tfr_xdoublt_xclust4/.object_stem_seurat_mean0.01_pct15.rds"
npcs = 20
cclust = "RNA_snn_res.0.4"
selectss = NULL
hvgf = NULL

### Reading
edata <- theObjectSavedIn(edataf)

### Operations
dname <- paste0(sub('clustering.*', 'trajectory', edataf), "/", basename(dirname(edataf)), "_", npcs, "PCs")
if(!dir.exists(dname)) dir.create(dname)
setwd(dname)
cat("Working at:", getwd(), "\n")
list.files()
scells = getsubset(selectss, edata@meta.data, v = TRUE)
edata <- edata[, scells]
if(!is.null(hvgf)){
  hvgdat <- read.csv(hvgf, stringsAsFactors = FALSE, row.names = 1)
  sum(hvgdat[, 'variable']); tvar <- rownames(hvgdat[hvgdat[, 'variable'], ])
  cat(commas(tvar), '\n')
  VariableFeatures(edata) <- tvar
}

## The first step in working with Monocle 3 is to load up your data into Monocle 3's main class, cell_data_set:
gene_annotation <- edata@assays$RNA@meta.features
gene_annotation <- cbind(gene_short_name = rownames(gene_annotation), gene_annotation)
cds <- new_cell_data_set(
  expression_data = edata@assays$RNA@counts,
  cell_metadata = edata@meta.data,
  gene_metadata = gene_annotation
)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(
  cds = cds,
  method = "PCA",
  num_dim = npcs,
  norm_method = "log",
  use_genes = VariableFeatures(edata),
  residual_model_formula_str = "~nCount_RNA+percent.mt",
  alignment_group = NULL,
  pseudo_count = NULL,
  scaling = TRUE,
  verbose = TRUE
)
cds[[cclust]] <- factormix(cds[[cclust]])

p <- plot_pc_variance_explained(cds)
pdf("1_variance_explained_elbow.pdf")
print(p)
dev.off()
cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings
p <- plot_pc_variance_explained(cds)
pdf("1_variance_explained_elbow_seurat.pdf")
print(p)
dev.off()
cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings[, 1:npcs]

## Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "batch")
cds <- align_cds(
  cds = cds,
  preprocess_method = "PCA",
  alignment_group = NULL,
  alignment_k = 20,
  residual_model_formula_str = "~nCount_RNA+percent.mt",
  verbose = TRUE
)

## Step 3: Reduce the dimensions using UMAP
cds@reducedDims@listData[["UMAP"]] <- edata@reductions[["umap"]]@cell.embeddings
# cds <- reduce_dimension(
#   cds = cds,
#   max_components = 2,
#   reduction_method = redu,
#   preprocess_method = NULL,
#   umap.metric = "cosine",
#   umap.min_dist = 0.1,
#   umap.n_neighbors = 15L,
#   umap.fast_sgd = TRUE,
#   umap.nn_method = "annoy",
#   cores = 1,
#   verbose = TRUE
# )

## Step 4: Cluster the cells
# Monocle is able to learn when cells should be placed in the same trajectory as
# opposed to separate trajectories through its clustering procedure. Recall that
# we run cluster_cells(), each cell is assigned not only to a cluster but also to
# a partition. When you are learning trajectories, each partition will eventually
# become a separate trajectory
cds <- cluster_cells(
  cds = cds,
  reduction_method = redu,
  k = 20,
  cluster_method = c("leiden", "louvain")[1],
  num_iter = 2,
  partition_qval = 0.05,
  weight = FALSE,
  resolution = NULL,
  random_seed = NULL,
  verbose = TRUE
)

## Step 5: Learn a graph
cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  close_loop = TRUE,
  learn_graph_control = NULL,
  verbose = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
pdf("2_partitions.pdf")
print(p)
dev.off()
table(cds@clusters@listData$UMAP$partitions)
if(!'partitions_bk' %in%  names(cds@clusters@listData$UMAP))
  cds@clusters@listData$UMAP$partitions_bk <- cds@clusters@listData$UMAP$partitions
tvar <- factor(rep('1', ncol(cds))); names(tvar) <- names(cds@clusters@listData$UMAP$partitions)
cds@clusters@listData$UMAP$partitions <- tvar

## Step 5: Learn a graph - one partition
cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  close_loop = TRUE,
  learn_graph_control = NULL,
  verbose = TRUE
)
p <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
pdf("2_partitions2.pdf")
print(p)
dev.off()
table(cds@clusters@listData$UMAP$partitions)

pdf("3_select_root.pdf", width = 9, height = 7)
plot_cells(
  cds = cds,
  reduction_method = redu,
  color_cells_by = cclust,
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5
)
dev.off()

## Step 6: Order cells
# root_clust <- names(which.max(table(cds[[cclust]])))
# root_cell <- rownames(colData(cds)[as.character(colData(cds)[[cclust]]) == root_clust, ])[1]
cds <- order_cells(
  cds = cds,
  reduction_method = "UMAP",
  root_pr_nodes = get_earliest_principal_node(cds),
  root_cells = NULL,#root_cell,
  verbose = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
)
pdf("4_pseudotime.pdf")
print(p)
dev.off()

cnames <- grep(paste0(cclust, "|orig"), colnames(colData(cds)), value = TRUE)
cnames <- cnames[sapply(colData(cds)[, cnames], class ) %in% c("character", "factor")]
cnames <- cnames[sapply(colData(cds)[, cnames], function(x) length(table(x)) ) > 1]
for(cname in cnames){
  p <- plot_cells(
    cds = cds,
    reduction_method = redu,
    color_cells_by = cname,
    group_cells_by = 'cluster',
    label_cell_groups = !TRUE,
    label_groups_by_cluster = !TRUE,
    # group_label_size = 6,
    # labels_per_group = 1,
    label_branch_points = TRUE, # black circles
    label_roots = TRUE, # white circles
    label_leaves = !TRUE, # gray circles
    graph_label_size = 3
  )
  pdf(paste0("5_", casefold(redu), "_", cname, ".pdf"), width = 7.5)
  print(p)
  dev.off()
}
graphics.off()

save(cds, file = "object.rdata")
