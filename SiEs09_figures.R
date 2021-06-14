#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2020-12-14
# ---

# This script will create plots for any figure from our Murine data
# Each figure's code will then probably be moved to its correct file

fig_dir <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/figures/SiEs09'
dir.create(fig_dir)
setwd(fig_dir)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
})
theme_set(theme_cowplot())

source('/home/ciro/scripts/handy_functions/devel/file_reading.R')
source('/home/ciro/scripts/handy_functions/devel/filters.R')
source('/home/ciro/scripts/handy_functions/devel/utilities.R')
source('/home/ciro/scripts/handy_functions/devel/plots.R')

### Global objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colours_f = "/home/ciro/scripts/handy_functions/data/colours.csv"
sc_sies09_f = "/home/ciro/large/simon/results/clustering/SiEs09_tfr_xdoublt_xclust4/.object_stem_seurat_mean0.01_pct15.rds"
sc_sies09_clust = 'RNA_snn_res.0.4'
colours <- readfile(colours_f, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
colours <- colours[colours[,1]!="", ]; colours["TFR", 1] <- "#ffd92f"
sc_sies09 <- theObjectSavedIn(sc_sies09_f)

sc_sies09_tcr_f = "/home/ciro/large/simon/results/clustering/SiEs09_tfr_xdoublt_xclust4/.object_stem_seurat_mean0.01_pct15_WithTCRTags_ChangedOn_2020-12-22.RDS"
# sc_sies09_tcr_f = "/home/ciro/large/simon/results/tcr/SiEs09_mod/tcr_metadata_2021-01-08.rds"
sc_sies09_tcr <- readfile(sc_sies09_tcr_f)
# sc_sies09_tcr_vf <- sc_sies09_tcr@meta.data
# all(sc_sies09_tcr_vf$clonotype.tag == sc_sies09_tcr$clonotype.tag, na.rm=T)
sc_sies09_tcr <- sc_sies09_tcr@meta.data

### Violin ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_configs = list(
  list(
    name = "fxx_dots_chilli",
    edataf = NULL,
    features = c("Il1r2", "Ctla4", "Ccr8"),
    type = "violin",
    axis_x = list(name = "RNA_snn_res.0.4"),
    axis_y = 'Seurat Normalized',
    individually = TRUE,
    ncol = 2, size = c(5, 5)
  )
)

fig_data = list(
  metadataf = sc_sies09_f,
  metadatafbk = sc_sies09_f,
  edataf = sc_sies09_f,
  edatafbk = sc_sies09_f,
  odata = sc_sies09
)
for(fig_config in fig_configs){
  # Set data # -----------------------------------------------------------------
  fig_data <- fig_set_data(
    config = fig_replace(config = fig_data, fig_config)
  )
  # str(fig_data, max.level = 1)

  # Plot # ---------------------------------------------------------------------
  pp <- lapply(X = fig_data$features, FUN = function(g){
    cat(g, "\n")

    # Create plot
    p <- violin(
      dat = fig_data$pdata,
      xax = 'Identity',
      yax = g,
      dots = grepl("dots", fig_data$name),
      colour_by = "pct",
      chilli = grepl("chilli", fig_data$name)
    )
    if(isTRUE(fig_data$individually)){
      pdf(paste0(fig_data$name, g, ".pdf"), width = fig_data$size[1], height = fig_data$size[2])
      print(p)
      dev.off()
      pdf(paste0(fig_data$name, g, "_blank.pdf"), width = fig_data$size[1], height = fig_data$size[2])
      print(plot_blank(p))
      dev.off()
    }
    return(p)
  })
}

### Volcano ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/home/ciro/scripts/handy_functions/devel/volcano.R')
source('/home/ciro/scripts/handy_functions/devel/overlap.R')
padjthr = 0.05
fcthr = 1
prefix = ""
volc_configs = list(
  list(
    trimmer = "c400",
    showgenes = c("Klf2", "S1pr1", "Sell", "Foxo1", "Tnfrsf9", "Pdcd1", "Tigit", "Tnfrsf4", "Tnfrsf1b", "Lag3", "Tnfrsf18", "Tox", "Tgfb1"),
    resname = "/home/ciro/large/simon/results/dgea/SiEs09_tfr_xdoublt_xclust4/comprs/global/2vsREST/results_2vsREST_mastlog2cpm.csv",
    fcflip = 1
  )
)

for(config_i in volc_configs){
  cat(basename(config_i$resname), "\n")
  res <- readfile(config_i$resname, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  str(res)
  dtype <- sub("Bmean", "", colnames(res)[grepl("Bmean", colnames(res))])
  grps <- unlist(strsplit(gsub("results_|_mast.*", "", basename(config_i$resname)), "vs"))
  tvar <- paste0(grps, "_mean", dtype)
  res_filt <- data.frame(res[!is.na(res$padj), ], stringsAsFactors = F, check.names = FALSE)
  res_filt$log2FoldChange <- res_filt$log2FoldChange * config_i$fcflip
  res_filt$gene <- gsub("'", "", res_filt$gene)
  res_filt$Mean <- round(log2(ifelse(res_filt$group == grps[1], res_filt[, tvar[1]], res_filt[, tvar[2]]) + 1), 1)
  res_filt$pcts <- res_filt$pct_diff
  genes2plot <- mysignames <- getDEGenes(res_filt, pv = padjthr, fc = fcthr, gene_name = "gene", further = NULL, v = TRUE)
  # tvar <- rownames(res_filt)[res_filt[, grep("minExp", colnames(res_filt), value = TRUE)]]
  # genes2plot <- genes2plot[genes2plot %in% tvar]
  res_filt$degs <- "Not_significant"
  res_filt[res_filt$gene %in% genes2plot, ]$degs <- "DEG"
  res_filt[!res_filt$gene %in% genes2plot, ]$Mean <- NA
  res_filt[!res_filt$gene %in% genes2plot, ]$pcts <- 0
  tvar <- cosmy(genes2plot, patties = "^rps|^rpl|^mt-|rp[0-9]{1,}-|^linc")
  showgenes <- if(is.null(config_i$showgenes) || any(config_i$showgenes %in% "add")){
    c(config_i$showgenes, bordering(res_filt[tvar, ], cnames = "log2FoldChange", n = 10))
  }else if(is.null(config_i$showgenes)){ FALSE }else{ config_i$showgenes }
  showgenes <- show_found(config_i$showgenes, rownames(res_filt), v = TRUE)
  tvar <- -log10(res_filt$padj); tvar[is.infinite(tvar)] <- max(tvar[is.finite(tvar)]); summary(tvar)
  trimit <- ifelse(isTRUE(max(tvar) < as.numeric(gsub("c", "", config_i$trimmer))), paste0("c", round(max(tvar))), config_i$trimmer)
  void <- volplot(
    res_filt,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = 'padj',
    lfctype = 'log2FoldChange',
    col_feature = "Mean",
    size_feature = "pcts",
    gene_name = 'gene',
    check_genes = list(text = parse_ens_name(config_i$showgenes)),
    titl = paste0("'", paste0(gsub("_", " ", grps), collapse = "' vs '"), "'"),
    return_plot = TRUE,
    clipp = trimit,
    v = TRUE
  ) + labs(size = "%", title = NULL)
  fname <- paste0(prefix, "volcano_", basename(dirname(config_i$resname)), "_fc", fcthr, "_trim", trimit)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(void)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
  print(plot_blank(void))
  dev.off()
}

#### Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_id = "proportion"
myvar = c("orig.donor", "orig.time_days")
myvar = c("orig.ht_id", sc_sies09_clust)
dfplot <- FetchData(
  object = sc_sies09,
  vars = c(colnames(sc_sies09@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Group <- factormix(dfplot[, myvar[2]])
table(dfplot[, myvar])
fname <- paste0(c(results_id, gsub("orig\\.|orig", "", myvar)), collapse = "_")

scells <- sample_grp(annot = dfplot, cname = myvar[1], v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Group)) +
  geom_point() +
  facet_wrap(facets = as.formula(paste0("~", myvar[1]))) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf(paste0(fname, "_umap.pdf"))
print(p)
dev.off()
pdf(paste0(fname, "_umap_blank.pdf"))
print(plot_blank(p))
dev.off()

pp <- plot_pct(
  x = dfplot, groups = myvar,
  normalise = myvar[2], return_table = TRUE,
  print_ptables = TRUE, v = T
)
propdf <- pp$table
write.csv(propdf, file = paste0(fname, ".csv"))

stupsize <- c(5, 5, 5, 5, 5, 0)
names(stupsize) <- c("BottleRocket1", "Rushmore1", "Darjeeling1", "FantasticFox1", "Zissou1", "custom")
for(i in 3){
  if(names(stupsize)[i] != "custom"){
    couls <- colorRampPalette(wesanderson::wes_palette(name = names(stupsize)[i], n = stupsize[[i]], type = "discrete"))(stupsize[[i]])
    if(names(stupsize)[i] %in% c("FantasticFox1", "Rushmore1", "Zissou1")) couls <- rev(couls)
  }else{
    couls <- v2cols(levels(pp$plot$data$fill_group), grcols)
  }
  ddf = reshape2::melt(propdf[grep("Percentage", propdf$Type), -2])
  ddf$value <- round(ddf$value * 100, 1)
  p2 <- ggplot(data = ddf, mapping = aes(x = variable, y = value, color = variable)) +
    geom_boxplot(size = 3) +
    geom_jitter(shape=16, position=position_jitter(0.2),color = "black", size = 3) +
    scale_color_manual(values = couls) + labs(x = NULL, y = "Percentage") +
    theme(legend.position = "none")
  pdf(paste0(fname, "_boxplot_", names(stupsize)[i], ".pdf"), width = 5, height = 5)
  print(p2)
  dev.off()
  pdf(paste0(fname, "_boxplot_", names(stupsize)[i], "_blank.pdf"), width = 5, height = 5)
  print(plot_blank(p2))
  dev.off()
}

pct_donut <- plot_pct(
  x = dfplot, groups = rev(myvar),
  normalise = FALSE,
  print_ptables = TRUE,
  type = "donut",
  return_table = TRUE
)
propdf <- pct_donut$table
write.csv(propdf, file = paste0(fname, ".csv"))
p <- pct_donut$plot +
  labs(fill = NULL) +
  theme_minimal() + theme(
    legend.position = "right",
    strip.text = element_text(face = 'bold', size = 10),
    axis.text.x = element_blank()
  )

pdf(paste0(fname, "_donut.pdf"), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0(fname, "_donut_blank.pdf"), width = 10, height = 10)
print(shut_it(p))
dev.off()

#### Gene stats ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/home/ciro/scripts/handy_functions/R/stats_summary_table.R')
results_id = "stats"
group_i = c("orig.time_days", "orig.donor")
genes = "Il1r2"
sc_sies09$group <- do.call(paste, c(sc_sies09[[group_i]], sep = "_"))
table(sc_sies09$group)
sstat <- data.frame(t(stats_summary_table(
  mat = sc_sies09@assays$RNA@data,
  groups = make_list(sc_sies09@meta.data, colname = "group", grouping = TRUE),
  rnames = genes,
  moments = "p",
  verbose = TRUE
)))
propdf <- data.frame(
  cbind(t(sapply(stringr::str_split(rownames(sstat), "_"), c)), sstat),
  row.names = NULL
)
colnames(propdf)[1:3] <- c("Day", "Donor", "Moment")
fname <- paste0(c(results_id, gsub("orig.", "", group_i)), collapse = "_")
write.csv(propdf, file = paste0(fname, ".csv"))
for(gene in genes){}
p <- ggplot(data = propdf, mapping = aes_string(x = "Day", y = gene, color = "Day")) +
  geom_boxplot(size = 1) +
  geom_jitter(shape=16, position=position_jitter(0.2),color = "black", size = 3) +
  # scale_color_manual(values = couls) +
  labs(x = NULL, y = paste0("Percentage:", gene)) +
  theme(legend.position = "none")
pdf(paste0(fname, "_", gene, "_boxplot.pdf"), width = 5, height = 5)
print(p)
dev.off()

#### GSEA #### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "signatures_sc_"
signatures_files = "/mnt/BioAdHoc/Groups/vd-vijay/Ariel/Simon/bulk/2019-06-06--heatmap/Summary-padj0.05-log2FC1/Supp_table.csv"
globalsign = list()

tvar <- read.csv(signatures_files, stringsAsFactors = FALSE)
tvar <- make_list(tvar, colname = "group", col_objects = "X")
tvar <- sapply(tvar, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
tvar <- tvar[sapply(tvar, length) > 0]
str(tvar)
signatures_list <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
str(signatures_list)
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/R/clustering_utilities.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
signatures_subset <- clean_feature_list(
 mat = sc_sies09@assays$RNA@data, features = signatures_list, filterby = "p~2", verbose = TRUE
)

# sc_sies09 <- signature_scoring(
#  object = sc_sies09,
#  prefix = paste0(result_id, "modulescore/"),
#  lsignatures = signatures_subset,
#  confounders = sc_sies09_clust,
#  verbose = TRUE
# )

gsea_results <- gsea_matrix(
  mat = expm1(sc_sies09@assays$RNA@data),
  groups = sc_sies09_clust,
  metadata = sc_sies09@meta.data,
  metric = "Signal2Noise",
  gsea_list = signatures_subset,
  method = "fgsea",
  path = paste0(result_id, "gsea/"),
  plot_it = TRUE,
  classical_plot = !TRUE,
  verbose = TRUE
)
x <- gsea_plot_summary(
  tests_list = gsea_results,
  # pathways = grep(pattern = "tcr", x = names(signatures_subset), value = TRUE, invert = TRUE),
  path = paste0(result_id,"gsea/"),
  padjthr = 0.05,
  nesthr = 0
)

### Dot-plot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "dotplot"
genes = c(
  "Klf2", "S1pr1", "Sell", "Foxo1", "Tcf7", "Mki67", "Top2a", "Tnfrsf9", "Pdcd1",
  "Tigit", "Tnfrsf4", "Tnfrsf1b", "Lag3", "Tnfrsf18", "Tox", "Tgfb1"
)

mygenes <- show_found(genes, rownames(sc_sies09), v = TRUE)
sc_sies09@meta.data$Identity <- factor(sc_sies09@meta.data[, sc_sies09_clust], c("0", "1", "2", "3", "5"))

dohclust <- grepl("hclust", result_id)
if(dohclust){
  suppressPackageStartupMessages({
    library("ggplot2")
    library("ggdendro")
  })
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  summarydf <- stats_summary_table(
    mat = sc_sies09@assays$RNA@data,
    group = make_list(sc_sies09@meta.data, colname = "Identity", grouping = TRUE),
    rnames = mygenes,
    moments = "mn"
  )
  hc <- hclust(dist(scale(summarydf)))
  mygenes <- hc$labels[hc$order]
  # dd <- reorder(as.dendrogram(hc), rev(hc$order))
  hcp <- ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE) +
    scale_y_reverse() + scale_x_reverse() + theme_void()
}

p <- DotPlot(
  object = sc_sies09,
  features = mygenes,
  group.by = "Identity",
  cols = c('#fff4ba', '#ff0000'),
  col.min = -1.5, col.max = 1.5
) + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)

fname <- paste0(result_id, "")
pdf(paste0(fname, ".pdf"), width = ifelse(dohclust, 15, 10), height = 10)
if(dohclust){ print(hcp | p) }else{ print(p) }
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 8, height = 10)
print(plot_blank(p))
graphics.off()

### Disease dimentional reduction plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
redu = list(umap = c("UMAP_1", "UMAP_2"))
dir.create("umap")
result_id = "umap/"
sselect = list(c("orig.time_days", "d11"))
sselect = list(c("orig.time_days", "d18"))
sselect = list(c("orig.time_days", "d18", "d11"))

p <- DimPlot(
  object = sc_sies09,
  cells = filters_subset_df(sselect, sc_sies09@meta.data, v = TRUE),
  reduction = names(redu),
  group.by = sc_sies09_clust,
  pt.size = 1.6
) + labs(x = redu[[1]][1], y = redu[[1]][2])

fname <- paste0(result_id, sub("orig.", "", filters_summary(sselect)))
pdf(paste0(fname, ".pdf"))
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"))
print(plot_blank(p))
graphics.off()

#### Cluster markers ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
result_id = "heatmap_markers_shared_mean"
ident_order = c("0", "1", "2", "3", "5")
identnames = paste0("CL", ident_order)
names(identnames) <- ident_order
ntop = Inf
marknames = "/home/ciro/large/simon/results/clustering/SiEs09_tfr_xdoublt_xclust4/seurat_mean0.01_pct15_pc20_res0.4/dgea_MAST_fc0.25_padj0.05_summary_stats.csv"
fc = 1

sc_sies09$Name = factor(identnames[as.character(sc_sies09@meta.data[, sc_sies09_clust])], levels = unique(identnames))
sc_sies09$Cluster = factor(sc_sies09@meta.data[, sc_sies09_clust], levels = names(identnames))
mygenes <- readfile(marknames, stringsAsFactors = FALSE, row.names = 1)

# upset plot
library(UpSetR)
mygenesl <- make_list(mygenes, colname = "cluster", col_objects = "gene")
str(mygenesl)
# overlap_calc(mygenesl[c(1, 3, 9)])
uddf <- as.data.frame.matrix(table(mygenes[, c("gene", "cluster")]))
uddf <- uddf[, names(identnames)]
colnames(uddf) <- paste0("C", colnames(uddf), " ", identnames[colnames(uddf)])
str(uddf)
pdf(paste0(result_id, "_upset.pdf"), onefile = TRUE)
upset(data = uddf, sets = colnames(uddf))
dev.off()

mygenes$gene <- gsub("'", "", mygenes$gene_name)
tvar <- mygenes$Dpct > .0; table(tvar)
mygenes <- mygenes[tvar, ]
tvar <- mygenes$avg_logFC > fc; table(tvar)
mygenes <- mygenes[tvar, ]
# tvar <- !grepl("&", mygenes$sCluster); table(tvar)
# mygenes <- mygenes[tvar, ]
mygenes$cluster <- factor(mygenes$cluster, names(identnames))
mygenes <- mygenes[order(mygenes$cluster), ]
# mygenes$nclust <- stringr::str_count(mygenes$sCluster, "&")
topgenes <- get_top_n(x = mygenes, n = ntop)
fname <- paste0(result_id, ifelse(nrow(topgenes) != nrow(mygenes), paste0("_top", ntop), ""))
# topgenes <- get_top_n(x = topgenes, n = ntop, orderby = 'nclust'); fname <- paste0(fname, "_nclustOrder")
genes <- gsub("'", "", topgenes$gene_name)
genes <- show_found(genes, rownames(sc_sies09), v = TRUE)
genesl <- make_list(x = topgenes, colname = "cluster", col_objects = "gene_name")

pdf(paste0(fname, "_fc", fc,".pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = sc_sies09,
  rnames = genes,
  orderby = "Cluster",
  use_mean = "Cluster",
  sample_it = c(cname = "Cluster", maxln = "-1000"),
  scale_row = TRUE,
  categorical_col = c("Cluster", "Name"),
  feature_order = TRUE,
  couls = NULL,
  hcouls = rev(c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c")),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE
)
graphics.off()

### Trajecroy figures ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "trajectory_"
sc_sies09_trj_f = "/home/ciro/large/simon/results/trajectory/SiEs09_tfr_xdoublt_xclust4_20PCs/object.rdata"
sc_sies09_trj <- theObjectSavedIn(sc_sies09_trj_f)

time_bin = 1
time_bin_i <- names(which.max(table(sc_sies09_trj[[time_bin]])))
cell_ids <- which(colData(sc_sies09_trj)[, time_bin] == time_bin_i)
closest_vertex <- sc_sies09_trj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(sc_sies09_trj), ])
tvar <- as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
root_pr <- igraph::V(principal_graph(sc_sies09_trj)[["UMAP"]])$name[tvar]

cds <- order_cells(
  cds = sc_sies09_trj,
  reduction_method = "UMAP",
  root_pr_nodes = root_pr,
  root_cells = NULL,#root_cell,
  verbose = TRUE
)
# Each leaf, denoted by light gray circles, corresponds to a different outcome
# (i.e. cell fate) of the trajectory. Black circles indicate branch nodes,
# in which cells can travel to one of several outcomes.
for(var2plot in c(sc_sies09_clust, "pseudotime")){
  p <- plot_cells(
    cds = cds,
    reduction_method = "UMAP",
    color_cells_by = var2plot,
    group_cells_by = "cluster",
    label_cell_groups = !TRUE,
    label_groups_by_cluster = !TRUE,
    cell_size = 1.6, cell_stroke = 0,
    show_trajectory_graph = var2plot == sc_sies09_clust,
    label_branch_points = var2plot == sc_sies09_clust, # black circles
    label_roots = var2plot == sc_sies09_clust, # white circles
    label_leaves = !TRUE, # gray circles
    graph_label_size = 3
  ) + theme(legend.title=element_blank())

  fname <- paste0(result_id, "umap_", var2plot)
  pdf(paste0(fname, ".pdf"))
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"))
  print(plot_blank(p + theme_cowplot()))
  dev.off()
}

### Clone overlap ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/home/ciro/scripts/handy_functions/devel/overlap.R')
library(UpSetR)
library(eulerr)
result_id = "clone_overlap"

# upset plot
sc_sies09_tcr$Clusters <- ifelse(sc_sies09_tcr[, sc_sies09_clust] == "2", "TFR", "REST")
clones_list <- make_list(sc_sies09_tcr, colname = "Clusters", col_objects = "clonotype.tag")
clones_list <- lapply(X = clones_list, FUN = function(x) x[!is.na(x)] )
clones_list_ov <- overlap_list(clones_list)
str(clones_list_ov)
euler_i <- clones_list_ov[sapply(clones_list_ov, length) > 0]
names(euler_i) <- gsub("n", "&", names(euler_i))
euler_i <- sapply(euler_i, length)

fit2 <- euler(uddf)
pdf(paste0(result_id, "_euler.pdf"))
plot(
  fit2,
  quantities = TRUE,
  # fill = "transparent",
  lty = 1:3,
  labels = list(font = 4)
)
dev.off()

uddf <- as.data.frame.matrix(table(sc_sies09_tcr[, c("clonotype.tag", sc_sies09_clust)]))
uddf <- as.data.frame.matrix(table(sc_sies09_tcr[, c("clonotype.tag", "Clusters")]))
colnames(uddf) <- paste0("C", colnames(uddf))
str(uddf)
head(uddf)
uddf[uddf > 0] <- 1

pdf(paste0(result_id, "_upset.pdf"))
upset(data = uddf, sets = colnames(uddf))
dev.off()

### Clone overlap ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "clone_size"
couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')

dfplot <- FetchData(
  object = sc_sies09,
  vars = c(sc_sies09_clust, c('UMAP_1', 'UMAP_2'))
)
dfplot$`Clone size` <- sc_sies09_tcr[rownames(dfplot), ]$clon.size.tag
hb_limit <- quantile(dfplot$`Clone size`, probs = .95, na.rm = T)

p <- ggplot(
  data = dfplot,
  aes(x = UMAP_1, y = UMAP_2, color = `Clone size`)
) +
  geom_point(alpha = 0.7) +
  scale_colour_gradientn(colours = couls, limits = c(1, hb_limit))

pdf(paste0(result_id, "_umap.pdf"))
print(p)
dev.off()
pdf(paste0(result_id, "_umap_blank.pdf"))
print(plot_blank(p))
dev.off()

dfplot$clone_size <- dfplot$`Clone size`
dfplot$clone_size[dfplot$clone_size > hb_limit[[1]]] <- hb_limit[[1]]
dfplot$cluster <- ifelse(dfplot[, sc_sies09_clust] == "2", "TFR", "TREG")
p <- ggplot(
    data = dfplot[which(dfplot$`Clone size` > 1), ],
    mapping = aes_string(x = "cluster", y = "clone_size", fill = "cluster")
  ) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  geom_jitter(width = 0.3, size = 0.4, color = "#595959", alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(values = v2cols(dfplot[, "cluster"], colours)) +
  labs(x = "Clusters", y = "Clone Size") +#expression("Log"[2]*"(Clone Size + 1)")) +
  theme(legend.position = "none")

pdf(paste0(result_id, "_boxplot.pdf"))
print(p)
dev.off()

### TCR overlap ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/home/ciro/scripts/handy_functions/devel/overlap.R')
library(UpSetR)
library(eulerr)
result_id = "clone_overlap_guo_"

metadata <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_annotation_tracer.RData')
metadata <- remove.factors(metadata)
tpm <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_tpm.RData')
str(metadata)
summary(metadata[metadata$BCL6.expr == "BCL6.neg", "BCL6.tpm"])
metadata$celltype <- "CD4"
summary(metadata$CD4.tpm)
table(metadata$orig.celltype)
tvar <- metadata$FOXP3.tpm>10
metadata$celltype[tvar] <- "TREG"
tvar <- (metadata$CXCR5.tpm>10 | metadata$BCL6.tpm>10) & metadata$FOXP3.tpm>10
metadata$celltype[tvar] <- "TFR"
tvar <- (metadata$CXCR5.tpm>10 | metadata$BCL6.tpm>10) & metadata$FOXP3.tpm<=10
metadata$celltype[tvar] <- "TFH"
table(metadata$celltype)

# upset plot
clones_df <- metadata[metadata$GreaterTHAN1 > 1 & metadata$celltype %in% c("TFR", "TFH"), ]
clones_df <- metadata[metadata$GreaterTHAN1 > 1, ]
clones_list <- make_list(clones_df, colname = "celltype", col_objects = "Name.TRA.and.TRB")
clones_list <- lapply(X = clones_list, FUN = function(x) x[!is.na(x)] )
suffix <- paste0(names(clones_list), collapse = "-")
str(clones_list)
euler_i <- overlap_list(clones_list)
str(euler_i)
euler_i <- euler_i[sapply(euler_i, length) > 0]
names(euler_i) <- gsub("n", "&", names(euler_i))
euler_i <- sapply(euler_i, length)

fit2 <- euler(euler_i)
pdf(paste0(result_id, suffix, "_euler.pdf"))
plot(
  fit2,
  quantities = TRUE,
  # fill = "transparent",
  lty = 1:3,
  labels = list(font = 4)
)
dev.off()

fit3 <- euler(euler_i, shape = "ellipse")
plot(fit3)
pdf(paste0(result_id, suffix, "_euler_ellipse.pdf"))
plot(fit3, quantities = TRUE, lty = 1:3, labels = list(font = 4))
dev.off()

uddf <- as.data.frame.matrix(table(clones_df[, c("Name.TRA.and.TRB", "celltype")]))
colnames(uddf) <- paste0("C", colnames(uddf))
str(uddf)
head(uddf)
uddf[uddf > 0] <- 1

pdf(paste0(result_id, suffix, "_upset.pdf"))
upset(data = uddf, sets = colnames(uddf))
dev.off()
