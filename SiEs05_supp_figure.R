mysample#!/bin/R5
#########################
# Supplementary figure ##
#########################

# This script creates the necessary plots for the TFR's supplementary figure
# using SiEs05 libraries


.libPaths('~/R/newer_packs_library/3.5/')
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
library(Seurat)

fname <- "/mnt/BioHome/ciro/simon/info/groupColours_art_ct.csv"
grcols <- read.csv(fname, stringsAsFactors = FALSE, row.names = 1)
grcols["6", ] <- grcols["TFR", ] <- "#ffd92f"
grcols["1", ] <- grcols["TREG", ] <- "#66c2a5"
grcols["common", ] <- "#810f7c"
grcols["follicular", ] <- "orangered2"
grcols["regulatory", ] <- "#6b6b6b"
grcols["Undef", ] <- "#FFFFFF"
grcols["Cluster1_TNFRSF9.neg", ] <- grcols["TREG_TNFRSF9.neg", ]
grcols["Cluster1_TNFRSF9.pos", ] <- grcols["TREG_TNFRSF9.pos", ]
grcols["Cluster6", ] <- grcols["TFR", ]

root = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/clustering_seurat/SiEs05_notcr'
npcs = 17
resolut = 'RNA_snn_res.0.6'
celltypes = c("1" = "TREG", "6" = "TFR", "2" = "CD8", "3" = "CD8", "8" = "CD8", "9" = "CD8")

setwd(paste0(root, '/summary/supp'))
clustname <- paste0(root, '/clustering/zetInfo/clustCells', npcs, 'PCs_30Ks_0.06667JD.RData')
mycells <- theObjectSavedIn(clustname)
mycells@reductions$umap <- mycells@reductions[[ grep("umap", reductions(mycells))[1] ]]

# Renaming
mycells@meta.data$celltypes <- celltypes[as.character(mycells@meta.data[, resolut])]
tvar <- as.character(mycells@meta.data[is.na(mycells@meta.data$celltypes), resolut])
mycells@meta.data[is.na(mycells@meta.data$celltypes), ]$celltypes <- tvar
mycells@meta.data$celltypes <- factormix(as.character(mycells@meta.data$celltypes))
table(mycells@meta.data[, resolut], mycells@meta.data$celltypes)

annot <- mycells@meta.data
void <- add_gene_tag(
  lgenes = 'TNFRSF9',
  annot = annot,
  mat = expm1(mycells@assays$RNA@data) * 100,
  thresh = 1, v = TRUE
)
# annot$tnfrsf9 <- paste0("Cluster", annot$RNA_snn_res.0.6, "_", void[rownames(annot), 'tag_TNFRSF9'])
# annot$tnfrsf9 <- sub("6_.*", "6", annot$tnfrsf9)
annot$tnfrsf9 <- paste0(annot$celltypes, "_", void[rownames(annot), 'tag_TNFRSF9'])
annot$tnfrsf9 <- sub("TFR_.*", "TFR", annot$tnfrsf9)
annot$tnfrsf9 <- sub("\\+", ".pos", annot$tnfrsf9)
annot$tnfrsf9 <- sub("\\-", ".neg", annot$tnfrsf9)
table(annot$tnfrsf9)
mycells@meta.data <- annot

mycols <- v2cols(select = mycells@meta.data$celltypes, sour = grcols, v = TRUE)
p <- DimPlot(object = mycells, reduction = "umap", group.by = "celltypes", cols = mycols, label = TRUE)
foxp3data <- p$data[p$data[, 3] %in% c("TREG", "TFR"), ]
foxp3data$celltypes <- "FOXP3"
p <- p + stat_ellipse(
  data = foxp3data, mapping = aes(x = UMAP_1, y =UMAP_2),
  type = "norm", show.legend = FALSE
)
pdf("clustering.pdf")
print(p)
dev.off()
p <- DimPlot(object = mycells, reduction = "umap", group.by = resolut, label = TRUE)
pdf("clustering_raw_clusters_def_colours.pdf")
print(p)
dev.off()
mycols <- v2cols(select = names(table(mycells@meta.data[, resolut])), v = TRUE)
mycols[names(mycols) %in% rownames(grcols)] <- grcols[names(mycols)[names(mycols) %in% rownames(grcols)], ]
p <- DimPlot(object = mycells, reduction = "umap", group.by = resolut, cols = mycols, label = TRUE)
pdf("clustering_raw_clusters.pdf")
print(p)
dev.off()

dir.create("markers")
genes <- c("FOXP3", "IL1R2", "TNFRSF9", "MKI67", "CD4", "CD8B", "CD8A")
genes <- getfound(genes, rownames(mycells), v = TRUE)
cofs <- range(FetchData(mycells, vars = genes))
void <- lapply(genes, function(g){
  p <- FeaturePlot(object = mycells, features = g,
      min.cutoff = cofs[1], max.cutoff = cofs[2]) +
    NoLegend() + NoAxes()
  # p <- p + stat_ellipse(
  #   data = foxp3data, mapping = aes(x = UMAP_1, y =UMAP_2),
  #   type = "norm", show.legend = FALSE
  # )
  pdf(paste0("markers/marker_gene_", g, ".pdf"))
  print(p)
  dev.off()
  return(p)
})

pdf("markers/marker_grid.pdf")
CombinePlots(void)
dev.off()

#### Heatmap ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- "~/large/simon/results/scdea/SiEs05_notcr/comprs/clusters/6vs1/results_6vs1_mastlog2cpm.csv"
edata <- expm1(mycells@assays$RNA@data)
res <- read.csv(fname, stringsAsFactors = FALSE, row.names = 1)
rownames(res) <- sub("'", "", res$gene_name)
res <- res[order(-res$log2FoldChange), ]
head(res)
pth <- 0.05; fth <- 0.25
genes <- getDEGenes(res, pv = pth, fc = fth, v = TRUE)
annot <- remove.factors(mycells@meta.data); annot <- annot[rev(order(annot$celltypes)), ]
cellnames <- rev(getsubset(c("celltypes", "TREG", "TFR"), annot, v = TRUE))
suffix = "all_cells"

cellnames <- rev(sample_grp(annot[cellnames, ], cname = "celltypes", v = TRUE))
suffix = "sampled_cells"

# using cells as columns
matex <- as.matrix(edata[genes, cellnames] * 100)
matcolms <- data.frame(celltype = annot[cellnames, "celltypes"])

# Using stats as columns
stattype <- "_mean"
suffix <- paste0(suffix, stattype)
tvar <- make_list(annot[cellnames, ], "celltypes", grouping = TRUE)
void <- get_stat_report(
  mat = edata, rnames = genes,
  groups = tvar,
  moments = c("mn", "p"), v = TRUE
)
matex <- as.matrix(void[, grep(stattype, colnames(void))] * 100)
colnames(matex) <- sub(stattype, "", colnames(matex))
matex <- matex[, c('TREG', 'TFR')]
matcolms <- data.frame(celltype = colnames(matex))

headmat(matex)
rownames(matcolms) <- colnames(matex) # column names and type
matrows <- data.frame(celltype = celltypes[as.character(res[rownames(matex), "group"])])
rownames(matrows) <- rownames(matex)

tfhtab <- readfile('/mnt/BioHome/ciro/divya/info/TFH_LIST_FROM_MICHELA_SHANE.txt', stringsAsFactors = FALSE, header = F)
fname <- "~/vdv/Ariel/Simon/bulk/2019-06-06--heatmap/Summary-padj0.05-log2FC1/Supp_table.csv"
ressupp <- readfile(fname, stringsAsFactors = FALSE, row.names = 1)
rownames(ressupp) <- casefold(rownames(ressupp), upper = TRUE)
# head(ressupp)
res$mouse_tag <- ressupp[rownames(res), "group"]
write.csv(res, file = "mast_results.csv")
matrows$box <- sub("\\..*", "", ressupp[rownames(matrows), "group"])
matrows$box <- sub("folicular", "follicular", matrows$box)
matrows$box[is.na(matrows$box)] <- "Undef"
matrows$TFHsig <- "Undef"
matrows[rownames(matrows) %in% tfhtab$V1, ]$TFHsig <- "TFH"
head(matrows); tail(matrows)

# order by
# source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
# gorder <- orderthis(annot = matrows, order_by = "pca", cname = "celltype", mat = t(matex), v = TRUE)
res$celltype <- celltypes[as.character(res$group)]
gorder <- orderthis(annot = matrows, order_by = "log2FoldChange", cname = "celltype", mat = res, v = TRUE)
suffix <- paste0(suffix, "_top50byFC_")
head(res[gorder[[1]], ])
head(res[rev(gorder[[2]]), ])
gorder[[2]] <- head(rev(gorder[[2]]), 50)
gorder[[1]] <- head(gorder[[1]], 50)
table(matrows[gorder[[1]], 1])
matrows <- matrows[unlist(gorder), ]
matex <- matex[unlist(gorder), ]

# matcouls <- list(celltype = v2cols(celltypes, grcols)) # colours for type
matcouls <- lapply(matrows, function(x) v2cols(x, grcols) )
# matcouls[[1]] <- rev(matcouls[[1]])
palettebreaks <- seq(-2, 2, 0.1) # Heatmap setting
mypalette <- colorRampPalette(c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"), space = 'Lab')
# grgenessep <- cumsum(table(topgenes[, "cluster"])) # to separate gene groups

matexz <- t(scale(t(matex)))
matexz[matexz < (-2)] <- -2
matexz[matexz > (2)] <- 2
library(pheatmap)
x <- pheatmap(
  mat               = matexz,
  cluster_rows      = 0,
  cluster_cols      = F,
  color             = mypalette(length(palettebreaks)-1),
  breaks            = palettebreaks,
  scale             = 'none',
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  main              = paste(length(genes), "differetially expressed genes - TFR vs TREG clusters"),
  fontsize          = 10,
  fontsize_col      = 5,
  fontsize_number   = 5,
  annotation_col    = matcolms,
  annotation_row    = matrows,
  annotation_colors = matcouls,
  annotation_legend = T,
  annotation_names_col = F,
  annotation_names_row = F,
  drop_levels       = TRUE,
  gaps_col          = table(matcolms), # gaps per type
  # gaps_row          = grgenessep,
  filename          = paste0('heatmap_', suffix, '_padj', pth, 'fc', fth, '.pdf'),
  width = 10, height = 15
)
try(dev.off(), silent = T)

#### violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- remove.factors(mycells@meta.data)
annot <- annot[getsubset(c("celltypes", "TFR", "TREG"), annot, v = TRUE), ]
expdata <- expm1(mycells@assays$RNA@data) * 100
lgenes <- list(
  TFR = c("CCR8", "IL1R2", "TIGIT", "TNFRSF9", "PDCD1", "LAG3", "TOX", "BATF"),
  TFH_TFR = c("TOX2", "ASCL2", "CXCR5", "BCL6", "IL21", "MAF"),
  TREG_TFR = c("TNFRSF4", "IKZF4", "FOXP3", "CTLA4", "IK2F2", "TNFRSF18", "PRDM1", "IL10"),
  TEFF_TFH_TFR = c("KLF2", "S1PR1", "S1PR4"),
  TEFF_TREG = c("IL7R", "VIM", "CD40LG"),
  defined = c("IL1R2", "TNFRSF9", "MKI67", "TOP2A", "CCR8", "CTLA4", "TNFRSF18",
    "TOX", "PDCD1", "TIGIT", "BATF", "IL10", "TGFB1", "TCF7")
)[6]
# lgenes <- lapply(lgenes, function(x) x[x %in% genes] )
lgenes <- lapply(lgenes, function(x) getfound(x, rownames(expdata), v = TRUE) )
lgenes <- lgenes[sapply(lgenes, length) > 0]
ressupp <- data.frame(gene = unlist(lgenes), stringsAsFactors = FALSE)
ressupp$group <- gsub("[0-9]", "", rownames(ressupp)); rownames(ressupp) <- ressupp$gene
# ressupp <- ressupp[ressupp$gene %in% genes, ]

brkies <- NULL #list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50))
colies <- NULL #c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
extramods <- NULL
mergegroups <- NULL #list(ACT = c("ACT1", "ACT2", "ACT3"))
nc <- 2
chilli <- FALSE
indiv <- TRUE
dir.create("violins")
for(gr in names(lgenes)){
  cat(gr, "\n")
  # system(paste("rm -r", gr))
  dname <- paste0("violins/", gr, ifelse(chilli, "_chilli", ""), ifelse(length(mergegroups), "_2g_", "_"))
  genes <- if(indiv){
    dname <- sub("_$", "/", dname)
    lgenes[[gr]]
  }else{
    lgenes[gr]
  }; if(!dir.exists(dname)) dir.create(dname, recursive = TRUE)
  if(length(genes) == 0) next
  # dir.create(dname)
  dname <- paste0(dname, ifelse(is.list(genes), "grid_", ""))
  if(!is.null(mergegroups)){
    annot$tmp <- ifelse(annot$celltypes %in% mergegroups[[1]], names(mergegroups), annot$celltypes)
    annot$tmp <- factor(annot$tmp, c(gorder[!gorder %in% mergegroups[[1]]], names(mergegroups)))
  }else{
    annot$tmp <- annot$celltypes
  }; annot$tmp <- factor(annot$tmp, level = unique(c(celltypes[celltypes %in% annot$tmp], unique(annot$tmp))))
  for(g in genes){
    cat(" ", g, "\n")
    g <- getfound(g, rownames(expdata), v = TRUE); if(length(g) == 0) next
    fname <- paste0(dname, ifelse(length(g) < 2, g, paste0(length(g), "genes")))
    fname <- paste0(fname, ifelse(is.null(brkies), "_pct100", tail(brkies[[2]], 1)))
    fname <- paste0(fname, ifelse(nc > 1 && length(g) > 1, paste0("_", nc, "col"), ""))
    # fname <- paste0(fname, "_manual")
    # source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
    p <- vlnplot(
      cpmdata = expdata,
      metadata = annot,
      gg = g,
      orderby = "tmp",
      noncero = chilli,
      # couls = c("#fffeee", "#ffe080", "#ffc100", "#ff9a00", "#ff7400", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000"),
      couls = colies,
      brks = brkies,
      vsetting = list(trim = FALSE),
      datatype = expression(bold('Log'[2]*'(CPM + 1)')),
      ncolp = 2,
      log2t = TRUE,
      legendt = 'posiv',
      cuof = 0,
      tags = FALSE,
      rotx = 45,
      return_plot = TRUE,
      v = TRUE
    )
    lenthy <- fitgrid(length(g), nCol = nc)[1] * 4 + 5
    pdf(paste0(fname, ".pdf"), 7, lenthy)
    print(p) #draw_vline(p)
    dev.off()
    pdf(paste0(fname, "blank.pdf"), 7, lenthy)
    print(p + shut_up)
    dev.off()
    if(!is.null(extramods)){
      pdf(paste0(fname, "_extramods.pdf"), 7, lenthy)
      print(p + extramods)
      dev.off()
      pdf(paste0(fname, "_extramodsblank.pdf"), 7, lenthy)
      print(p + extramods + shut_up)
      dev.off()
    }
  }
}

#### GSEA ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- "~/large/simon/results/scdea/SiEs05_notcr/comprs/clusters/6vs1/results_6vs1_mastlog2cpm.csv"
res <- readfile(fname, stringsAsFactors = FALSE)
mgenes <- 'CXCR5|BCL6|DUSP4|CCR8|CTLA4|PDCD1|TIGIT|BATF|IL10|TGFB1|TCF7|TNFRSF9|TNFRSF18|TOX'
rownames(res) <- gsub("'", "", res$gene_name)
res[grepl(mgenes, rownames(res)), 1:5]
# res$log2FoldChange <- -res$log2FoldChange

annot <- mycells@meta.data
annot <- annot[getsubset(c("celltypes", "TFR", "TREG"), annot, v = TRUE), ]
expdata <- expm1(mycells@assays$RNA@data) * 100
tvar <- rev(make_list(annot, "celltypes", grouping = T))
# res <- data.frame(gene_name = rownames(expdata), dummy = expdata[, 1])
# rownames(res) <- as.character(res$gene_name)
source('/mnt/BioHome/ciro/scripts/gsea/gsea_liger.R')
res$s2n <- gsea_metric(
  groups = tvar,
  mat = as.matrix(expdata)[, names(tvar)],
  metric = 'Signal2Noise',
  rnames = rownames(res),
  v = TRUE
)
head(res$s2n)

mymetric <- c('s2n', 'log2FoldChange')[1]
# glist <- readfile("~/simon/info/gsea_lists.csv", stringsAsFactors = F, header = TRUE)
# tmp <- readfile('/mnt/BioHome/ciro/simon/info/Guo_TFR_vs_TREG_DEGs.csv', stringsAsFactors = FALSE, row.names = 1)
# glist$tfr_signature_guo <- ""
# glist[1:nrow(tmp), ]$tfr_signature_guo <- tmp[, 1]
# write.csv(glist, "~/simon/info/gsea_lists_extended.csv", row.names = F, quote = FALSE)
source('/mnt/BioHome/ciro/scripts/gsea/gsea_liger.R')
void <- gsea_liger(
  res = res,
  gene_name = "gene_name",
  lfc.type = mymetric,
  gsea_file = "~/simon/info/gsea_lists_extended.csv",
  method = c("liger", "fgsea")[2],
  myseed = 27,
  path = paste0('gsea_', npcs, 'PCs_', resolut, "/"),
  plot_all = TRUE,
  v = TRUE
)
# gene ontology
dir.create("gene_ontology")
cl <- 6
reson <- read.csv(sub("results_.*", paste0("a1_GOA_mastlog2cpm/", cl, ".csv"), fname), stringsAsFactors = F)
write.table(reson, file = paste0("gene_ontology/cluster", cl, ".txt"), quote = F, sep = "\t")

# how do the matrics correlate
ddf <- res[, c('s2n', 'log2FoldChange')]
head(ddf)
glist <- readfile("~/simon/info/gsea_lists_extended.csv", stringsAsFactors = F, header = TRUE)
fname <- paste0('gsea_', npcs, 'PCs_', resolut, "/", paste0(colnames(ddf)[1:2], collapse = '_'), '_correlation.pdf')
pdf(fname)
for(x in colnames(glist)){
  genes <- glist[[x]]; genes <- genes[!is.na(genes)]; genes <- genes[genes != ""]
  ddf[, x] <- rownames(ddf) %in% genes
  p <- ggplot(ddf, aes_string(x = colnames(ddf)[1], y = colnames(ddf)[2], colour = x)) +
    geom_point(alpha = 0.7) + geom_point(data = ddf[ddf[, x], ], alpha = 0.7) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  print(p)
}; dev.off()

#### Trajectory analysis ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(monocle)
fname <- paste0("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/monocle/clusters1n6/data/DDRTree_mono_15PC_counts_data.RData")
cds <- theObjectSavedIn(fname)
head(cds@phenoData@data)

genes <- c("TOP2A", "MKI67", "IL1R2","FOXP3", "TNFRSF9")
genes <- c("IL1R2", "TNFRSF9", "MKI67", "TOP2A", "CCR8", "CTLA4", "TNFRSF18",
  "TOX", "PDCD1", "TIGIT", "BATF", "IL10", "TGFB1", "TCF7")
genes <- getfound(genes, rownames(cds), v = TRUE)
cds@phenoData@data <- cds@phenoData@data[, !colnames(cds@phenoData@data) %in% genes]
cds@phenoData@data$celltypes <- as.character(mycells@meta.data[rownames(cds@phenoData@data), "celltypes"])
cds@phenoData@data$tnfrsf9 <- as.character(mycells@meta.data[rownames(cds@phenoData@data), "tnfrsf9"])
mycells_wtcr = readRDS('vfinput/mycells_cl1n6_WithTCRTags_Modified.RDS')
percell_clones <- remove.factors(mycells_wtcr@meta.data)
cds@phenoData@data$clonal <- as.character(percell_clones[rownames(cds@phenoData@data), "TCR.tag"])
# tvar <- as.matrix(t(exprs(cds[genes, rownames(cds@phenoData@data)])))
expdata <- expm1(mycells@assays$RNA@data) * 100
edata_genes <- log2(as.matrix(t(expdata[genes, rownames(cds@phenoData@data)])) + 1)
# apply(edata_genes, 2, max)
cds@phenoData@data <- cbind(cds@phenoData@data, edata_genes)
cds@phenoData@varMetadata[colnames(cds@phenoData@data), ]$labelDescription <- NA

mycols <- v2cols(select = cds@phenoData@data$RNA_snn_res.0.6, sour = grcols, v = TRUE)
pdf(paste0('trajectory_clusters1n6.pdf'), height = 12, width = 12)
plot_cell_trajectory(cds, color_by = "RNA_snn_res.0.6", cell_size = 2) + scale_color_manual(values = mycols)
graphics.off()

mycols <- v2cols(select = cds@phenoData@data$tnfrsf9, sour = grcols, v = TRUE)
pdf(paste0('trajectory_celltype_tnfrsf9.pdf'), height = 12, width = 12)
plot_cell_trajectory(cds, color_by = "tnfrsf9", cell_size = 2) + scale_color_manual(values = mycols)
graphics.off()

mycols <- c("Non-expanded" = 'darkslategray3', "Expanded" = 'red')
pdf(paste0('trajectory_celltype_clonal.pdf'), height = 12, width = 12)
plot_cell_trajectory(cds, color_by = "clonal", cell_size = 2) + scale_color_manual(values = mycols)
graphics.off()

dir.create("trajectory_genes")
pp <- lapply(genes, function(x){
  p <- plot_cell_trajectory(cds, color_by = x, cell_size = 1) +
    scale_color_gradient(low = "#BEBEBE", high = "blue") +
    NoAxes() #+ NoLegend()
  pdf(paste0("trajectory_genes/", x, ".pdf"), 8, 8)
  print(p)
  dev.off()
})
pdf(paste0('trajectory_genes.pdf'), height = 16, width = 19)
CombinePlots(pp)
graphics.off()

### Now the tree
library(RColorBrewer)
library(dendextend)
library(scales)

thesecells <- rownames(cds@phenoData@data)

pbmc <- mycells[, thesecells]
annot <- pbmc@meta.data
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 500)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCount_RNA"))
pbmc <- RunPCA(pbmc, npcs = 30)

presuffix <- '_cluster1n6'
pdf(paste0("elbow_plot", suffix,".pdf"))
print(ElbowPlot(pbmc, ndims = 30))
dev.off()

elbows <- get_elbow(x = pbmc@reductions$pca@stdev, threshold = seq(.7, .9, .02), decide = TRUE)
npcs <- 20
suffix <- paste0(presuffix, "_", npcs, "pcs")
## Annotation
pc <- Embeddings(object = pbmc, reduction = "pca")
embed <- pc[, paste0("PC_", 1:npcs)]

met<- "average"
sampleTree2 = hclust(dist((embed)), method = met)
dend2<- as.dendrogram(sampleTree2)

## Dendrogram
dend2<- dend2 %>% set("branches_lwd", 0.03)

print(suffix)
for ( var.de in names(annot)[12] ){
  print(var.de)
  vec.use <- v2cols(annot[,var.de], grcols)
  color.vec <- vec.use[annot[,var.de]]
  pdf(paste0("Dendrogram_", var.de, suffix,".pdf"),14, 5)
  par(mar = c(1.5,3,0.3,0.5), xaxs="i")
  plot(dend2, leaflab="none", axes=F, yaxt= "n")
  axis( side =2, pos=-3.5,labels=TRUE, lwd=1.5)
  colored_bars(colors=color.vec, dend2, rowLabels = c(""), y_scale=2, y_shift=-0.5,bord.col="black", bord.width=2)
  dev.off()
  ## Plotting legend
  pdf( paste0( "Dendrogram_Legend_", var.de, ".pdf"))
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend = names(vec.use), pch=16, pt.cex=3, cex=1.5, bty='n', col = vec.use)
  mtext(var.de, at=0.2, cex=2)
  dev.off()
}


#### TCR labels ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cr_outs <- c(
  "~/large/hayley/raw/NV023/COUNTS/SiEs05_LN_TCR/outs",
  "~/large/hayley/raw/NV023/COUNTS/SiEs05_T_TCR/outs"
)

## Vicente's input ##
dir.create('vfoutput')
dir.create('vfinput')
# aggr table
aggr_tab <- data.frame(
  library_id = basename(dirname(cr_outs)),
  clonotypes = list.files(cr_outs, pattern = 'clonotypes', full.names = TRUE),
  annotations = list.files(cr_outs, pattern = 'filtered_contig_annotations', full.names = TRUE),
  row.names = basename(dirname(cr_outs)),
  stringsAsFactors = FALSE
)
aggr_tab <- aggr_tab[c('SiEs05_T_TCR', 'SiEs05_LN_TCR'), ]
write.csv(aggr_tab, file = "vfinput/vdj_aggr_table.csv", quote = FALSE, row.names = FALSE)
saveRDS(mycells, "vfinput/mycells.RDS")
thesecells <- getsubset(c("celltypes", "TREG", "TFR"), mycells@meta.data, v = TRUE)
mycells_ss <- mycells[, thesecells]
saveRDS(mycells_ss, "vfinput/mycells_cl1n6.RDS")

# Checking what frequency and proportion exactly are
fname <- read.csv("vfinput/vdj_aggr_table.csv", stringsAsFactors = FALSE)[2, 2]
tcrtab <- read.csv(fname, stringsAsFactors = FALSE)
head(tcrtab)
tcrtab[1, 2] / sum(tcrtab[, 2]) # fraction of barcodes with the clonotype

packs <- c('Seurat', 'optparse', 'ggplot2', 'reshape2', 'english', 'VennDiagram', 'stringr')
sapply(packs, require, character.only = TRUE)

command <- paste0("/share/apps/R/3.5/bin/Rscript ",
  "/home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/preliminary_TCR_data_analysis.1.5.1.R ",
  # "/home/ciro/scripts/functions/preliminary_TCR_data_analysis.1.5._snap_GEB72019_3_15PM_CDMXtime.R ",
  "--ReportsPath ", getwd(), "/vfoutput ",
  "--TCRContigs ", getwd(), "/vfinput/filtered_contig_annotations_aggr.csv ",
  "--TCRClonotypes ", getwd(), "/vfinput/clonotypes_aggr.csv ",
  "--SeuratObj ", getwd(), "/vfinput/mycells_cl1n6.RDS ",
  "--Tags \"c('RNA_snn_res.0.6', 'celltypes', 'origlib')\""
)
write.table(command, file = "vfinput/_command", quote = FALSE, row.names = FALSE, col.names = FALSE)
system('sh vfinput/_command')

# setwd('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/clustering_seurat/SiEs05_notcr/summary/supp')
mycells_wtcr = readRDS('vfinput/mycells_cl1n6_WithTCRTags_Modified.RDS')
percell_clones <- remove.factors(mycells_wtcr@meta.data)
sum(table(percell_clones$clonotype.tag) > 1)
table(percell_clones$TCR.tag)
write.csv(percell_clones, file = 'vfoutput/percell_clones.csv')
cellclonotype_tab = read.csv('vfinput/filtered_contig_annotations_aggr.csv', stringsAsFactors = FALSE)
clonotype_tab = read.csv('vfinput/clonotypes_aggr.csv', stringsAsFactors = FALSE)
rownames(clonotype_tab) <- clonotype_tab[, 1]
colnames(clonotype_tab) <- sub('frequency', 'freq_total', colnames(clonotype_tab))
colnames(clonotype_tab) <- sub('proportion', 'prop_total', colnames(clonotype_tab))
head(clonotype_tab)
head(clonotype_tab[, 2]) / sum(clonotype_tab[, 2]) # it's a bit off
subset_clones_tab <- clonotype_tab[clonotype_tab$clonotype_id %in% percell_clones$clonotype.tag, ]
head(subset_clones_tab)
# clone size, tissue (LN, T, both), cluster (each, or both)
percell_clones$clusters <- paste0("Cluster_", percell_clones$RNA_snn_res.0.6)
cnames <- c('orig.tissue', 'celltypes', 'clusters')
final_clonetype_tab <- apply(subset_clones_tab, 1, function(x){
  props <- lapply(percell_clones[which(percell_clones$clonotype.tag %in% x[1]), cnames], table)
  propsdf <- data.frame(t(melt(props)), stringsAsFactors = FALSE)
  colnames(propsdf) <- unlist(lapply(props, names))
  propsdf <- propsdf[2, ]; rownames(propsdf) <- x[1]
  propsdf$clone_size <- sum(percell_clones$clonotype.tag %in% x[1])
  return(propsdf)
})
final_clonetype_tab_sum <- data.frame(rbindlist(final_clonetype_tab, fill = TRUE), stringsAsFactors = FALSE)
rownames(final_clonetype_tab_sum) <- subset_clones_tab[, 1]
final_clonetype_tab_sum[is.na(final_clonetype_tab_sum)] <- 0
head(final_clonetype_tab_sum)
final_clonetype_tab_sum$proportion <- final_clonetype_tab_sum$clone_size / sum(final_clonetype_tab_sum$clone_size)
head(final_clonetype_tab_sum$clone_size) / nrow(final_clonetype_tab_sum)
overlaptreg_tfr <- rowSums(final_clonetype_tab_sum[, 3:4] > 0) > 1
sum(overlaptreg_tfr)
head(final_clonetype_tab_sum[overlaptreg_tfr, ])
final_clonetype_tab_sum <- cbind(final_clonetype_tab_sum, subset_clones_tab[rownames(final_clonetype_tab_sum), -1])
headmat(final_clonetype_tab_sum, 10)
headmat(clonotype_tab[rownames(final_clonetype_tab_sum), ])
write.csv(final_clonetype_tab_sum, file = 'vfoutput/clonotype_table.csv')

#### TRACER plots ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/home/ciro/simon/scripts/tumor_tfr/tracer/run-tracer-ciro.R')

#### Overlaps plots ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(VennDiagram)
mycells_wtcr = readRDS('vfinput/mycells_cl1n6_WithTCRTags_Modified.RDS')
percell_clones <- remove.factors(mycells_wtcr@meta.data)
percell_clones <- percell_clones[!is.na(percell_clones$clonotype.tag), ]
table(percell_clones$TCR.tag)
# fname <- read.csv("vfinput/filtered_contig_annotations_aggr.csv", stringsAsFactors = FALSE)
tagclonlist <- make_list(x = percell_clones, colname = resolut, col_objects = 'clonotype.tag')
tagclonlist <- make_list(x = percell_clones[which(percell_clones$TCR.tag == 'Expanded'), ], colname = resolut, col_objects = 'clonotype.tag')
fname <- paste0('vfoutput/tag_specific_analysis/RNA_snn_res.0.6/clonotypes_sharing/TagByClonotypesSharingVennDiagram.pdf')
pdf(fname)
grid::grid.draw(venn.diagram(
  x = tagclonlist, file = NULL,
  main = expression('Clonotypes sharing between FOXP3'^'+ '*'clusters'),
  fill = grcols[names(tagclonlist), ],
  alpha = rep(0.5, length(tagclonlist)),
  force.unique = TRUE
))
dev.off()

#### TFR integration ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adjust TFR integration colours
fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/integration/seurat_6sets_tfr_500hvg/integrated.RData'
idata <- theObjectSavedIn(fname)
mycols <- v2cols(select = idata@meta.data$orig.set, sour = grcols, v = TRUE)
p <- DimPlot(object = idata, reduction = "umap", group.by = "orig.set", cols = mycols, label = TRUE)
pdf("tfr_integration_colour_set.pdf")
print(p)
dev.off()

#### GEO submission #### -------------------------------------------------------
### 10x data ###
source("/home/ciro/scripts/functions/geo_functions.R")

setwd('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/clustering_seurat/SiEs05_notcr/summary/supp/geo')
sourcedir <- c('~/large/hayley/raw/NV020', '~/large/hayley/raw/NV023')
fnames <- list.files(sourcedir, pattern = 'SiEs05', full.names = TRUE)
fnames <- list.files(fnames, pattern = 'fastq', recursive = TRUE, full.names = TRUE)
fnames <- get_unique_names(fnames)
commands <- merge_fastqs(fnames)
void <- sapply(commands, function(x){
  cat(x, '\n'); system(x)
})

system("md5sum *fastq* >> checksum_processed.txt")

prefix <- 'eschweiler'
# Seurat objects
fnames <- c(
  "~/large/simon/results/clustering_seurat/SiEs05_notcr/clustering/zetInfo/clustCells17PCs_30Ks_0.06667JD.RData"
)
loadRData <- function(fileName) load(fileName); get(ls()[ls() != "fileName"])
fname <- fnames
# for(fname in fnames){
cat(basename(fname), '\n')
mycells <- theObjectSavedIn(fname)
edata <- mycells@assays$RNA@counts
metadata <- remove.factors(mycells@meta.data)
# }
umi <- data.frame(edata, stringsAsFactors = FALSE)
umi <- cbind(gene_name = rownames(umi), umi)
umi <- remove.factors(umi)
write.table(edata, file = paste0(prefix, "_umi_data.txt"), sep = '\t', quote = FALSE)
write.table(metadata, file = paste0(prefix, "_annotation.txt"), sep = '\t', quote = FALSE)
system(paste0('gzip ', prefix, "_umi_data.txt"))
system(paste0('gzip ', prefix, "_annotation.txt"))

mycells_wtcr = readRDS('../vfinput/mycells_cl1n6_WithTCRTags_Modified.RDS')
percell_clones <- remove.factors(mycells_wtcr@meta.data)
write.table(percell_clones, file = paste0(prefix, "_tcr_data.txt"), sep = '\t', quote = FALSE)
system(paste0('gzip ', prefix, "_tcr_data.txt"))

system("md5sum *txt.gz >> checksum_processed.txt")
system("more checksum_processed.txt")

cat("Instrument information\n")
runs <- c(
  '/mnt/NovaSeq/191113_A00475_0141_BHVGMMDSXX_NV020',
  '/mnt/NovaSeq/191216_A00475_0163_AHVGHHDSXX_NV023'
)
ilines <- c(
  "Instrument", "Paired", "Forward", "Reverse", "genome", # core mappping
  "Flowcell", "Read", "Application" # run reports
)
elines <- c(
  "Version", "Health", "Custom", "Planned", "IndexRead", "IsIndexedRead" # run reports
)
ilines <- paste0("-P '", paste0(ilines, collapse = "|"), "'")
elines <- paste0("-Pv '", paste0(elines, collapse = "|"), "'")
system("echo META INFO > instrument.txt")
for(run in runs){
  # SampleSheet.csv
  fname <- paste0(run, "/", c("InputMetadata.csv", "RunParameters.xml", "runParameters.xml"))
  fname <- fname[file.exists(fname)]
  system(paste0("echo '---- ", basename(run), " ----' >> instrument.txt"))
  system(paste0("echo 'Parent folder: ", basename(dirname(run)), "' >> instrument.txt"))
  void <- sapply(fname, function(x){
    system(paste0("echo ", basename(x), " >> instrument.txt"))
    if(!file.exists(x)){ system(paste0("echo 'No Metadata file(s)' >> instrument.txt")); return(NULL) }
    system(paste0("grep ", ilines, " ", x, " | grep ", elines, " | sort -u >> instrument.txt"))
    system(paste0("echo >> instrument.txt"))
    NULL
  })
}
system("cat instrument.txt")
