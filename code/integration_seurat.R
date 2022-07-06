#!/usr/bin/R5

######################
# Seurat integration #
######################

# This is going to be a code to analyze single cell data with Seurat alignment

### Installing ### ---
# devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
.libPaths('~/R/newer_packs_library/3.5/')
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

deps <- c('Seurat', 'ggplot2', 'cowplot')
load_packs(deps, v = T)

root <- '~/large/simon/results/integration'
setwdc(root)

annot_list_tags <- theObjectSavedIn(paste0(root, '/data/annot_list_tags_10TPM.RData'))
edata_list <- theObjectSavedIn(paste0(root, '/data/edata_list.RData'))
names(edata_list) # we will take Lambrechts out
edata_list <- edata_list[-5]
annot_list_tags <- annot_list_tags[-5]

# Create Seurat object per data set
pancreas.list <- lapply(names(annot_list_tags), function(x){
  CreateSeuratObject(counts = edata_list[[x]], meta.data = annot_list_tags[[x]])
}); names(pancreas.list) <- names(annot_list_tags)

sets <- names(pancreas.list)#[c(1:4)]
setwdc(paste0(root, '/seurat_', length(sets), "sets"))

dir.create('qcs')

for (i in 1:length(x = pancreas.list)) {
  cat(names(pancreas.list)[i], '\n')
  pancreas.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pancreas.list[[i]], pattern = "^MT-")
  plot1 <- FeatureScatter(pancreas.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pancreas.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pdf(paste0("qcs/", names(pancreas.list)[i], ".pdf"), 12, 7)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()

  # thesecells <- rownames(pancreas.list[[i]]@meta.data[pancreas.list[[i]]@meta.data[, "nFeature_RNA"] < tvar, ])
  # plot1 <- FeatureScatter(pancreas.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt", cells = thesecells)
  # plot2 <- FeatureScatter(pancreas.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = thesecells)
  # pdf(paste0("qcs/", names(pancreas.list)[i], "_filtered.pdf"), 12, 7)
  # print(CombinePlots(plots = list(plot1, plot2)))
  # dev.off()
  # pancreas.list[[i]] <- subset(pancreas.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
}

summ_filters <- data.frame(rbindlist(lapply(pancreas.list, function(x){
  mytab <- t(data.frame("X" = c(range(x@meta.data[, "percent.mt"]),
  mean(x@meta.data[, "percent.mt"]),
  range(x@meta.data[, "nFeature_RNA"]),
  quantile(x@meta.data[, "nFeature_RNA"], prob = 0.998))))
  colnames(mytab) <- c('MinMTpct', 'MaxMTpct','MeanMTpct' , 'MinFeat', 'MaxFeat', 'Q99.8%')
  data.frame(mytab)
})))
rownames(summ_filters) <- names(pancreas.list)
summ_filters

# columns for visualisation
orignames <- c("orig.set", "orig.majorCluster", 'tag_FOXP3', 'tag_ct')
npcs <- 30

if(file.exists('integrated.RData')) cat('Go to file creating line\n')

for (i in 1:length(x = pancreas.list)) {
  cat(names(pancreas.list)[i], '\n')
  pancreas.list[[i]] <- NormalizeData(object = pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(object = pancreas.list[[i]],
      selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  cat(commas(VariableFeatures(object = pancreas.list[[i]]), 10), '\n')
}
# Check most variable genes overlaps
myvargenes <- unique(unlist(lapply(pancreas.list, VariableFeatures)))
ogenes <- sapply(pancreas.list, function(x) myvargenes %in% VariableFeatures(x) )
rownames(ogenes) <- myvargenes
head(ogenes, 20)
myvargenes <- myvargenes[apply(ogenes, 1, all)]
length(myvargenes)

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list[sets], dims = 1:npcs)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:npcs)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(object = pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = npcs, verbose = FALSE)

pdf(paste0('sdevPCs_', npcs,'PCs.pdf'), width = 10, height = 8)
ElbowPlot(object = pancreas.integrated, ndims = npcs)
graphics.off()

pc_sdev <- pancreas.integrated@reductions$pca@stdev
get_elbow(1:length(pc_sdev), pc_sdev, seq(95, 70, by = -5)/100)

chnpcs <- 15
nres <- 0.2
redu <- "umap"

setwdc(paste0('~/large/simon/results/integration/seurat_', length(sets), "sets/PC", chnpcs, 'R', nres))

if(redu == "umap"){
  cat("Runnning UMAP\n")
  pancreas.integrated <- RunUMAP(object = pancreas.integrated, reduction = "pca", dims = 1:chnpcs, min.dist = 0.05, spread = 2)
}else{
  cat("Runnning t-SNE\n")
  redu <- "tsne"
  pancreas.integrated <- RunTSNE(object = pancreas.integrated, reduction = "pca", dims = 1:chnpcs, check_duplicates = FALSE,
    tsne.method = "FIt-SNE", fast_tsne_path = '/mnt/BioHome/ciro/bin/FIt-SNE2/bin/fast_tsne')
}
for(orig in orignames){
  tvar <- length(unique(pancreas.integrated@meta.data[, orig]))
  pdf(paste0('integrated_', sub('orig.', '', orig), '.pdf'), height = 10, width = ifelse(tvar > 15, 14, 10))
  print(DimPlot(object = pancreas.integrated, reduction = redu, group.by = orig))
  graphics.off()
}
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:chnpcs)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = nres)
tailmat(pancreas.integrated[[]], 10)
gby <- paste0("integrated_snn_res.", nres)
pancreas.integrated@meta.data[, gby] <- as.character(pancreas.integrated@meta.data[, gby])

pdf(paste0('clusters_', gby, '.pdf'), height = 8, width = 8)
DimPlot(object = pancreas.integrated, reduction = redu, group.by = gby)
graphics.off()

freq_tablep(metadata = pancreas.integrated@meta.data, cnames = c(gby, 'orig.set'),
  pnames = c('Clusters in sets', 'Sets in clusters'), dowrite = TRUE)

markers <- read.csv('/mnt/BioHome/ciro/simon/info/markers.csv', stringsAsFactors = F)
markers
mymarkers <- unique(c('FOXP3', markers[, 1], 'IL2RA', 'TNFRSF9', 'TNFRSF18', 'DUSP4',
  'CCR8', 'IL1R2', 'IKZF2', 'ENTPD1', 'LAG3', 'TIGIT', 'CTLA4', 'PDCD1','TOX'))[1]
mymarkers <- getfound(mymarkers, rownames(pancreas.integrated@assays$RNA), v = T)

fname <- paste0('markers_', redu, '_vln.pdf')
pdf(fname, height = 7.5, width = 15)
for(i in 1:length(mymarkers)){
  print(plot_grid(FeaturePlot(pancreas.integrated, features = mymarkers[i], min.cutoff = 0),
    VlnPlot(pancreas.integrated, mymarkers[i], group.by = gby, assay = "RNA") + NoLegend()))
}
# vlnplot(pancreas.integrated, gg = mymarkers[i], orderby = gby, plotdots = T, noncero = T, v = T)
graphics.off()

pdf(paste0('markers_', redu, '.pdf'), height = 8, width = 8)
for(i in 1:length(mymarkers)){
  print(FeaturePlot(pancreas.integrated, features = mymarkers[i], min.cutoff = 0))
}
graphics.off()

pdf(paste0('markers_', redu, '_split.pdf'), height = 5, width = 25)
for(i in 1:length(mymarkers)){
  print(FeaturePlot(pancreas.integrated, features = mymarkers[i], min.cutoff = 0, split.by = 'orig.set'))
}
graphics.off()

DefaultAssay(object = pancreas.integrated)
dim(pancreas.integrated@assays$integrated@counts)
dim(pancreas.integrated@assays$integrated@data)

save(pancreas.integrated, file = '../integrated.RData')
load('integrated.RData') #### --------------------------------------------------

# pdf('combined_genes.pdf', 16, 5)
# FeaturePlot(pancreas.integrated, features = c('BCL6', 'CXCR5'), blend = T)
# graphics.off()

## BCL6+, CXCR5+, BCL6+CXCR5+ cells within the 4 FoxP3+ clusters
metadata <- pancreas.integrated[[]]#[, c('orig.fcmarkers', gby)]
metadata <- remove.factors(metadata)
colnames(metadata) <- sub(gby, 'Cluster', colnames(metadata))
# void <- theObjectSavedIn('../../data/metadata.RData')
# metadata <- cbind_repcol(void[getfound(rownames(metadata), rownames(void), v = T), ], metadata)
# metadata <- metadata[getsubset(c('Cluster', '1', '4', '5', '7'), metadata, v = T), ]
metadata <- metadata[, sapply(metadata, function(x) all(!is.na(x)) && length(table(x)) < 100 ) ]
metadata <- metadata[, getpats(colnames(metadata), c('tag', 'orig.fc', 'Cluster'), 'major')]
head(metadata)
sapply(metadata, table)
tvar <- sapply(head(colnames(metadata), -1), function(x) table(metadata[, 'Cluster'], metadata[, x]) )
tvar <- t(do.call(cbind, tvar))
tvar
write.csv(tvar, file = 'clusters_markers.csv')

mymat <- pancreas.integrated@assays$RNA #integrated
mymat <- as.matrix(mymat[getfound(mymarkers, rownames(mymat), v = T), ])
tvar <- make_list(remove.factors(pancreas.integrated@meta.data), gby, grouping = T)
tvar <- mixedsort(tvar)
void <- get_stat_report(mymat[, names(tvar)], groups = tvar, moments = c('bm', 'mn', 'p'), v = T)
rownames(void) <- paste0("'", rownames(void))
head(void)
write.csv(void, file = 'genes_stats_merged.csv')

freq_tablep(metadata = metadata, cnames = c('tag_FOXP3', 'Cluster'))

#### Differential Expression #### ------
myidents <- c(1, 4, 6, 8)[-4]
prefix <- 'dea_global'
sset <- c('orig.set', 'guo', 'zheng', 'zhang')
prefix <- 'dea_foxp3'
sset <- list(c('orig.set', 'guo', 'zheng', 'zhang'), c(gby, myidents))
prefix <- 'dea_foxp3_gs'

dir.create(prefix)
pancreas.subset <- SubsetData(pancreas.integrated, cells = getsubset(sset, pancreas.integrated[[]], v = T))
table(pancreas.subset@meta.data[, c('orig.set', gby)])
idents <- matrix(myidents)#unique(pancreas.integrated[, gby])
idents <- combinations(nrow(idents), r = 2, v = idents[, 1], set = TRUE, repeats.allowed = FALSE)
for(i in 1:nrow(idents)){
  identy <- idents[i, 1]
  if(ncol(idents) > 1) identy2 <- idents[i, 2] else identy2 <- NULL
  cat('Group(s)', commas(idents[i, ]), '\n')
  fname <- paste0(c(paste0(prefix, '/fdiffExp'), identy, identy2, sset[[1]][-1], '.csv'), collapse = "_")
  if(!file.exists(fname)){
    cmarkers <- FindConservedMarkers(object = pancreas.subset, ident.1 = identy, ident.2 = identy2, grouping.var = "orig.set", logfc.threshold = 0.1)
    cmarkers$min_avg_logFC <- apply(cmarkers[, getpats(colnames(cmarkers), 'avg_logFC')], 1, function(x){
      ifelse(all(min(x) * x > 0), min(x), 0)
    })
    write.csv(cmarkers, file = fname)
  }else cmarkers <- read.csv(fname, stringsAsFactors = F, check.names = F, row.names = 1)
  head(cmarkers); dim(cmarkers)
  # degs <- cmarkers$max_pval < 0.05
  # [rowSums(cmarkers[, getpats(colnames(cmarkers), '_pct.')] > 0.1) == 4, ]
  degs <- getDEGenes(cmarkers, pv = 0.05, upreg = T, pvtype = 'minimump_p_val', lfc.type = 'min_avg_logFC')
  degs <- rownames(cmarkers)[rownames(cmarkers) %in% degs]
  head(cmarkers[degs, ], 30)
  summary(abs(cmarkers[degs, ]$zheng_avg_logFC))

  # fname <- sub(paste0(prefix, "/f"), paste0(prefix, "/"), fname)
  # if(!file.exists(fname)){
  #   void <- DoHeatmap(pancreas.subset, features = degs, group.bar = T) + theme(axis.text.y = element_text(size = 4))
  #   thesegenes <- getfound(mymarkers, rownames(cmarkers[degs, ]), v = T)
  #   if(length(thesegenes) > 0) plots <- VlnPlot(object = pancreas.subset, features = thesegenes, split.by = 'orig.set', pt.size = 0, combine = FALSE)
  #   pdf(sub('\\.csv', '.pdf', fname), width = 16, height = 8)
  #   print(void)
  #   if(length(thesegenes) > 0) print(CombinePlots(plots = plots, ncol = fitgrid(thesegenes)[2]))
  #   graphics.off()
  # }
  cat('Done!\n')
}

mymarkersf <- list.files(prefix, pattern = 'fdiffExp.*csv', full.names = T)
cmarkers <- lapply(mymarkersf, read.csv, stringsAsFactors = F, check.names = F, row.names = 1)
tvar <- gsub(paste0(c(prefix, "/fdiffExp", "_", ".csv", sset[[1]][-1]), collapse = "|"), "", mymarkersf)
tvar <- sapply(strsplit(tvar, ""), paste, collapse =  "vs")
names(cmarkers) <- tvar
head(cmarkers[[1]])
degs <- unique(unlist(lapply(cmarkers, function(x){
  cnames <- getpats(colnames(x), '_pct.')
  getDEGenes(x[rowSums(x[, cnames] > 0.1) == length(cnames), ], pv = 0.05, upreg = T,
    pvtype = 'minimump_p_val', lfc.type = 'min_avg_logFC')
})))
length(degs)

fname <- paste0(c(prefix, 'degs', sset[[1]][-1], '.csv'), collapse = "_")
write.csv(degs, file = sub("degs", "degs_list", fname), row.names = F)
void <- DoHeatmap(pancreas.subset, features = degs, group.bar = T) + theme(axis.text.y = element_text(size = 1))
pdf(sub('\\.csv', '.pdf', fname), width = 16, height = 8)
print(void)
graphics.off()

cnames <- c('min_pct.1', 'min_pct.2', 'max_pval', 'minimump_p_val', 'min_avg_logFC')
cmarkerscombine <- rbindlist(lapply(names(cmarkers), function(y){
  x <- cmarkers[[y]]
  pct_names <- getpats(colnames(x), '_pct.')
  x <- x[getDEGenes(x, pv = 0.05, upreg = T, pvtype = 'minimump_p_val', , lfc.type = 'min_avg_logFC'), ]
  x$min_pct.1 <- apply(x[, getpats(colnames(x), 'pct.1')], 1, min)
  x$min_pct.2 <- apply(x[, getpats(colnames(x), 'pct.2')], 1, min)
  x <- x[abs(x$min_pct.1 - x$min_pct.2) > 0.01, ]
  cbind(gene_name = paste0("'", rownames(x)), cluster = y, x[, cnames])
}))
cmarkerscombine
write.csv(cmarkerscombine, file = fname, row.names = F)

ntopg <- 12
topgenes <- as.data.frame(rbindlist(lapply(levels(cmarkerscombine$cluster), function(x){
  dat <- cmarkerscombine[as.character(cmarkerscombine$cluster) == x, ]
  setorder(dat, minimump_p_val)
  head(dat, ntopg)
})))
mymarkers <- thesegenes <- unique(sub("'", "", as.character(topgenes$gene_name)))
void <- DotPlot(pancreas.subset, thesegenes, group.by = gby) +
  theme(axis.text.x = element_text(angle = 45, face = "bold", hjust = 1))
pdf(sub('\\.csv', paste0('top', ntopg, '.pdf'), fname), width = 16, height = 8)
print(void)
graphics.off()
fname <- sub('\\.csv', paste0('top', ntopg, '_vln.pdf'), fname) # back to tsne and vlnplots

mymat <- as.matrix(GetAssayData(object = pancreas.subset, slot = "data"))
# mymat <- as.matrix(GetAssayData(object = pancreas.subset, assay = "RNA"))
tmp <- getfound(degs, rownames(mymat), v = T)
metadata <- pancreas.subset[[]]
metadata <- remove.factors(metadata)
colnames(metadata) <- sub(gby, 'Cluster', colnames(metadata))
headmat(metadata); tailmat(metadata)

source('/mnt/BioHome/ciro/scripts/functions/group_specificity.R')
group_spec <- g_sp(
  cmarkers, # comparisons stats
  this_degs = NULL, # list of DEGs vectors per comparison
  fpvtype = 'minimump_p_val', # significance
  ffctype = 'min_avg_logFC', # fold-change
  padjthr = 0.05, # significance threshold
  fcthr = 0, # fc threshold
  methd = 'suas', # method
  gglobal = FALSE, # if data struture is for global
  gref = NULL, # reference group for activation
  sharedmax = 2, # maximum number of groups sharing a gene
  groups_cols = NULL, # groups colours data.frame
  expr_mat = mymat, # matrix for visualisation
  expr_mattype = 'SeuratIntegrated', # chosen matrix for visualisation
  datatype = 'sc', # to choose the visualisation
  vs = 'vs', # string splitting comparison names
  gtf = NULL, # extra info for genes, data.frame
  gtfadd = NULL, # columns to add
  path_plot = 'dea_foxp3_gsa', # path to plot
  annotation = metadata, # annotation for samples
  cname = 'Cluster', # name id for files
  hmg_order = NULL, # order of groups data.frame
  ngenes = 20, # number of genes to plot
  sufix = 'mean', # sufix to add to file names
  order_by = 'column_name', # order samples in heatmap
  hm_order = 'minFC_p', # gene order per group
  sepchar = 'n',
  log_norm = FALSE,
  coulrange = c('blue', 'black', 'yellow'), # colours to use
  groupsamp = FALSE, # sample samples in group to plot
  verbose = TRUE, # Print progress
  myseed = 27 # seed for determinism
)
