``#!/usr/bin/R5

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

# columns for visualisation
orignames <- c("orig.set", "orig.majorCluster", 'tag_FOXP3', 'tag_ct')
npcs <- 30

## Select genes
anfeats <- 500
anfeats <- anfeats[rowSums(sapply(edata_list, function(x) anfeats %in% rownames(x) )) == 9]

root <- '~/large/simon/results/integration'
setwd(root)

annot_list_tags <- theObjectSavedIn(paste0(root, '/data/annot_list_tags_10TPM.RData'))
cellnames <- lapply(annot_list_tags, function(x) getsubset(c("tag_ct", "TFR"), x, v = TRUE) )
cellnames <- cellnames[sapply(cellnames, length)  > 30]
edata_list <- theObjectSavedIn(paste0(root, '/data/edata_list.RData'))
studname <- names(cellnames)
annot_list_tags <- lapply(studname, function(x) annot_list_tags[[x]][cellnames[[x]], ] )
edata_list <- lapply(studname, function(x) edata_list[[x]][, cellnames[[x]]] )
names(edata_list) <- studname
names(annot_list_tags) <- studname
gc()

# Create Seurat object per data set
pancreas.list <- lapply(names(annot_list_tags), function(x){
  CreateSeuratObject(counts = edata_list[[x]], meta.data = annot_list_tags[[x]])
}); names(pancreas.list) <- names(annot_list_tags)

sets <- names(pancreas.list)
selected <- "sets_tfr_"
setwdc(paste0(root, '/seurat_', length(sets), selected, ifelse(is.numeric(anfeats), anfeats, length(anfeats)), "hvg"))

if(file.exists('integrated.RData')) cat('Go to file creating line\n')

for (i in 1:length(x = pancreas.list)) {
  cat(names(pancreas.list)[i], '\n')
  pancreas.list[[i]] <- NormalizeData(object = pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(object = pancreas.list[[i]],
      selection.method = "vst", nfeatures = ifelse(is.numeric(anfeats), anfeats, length(anfeats)), verbose = FALSE)
  cat(commas(VariableFeatures(object = pancreas.list[[i]]), 10), '\n')
}
# Check most variable genes overlaps
myvargenes <- unique(unlist(lapply(pancreas.list, VariableFeatures)))
ogenes <- sapply(pancreas.list, function(x) myvargenes %in% VariableFeatures(x) )
rownames(ogenes) <- myvargenes
head(ogenes, 20)
myvargenes <- myvargenes[apply(ogenes, 1, all)]
length(myvargenes)

npcs <- 20
kfilter <- min(c(200, sapply(pancreas.list, ncol)))
sets <- names(sort(sapply(pancreas.list, ncol)))
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list[sets], anchor.features = anfeats, dims = 1:npcs, k.filter = kfilter)#, k.anchor = 5)#, k.score = npcs)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:npcs)#, k.weight = npcs)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(object = pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = npcs, verbose = FALSE)

pdf(paste0('sdevPCs_', npcs,'PCs.pdf'), width = 10, height = 8)
ElbowPlot(object = pancreas.integrated, ndims = npcs)
graphics.off()

elbows <- get_elbow(x = pancreas.integrated@reductions$pca@stdev, threshold = seq(.7, .9, .02), decide = TRUE)
print(elbows);

chnpcs <- 9
nres <- 0.2
redu <- "umap"

setwdc(paste0(root, '/seurat_', length(sets), selected, ifelse(is.numeric(anfeats), anfeats, length(anfeats)), "hvg/PC", chnpcs, 'R', nres))

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
  pdf(paste0('integrated_', sub('orig.', '', orig), "_", redu, '.pdf'), height = 10, width = ifelse(tvar > 15, 14, 10))
  print(DimPlot(object = pancreas.integrated, reduction = redu, group.by = orig))
  graphics.off()
}
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:chnpcs)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = nres)
tailmat(pancreas.integrated[[]], 10)
gby <- paste0("integrated_snn_res.", nres)
pancreas.integrated@meta.data[, gby] <- as.character(pancreas.integrated@meta.data[, gby])

pdf(paste0('clusters_', gby, "_", redu, '.pdf'), height = 8, width = 8)
DimPlot(object = pancreas.integrated, reduction = redu, group.by = gby)
graphics.off()

freq_tablep(metadata = pancreas.integrated@meta.data, cnames = c(gby, 'orig.set'),
  pnames = c('Clusters in sets', 'Sets in clusters'), dowrite = TRUE)
freq_tablep(metadata = pancreas.integrated@meta.data, cnames = c(gby, 'tag_ct'),
  pnames = c('Clusters in cell types', 'Cell types in clusters'), dowrite = TRUE)

markers <- read.csv('/mnt/BioHome/ciro/simon/info/markers.csv', stringsAsFactors = F)
markers
mymarkers <- unique(c('FOXP3', markers[, 1], 'IL2RA', 'TNFRSF9', 'TNFRSF18', 'DUSP4',
  'CCR8', 'IL1R2', 'IKZF2', 'ENTPD1', 'LAG3', 'TIGIT', 'CTLA4', 'PDCD1','TOX'))[1]
mymarkers <- getfound(mymarkers, rownames(pancreas.integrated@assays$RNA), v = T)

DefaultAssay(object = pancreas.integrated) <- "RNA"

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

DefaultAssay(object = pancreas.integrated) <- "integrated"

DefaultAssay(object = pancreas.integrated)
dim(pancreas.integrated@assays$integrated@counts)
dim(pancreas.integrated@assays$integrated@data)

save(pancreas.integrated, file = '../integrated.RData')
load('../integrated.RData') #### --------------------------------------------------
