#!/usr/bin/R

# This script is to create visualise the results from monocle and signature
# analyses with analysis of data from Guo et al

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

#### loading files #### --------------------------------------------------------
# Colours
colfile <- '/mnt/BioHome/ciro/simon/info/groupColours_art_ct.csv'
gr.cols <- read.csv(colfile, row.names = 1, stringsAsFactors = F)
# Classifed cells annotation
# fname = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/seurat/guo_lung/requests/dens_props/tumour_cd4.csv'
fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/treg_classification.csv'
annot <- read.csv(fname, stringsAsFactors = F, check.names = F, row.names = 1)
# grps <- c('cTREGoFOXP3N', 'TREG_TNFRSF9-', 'TREG_TNFRSF9+', 'cTREGoTFR', 'cTFR')[-1]
grps <- c('TREG_TNFRSF9.neg', 'TREG_TNFRSF9.pos', 'cTFR', 'cTREGoTFR')
tannot <- annot[getsubset(c('orig.ctact', grps), annot, v = T), ]
tannot$cellname <- tannot$orig.ctact
# Expression data
tpm <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_tpm.RData')
# log2tpm <- log2(tpm[, rownames(tannot)] + 1)
# signature
# sigfile <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/seurat/guo_lung/requests/cTREGoTFR_19g/cTREGoTFR_19g_signature.csv'
# sigfile <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/TFR_94g/TFR_94g_signature.csv'
sigfile <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/TFR_combined_34g/TFR_combined_34g_signature.csv'
signature_tab <- read.csv(sigfile, stringsAsFactors = F, row.names = 1)
# Monocle components
# comps <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/monocle/treg_tfr2/data/DDRTree_gs_14PC_tpm_components.RData')
comps <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/monocle/guo_lung/data/DDRTree_seu_14PC_tpm_components.RData')
setwd('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/')

### ploting against monocle result ### -----------------------------------------
df <- cbind_repcol(comps, tannot[, 'cellname', drop = F])
df <- df[getorder(grps, df, 'cellname'), ]
df <- cbind_repcol(df, signature_tab)
head(df)
df$cellname <- factor(as.character(df$cellname), levels = grps)
levels(df$cellname)
unique(df$cellname)
mycols <- gr.cols[grps, ]
p <- ggplot(df, aes(-MC1, TFR_combined_34g1, color = cellname)) + geom_point(size = 2.3) +
  scale_color_manual(values = mycols) + geom_smooth(se = FALSE, color = 'gray15', size = 0.5) + # LOESS is used for the line by default if < 1000 obs
  labs(title = 'Component vs Signature', y = expression(T[F*R]*" signature"), x = "Component 1")
  # "TFR signature"
pdf('monocle_signature14pc_TFR_combined_34g.pdf', 10, 8); p; dev.off()

### Plotting Monocle ### -------------------------------------------------------
cds <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/monocle/guo_lung/data/DDRTree_seu_14PC_tpm_data.RData')
# cds@reducedDimS[2, ] <- -1 * cds@reducedDimS[2, ] # this disrupts the tree plotted behind
pdfsize <- round(log10(nrow(tannot)) * 3.5, 1)
catg <- 'orig.ctact'
fsufix <- 'monocle_'
cell_type_color <- lapply(tannot[, catg, drop = FALSE], function(x){
  v2cols(select = x, sour = gr.cols, fw = 'gg')
})
library(monocle)
void <- lapply(catg, function(catgy){
  plot_cell_trajectory(cds, color_by = catgy, cell_size = 3) + scale_color_manual(values = cell_type_color[[catgy]])
})
lablsx <- c(3, 0, -3, -6); names(lablsx) <- lablsx * -1
lablsy <- c(0, -2, -4); names(lablsy) <- lablsy * -1
pdf(paste0(fsufix, 'trajectory', pdfsize, '_flipped.pdf'), height = pdfsize, width = pdfsize)
void[[1]] + scale_y_reverse(breaks = as.numeric(names(lablsy)), labels = lablsy) + scale_x_reverse(breaks = as.numeric(names(lablsx)), labels = lablsx)
graphics.off()

## Violin IL2 transcript
edata_list <- theObjectSavedIn('../data/edata_list.RData')
annot_list_tags <- theObjectSavedIn('../data/annot_list_tags_10TPM.RData')
g2plot <- c("IL2", "IL10", "TCF7")
c2see <- 'tag_ct'
fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/tfr_signature_genes.csv'
mytab <- readfile(fname, row.names = 1, stringsAsFactors = F)
gens <- rownames(mytab)
setwdc('~/large/simon/results/integration/gene_expression')
voidz <- list()
for (i in names(x = annot_list_tags)) {
  print(i)
  annot <- annot_list_tags[[i]]
  mat <- edata_list[[i]]
  tvar <- annot#[getsubset(c(c2see, names(table(annot[, c2see]))), annot, v = T), ]
  fname <- paste0('violins_', i, '.pdf')
  if(!file.exists(fname)){
    pdf(fname)
    vlnplot(mat, tvar, gg = g2plot, orderby = c2see, log2t = TRUE, v = TRUE, rotx = 90, plotdots = T)
    vlnplot(mat, tvar, gg = g2plot, orderby = c2see, log2t = TRUE, v = FALSE, rotx = 90, legendt = 'means', plotdots = T)
    vlnplot(mat, tvar, gg = g2plot, orderby = c2see, log2t = TRUE, v = TRUE, rotx = 90, noncero = FALSE, plotdots = T)
    graphics.off()
  }
  tmp <- make_list(tvar, colname = c2see, grouping = TRUE)
  mymarkers <- getfound(c(g2plot, gens), rownames(mat), v = T)
  void <- get_stat_report(mat[, names(tmp)], groups = tmp, rnames = mymarkers, v = T)
  head(void)
  voidz[[i]] <- void[g2plot, ]
  # fname <- paste0('markers_', i ,'.csv')
  # if(!file.exists(fname)){
  #   write.csv(void, file = fname)
  # }
}
il2_across_sets <- data.frame(rbindlist(voidz))
rownames(il2_across_sets) <- paste0(rep(g2plot, each = length(annot_list_tags)), "_", names(voidz))
write.csv(il2_across_sets, file = 'markers.csv')

load('../seurat_9sets/integrated.RData')
mymarkers <- g2plot
pdf(paste0('markers_split', paste0(g2plot, collapse = "_"), '.pdf'), height = 15, width = 25)
print(FeaturePlot(pancreas.integrated, features = mymarkers, min.cutoff = 0, split.by = 'orig.set'))
for(i in 1:length(mymarkers)){
}
graphics.off()
dim(pancreas.integrated@assays$integrated@data)
mat <- pancreas.integrated@assays$integrated@data
mat <- do.call(cbind, lapply(edata_list, function(x) x[g2plot, ] ))
fname <- paste0('violins_all_raw.pdf')
if(!file.exists(fname)){
  pdf(fname, 8, 24)
  vlnplot(mat, pancreas.integrated@meta.data, gg = g2plot, orderby = 'orig.set', covi = c2see, log2t = grepl("raw", fname), v = TRUE, rotx = 90, plotdots = T)
  vlnplot(mat, pancreas.integrated@meta.data, gg = g2plot, orderby = 'orig.set', covi = c2see, log2t = grepl("raw", fname), v = FALSE, rotx = 90, legendt = 'means', plotdots = T)
  vlnplot(mat, pancreas.integrated@meta.data, gg = g2plot, orderby = 'orig.set', covi = c2see, log2t = grepl("raw", fname), v = FALSE, rotx = 90, noncero = TRUE, plotdots = T)
  graphics.off()
}
