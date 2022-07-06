#!/usr/bin/R

# This script is to classify TREG cell types from data from Guo et al

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
### Loading data ### -----------------------------------------------------------
# metadata <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_annotation.RData')
# ln -s /mnt/BioAdHoc/Groups/vd-vijay/Ariel/SingleCell/Simon/guo-2018-reanalysis/tracer/2019-06-03-analysis-cutoff10/tracer-annotation.RData \
#   /mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_annotation_tracer.RData
metadata <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_annotation_tracer.RData')
metadata <- remove.factors(metadata)
str(metadata)
tpm <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_tpm.RData')
setwdc('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/')

### Setting relevant columns ###
metadata$comp_cloneID <- paste0(metadata$Name.TRA.and.TRB)
table(metadata$GreaterTHAN1, metadata$orig.cloneStatus)
metadata$orig.clonexp <- ifelse(metadata$GreaterTHAN1 > 1, 'Clonal', 'NoClonal')
table(metadata$GreaterTHAN1, metadata$orig.clonexp)
# celltypes
table(metadata$group)
# FOXP3.neg              TFR TREG_TNFRSF9.neg TREG_TNFRSF9.pos
#      1675              289              471              348
metadata$orig.fctype <- gsub('(^....).*', '\\1', metadata$group)
metadata$orig.fctype <- gsub('FOXP', 'FOXP3-', metadata$orig.fctype)
table(metadata$orig.fctype)
# FOXP  TFR TREG
# 1675  289  819

# NOTE: At this point, everything is clonally expanded <<<<<<<<<<<<<<<<<<<<<<<<<<
cells <- list(c('orig.tissue', 'T'), c('orig.celltype','CD4'), c('orig.clonexp', 'Clonal')) # taking only tumour CD4 cells
# cells[[length(cells) + 1]] <- c('orig.Patient', 'P0706', 'P0913', 'P1120', 'P1010')
tmetadata <- metadata[getsubset(cells, metadata, v = T), ]

# Proportion of Tregs sharing with Tfrs and not Tfrs DE comparison
table(tmetadata$orig.fctype)
# make table with sharing cells
tmp <- vlist2df(count_overlap(c(c('orig.fctype', 'TREG', 'TFR', 'FOXP3-', 'comp_cloneID')), tmetadata, combs = 'r'))
sapply(tmp, function(x) sum(!is.na(x))); head(tmp)
# split cells
tregs <- tmetadata[getsubset(c('orig.fctype', 'TREG'), tmetadata, v = T), ]
tfrs <- tmetadata[getsubset(c('orig.fctype', 'TFR'), tmetadata, v = T), ]
nfoxp3 <- tmetadata[getsubset(c('orig.fctype', 'FOXP3-'), tmetadata, v = T), ]
treg_tfr <- tregs[tregs[, 'comp_cloneID'] %in% tfrs[, 'comp_cloneID'], ] # TREG clones in TFR
dim(treg_tfr)
treg_nfoxp3 <- tregs[tregs[, 'comp_cloneID'] %in% nfoxp3[, 'comp_cloneID'], ] # TREG clones in FOXP3- cells
dim(treg_nfoxp3)
# myfreqs <- share_counter(selection = c('orig.fctype', 'TREG', 'TFR', 'FOXP3-', 'comp_cloneID'), doby = 'orig.Patient', annot = tmetadata)
# myfreqs <- share_counter(selection = c('orig.fctype', 'TFR', 'TREG', 'FOXP3-', 'comp_cloneID'), doby = 'orig.Patient', annot = tmetadata)
# myfreqs <- share_counter(selection = c('orig.fctype', 'TFR', 'TREG', 'FOXP3-', 'comp_cloneID'), doby = 'orig.Patient', annot = tmetadata, combs = TRUE)

# testing if cell names are mutually excluded
tvar <- vlist2df(overlap_calc(list(nfoxp3 = rownames(nfoxp3), tfrs = rownames(tfrs), tregs = rownames(tregs))))
head(tvar); sapply(tvar, function(x) sum(!is.na(x)) )
mytabl <- list(
  TREGoFOXP3N = setdiff(rownames(treg_nfoxp3), rownames(treg_tfr)),
  TREGoTFR = setdiff(rownames(treg_tfr), rownames(treg_nfoxp3)),
  TFR = rownames(tfrs),
  FOXP3N = rownames(nfoxp3)
)
headmat(tmetadata[intersect(rownames(treg_tfr), rownames(treg_nfoxp3)), ], 20)
mytabl$TREG <- setdiff(rownames(tregs), unlist(mytabl))
names(mytabl) <- paste0('c', names(mytabl)) # tagging them as clonal
str(mytabl)

# Now taking from all the cells only the no clonally expanded
cells <- list(c('orig.tissue','T'), c('orig.celltype','CD4'), c('orig.clonexp', 'NoClonal'))
nc_metadata <- metadata[getsubset(cells, metadata, v = T), ]
tvar <- make_list(nc_metadata, 'orig.fctype')
str(c(mytabl, tvar))
mytab <- melt(c(mytabl, tvar))
head(mytab)
rownames(mytab) <- mytab[, 1]; mytab <- mytab[, -1, drop = FALSE]
colnames(mytab) <- 'orig.cell_set'
mytab <- cbind_repcol(mytab, metadata) # repcol removes the repeated columns and binds matching the rownames
headmat(mytab); table(mytab[, 1])
cat('Total:', nrow(mytab), '\nAll metas Clonal or NoClonal:', nrow(nc_metadata) + nrow(tmetadata), '\n')
table(mytab$orig.cell_set, mytab$group)
mytab$orig.ctact <- mytab$orig.cell_set
tvar <- mytab$orig.cell_set == "TREG"
mytab[tvar, ]$orig.ctact <- mytab$group[tvar]
table(mytab$orig.ctact)
# write.csv(mytab, file = 'treg_classification.csv')
# save(mytab, file = 'treg_classification.RData')
# thistab = read.csv(fname, row.names = 1)
