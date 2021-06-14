#!/usr/bin/R

##########################
# Fetch single-cell data #
##########################

# This script is to download and gather the public single-cell data

ROOTDIR=/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw
cd ${ROOTDIR}
mkdir ${ROOTDIR}/GSEs

#### Lambrechts - Lung ##### ---------------------------------------------------
mkdir ${ROOTDIR}/lambrechts
cd ${ROOTDIR}/lambrechts
## --------- deprecated ---------
# # manually downloaded files from: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6149/files/
# unzip E-MTAB-6149.processed.4.zip
# # go to R
# load('T_cell.Cellview.Rds')
# head(tsne.data)
# annot <- tsne.data[, -3]
# colnames(annot) <- c('tSNE_1', 'tSNE_2', 'cluster_number')
# head(annot)
# annot$cluster_number <- paste0("G", annot$cluster_number)
# table(annot$cluster_number)
# ident <- c('C1', 'C2', 'C3', 'C4', 'C5', 'NK', 'Treg', 'Tprolif', 'C9')
# names(ident) = paste0("G", 1:length(ident))
# annot$orig.majorCluster <- ident[annot$cluster_number]
# table(annot$orig.majorCluster)
# annot$orig.celltype <- annot$orig.majorCluster
# annot[annot$orig.majorCluster %in% c('Treg', 'C1', 'C3', 'C9'), 'orig.celltype'] <- 'CD4'
# annot[annot$orig.majorCluster %in% c('C2', 'C4', 'C5', 'Tprolif'), 'orig.celltype'] <- 'CD8'
# table(annot$orig.celltype)
# annot <- cbind(UniqueCell_ID = rownames(annot), annot)
# save(annot, file = '../lambrechts_EMTAB6149Tcells_annotation.RData')
# head(log2cpm[, 1:10])
# all(colnames(log2cpm) == rownames(annot))
# all(rownames(log2cpm) == rownames(featuredata))
# sum(rownames(log2cpm) %in% rownames(featuredata))
# featuredata <- featuredata[rownames(log2cpm), ]
# head(featuredata)
# tvar <- duplicated(featuredata[, 'Associated.Gene.Name'])
# sum(tvar)
# rownames(log2cpm) <- featuredata[, 'Associated.Gene.Name']
# log2cpm <- as.matrix(log2cpm)
# save(log2cpm, file = '../lambrechts_EMTAB6149Tcells_log2cpm.RData')
## --------- deprecated ---------

# from https://gbiomed.kuleuven.be/english/research/50000622/laboratories/54213024/scRNAseq-NSCLC
# download the T cells
# go to R - this was done with Ariel's session bc I had issues with loomR
library(loomR)
fname <- list.files(pattern = 'loom')[2]
lfile <- connect(filename = fname, mode = "r+")
lfile
umi <- t(lfile$matrix[, ])
rownames(umi) <- lfile[["row_attrs/Gene"]][]
colnames(umi) <- lfile[["col_attrs/CellID"]][]
head(umi[, 1:8])
dim(umi)
attrs <- head(names(lfile[['col_attrs']]), -1)
annot <- lfile$get.attribute.df(MARGIN = 2, col.names = "CellID", attribute.names = attrs)[, attrs]
dim(x = annot)
head(x = annot)
table(annot$ClusterName)
table(annot$ClusterName, annot$ClusterID)
annot$orig.celltype <- sub(".*\\((...).*", "\\1", annot$ClusterName)
annot$orig.celltype <- sub("nat", "NK", annot$orig.celltype)
annot$orig.celltype <- sub("reg", "CD4", annot$orig.celltype)
table(annot$orig.celltype)
annot$orig.majorCluster <- sub(".*([0-9] \\(.*)cells.*", "\\1", annot$ClusterName)
annot$orig.majorCluster <- gsub("\\+| |-", "", annot$orig.majorCluster)
annot$orig.majorCluster <- gsub("\\(", "_", annot$orig.majorCluster)
annot$orig.majorCluster <- gsub("naturalkiller", "NK", annot$orig.majorCluster)
annot$orig.majorCluster <- gsub("regulatoryT", "Treg", annot$orig.majorCluster)
table(annot$orig.majorCluster)
annot <- annot[, c('CellID', 'orig.majorCluster', 'orig.celltype', 'Embeddings_X', 'Embeddings_Y')]
colnames(annot) <- sub('CellID', 'UniqueCell_ID', colnames(annot))
annot$orig.tissue <- 'T'

markers <- c('FOXP3', 'CXCR5', 'BCL6', 'CD4')
cts <- sweep(umi, 2, colSums(umi), '/') * 1e+4
tvar <- cbind(annot, t(log2(cts[markers, rownames(annot)] + 1)))
g <- lapply(markers, function(g){
  ggplot(tvar, aes_string('Embeddings_Y', 'Embeddings_X', color = g)) +
    geom_point(size = 0.2) + scale_color_gradient(low = 'yellow', high = 'red')
})
g[[length(g) + 1]] <- ggplot(tvar, aes(Embeddings_Y, Embeddings_X, color = orig.majorCluster)) + geom_point(size = 0.2)
g[[length(g) + 1]] <- ggplot(tvar, aes(Embeddings_Y, Embeddings_X, color = orig.celltype)) + geom_point(size = 0.2)
pdf('~/large/simon/raw/lambrechts/tSNE.pdf', 18, 10)
plot_grid(plotlist = g)
dev.off()

save(annot, file = '../lambrechts_LOOMTCELLS_annotation.RData')
save(umi, file = '../lambrechts_LOOMTCELLS_umi.RData')
# from Ciro's terminal
cp /mnt/BioHome/amadrigal/lambrechts_loomTcells_annotation.RData ../
cp /mnt/BioHome/amadrigal/lambrechts_loomTcells_umi.RData ../
# back to Ariel's tarminal
rm ~/lambrechts_loomTcells_annotation.RData
rm ~/lambrechts_loomTcells_umi.RData

# done with created RData objects
load('../lambrechts_LOOMTCELLS_annotation.RData')
load('../lambrechts_LOOMTCELLS_umi.RData')
mycells <- CreateSeuratObject(raw.data = umi, meta.data = annot)
save(mycells, file = '../lambrechts_loomTcells_seurat.Robj')

#### Tirosh - oligodendroglioma ##### ------------------------------------------
cd ${ROOTDIR}/GSEs
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70630/soft/GSE70630_family.soft.gz
gunzip GSE70630_family.soft.gz
less GSE70630_family.soft | grep TPM | head | grep TPM # to confirm what type of counts they are
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70630/suppl/GSE70630_OG_processed_data_v2.txt.gz
gunzip GSE70630_OG_processed_data_v2.txt.gz
# go to R
edata <- read.table('GSE70630_OG_processed_data_v2.txt', row.names = 1, header = T, check.names = F, stringsAsFactors = F)
tpm = as.matrix(edata)
colSums(((2**tpm[, 1:20]) - 1))# * 10)
tpm <- ((2**tpm) - 1)# * 10
summary(tpm[, 1:10]); dim(tpm)
cnames <- colnames(tpm); head(cnames)
cnames[!grepl("MGH", cnames)] <- paste0("MGH", cnames[!grepl("MGH", cnames)])
colnames(tpm) <- cnames
annot <- data.frame(UniqueCell_ID = cnames,
  orig.tumour = sub("(.*)_.*_.*", "\\1", cnames),
  orig.pnumber = sub(".*_(.*)_.*", "\\1", cnames),
  orig.tag_number = sub(".*_.*_(.*)", "\\1", cnames),
  row.names = cnames,
stringsAsFactors = FALSE)
annot$orig.tissue <- 'T'
sapply(annot[, -1], table)

newann <- read.table('GSE70630_cell_type_assignment.txt', sep = '\t', header = 1, stringsAsFactors = FALSE)
rownames(newann) <- newann[, 1]; newann <- newann[-1, -1]; head(newann)
cnames <- rownames(newann); head(cnames)
cnames[!grepl("MGH", cnames)] <- paste0("MGH", cnames[!grepl("MGH", cnames)])
rownames(newann) <- cnames
table(newann[, 2]) == table(newann[, 1]) # same number
annot$orig.celltype <- newann[rownames(annot), 1]
annot$orig.celltype[annot$orig.celltype == 'Microglia/Macrophage'] <- 'microglia_macrophage'
annot$orig.celltype[grepl('Oligodendrocytes', annot$orig.celltype)] <- 'non_malignant'
table(annot$orig.celltype)
head(annot)

rownames(tpm)[grepl('^IG', rownames(tpm))]
# Markers from: https://satijalab.org/seurat/pbmc3k_tutorial_1_4.html and wikipedia
markers <- read.csv('/mnt/BioHome/ciro/simon/info/markers.csv', stringsAsFactors = F)
markers
log2tpm <- log2(tpm[markers[, 1], ] + 1)

source('~/scripts/functions/handy_functions.R')
pdf(paste0('tirosh_expression_', sub("orig.", "", "groups"),'.pdf'), 10, 10)
for(tvar in c('orig.tumour', 'orig.celltype', 'orig.majorCluster')){
  if(!tvar %in% colnames(annot)) next
  if(length(unique(annot[, tvar])) <= 1) next
  vlnplot(cpmdata = log2tpm, metadata = annot, orderby = tvar, v = T,
    gg = getfound(markers[, 1], rownames(log2tpm), v = T), rotx = 45, datatype = 'log2(tpm + 1)')
}
dev.off()

save(annot, file = '../tirosh_GSE70630_all_annotation.RData')
save(tpm, file = '../tirosh_GSE70630_all_tpm.RData')

#### Puram - head_neck ##### ---------------------------------------------------
cd ${ROOTDIR}/GSEs
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz
gunzip GSE103322_HNSCC_all_data.txt.gz
# go to R
edata <- data.table::fread('GSE103322_HNSCC_all_data.txt')
edata[1:20, 1:5]
preannot <- data.frame(t(edata[1:5, ]), stringsAsFactors = FALSE)[-1, ]
head(preannot)
# colnames(preannot) <- c('enzyme', 'LN', 'orig.tissue', 'noncancer', 'orig.majorCluster')
# annot <- data.frame(enzyme = ifelse(preannot$enzyme == 1, 'max', 'other'),
#   orig.tumour = ifelse(preannot$orig.tissue == 1, 'T', 'NT'),
#   orig.tissue = "T",
#   orig.majorCluster = sub("-| ", "", preannot$orig.majorCluster),
#   row.names = colnames(edata)[-1],
# stringsAsFactors = FALSE)
lapply(preannot, table, useNA = 'always')
colnames(preannot) <- c('enzyme', 'LN', 'cancer', 'noncancer', 'majorCluster')
annot <- data.frame(enzyme = ifelse(preannot$enzyme == 1, 'maxima', 'other'),
  orig.tissue = ifelse(preannot$LN == 1, 'LN', 'T'),
  orig.cancer = ifelse(preannot$cancer == 1, 'C', 'NC'),
  orig.noncancer = ifelse(preannot$noncancer == 1, 'NC', 'C'),
  orig.majorCluster = sub("-| ", "", preannot$majorCluster),
  row.names = colnames(edata)[-1],
stringsAsFactors = FALSE)
sapply(annot, table)
table(annot$orig.cancer, annot$orig.tissue)
table(annot$orig.cancer, annot$orig.majorCluster)
table(annot$orig.noncancer, annot$orig.majorCluster)
table(annot$orig.tissue, annot$orig.majorCluster)
annot$orig.celltype <- annot$orig.majorCluster
# annot[annot$orig.majorCluster %in% c('Dendritic', 'Tcell', 'Macrophage'), 'orig.celltype'] <- 'CD4'
annot[annot$orig.majorCluster %in% 'Tcell', 'orig.celltype'] <- 'CD4'
table(annot$orig.celltype)
annot <- cbind(UniqueCell_ID = rownames(annot), annot)
annot$UniqueCell_ID <- as.character(annot$UniqueCell_ID)
head(annot)
save(annot, file = '../puram_GSE103322HNSCC_all_annotation.RData')
tpm <- data.frame(edata[-c(1:5), ], stringsAsFactors = FALSE)[, -1]
tpm <- apply(tpm, 2, as.numeric)
rownames(tpm) <- gsub("'", "", unlist(edata[-c(1:5), 1]))
summary(tpm[, 1:5])
colSums(((2**tpm[, 1:20]) - 1))
tpm <- ((2**tpm) - 1)
tpm <- tpm * 10 # it is in the lower range for a TPM
colSums(tpm[, 1:20])
save(tpm, file = '../puram_GSE103322HNSCC_all_tpm.RData')
# for Divya and Oliver
library(Seurat)
mycells <- CreateSeuratObject(raw.data = tpm, meta.data = annot)
save(mycells, file = '../puram_GSE103322HNSCC_seurat.Robj')

#### Arnon - head_neck ##### ---------------------------------------------------
mkdir ${ROOTDIR}/arnon
cd ${ROOTDIR}/arnon
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_cell.annotations.csv.gz
gunzip GSE115978_cell.annotations.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_tpm.csv.gz
gunzip GSE115978_tpm.csv.gz
# go to R
preannot <- read.csv('GSE115978_cell.annotations.csv', row.names = 1, header = 2, check.names = F, stringsAsFactors = F)
head(preannot)
annot <- preannot[, 1:4]; colnames(annot) <- c('orig.tumour', 'orig.majorCluster', 'orig.treatment', 'orig.cohort')
sapply(annot, table)
annot$orig.treatment <- ifelse(annot$orig.treatment == "treatment.naive", "Pre", "Post") # following Sade's name
annot$orig.BL_TRT <- annot$orig.treatment
table(annot$orig.BL_TRT, annot$orig.treatment)
annot$orig.majorCluster[annot$orig.majorCluster == "?"] <- "Unknown"
annot$orig.celltype <- annot$orig.majorCluster
annot[annot$orig.majorCluster %in% c('T.CD4', 'T.cell'), 'orig.celltype'] <- 'CD4'
annot[annot$orig.majorCluster %in% c('T.CD8'), 'orig.celltype'] <- 'CD8'
table(annot$orig.celltype)
annot$orig.tissue <- "T"
annot <- cbind(UniqueCell_ID = rownames(annot), annot)
annot$UniqueCell_ID <- as.character(annot$UniqueCell_ID)
save(annot, file = '../arnon_GSE115978_annotation.RData')

edata <- data.table::fread('GSE115978_tpm.csv')
head(edata[, 1:10])
summary(edata[, 2:10])
tpm <- as.matrix(edata[, -1])
rownames(tpm) <- unlist(edata[, 1])
str(tpm)
colSums(((2**tpm[, 1:20]) - 1))
tpm <- ((2**tpm) - 1)
summary(tpm[, 1:10])
colSums(tpm[, 1:10])
tpm <- tpm * 10 # in the paper they said they divided by 10
save(tpm, file = '../arnon_GSE115978_tpm.RData')
# # Download from: https://portals.broadinstitute.org/single_cell/study/SCP109/melanoma-immunotherapy-resistance
# # go to R
# preannot <- read.table('metadata_all.txt', row.names = 1, header = 2, check.names = F, stringsAsFactors = F)
# head(preannot)
# annot <- preannot[-1, 1:2]; colnames(annot) <- c('orig.clinic', 'orig.majorCluster')
# head(annot)
# sapply(annot, table)
# preannot <- read.table('tumors.nonmal_tsne_anno.txt', row.names = 1, header = 2, check.names = F, stringsAsFactors = F)
# preannot <- read.table('tumors.mal_tsne_anno.txt', row.names = 1, header = 2, check.names = F, stringsAsFactors = F)

#### Li - Melanoma ##### -------------------------------------------------------
mkdir ${ROOTDIR}/data_Li2018/matrix
cd ${ROOTDIR}/data_Li2018/matrix
# download from Quick Links: https://bitbucket.org/tanaylab/li_et_al_cell_2018_melanoma_scrna/src/default/
wget http://www.wisdom.weizmann.ac.il/~lubling/Li2018/Li2018_umi_counts.tar.gz
tar -xvzf Li2018_umi_counts.tar.gz
mv Li2018_cells.tsv barcodes.tsv
mv Li2018_matrix.mtx matrix.mtx
paste Li2018_genes.tsv Li2018_genes.tsv > genes.tsv
# go to R
mat <- Seurat::Read10X(data.dir = './')
umi <- as.matrix(mat)
dim(umi)
head(umi[, 1:10])
annot <- read.table('Li2018_metadata.tsv', sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
head(annot)

annoth <- read.table('dep/Li2018_metadata_header.tsv', sep = "\t", stringsAsFactors = FALSE)
annotd <- read.table('dep/Li2018_metadata.tsv', skip = 1, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
colnames(annotd) <- gsub(" ", "_", annoth[1, ])

data.frame(colnames(annotd), new =  c(colnames(annot)[1:2], rep("================", 2),
  colnames(annot)[3:25], "================", colnames(annot)[26:35]))
all(annot$non_t_nk_mc_group == annotd$non_t_nk_mc_group, na.rm = T)
all(annot$TCRseq.available == annotd$TCRseq_available, na.rm = T)
all(annot$Median.Net.gene.UMIs == annotd$Median_Net_gene_UMIs, na.rm = T)
# Using not updated one, given that his update eliminate some columns
annot <- annotd

sapply(annot[, grepl("_mc", colnames(annot))], table, useNA = 'always')
sapply(annot[, grepl("group", colnames(annot))], table, useNA = 'always')
table(annot$all_mc_group, annot$t_nk_mc_group)
table(annot$non_t_nk_mc_group, annot$t_nk_mc_group) # no overlap
table(annot$all_mc_group, annot$non_t_nk_mc_group)
annot$orig.majorCluster <- annot$t_nk_mc_group # merging the sub-groups
annot$orig.majorCluster[is.na(annot$orig.majorCluster)] <- annot$non_t_nk_mc_group[is.na(annot$orig.majorCluster)]
annot$orig.majorCluster[is.na(annot$orig.majorCluster)] <- annot$all_mc_group[is.na(annot$orig.majorCluster)]
annot$orig.majorCluster[is.na(annot$orig.majorCluster)] <- "Unknown"
table(annot$orig.majorCluster, useNA = 'always')
table(annot$orig.majorCluster, annot$Cell_type, useNA = 'always') # to check how they are subdivided
table(annot$Cell_type)
annot$orig.celltype <- casefold(sub("\\(", "_", gsub(" |cells|\\+|\\)", "", annot$Cell_type)), upper = T)
annot[grepl('Treg|Tfh|cd4', annot$orig.majorCluster), ]$orig.celltype <- "CD4"
table(annot$orig.celltype)
table(annot$orig.majorCluster, annot$orig.celltype)
colnames(annot) <- sub("^PatientID$", "orig.Patient_ID", colnames(annot))
colnames(annot) <- sub("^location$", "orig.location", colnames(annot))
colnames(annot) <- sub("^treatment$", "orig.treatment", colnames(annot))
colnames(annot) <- sub("Batch.Set.ID", "orig.tissue", colnames(annot))
colnames(annot) <- sub("^LibPrepBy$", "orig.LibPrepBy", colnames(annot))
annot$orig.BL_TRT <- ifelse(annot$orig.treatment == "N", "Pre", "Post") # following Sade's name
annot$orig.tissue <- sub("Tumor", "T", annot$orig.tissue)
table(annot$orig.LibPrepBy, useNA = "always")
annot <- cbind(UniqueCell_ID = rownames(annot), annot)
annot$UniqueCell_ID <- as.character(annot$UniqueCell_ID)
str(annot)
sapply(annot[, grepl("orig", colnames(annot))], table)

save(annot, file = '../../li_MELANOMA_all_annotation.RData')
save(umi, file = '../../li_MELANOMA_all_umi.RData')

library(Seurat)
mycells <- CreateSeuratObject(raw.data = umi, meta.data = annot)
save(mycells, file = '../../li_MELANOMA_seurat.Robj')

# # download the repository: from https://bitbucket.org/tanaylab/li_et_al_cell_2018_melanoma_scrna/downloads/
# mv tanaylab-li_et_al_cell_2018_melanoma_scrna-db49e8dfcb3d/* ./
# rm -r tanaylab-li_et_al_cell_2018_melanoma_scrna-db49e8dfcb3d/
# wget http://www.wisdom.weizmann.ac.il/~lubling/Li2018/data_Li2018.tar.gz
# tar -xvzf data_Li2018.tar.gz
# wget http://www.wisdom.weizmann.ac.il/~lubling/Li2018/scrna_db_Li2018.tar.gz
# tar -xvzf scrna_db_Li2018.tar.gz
# # go to R5
# # Loading code and downloading required data files
# source("pipe.r")
#
# # Building the metacells
# build_metacells()
#
# # Generate figures
# generate_figs()
#
# # Reproduce the Guo et al 2018 (lung cancer scRNA data) analysis
# source("guo2018.r")

# cd ${ROOTDIR}/GSEs
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/suppl/GSE123139_T_cells_tcrb_v2.txt.gz
# gunzip GSE123139_T_cells_tcrb_v2.txt.gz
# # go to R
# edata <- data.table::fread('GSE123139_T_cells_tcrb_v2.txt')
# head(edata[, 1:40])
# cd ${ROOTDIR}/data_Li2018
# # go to R
# library(data.table)
# fnames <- paste0('umi.tab/', list.files('umi.tab'))
# edatal <- lapply(fnames, function(x){
#   tvar <- t(read.table(x))
#   data.frame(Cell = rownames(tvar), tvar, check.names = FALSE, stringsAsFactors = FALSE)
# })
# sum(sapply(edatal, ncol))
# sapply(edatal[c(1:3, 292:295)], function(x) tail(colnames(x)) )
# sapply(edatal[c(1:3, 292:295)], function(x) head(colnames(x)) )
# edata <- rbindlist(edatal)
# head(edata[, 1:15])
# unique(unlist(edata[, 1]))

#### Guo - Lung cancer #### ----------------------------------------------------
# Go to R
setwd('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/GSEs')
list.files()
### first annotation ###
annot <- read.table('GSE99254_family.soft.txt', header = 1)
annot <- taRifx::remove.factors(annot)
rownames(annot) = annot$UniqueCell_ID
annot$celltype = sub("\\.", "", sub("(^...).*", "\\1", annot$majorCluster))
annot$tissue = sub("(^.).*", "\\1", annot$sampleType)
table(annot$celltype)
table(annot$tissue)
head(annot)
table(annot$majorCluster)
colnames(annot)[-1] <- paste0('orig.', colnames(annot)[-1])
head(annot)
save(annot, file = '../guo_GSE99254_all_annotation.RData')

### now count data ###
# mat <- fread('GSE99254_NSCLC.TCell.S12346.TPM.txt')
mat <- fread('GSE99254_NSCLC.TCell.S12346.count.txt')
# fixing NAs
sum(is.na(mat$symbol))
mat$geneID[is.na(mat$symbol)]
head(mat[is.na(mat$symbol), 1:10])
mat <- mat[!is.na(mat$symbol), ]

head(mat[, 1:10])
tpm <- as.matrix(mat[, -c(1:2)])
rownames(tpm) <- mat$symbol
str(tpm)
head(tpm[, 1:10])
save(tpm, file = '../guo_GSE99254_all_tpm.RData')
all(rownames(annot) %in% colnames(tpm))

library("XLConnect")
clondata <- readWorksheetFromFile("/mnt/BioHome/ciro/simon/info/guo_GSE99254_tcr.xlsx", sheet=1, startRow = 2, endRow = 10204)

# trying to find missing cells
cat("In annotation:", sum(clondata$Cell.Name %in% rownames(annot)),
  "\nIn this extra annotation:", nrow(clondata),
  "\nTotal cells", nrow(annot), "\n")
clondata$Cell.Name[!clondata$Cell.Name %in% rownames(annot)]
rownames(annot)[grepl("20171208", rownames(annot))] # looks like some cells have T-CD4 pattern instead of sampleType
rownames(annot)[grepl("T\\-CD4", rownames(annot))]
head(annot[grepl("T\\-CD4", rownames(annot)), ])
head(annot) # they probably should have the sampleType in their name
clondata$Cell.Name[grepl("T\\-CD4", clondata$Cell.Name)] # no cells with the same pattern
head(clondata[!clondata$Cell.Name %in% rownames(annot), ])
# maybe the are in the data from Li
fname <- "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/data_Li2018/data/Guo2018_cell_clust_md.txt"
annot_from_li <- read.table(fname, stringsAsFactors = FALSE, header = 1, row.names = 1)
head(annot_from_li)
clondata$Cell.Name[!clondata$Cell.Name %in% rownames(annot_from_li)] # they are not!!
head(clondata[!clondata$Cell.Name %in% rownames(annot_from_li), ]) # they are not!!
table(clondata[!clondata$Cell.Name %in% rownames(annot_from_li), ]$Clone.Status) # they are not!!
# NOTE: giving an inferred name
tvar <- grepl("T\\-CD4", rownames(annot))
head(annot[tvar, ])
cellnames <- rownames(annot)[tvar]
cellnames <- gsub("T\\-CD4-", "", cellnames)
cellnames <- paste0(annot[tvar, "orig.sampleType"], cellnames)
sum(cellnames %in% clondata$Cell.Name) # Found them...
rownames(annot)[tvar] <- cellnames
all(rownames(annot) %in% colnames(tpm))
tpm <- tpm[, annot$UniqueCell_ID]
all(annot$UniqueCell_ID == colnames(tpm))
annot$UniqueCell_ID <- rownames(annot)
colnames(tpm) <- rownames(annot)

rownames(clondata) <- clondata$Cell.Name
clondata <- clondata[rownames(annot), ]
str(clondata)
head(clondata)
rownames(clondata) <- rownames(annot)
# myinterest <- c('Clone.ID', 'Clone.Status', 'Clone.Frequency', 'Identifier.Alpha1.', 'Identifier.Beta1.')
myinterest <- colnames(clondata)[-c(1:2)]
clones <- clondata[, myinterest]
clones[, 'alphav'] <- gsub("_.*", "", clones[, 'Identifier.Alpha1.'])
clones[, 'alphaj'] <- gsub(".*_", "", clones[, 'Identifier.Alpha1.'])
clones[, 'betav'] <- gsub("_.*", "", clones[, 'Identifier.Beta1.'])
clones[, 'betaj'] <- gsub(".*_", "", clones[, 'Identifier.Beta1.'])
head(clones)
clones[is.na(clones)] <- 'void'
void <- apply(clones, 2, unique)
orignames <- sapply(void, length) < 36
head(clones[, orignames])
colnames(clones)[orignames] <- paste0('orig.', colnames(clones)[orignames])
head(clones)
str(clones)
annot <- cbind(annot, clones[rownames(annot), ])
colnames(annot) = gsub("Clone.Frequency", "cloneFreq", colnames(annot))
colnames(annot) = gsub("Clone.Status", "cloneStatus", colnames(annot))
annot$orig.cloneFreq = paste0("F", annot$orig.cloneFreq)
annot$orig.cloneFreq = gsub("Fvoid", "void", annot$orig.cloneFreq)
head(annot)
save(annot, file = '../guo_GSE99254_all_annotation.RData')
save(tpm, file = '../guo_GSE99254_all_tpm.RData')

### Classifaying TREG comparment and activated Tregs ### -----------------------
# A) FOXP3+ CXCR5 and/or BCL6+
# B) FOXP3+ CXCR5- and BCL6-
# C) FOXP3-
source('~/scripts/functions/handy_functions.R')
genes <- getfound(c('FOXP3', 'CXCR5', 'BCL6'), rownames(tpm), v = T)
void <- add_gene_tag(lgenes = genes, annot = annot, mat = tpm, thresh = c(10, 10, 10), tag = c('tag', '+', '-'), prefix = FALSE, v = F)
head(void)
sapply(void, table)
void$orig.fcmarkers <- apply(void, 1, paste, collapse="")
void$orig.fctype <- void$tag_FOXP3
void[void$orig.fcmarkers == "FOXP3+CXCR5-BCL6-", 'orig.fctype'] <- 'TREG'
void[void$orig.fcmarkers %in% c("FOXP3+CXCR5-BCL6+", "FOXP3+CXCR5+BCL6-", "FOXP3+CXCR5+BCL6+"), 'orig.fctype'] <- 'TFR'
table(void$orig.fctype)
table(void$orig.fcmarkers)
void <- cbind(void, add_gene_tag(lgenes = 'TNFRSF9', annot = annot, mat = tpm, thresh = 5, tag = c('tag', '+', '-'), prefix = FALSE, v = F))
void$orig.fca <- paste0(void$orig.fctype, '_', void$tag_TNFRSF9)
table(void$orig.fctype, void$tag_TNFRSF9)
table(void$orig.fca)
annot <- cbind(annot, void[rownames(annot), grepl("orig", colnames(void))])
head(annot)

#### Sade - Melanoma #### ------------------------------------------------------
cd ${ROOTDIR}/GSEs
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz
gunzip GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz
gunzip GSE120575_patient_ID_single_cells.txt.gz
grep "Sample " GSE120575_patient_ID_single_cells.txt | tail
# go to R
# annot <- read.table('~/large/divya/raw/GSE120575_patient_ID_single_cells.txt', header = 1, sep = '\t', check.names = F)
annot <- read.table('GSE120575_patient_ID_single_cells.txt', header = 1, skip = 19, sep = '\t', check.names = F, nrows = 16291, stringsAsFactors = FALSE)
rownames(annot) = annot$title
annot <- annot[, !apply(annot, 2, function(x) sum(is.na(x))) > nrow(annot) / 20]
colnames(annot) <- gsub("characteristics: ", "", colnames(annot))
colnames(annot) <- gsub("patinet ID \\(Pre=baseline; Post= on treatment\\)", "PID_BL_TRT", colnames(annot))
colnames(annot) <- gsub("title", "UniqueCell_ID", colnames(annot))
annot$BL_TRT <- gsub("_P.*$", "", annot$PID_BL_TRT)
annot$Patient <- gsub("^.*_P", "P", annot$PID_BL_TRT)
head(annot)

apply(annot[, -c(1, 2)], 2, function(x) unique(x))
tvar <- apply(annot[, -c(1, 2)], 2, function(x) length(unique(x)))
tmp <- names(tvar[tvar > 1])
colnames(annot) <- gsub(paste0("(", paste0(tmp, collapse = "|"), ")"), "orig.\\1", colnames(annot))

# adding cell cluster identity
clusters <- read.csv('/mnt/BioHome/ciro/simon/info/sade_cluster_cell_name_mapping.csv', row.names = 1, stringsAsFactors = F)
all(rownames(clusters) %in% rownames(annot))
annot$cluster_number <- paste0("G", clusters[rownames(annot), 1])
table(annot$cluster_number)
ident <- c('b_cells', 'plasma', 'mono_macro', 'dendritic', 'lymphocytes', 'exCD8', 'Treg', 'cytotoxic_lymph', 'ex_hsCD8', 'Tmem', 'exlymphocytes_CC')
names(ident) = paste0("G", 1:length(ident))
annot$orig.majorCluster <- ident[annot$cluster_number]
table(annot$orig.majorCluster)
annot$orig.celltype <- annot$orig.majorCluster
annot[annot$orig.majorCluster %in% c('Treg', 'Tmem', 'mono_macro', 'dendritic')[1:2], 'orig.celltype'] <- 'CD4'
annot[annot$orig.majorCluster %in% c('exCD8', 'ex_hsCD8'), 'orig.celltype'] <- 'CD8'
table(annot$orig.celltype, annot$orig.majorCluster)
annot$orig.tissue <- "T"
head(annot)
save(annot, file = '../sade_GSE120575_all_annotation.RData')

### now count data ###
# mat2 <- data.table::fread('~/large/divya/raw/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_chopped_2nd_line.txt')
# mathead2 <- as.matrix(read.table('~/large/divya/raw/GSE120575_header.txt'))[1, ]
mat <- data.table::fread('GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt')
mathead <- scan('GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt', what = character(), nline = 1, sep = "\t")
mathead[1] <- "X"
head(mathead)
# all(mathead == mathead2)
head(mat[, 1:10])
head(mat[, (ncol(mat) - 10):ncol(mat)])
mat <- mat[,V16293:=NULL]
length(mathead)
dim(mat)
colnames(mat) <- mathead
# fixing NAs
sum(is.na(mat$X))
mat$geneID[is.na(mat$X)]
head(mat[is.na(mat$X), 1:10])
mat <- mat[!is.na(mat$X), ]

head(mat[, 1:10])
tpm <- as.matrix(mat[, -c(1)])
rownames(tpm) <- mat$X
str(tpm)
head(tpm[, 1:10])
colSums(tpm[, 1:10])
summary(tpm[, 1:10])
# all(colSums(tpm1[, 1:300]) == colSums(tpm[, 1:300]))
save(tpm, file = '../sade_GSE120575_all_tpm.RData')

#### Savas - Breast #### -------------------------------------------------------
mkdir ${ROOTDIR}/savas_GSE110686
cd ${ROOTDIR}/savas_GSE110686
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110686/suppl/GSE110686_tils20%2B32_barcodes.tsv.gz
gunzip GSE110686_tils+32_barcodes.tsv.gz; mv GSE110686_tils20+32_barcodes.tsv barcodes.tsv
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110686/suppl/GSE110686_tils20%2B32_genes.tsv.gz
gunzip GSE110686_tils+32_genes.tsv.gz; mv GSE110686_tils20+32_genes.tsv genes.tsv
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110686/suppl/GSE110686_tils20%2B32_matrix.mtx.gz
gunzip GSE110686_tils+32_matrix.mtx.gz; mv GSE110686_tils20+32_matrix.mtx matrix.mtx
# go to R
umi <- Seurat::Read10X(data.dir = './')
umi <- as.matrix(umi)
dim(umi)
head(umi[, 1:10])
summary(umi[, 1:10])
tvar <- sub('[A-z]+\\-', '', colnames(umi), perl = TRUE)
table(tvar); sum(table(tvar)[1:2])
# we see that donor 1 is the first two libraries and the 2nd donor is the 3rd
# as described in the GEO
annot <- data.frame(row.names = colnames(umi),
  orig.Patient_ID = rep(paste0("D", c(1:2,2)), table(tvar)),
  orig.lib = paste0("L", tvar), stringsAsFactors = FALSE)
table(annot)
head(annot)
clusters = readRDS('large/simon/raw/savas_GSE110686/cluster.rds')
head(clusters); table(clusters)
annot$orig.majorCluster <- as.character(clusters[rownames(annot)])
annot$orig.majorCluster[is.na(annot$orig.majorCluster)] <- 'filtered'
table(annot$orig.majorCluster)
annot$orig.celltype <- sub("(^...).*", "\\1", annot$orig.majorCluster)
table(annot$orig.celltype)
annot$orig.tissue <- "T"
annot <- cbind(UniqueCell_ID = rownames(annot), annot)

save(umi, file = '../savas_GSE110686_all_umi.RData')
save(annot, file = '../savas_GSE110686_all_annotation.RData')

#### Yosy - Carcinoma ##### ----------------------------------------------------
cd /mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw
if [ ! -d yost ]; then mkdir yost; fi
cd yost
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fbcc%5FscRNA%5Fcounts%2Etxt%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fbcc%5Ftcr%2Etxt%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fscc%5FscRNA%5Fcounts%2Etxt%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fscc%5Ftcr%2Etxt%2Egz
# Extra data sent by Yost herself through email

for fname in `ls`; do
  echo ${fname}
  gunzip ${fname}
done

head -n 1 GSE123813_bcc_scRNA_counts.txt > GSE123813_bcc_scRNA_counts_header.txt

## Go to R
mat <- data.table::fread('GSE123813_bcc_scRNA_counts.txt')
math <- data.table::fread('GSE123813_bcc_scRNA_counts_header.txt')
head(colnames(math))
dim(mat)
mat[1:20, 1:5]
umi <- as.matrix(mat[, -1])
rownames(umi) <- unlist(mat[, 1])
colnames(umi) <- colnames(math)
umi[1:20, 1:2]
str(umi)

preannot <- read.csv('bcc_all_metadata.txt', row.names = 1, header = 1, sep = "\t", check.names = F, stringsAsFactors = F)
# preannot <- read.csv('scc_metadata.txt', row.names = 1, header = 1, sep = "\t", check.names = F, stringsAsFactors = F)
head(preannot)
sapply(preannot, table)
table(preannot$cluster, preannot$sort)
table(preannot$patient, preannot$sort)
annot <- preannot
colnames(annot) <- sub("cluster", "orig.majorCluster", colnames(annot))
colnames(annot) <- sub("treatment", "orig.BL_TRT", colnames(annot))
colnames(annot) <- sub("patient", "orig.Patient", colnames(annot))
annot$orig.BL_TRT <- ifelse(annot$orig.BL_TRT == "pre", "Pre", "Post") # following Sade's name

annot$orig.celltype <- annot$orig.majorCluster
annot[annot$orig.majorCluster %in% c('CD4_T_cells', 'Tregs'), "orig.celltype"] <- 'CD4'
table(annot$orig.celltype)
sapply(annot, table)

save(umi, file = '../yost_GSE123813_bcc_umi.RData')
save(annot, file = '../yost_GSE123813_bcc_annotation.RData')
