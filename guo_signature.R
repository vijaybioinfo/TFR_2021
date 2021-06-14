#!/usr/bin/R

# This script is to create visualise the results from monocle and signature
# analyses with analysis of data from Guo et al

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

#### loading files #### --------------------------------------------------------
# Colours
degs_path <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/gtumour/comprs/treg_sharing/'
colfile <- '/mnt/BioHome/ciro/simon/info/groupColours_art_ct.csv'
gr.cols <- read.csv(colfile, row.names = 1, stringsAsFactors = F)
fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/treg_classification.csv'
annot <- read.csv(fname, stringsAsFactors = F, check.names = F, row.names = 1)
grps <- c('cTFR', 'cTREGoTFR', 'TREG')
tannot <- annot[getsubset(c('orig.cell_set', grps), annot, v = T), ]
tannot <- tannot[getorder(grps, tannot, 'orig.cell_set'), ]
table(tannot$orig.ctact)
tpm <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/raw/guo_GSE99254_all_tpm.RData')

### Getting the signature ### -----------
# get tannot from pca code
setwd('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/')
fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/tfr_signature_genes.csv'
mytab <- readfile(fname, row.names = 1, stringsAsFactors = F)
myname <- "TFR"

# fname <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon/results/mast/guo_lung/combined_dea_up_down.csv'
# mytab <- readfile(fname, row.names = 1, stringsAsFactors = F)
# mytab <- mytab[mytab$up, ]
# myname <- "TFR_combined"

gens <- rownames(mytab)
cname <- 'cellname'
tannot <- creplace(tannot, cin = 'orig.cell_set', cfrom = 'orig.ctact', subs = 'TREG', newname = cname)
grps <- names(table(tannot$cellname))

source('/mnt/BioHome/ciro/scripts/functions/signature.R')
void <- get_signature(
  mat = tpm,
  genes = gens,
  annot = tannot,
  catg = cname,
  sname = myname,
  smethods = 1,
  feat2see = c(cname, rev(grps)),
  grpct = NULL,
  filtercells = FALSE,
  gcouls = gr.cols,
  path = './',
  reversing = FALSE,
  v = TRUE
)
