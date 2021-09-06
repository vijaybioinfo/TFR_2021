#!/bin/bash
######################
# MONOCLE 3 WRAPPER ##
######################

## DESCRIPTION ##
# Wrapper for Monocle 3 analysis

COMDIR=/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon
#### parameters #### -----------------------------------------------------------
MYCELLSF=${COMDIR}/results/seurat/guo_lung/requests/dens_props/treg_sharing.csv # metadata, CSV or Seurat object
RAWDATA=${COMDIR}/raw/guo_GSE99254_all_tpm.RData # 10X folder PATH (before outs), CSV or Seurat object (same as MYCELLSF)
OUTDIR=${COMDIR}/results/monocle/treg_sharing
CATG=orig.cell_set~orig.treg_fca~orig.fca~orig.majorCluster~orig.sampleType~orig.fctype~orig.fcmarkers # colunm names to visualize in meta data or annotation file
GCOLOUR=/mnt/BioHome/ciro/simon/info/groupColours_art_ct.csv # header: group,colour
# specific paramaters #
REDTYPE=DDRTree # type of dimReduction, options: "DDRTree", "ICA", "tSNE", "UMAP"
SEUDISPS=FALSE # use seurat HVGs
PCS_COMP=20 # number of PCs to use
GENESETF=~/simon/large/results/mast/gtumour/comprs/treg_sharing/cTREGoTFRvscTREGoFOXP3N/degs_cTREGoTFR.RData
CTSTYPE=tpm # type of data
FILTCELLS=TRUE # filter cells according to distribution
#### parameters #### -----------------------------------------------------------

# ## Now on Treg, Tfrs cells that are clonaly expanded
# MYCELLSF=${COMDIR}/results/seurat/guo_lung/requests/dens_props/tumour_cd4.csv
# OUTDIR=${COMDIR}/results/monocle/treg_tfr_clonal
# CATGFILT=orig.cell_set~cTFR~cTREGoTFR~cTREGoFOXP3N

# ## Now we add the not clonaly expanded
# OUTDIR=${COMDIR}/results/monocle/treg_tfr
# CATGFILT=orig.cell_set~cTFR~cTREGoTFR~cTREGoFOXP3N~TREG

## Now we add the not clonaly expanded but no cTREGoFOXP3N
GENESETF=~/simon/large/results/mast/gtumour/comprs/treg_sharing/cTREGoTFRvscTREGoFOXP3N/degs.RData
OUTDIR=${COMDIR}/results/monocle/treg_tfr2
CATGFILT=orig.cell_set~cTFR~cTREGoTFR~TREG # filtering for these cells
PCS_COMP=14

. jobConstMon.sh
