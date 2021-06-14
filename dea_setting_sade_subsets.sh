#!/bin/bash
###############################
# Single-Cell DEA: SCDE/MAST ##
###############################

# Wrapper for script scde/mast.R
# run this script as: sh jobConstDEA.sh

# Call this script as:
# . /path/to/script/thisScript.sh

COMDIR=/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon
##### Default parameters ####
# directories need / in the end except OUTDIR
DEA=scde # type of differental expression analysis
DEA=mastlog2 # type of differental expression analysis
MYCELLSF=${COMDIR}/results/integration/seurat_9sets/sade_proportions10/pid_therapy_categories_annotation.RData # metadata, CSV or Seurat object
RAWDIR=${COMDIR}/raw/sade_GSE120575_all_tpm_bk.RData # 10X folder PATH (before outs), CSV or Seurat object (same as MYCELLSF)
OUTDIR=${COMDIR}/results/mast/sade # dont put "/" at the end
GROUP1=Post; GROUP2=Pre
CATGS=orig.BL_TRT # colunm name in meta data or annotation file, example: "res.0.8", ignored if COMPRS is a file
CPMCONV=log2tpm # log2cpm calculates CPM for the samples and logs it; log2 just logs them
#####

# ### CD4 ###
# SUBNAME=CD4 # significative name for comparison
# COLOURSF=/mnt/BioHome/ciro/simon/info/groupColours_art_ct.csv # header: group,colour

# ### FOXP3 ###
# CATGFILT="c('tag_FOXP3', 'FOXP3+')"
# SUBNAME=FOXP3

# ### Treg ###
# CATGFILT="c('tag_ct', 'TREG')"
# SUBNAME=TREG

### TFR ###
SUBNAME=TFR
CATGFILT="c('tag_ct', 'TFR')"

### filter for PD1 therapy
SUBNAME=${SUBNAME}_pd1_10TPM
CATGFILT="list(c('tag_ct', 'TFR'), c('orig.therapy', 'anti-PD1'), c('orig.Patient', '-P1', '-P6'))"

OVR=IGN # ignore defaults for job parameters
MEM=8gb
PPN=1

. /mnt/BioHome/ciro/scripts/scdea/jobConstDEA.sh
