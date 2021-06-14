# TFR 2021
This repository contains scripts utilised to create figures and analyse
our TFR paper related data.

### downloads.R
Downloading all the data sets. You can find how the expression matrices and annotations were set.

### guo_treg_classification.R
Further classification of cells to find the clonotype sharing between cell types, i.e., TREGoFOXP3N or TREGoTFR

### dea_setting_sade_subsets.sh
Wrapper for comparison between pre and post PD1 therapy from the TFR cells identified in Sade-Feldman data.
This calls the script that performs the differential expression /mnt/BioHome/ciro/scripts/scdea/scDEA.R
which selects the tool from /mnt/BioHome/ciro/scripts/scdea/scDEA_methods.R.

### integration_preprocesing.R
Filters applied to the 9 studies. Just after the function you can find the threshold that where used.

### integration_seurat.R
Integrated analysis of the 9 datasets following the instructions given in the
online [tutorial](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html).

### jobConstMon_setting_treg_sharing.sh
Wrapper for /mnt/BioHome/ciro/scripts/monocle/monocle3_v1.R that performs Monocle 3 analysis. This was done with the genes specified by the parameter GENESETF.

### Signature vs monocle's 1st component
guo_monocleVSsignature_plotting.R
guo_signature.R

### integration_seurat_tfr.R
Integrated analysis was filtered for SmartSeq2 studies and only TFR cells were reintegrated

### SiEs05_supp_figure.R
It includes figures from the T/LN and TCR libraries

### SiEs09_workflow.rmd
Workflow used for the murine data. In here, other SiEs09 YAML/R files are covered.

