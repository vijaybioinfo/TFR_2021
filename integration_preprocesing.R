#!/usr/bin/R5

####################
# Treg integration #
####################

# this script is to process and set cutoffs for data sets involved in the
# integration of T cells populations

.libPaths('~/R/newer_packs_library/3.5/')
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/simon'
datadir <- paste0(root, '/raw')

#### Functions #### ------------------------------------------------------------
source('/mnt/BioHome/ciro/simon/scripts/integration_preprocesing_functions.R')

#### ######### #### ------------------------------------------------------------
sets_info <- list(
  guo = list(ctype = 'lung', dtype = 'tpm', norm = TRUE, gthres = c(10, 10, 10)),
  zheng = list(ctype = 'liver', dtype = 'tpm', norm = TRUE, gthres = c(10, 10, 10)),
  zhang = list(ctype = 'colorectal', dtype = 'tpm', norm = TRUE, gthres = c(10, 10, 10)),
  sade = list(ctype = 'melanoma', dtype = 'tpm', norm = TRUE, gthres = c(10, 10, 10)),
  lambrechts = list(ctype = 'lung', dtype = 'umi', norm = TRUE, gthres = c(1, 1, 1)),
  puram = list(ctype = 'head_neck', dtype = 'tpm', norm = TRUE, gthres = c(10, 10, 10)),
  arnon = list(ctype = 'melanoma', dtype = 'tpm', norm = TRUE, gthres = c(1, 1, 1)),
  savas = list(ctype = 'breast', dtype = 'umi', norm = TRUE, gthres = c(1, 1, 1)),
  li = list(ctype = 'melanoma', dtype = 'umi', norm = TRUE, gthres = c(1, 1, 1)),
  tirosh = list(ctype = 'oligodendroglioma', dtype = 'tpm', norm = TRUE, gthres = c(4, 0, 0))
)[1:9]
submark = c('CD8B', 'CD8A', 'CD4') # genes to filter cells by
ct_markers = c('BCL6', 'CXCR5', 'FOXP3') # genes to classify cells by
# sets where filtering will be applied
cleaning = list(sade = c('all', 10), lambrechts = c('all', 1), savas = c('CD4+ RGCC+', 'CD4+ others', 1), puram = c('Tcell', 10))
cname = c('orig.celltype', 'orig.majorCluster') # major groups where to look
selected = c('sade', 'puram', 'arnon', 'li', 'lambrechts')[-c(1:5)] # sets in which the next groups will be extracted from
only_them = c('orig.majorCluster', 'Tcell', 'T.CD4', 'Treg', 'Tfh', 'em-cd4', 'dysf-cd4', '7_Treg')

summ <- t(sapply(sets_info, function(x){
  names(x[['gthres']]) <- c('BCL6', 'CXCR5', 'FOXP3')
  # round(x[['gthres']], 3)
  x[['gthres']]
}))
summ <- cbind(set = names(sets_info),
  rbindlist(lapply(sets_info, function(x) data.frame(x[names(x) != 'gthres']) )),
  summ
)
summ <- remove.factors(data.frame(summ))
summ

### Getting datasets annotations and expression data ### -----------------------
annotf <- list.files(datadir); annotf <- annotf[grepl('_annot', annotf)]
# exluding Guo without extra fieldss, lambrechts log2cpm and others
annotf <- annotf[!grepl('EMT|gomes|extended|dum', annotf)]
ctsf <- list.files(datadir); ctsf <- ctsf[grepl('tpm\\.|umi\\.', ctsf)]
# Naming the files
names(ctsf) <- sub('(^.*)_[A-Z].*', '\\1', ctsf)
names(annotf) <- sub('(^.*)_[A-Z].*', '\\1', annotf)
annotf <- annotf[names(sets_info)]
ctsf <- ctsf[names(sets_info)]

### Loading data ### -----------------------------------------------------------
cat('Loading', length(annotf), 'data sets:', commas(names(annotf), Inf), '\n')
annot_list <- lapply(paste0(datadir, '/', annotf), theObjectSavedIn)
names(annot_list) <- names(annotf)
# str(annot_list)
sapply(annot_list, function(x) colnames(x)[grepl('orig', colnames(x))] )
edata_listi <- lapply(paste0(datadir, '/', ctsf), theObjectSavedIn)
names(edata_listi) <- names(ctsf)
str(edata_listi)

# checking average count content
# Puram is in the lower range for a TPM - modified during download
# Arnon said they divided by 10 - modified during download
# edata_listi[['puram']] <- edata_listi[['puram']] * 10
# edata_listi[['arnon']] <- edata_listi[['arnon']] * 10
set.seed(27)
summ$mean_nCTS <- sapply(edata_listi, function(x){
  mean(colSums(x[, sample(1:ncol(x), 500)]) )
})

check_treament(listannot = annot_list)
# arnon has per tumour only pre or post but not both
# li different types of treament and a naive group - corrected in download
# annot_list[['li']][, 'orig.BL_TRT'] <- ifelse(annot_list[['li']][, 'orig.treatment'] == "N", "Pre", "Post")
# annot_list[['arnon']][, 'orig.BL_TRT'] <- ifelse(annot_list[['arnon']][, 'orig.treatment'] == "treatment.naive", "Pre", "Post")
check_treament(listannot = annot_list)
summ$treatment <- !sapply(check_treament(listannot = annot_list), is.null)
summ[, -c(4:7)]

### Preprocess data ### --------------------------------------------------------
maxchar <- max(nchar(sapply(names(edata_listi), head, 1)))
genenames <- unique(unlist(sapply(edata_listi, rownames)))
head(genenames); length(genenames)
ogenes <- sapply(edata_listi, function(x) genenames %in% rownames(x) )
apply(ogenes, 2, sum)
cgenes <- genenames[apply(ogenes, 1, all)]
head(cgenes); length(cgenes)
edata_list <- sapply(names(edata_listi), function(dname){
  cat(addspaces(dname, maxchar), '---------------------------------------\n')
  cts <- edata_listi[[dname]]
  if(sets_info[[dname]][['dtype']] == 'umi'){
    cts <- sweep(cts, 2, colSums(cts), '/') * 1e+6
    sets_info[[dname]][['dtype']] <<- 'cpm'
  }
  return(cts)
})
set.seed(27)
summ$mean_nCTSnorm <- sapply(edata_list, function(x){
  mean(colSums(x[, sample(1:ncol(x), 500)]) )
})

### Checking markers and groups ### --------------------------------------------
setwdc(paste0(root, '/results/integration/markers_preview/prefilter'), recursive = T)
markers <- read.csv('/mnt/BioHome/ciro/simon/info/markers.csv', stringsAsFactors = F)
markers
mymarkers <- markers[, 1]
mymarkers <- unique(c('IL2RA', 'TNFRSF9', 'TNFRSF18', 'DUSP4', 'TNFRSF14', mymarkers))
mymarkers <- c('BCL6', 'CXCR5', 'FOXP3', 'CD4', 'CD8A', 'CD8B', 'IL2RA')

tvar <- data.frame(dcast.data.table(setDT(melt(lapply(annot_list, function(x) table(x$orig.celltype) ))), Var1 ~ L1))
write.table(tvar, file = 'celltypes_per_set.csv', sep = ',', row.names = FALSE, col.names = TRUE)
sname <- 1:length(annot_list)
void <- explots(lsets = list(annot = annot_list[sname], edata = edata_list[sname], info = sets_info[sname]),
  mymarkers = mymarkers, submark = submark,
  cname = cname, cleaning = cleaning)
saveWorkbook(void, file = "a1_clusters_celltypes.xlsx" , overwrite = FALSE)

### Subsetting to CD4 tumour ### -----------------------------------------------
setwdc(paste0(root, '/results/integration/markers_preview/midfilter'))
cells <- list(c('orig.tissue', 'T'), c('orig.celltype', 'CD4'), c('orig.BL_TRT', 'Pre'))
sapply(annot_list, function(x) table(x[x$orig.celltype == 'CD4', 'orig.majorCluster']) ) # virtually the excel sheet
sapply(annot_list[selected], function(x) table(x[x$orig.celltype == 'CD4', 'orig.majorCluster']) )
annot_list_ss <- list()
for(dname in names(annotf)){
  cat(addspaces(dname, maxchar), '---------------------------------------\n'); dir.create(dname)
  annot <-  annot_list[[dname]]
  annot$orig.set <- dname
  cfilt <- cells[sapply(cells, head, 1) %in% colnames(annot)]
  if(dname %in% selected) cfilt[[length(cfilt) + 1]] <- only_them
  if(length(cfilt)){
    annot <- annot[getsubset(cfilt, annot, v = TRUE), ]
  }else{
    cat('Taking all cells in dataset\n')
  }
  fname <- paste0(dname, "/", sub("orig\\.", "", cname[1]), "_markers.pdf")
  if(!file.exists(fname)){
    cts <- edata_list[[dname]][, rownames(annot)]
    tvar <- ifelse(dname %in% names(cleaning), as.numeric(tail(cleaning[[dname]], 1)), 0)
    get_densities(mat = cts, log2t = sets_info[[dname]][['norm']], genes = submark, subname = fname, cuof = tvar)
  }
  annot_list_ss[[dname]] <- annot
}; names(annot_list_ss) <- names(annot_list)

### Filtering out to CD8+ ### --------------------------------------------------
# tannot_list_ss <- annot_list_ss
# annot_list_ss <- tannot_list_ss
for(dname in names(annot_list_ss)){
  cat(addspaces(dname, maxchar), '---------------------------------------\n')
  dir.create(dname)
  annot <-  annot_list_ss[[dname]]
  if(dname %in% names(cleaning)){
    i <- 1 # if(dname == 'lambrechts') i <- 1:2 else
    cat('Filtering ', commas(submark[i]), '+\n', sep = ''); cleaner <- cleaning[[dname]]
    cts <- edata_list[[dname]][, rownames(annot)]
    if('all' %in% cleaner) cleaner <- c(unique(annot[, 'orig.majorCluster']), tail(cleaner, 1))
    tvar <- unlist(lapply(head(cleaner, -1), function(x){
      cat('>', x, '\n')
      rnames <- initrnames <- getsubset(c('orig.majorCluster', x), annot)
      void <- add_gene_tag(lgenes = submark, annot = annot[rnames, ], mat = cts[, rnames], thresh = as.numeric(tail(cleaner, 1)))
      cat('Threshold:', tail(cleaner, 1)); print(table(void[, paste0('tag_', submark[i])]))
      cat("."); get_densities(mat = cts[, rnames], log2t = sets_info[[dname]][['norm']], genes = submark,
        subname = paste0(dname, "/", x, "_bf"), cuof = as.numeric(tail(cleaner, 1)))
      rnames <- getsubset(gettag(submark[i], "n_"), void)
      cat("|"); get_densities(mat = cts[, rnames], log2t = sets_info[[dname]][['norm']], genes = submark,
        subname = paste0(dname, "/", x), cuof = as.numeric(tail(cleaner, 1)))
      cat('\n'); return(getsubset(gettag(submark[i], "p_"), void, op = ifelse(dname == 'no_set', 'or', 'and')))
    }))
    tvar <- rownames(annot)[!rownames(annot) %in% tvar]
  }else tvar <- rownames(annot)
  annot_list_ss[[dname]] <- annot[tvar, ]
}
sapply(annot_list_ss, function(x) unique(x[, 'orig.tissue']) )
sapply(annot_list_ss, function(x) unique(x[, 'orig.celltype']) )
sapply(annot_list_ss, function(x) table(x[, 'orig.majorCluster']) )
sapply(tannot_list_ss, function(x) table(x[, 'orig.majorCluster']) )
t(sapply(tannot_list_ss, nrow ))
# guo zheng zhang sade lambrechts puram arnon savas   li
# 2783  1163  2473  608      11130  1237   856  3726 4904
t(sapply(annot_list_ss, nrow )) # percentage of cells kept
# guo zheng zhang sade lambrechts puram arnon savas   li
# 2783  1163  2473  599      10447   995   856  3635 4904
t(sapply(annot_list_ss, nrow )) / t(sapply(annot_list, nrow )) * 100
check_treament(listannot = annot_list_ss)
edata_list <- sapply(names(edata_list), function(x){
  cts <- edata_list[[x]]
  cts[, rownames(annot_list_ss[[x]])]
})
gc()

### Checking markers on filtered ones ### --------------------------------------
setwdc(paste0(root, '/results/integration/markers_preview/postfilter'))
void <- explots(lsets = list(annot = annot_list_ss, edata = edata_list, info = sets_info),
  mymarkers = mymarkers, submark = submark, cleaning = cleaning)
saveWorkbook(void, file = "a1_clusters_celltypes.xlsx" , overwrite = FALSE)

### Plotting thresholds ### ----------------------------------------------------
setwdc(paste0(root, '/results/integration/markers_treg'))
annot_list_tags <- list()
for(dname in names(edata_list)){ # dname = names(edata_list)[8]
  cat(addspaces(dname, maxchar), '---------------------------------------\n')
  cts <- edata_list[[dname]]; annot <- annot_list_ss[[dname]]
  get_densities(mat = cts, log2t = sets_info[[dname]][['norm']], genes = ct_markers, subname = dname, cuof = sets_info[[dname]][['gthres']])
  tvar <- c(ct_markers, 'CD4'); tmp <- c(sets_info[[dname]][['gthres']], 0)
  void <- add_gene_tag(lgenes = tvar, annot = annot, mat = cts, thresh = tmp, tag = c('tag', '+', '-'), prefix = FALSE)
  tcst <- cts[, rownames(void)[void$tag_FOXP3 == 'FOXP3+']]
  get_densities(mat = tcst, log2t = sets_info[[dname]][['norm']], genes = ct_markers, subname = paste0(dname, '_foxp3'), cuof = sets_info[[dname]][['gthres']])
  void$orig.fcmarkers <- apply(void[, 3:1], 1, paste, collapse="")
  mytab <- sapply((length(tvar) - 1):1, function(x){
    apply(void[, c(length(tvar), x)], 1, paste, collapse="")
  })
  void <- cbind(void, mytab)
  void$orig.allmarkers <- apply(void[, 4:1], 1, paste, collapse="")
  void$tag_ct[void$orig.fcmarkers == 'FOXP3+CXCR5-BCL6-'] <- 'TREG'
  void$tag_ct[void$orig.fcmarkers != 'FOXP3+CXCR5-BCL6-' & void$tag_FOXP3 == 'FOXP3+'] <- 'TFR'
  void$tag_ct[is.na(void$tag_ct)] <- 'NCT'
  annot <- cbind_repcol(void, annot)
  annot_list_tags[[dname]] <- annot
}

mytab <- data.frame(rbindlist(
  lapply(c(names(void), 'orig.celltype'), function(cname){
    allfacts <- unique(unlist(lapply(annot_list_tags, function(x)
      names(table(x[, cname]))
    )))
    tvar <- sapply(annot_list_tags, function(x)
      table(factor(x[, cname], allfacts), useNA = 'always')
    )
    data.frame(cbind(rownames(tvar), tvar), stringsAsFactors = FALSE)
  })
), stringsAsFactors = FALSE)
mytab <- mytab[complete.cases(mytab), ]
rownames(mytab) <- mytab[, 1]; mytab <- mytab[, -1]
head(mytab); tail(mytab)
x <- t(summ[, -1]); colnames(x) <- summ[, 1]
rownames(x)[rownames(x) %in% ct_markers] <- paste0('Threshold_', rownames(x)[rownames(x) %in% ct_markers])
tmp <- mat_names('Threshold_CD4', colnames(x)); tmp['Threshold_CD4', ] <- 0
tvar <- data.frame(rbindlist(lapply(summ$set, function(x){
  y <- ifelse(is.null(tail(cleaning[[x]], 1)), NA, tail(cleaning[[x]], 1))
  z <- head(cleaning[[x]], -1); if(is.null(z)) z <- NA
  if(is.na(y)) a <- NA else a <- submark[1]; #if(x == 'lambrechts') a <- submark[-3]
  data.frame(FiltOutThres = y, FiltGenes = commas(a), FiltGroups = commas(z, Inf), stringsAsFactors = FALSE)
})), stringsAsFactors = FALSE); rownames(tvar) <- summ$set
x <- rbind(x, tmp); tvar <- rbind(x, t(tvar))
mytab <- rbind(tvar, mytab)
head(mytab, 15); tail(mytab, 15)

annot <- remove.factors(data.frame(rbindlist(annot_list_tags, fill = TRUE), stringsAsFactors = F))
rownames(annot) <- annot[, 'UniqueCell_ID']
# annot <- annot[, sapply(annot, function(x) all(!is.na(x)) )]
head(annot[, head(1:ncol(annot), 20)])
table(annot$tag_ct) / nrow(annot) * 100
table(annot$tag_FOXP3) / nrow(annot) * 100

### Plotting proportions ### ---------------------------------------------------
setwdc(paste0(root, '/results/integration/labels_proportions5'))
# g <- list(get_props(metadata = annot, group_by = 'tag_FOXP3', resolution = 'orig.set', return_plot = TRUE, norm_type = 'n', v = T)[[1]])
# g[[2]] <- get_props(metadata = annot, group_by = 'orig.fcmarkers', resolution = 'orig.set', return_plot = TRUE, norm_type = 'n', v = T)[[1]]
g <- list(make_bplot(annot))
annot$gene_tag <- annot$tag_FOXP3
pdf(paste0('1a_FOXP3_across_sets.pdf'), height = 8, width = 16)
plot_grid(plotlist = g, 'Tumour CD4 T cells')
graphics.off()

annot_f3 <- annot[annot$tag_FOXP3 == 'FOXP3+', ]
g <- list(get_props(metadata = annot_f3, group_by = 'tag_ct', resolution = 'orig.set', return_plot = TRUE, norm_type = 'n', v = T)[[2]])
g[[2]] <- get_props(metadata = annot_f3, group_by = 'orig.fcmarkers', resolution = 'orig.set', return_plot = TRUE, norm_type = 'n', v = T)[[2]]
pdf(paste0('1a_FOXP3p_across_sets.pdf'), height = 8, width = 16)
plot_grid(plotlist = g)
graphics.off()

### Writing proportions ### ----------------------------------------------------
tvar <- cbind(table(annot$orig.set, annot$tag_FOXP3), table(annot_f3$orig.set, annot_f3$tag_ct))
write.csv(mytab, file = 'proportions_all.csv')
write.csv(tvar, file = 'proportions.csv')

### Prepearing for integration ### ---------------------------------------------
setwdc(paste0(root, '/results/integration/data'))
write.csv(mytab, file = 'summary.csv', row.names = F)
save(annot, file = 'metadata.RData')
save(annot_list_tags, file = 'annot_list_tags.RData')
rm(edata_list); gc()

# expression matrices
edata_list <- lapply(names(annot_list_tags), function(x){
   edata_listi[[x]][, rownames(annot_list_tags[[x]])]
}); names(edata_list) <- names(annot_list_tags)
save(edata_list, file = 'edata_list.RData')
