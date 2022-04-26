### functions from snpnet package that are used for snpboost and are not adapted

cat_or_zcat <- function(filename, configs=list(zstdcat.path='zstdcat', zcat.path='zcat')){
  if(stringr::str_ends(basename(filename), '.zst')){
    return(configs[['zstdcat.path']])
  }else if(stringr::str_ends(basename(filename), '.gz')){
    return(configs[['zcat.path']])
  }else{
    return('cat')
  }
}

readIDsFromPsam <- function(psam){
  FID <- IID <- NULL  # to deal with "no visible binding for global variable"
  df <- data.table::fread(psam) %>%
    dplyr::rename('FID' = '#FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='_'))
  df$ID
}


####try to import from snpnet_compact
readPheMaster <- function(phenotype.file, psam.ids, family, covariates, phenotype, status, split.col, configs){
  
  sort_order <- . <- ID <- NULL  # to deal with "no visible binding for global variable"
  
  if(!is.null(family) && family == 'cox'){
    selectCols <- c("FID", "IID", covariates, phenotype, status, split.col)
  } else{
    selectCols <- c("FID", "IID", covariates, phenotype, split.col)
  }
  
  phe.master.unsorted <- data.table::fread(
    cmd=paste(cat_or_zcat(phenotype.file, configs), phenotype.file, ' | sed -e "s/^#//g"'),
    colClasses = c("FID" = "character", "IID" = "character"), select = selectCols
  )
  phe.master.unsorted$ID <- paste(phe.master.unsorted$FID, phe.master.unsorted$IID, sep = "_")
  
  # make sure the phe.master has the same individual ordering as in the genotype data
  # so that we don't have error when opening pgen file with sample subset option.
  phe.master <- phe.master.unsorted %>%
    dplyr::left_join(
      data.frame(ID = psam.ids, stringsAsFactors=F) %>%
        dplyr::mutate(sort_order = 1:n()),
      by='ID'
    ) %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()
  rownames(phe.master) <- phe.master$ID
  
  for (name in c(covariates, phenotype)) {
    set(phe.master, i = which(phe.master[[name]] == -9), j = name, value = NA) # missing phenotypes are encoded with -9
  }
  
  # focus on individuals with complete covariates values
  if (is.null(covariates)) {
    phe.no.missing <- phe.master
  } else {
    phe.no.missing <- phe.master %>%
      dplyr::filter_at(dplyr::vars(covariates), dplyr::all_vars(!is.na(.)))
  }
  
  # focus on individuals with at least one observed phenotype values
  phe.no.missing <- phe.no.missing %>%
    dplyr::filter_at(dplyr::vars(phenotype), dplyr::any_vars(!is.na(.))) %>%
    dplyr::filter(ID %in% psam.ids) # check if we have genotype
  
  phe.no.missing.IDs <- phe.no.missing$ID
  
  if(!is.null(split.col)){
    # focus on individuals in training and validation set
    phe.no.missing.IDs <- intersect(
      phe.no.missing.IDs,
      phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]
    )
  }
  if(!is.null(configs[['keep']])){
    # focus on individuals in the specified keep file
    phe.no.missing.IDs <- intersect(phe.no.missing.IDs, readPlinkKeepFile(configs[['keep']]))
  }
  checkMissingPhenoWarning(phe.master, phe.no.missing.IDs)
  
  phe.master[ phe.master$ID %in% phe.no.missing.IDs, ]
}

checkMissingPhenoWarning <- function(phe.master, phe.no.missing.IDs){
  # Show warning message if there are individuals (in phe file)
  # that have (genotype or phenotype) missing values.
  phe.missing.IDs <- phe.master$ID[ ! phe.master$ID %in% phe.no.missing.IDs ]
  if(length(phe.missing.IDs) > 0){
    warning(sprintf(
      'We detected missing values for %d individuals (%s ...).\n',
      length(phe.missing.IDs),
      paste(utils::head(phe.missing.IDs, 5), collapse=", ")
    ))
  }
}

readPlinkKeepFile <- function(keep_file){
  ID <- NULL  # to deal with "no visible binding for global variable"
  keep_df <- data.table::fread(keep_file, colClasses='character', stringsAsFactors=F)
  keep_df$ID <- paste(keep_df$V1, keep_df$V2, sep = "_")
  keep_df %>% dplyr::pull(ID)
}

inferFamily <- function(phe, phenotype, status){
  if (all(unique(phe[[phenotype]] %in% c(0, 1, 2, -9)))) {
    family <- "binomial"
  } else if(!is.null(status) && (status %in% colnames(phe))) {
    family <- "cox"
  } else {
    family <- "gaussian"
  }
  family
}

updateConfigsWithFamily <- function(configs, family){
  out <- configs
  out[['family']] <- family
  if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
  out
}

readBinMat <- function(fhead, configs){
  # This is a helper function to read binary matrix file (from plink2 --variant-score zs bin)
  rows <- data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', fhead, '.vars.zst'), head=F)$V1
  cols <- data.table::fread(paste0(fhead, '.cols'), head=F)$V1
  bin.reader <- file(paste0(fhead, '.bin'), 'rb')
  M = matrix(
    readBin(bin.reader, 'double', n=length(rows)*length(cols), endian = configs[['endian']]),
    nrow=length(rows), ncol=length(cols), byrow = T
  )
  close(bin.reader)
  colnames(M) <- cols
  rownames(M) <- rows
  if (! configs[['save.computeProduct']]) system(paste(
    'rm', paste0(fhead, '.cols'), paste0(fhead, '.vars.zst'),
    paste0(fhead, '.bin'), sep=' '
  ), intern=F, wait=T)
  M
}

cleanUpIntermediateFiles <- function(configs){
  for(subdir in c(configs[["save.dir"]], configs[["meta.dir"]])){
    system(paste(
      'rm', '-rf', file.path(configs[['results.dir']], subdir), sep=' '
    ), intern=F, wait=T)
  }
}