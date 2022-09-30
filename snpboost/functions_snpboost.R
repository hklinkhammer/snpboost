### newly implemented functions

get_coefficients <- function(coeff_path,index=nrow(coeff_path),covariates=NULL){
  coeffs_df <- coeff_path[1:index,]
  coeffs_df <- coeffs_df[rev(!duplicated(rev(coeffs_df$name))),]
  if(length(covariates)>0){
    coeffs <- c(coeff_path$intercept[index],coeff_path[index,covariates],coeffs_df$variant[-1])
  }else{
    coeffs <- c(coeff_path$intercept[index],coeffs_df$variant[-1])
  }
  names(coeffs)=c("(Intercept)",covariates,coeffs_df$name[-1])
  coeffs <- unlist(coeffs)
  coeffs
}

build_coefficients_path <- function(coeff_path,variants=na.omit(unique(coeff_path$name))){
  coeff_df <-  coeff_path[coeff_path$name %in% variants,] %>%
    add_column(iteration=0:(nrow(.)-1),.before=0) %>%
    pivot_wider(values_from = variant) %>% 
    fill(3:ncol(.)) %>% mutate_all(~replace(.,is.na(.),0)) %>%
    pivot_longer(cols=-1)
  
  coeff_df
}

plot_coefficients_path <- function(coeff_path,variants=na.omit(unique(coeff_path$name)),covariates=NULL,plot_intercept=F,plot_covariates=F,legend=T){
  coeff_df <- build_coefficients_path(coeff_path,variants)
  
  if(!plot_intercept){coeff_df <- coeff_df %>% filter(name!="intercept")}
  if(!plot_covariates){coeff_df <- coeff_df %>% filter(!name %in% covariates)}
  
  ggplot(coeff_df,aes(x=iteration,y=value,color=name))+
    geom_line()+
    theme(legend.position = ifelse(legend,"bottom","none"))
}

predict_SNP_boost <- function(fit, new_genotype_file, new_phenotype_file, phenotype, subset=NULL,idx=NULL){
  if(is.null(idx)){
    if(!is.null(fit$metric.val)){
      idx <- which.min(fit$metric.val)
    }else{
      idx <- which.min(fit$metric.train)
    }
  }
  final_betas <- get_coefficients(fit$beta,idx,covariates=fit$configs[['covariates']])
  chosen_SNPs <- cbind(str_split(names(final_betas)[-(1:(length(fit$configs[['covariates']])+1))],"_",simplify = T),
                       final_betas[-(1:(length(fit$configs[['covariates']])+1))])
  rownames(chosen_SNPs)=rep("",nrow(chosen_SNPs))
  fwrite(chosen_SNPs,paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),sep="\t")
  
  phe <- fread(new_phenotype_file)
  if(!is.null(subset)){
    phe <- phe[phe[[fit$configs[['split.col']]]]==subset,]
  }
  
  if(nrow(chosen_SNPs)>0){
    plink2_cmd <- paste(fit$configs['plink2.path'],"--pfile",new_genotype_file,"vzs","--score",
                        paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),1,2,3,"header cols=maybefid,nallele,denom,dosagesum,scoreavgs,scoresums","--out",paste0(fit$configs['results.dir'],"/PRS"))
    system(plink2_cmd, intern=F, wait=T)
    
    PRS <- fread(paste0(fit$configs['results.dir'],"/PRS.sscore")) %>% rename(PRS = SCORE1_SUM)
    
    phe_PRS <- left_join(phe,PRS,by="IID")
  }else{
    phe_PRS <- phe %>% mutate(PRS=0) 
  }
  
  pred <- matrix(final_betas[1] + phe_PRS$PRS + if(!is.null(fit$configs[['covariates']])){as.matrix(phe_PRS%>%select(fit$configs[['covariates']])) %*% final_betas[fit$configs[['covariates']]]}else{0},ncol=1)
  
  rownames(pred)=phe_PRS$IID
  
  residuals <- computeResiduals(phe_PRS[[phenotype]], pred, fit$configs[['family']])
  rownames(residuals) <- phe_PRS$IID
  
  metric <- computeMetric(residuals, phe_PRS[[phenotype]], pred, fit$configs[['metric']])
  
  if(fit$configs[['family']]=='gaussian'){
    Rsquare <- cor(phe_PRS[[phenotype]],pred)^2
  }else if(fit$configs[['family']]=='binomial'){
    AUC <- as.numeric(pROC::auc(phe_PRS[[phenotype]],exp(pred)/(1+exp(pred))))
  }
  out <- list(prediction = pred, residuals=residuals, metric = metric,
              Rsquare=if(fit$configs[['family']]=='gaussian'){Rsquare}else{NULL}, AUC=if(fit$configs[['family']]=='binomial'){AUC}else{NULL})
}

computeResiduals <- function(phenotype, prediction, family) {
  if (family == 'gaussian') {
    residuals <- phenotype - prediction
  } else if (family == 'binomial'){
    p.hat <- exp(prediction)/(1+exp(prediction))
    residuals <- (phenotype - p.hat)
  }
  residuals
}

### so far only MSEP implemented, to be extended for further metrics
computeMetric <- function(residual, phenotype, prediction, metric.type) {
  if (metric.type == 'MSEP') {
    metric <- mean(residual^2)
  } else if (metric.type == 'log loss'){
    p.hat <- exp(prediction)/(1+exp(prediction))
    metric <- - mean(phenotype * log(p.hat) + (1-phenotype) * log( 1-p.hat), na.rm=T)
  }
  metric
}

### check outer stopping criterion, leant on checkEarlyStopping from snpnet
checkEarlyStoppingBatches <- function(metric.val, iter, b_stop, configs, m_batch_list){
  if(iter <= b_stop){
    earlyStop <- FALSE
  }else{
    max.valid.idx.lag <- 1+sum(m_batch_list[1:(iter-b_stop)])
    if(configs[['metric']]=='MSEP' || configs[['metric']]=='log loss'){
      min.val.1 <- min(metric.val[1:max.valid.idx.lag])
      min.val.2 <- min(metric.val[(max.valid.idx.lag+1):sum(!is.na(metric.val))])
      snpboostLogger(sprintf('batch=%g, stopping lag=%g, min.val.1=%g min.val.2=%g', iter, max.valid.idx.lag, min.val.1, min.val.2))
      if (
        (configs[['early.stopping']]>0) &&
        (min.val.1 < min.val.2)
      ) {
        snpboostLogger(sprintf(
          "Early stopped at iteration %d (step =%d ) with validation metric: %.14f.",
          iter, which.min(metric.val), min(metric.val, na.rm = T)
        ))
        earlyStop <- TRUE
      } else {
        earlyStop <- FALSE
      }
    }else{ ### so far no other metric is supported
      earlyStop <- FALSE
    }
  }
  earlyStop
}

### functions from snpnet package that are used for snpboost and are adapted

####try to import from snpnet_compact
#### setupConfigs: default value for standardize.variant = TRUE
setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, covariates, split.col, status.col, mem){
  out.args <- as.list(environment())
  defaults <- list(
    missing.rate = 0.1,
    MAF.thresh = 0.001,
    nCores = 1,
    glmnet.thresh = 1e-07,
    num.snps.batch = 1000,
    vzs=TRUE, # geno.pfile vzs
    increase.size = NULL,
    standardize.variant = TRUE,
    early.stopping = TRUE,
    stopping.lag = 2,
    niter = 50,
    keep = NULL,
    KKT.verbose = FALSE,
    save = FALSE,
    save.computeProduct = FALSE,
    prevIter = 0,
    results.dir = NULL,
    meta.dir = 'meta',
    save.dir = 'results',
    verbose = FALSE,
    KKT.check.aggressive.experimental = FALSE,
    gcount.basename.prefix = 'snpboost.train',
    gcount.full.prefix=NULL,
    endian="little",
    metric=NULL,
    plink2.path='plink2',
    zstdcat.path='zstdcat',
    zcat.path='zcat',
    rank = TRUE,
    excludeSNP = NULL
  )
  out <- defaults
  
  # store additional params
  for (name in setdiff(names(out.args), "configs")) {
    out[[name]] <- out.args[[name]]
  }
  
  # update the defaults with the specified parameters and keep redundant parameters from configs
  for (name in names(configs)) {
    out[[name]] <- configs[[name]]
  }
  
  # update settings
  out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
  if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2
  
  # configure the temp file locations
  #   We will write some intermediate files to meta.dir and save.dir.
  #   those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
  if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
  dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
  dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
  if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
    out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
  )
  
  out
}

### for gaussian we use MSEP
### for binomial we use log loss
setDefaultMetric <- function(family){
  if (family == "gaussian") {
    metric <- 'MSEP'
  } else if (family == "binomial") {
    metric <- 'log loss'
#  } else if (family == "cox") {
#    metric <- 'C'
  } else {
    stop(paste0('The specified family (', family, ') is not supported!'))
  }
  metric
}

### computeStats: calculation of means corrected
computeStats <- function(pfile, ids, configs) {
  ID <- original_ID <- NULL  # to deal with "no visible binding for global variable"
  ALT <- MISSING_CT <- OBS_CT <- HAP_ALT_CTS <- HET_REF_ALT_CTS <- TWO_ALT_GENO_CTS <- NULL
  stats_msts <- stats_means <- stats_pNAs <- stats_SDs <- NULL
  
  keep_f       <- paste0(configs[['gcount.full.prefix']], '.keep')
  gcount_tsv_f <- paste0(configs[['gcount.full.prefix']], '.gcount.tsv')
  
  dir.create(dirname(configs[['gcount.full.prefix']]), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
    gcount_df <- data.table::fread(gcount_tsv_f)
  } else {
    # To run plink2 --geno-counts, we write the list of IDs to a file
    data.frame(ID = ids) %>%
      tidyr::separate(ID, into=c('FID', 'IID'), sep='_') %>%
      data.table::fwrite(keep_f, sep='\t', col.names=F)
    
    # Run plink2 --geno-counts
    cmd_plink2 <- paste(
      configs[['plink2.path']],
      '--threads', configs[['nCores']],
      '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
      '--keep', keep_f,
      '--out', configs[['gcount.full.prefix']],
      '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs'
    )
    if (!is.null(configs[['mem']])) cmd_plink2 <- paste(cmd_plink2, '--memory', configs[['mem']])
    
    system(cmd_plink2, intern=F, wait=T)
    
    # read the gcount file
    ### means ans msts corrected (count HAP_ALT_CTS as two instead of one)
    gcount_df <-
      data.table::fread(paste0(configs[['gcount.full.prefix']], '.gcount')) %>%
      dplyr::rename(original_ID = ID) %>%
      dplyr::mutate(
        ID = paste0(original_ID, '_', ALT),
        stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
        stats_means = (2*HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_msts  = (4*HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_SDs   = (stats_msts - stats_means * stats_means) * OBS_CT / (OBS_CT - 1)
      )
  }
  
  out <- list()
  out[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
  out[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
  out[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)
  
  for(key in names(out)){
    names(out[[key]]) <- gcount_df %>% dplyr::pull(ID)
  }
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
  out[["excludeSNP"]] <- out[["excludeSNP"]][ ! is.na(out[["excludeSNP"]]) ]
  out[["excludeSNP"]] <- base::unique(c(configs[["excludeSNP"]], out[["excludeSNP"]]))
  
  if (configs[['save']]){
    gcount_df %>% data.table::fwrite(gcount_tsv_f, sep='\t')
    saveRDS(out[["excludeSNP"]], file = file.path(dirname(configs[['gcount.full.prefix']]), "excludeSNP.rda"))
  }
  
  out
}

### computeProduct: prod.full adapted to be the exact correlation
### maybe algorithm could be fastened by using single precision in plink (as in snpnet-compact)
computeProduct <- function(residual, pfile, vars, stats, configs, iter) {
  ID <- NULL  # to deal with "no visible binding for global variable"
  time.computeProduct.start <- Sys.time()
  snpboostLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)
  
  gc_res <- gc()
  if(configs[['KKT.verbose']]) print(gc_res)
  
  snpboostLogger('Start plink2 --variant-score', indent=3, log.time=time.computeProduct.start)
  dir.create(file.path(configs[['results.dir']], configs[["save.dir"]]), showWarnings = FALSE, recursive = T)
  
  residual_f <- file.path(configs[['results.dir']], configs[["save.dir"]], paste0("residuals_iter_", iter, ".tsv"))
  
  # write residuals to a file
  residual_df <- data.frame(residual)
  colnames(residual_df) <- paste0('lambda_idx_', colnames(residual))
  residual_df %>%
    tibble::rownames_to_column("ID") %>%
    tidyr::separate(ID, into=c('#FID', 'IID'), sep='_') %>%
    data.table::fwrite(residual_f, sep='\t', col.names=T)
  
  # Run plink2 --geno-counts
  cmd_plink2 <- paste(
    configs[['plink2.path']],
    '--threads', configs[['nCores']],
    '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
    '--read-freq', paste0(configs[['gcount.full.prefix']], '.gcount'),
    '--keep', residual_f,
    '--out', stringr::str_replace_all(residual_f, '.tsv$', ''),
    '--variant-score', residual_f, 'zs', 'bin' #'bin4', 'single-prec'
  )
  if (!is.null(configs[['mem']])) {
    cmd_plink2 <- paste(cmd_plink2, '--memory', max(as.integer(configs[['mem']]) - ceiling(sum(as.matrix(gc_res)[,2])),640))
  }
  
  system(cmd_plink2, intern=F, wait=T)
  
  prod.full <- readBinMat(stringr::str_replace_all(residual_f, '.tsv$', '.vscore'), configs)
  if (! configs[['save.computeProduct']] ) system(paste(
    'rm', residual_f, stringr::str_replace_all(residual_f, '.tsv$', '.log'), sep=' '
  ), intern=F, wait=T)
  
  snpboostLoggerTimeDiff('End plink2 --variant-score.', time.computeProduct.start, indent=4)
  
  rownames(prod.full) <- vars 
  prod.full[stats[["excludeSNP"]], ] <- NA
  if (configs[["standardize.variant"]]) {
    prod.full <- 1/(length(residual)-1)*(prod.full-length(residual)*stats[["means"]]*mean(residual))/(sqrt(stats[["sds"]])*sd(residual))
  }
  snpboostLoggerTimeDiff('End computeProduct().', time.computeProduct.start, indent=3)
  prod.full
}


prepareFeatures <- function(pgen, vars, names, stat) {
  buf <- pgenlibr::ReadList(pgen, match(names, vars), meanimpute=F)
  features.add <- as.data.table(buf)
  colnames(features.add) <- names
  for (j in 1:length(names)) {
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stat[["means"]][names[j]])
  }
  features.add
}

## logger functions

snpboostLogger <- function(message, log.time = NULL, indent=0, funcname='SNP_boost'){
  if (is.null(log.time)) log.time <- Sys.time()
  cat('[', as.character(log.time), ' ', funcname, '] ', rep(' ', indent * 2), message, '\n', sep='')
}

timeDiff <- function(start.time, end.time = NULL) {
  if (is.null(end.time)) end.time <- Sys.time()
  paste(round(end.time-start.time, 4), units(end.time-start.time))
}

snpboostLoggerTimeDiff <- function(message, start.time, end.time = NULL, indent=0){
  if (is.null(end.time)) end.time <- Sys.time()
  snpboostLogger(paste(message, "Time elapsed:", timeDiff(start.time, end.time), sep=' '), log.time=end.time, indent=indent)
}
