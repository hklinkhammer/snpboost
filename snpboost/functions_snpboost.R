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

split_last <- function(string, pattern){
  last_pos <- stringi::stri_locate_last(string, regex = pattern)
  
  if (!is.na(last_pos[1])) {
    # Extract the part before the last occurrence
    before_last <- substring(string, 1, last_pos[1] - 1)
    
    # Extract the part after the last occurrence
    after_last <- substring(string, last_pos[2] + 1)
    
    # Combine into a list
    result <- c(before_last, after_last)
  } else {
    # If the pattern is not found, return the original string in the first part
    result <- c(string, "")
  }
  
  result
}

predict_snpboost <- function(fit, new_genotype_file, new_phenotype_file, phenotype, subset=NULL,idx=NULL){
  if(is.null(idx)){
    if(!is.null(fit$metric.val)){
      idx <- which.min(fit$metric.val)
    }else{
      idx <- which.min(fit$metric.train)
    }
  }
  final_betas <- get_coefficients(fit$beta,idx,covariates=fit$configs[['covariates']])

  chosen_SNPs <- cbind(t(mapply(split_last, names(final_betas)[-(1:(length(fit$configs[['covariates']])+1))], MoreArgs = list(pattern = "_"))),
                       final_betas[-(1:(length(fit$configs[['covariates']])+1))])
  rownames(chosen_SNPs)=rep("",nrow(chosen_SNPs))
  fwrite(chosen_SNPs,paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),sep="\t")
  
  phe <- fread(new_phenotype_file)
  if(!is.null(subset)){
    phe <- phe[phe[[fit$configs[['split.col']]]]==subset,]
  }
  
  if(nrow(chosen_SNPs)>0){
    plink2_cmd <- paste(fit$configs['plink2.path'],"--pfile",new_genotype_file,"vzs","--score",
                        paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),1,2,3,"header cols=maybefid,denom,dosagesum,scoreavgs,scoresums","--out",paste0(fit$configs['results.dir'],"/PRS"))
    system(plink2_cmd, intern=F, wait=T)
    
    PRS <- fread(paste0(fit$configs['results.dir'],"/PRS.sscore")) %>% rename(PRS = SCORE1_SUM)
    
    phe_PRS <- left_join(phe,PRS,by="IID")
  }else{
    phe_PRS <- phe %>% mutate(PRS=0) 
  }
  
  pred <- matrix(final_betas[1] + phe_PRS$PRS + if(!is.null(fit$configs[['covariates']])){as.matrix(phe_PRS%>%select(fit$configs[['covariates']])) %*% final_betas[fit$configs[['covariates']]]}else{0},ncol=1)
  
  rownames(pred)=phe_PRS$IID
  if(fit$configs[['family']]=='cox'){
      surv <- survival::Surv(phe_PRS[[phenotype]], phe_PRS[[fit$configs[['status.col']]]])
      T_order <- order(surv[,1])
      IPCw <- IPCweights(surv)
      sum_help <- function(i){ifelse(i<length(surv[,1]),IPCw[T_order][i]^2*sum((surv[,1][T_order][(i+1):length(surv)]-surv[,1][T_order][i])!=0)+
                                       0.5*IPCw[T_order][i]^2*(sum((surv[,1][T_order]-surv[,1][T_order][i])==0)-1),0.5*IPCw[T_order][i]^2*(sum((surv[,1][T_order]-surv[,1][T_order][i])==0)-1))}
      
      w_denom <- sum(mapply(sum_help,1:length(surv)))
      wweights <- IPCw^2/w_denom
      
      residuals <- computeResiduals(surv, pred, fit$configs, IPCw, w_denom, T_order, wweights, fit$sigma)
      rownames(residuals) <- phe_PRS$IID
      
      metric <- computeMetric(residuals, surv, pred, fit$configs, surv, IPCw, fit$sigma, T_order)
    }else{
      residuals <- computeResiduals(phe_PRS[[phenotype]], pred, fit$configs)
      rownames(residuals) <- phe_PRS$IID
    
      metric <- computeMetric(residuals, phe_PRS[[phenotype]], pred, fit$configs, sigma = fit$sigma)
  }
  
  if(fit$configs[['family']]=='gaussian'){
    Rsquare <- cor(phe_PRS[[phenotype]],pred)^2
  }else if(fit$configs[['family']]=='binomial'){
    AUC <- as.numeric(pROC::auc(phe_PRS[[phenotype]],exp(pred)/(1+exp(pred))))
  }
  out <- list(prediction = pred, residuals=residuals, metric = metric,
              Rsquare=if(fit$configs[['family']]=='gaussian'){Rsquare}else{NULL}, AUC=if(fit$configs[['family']]=='binomial'){AUC}else{NULL})
}

computeResiduals <- function(phenotype, prediction, configs, IPCw=NULL, w_denom = NULL, T_order = NULL, wweights = NULL, sigma=NULL) {
  if (configs[['metric']] == 'MSEP') {
    residuals <- phenotype - prediction
  } else if(configs[['metric']] == 'quantilereg'){
    residuals <- -(1-configs[['quantile']])*(prediction > phenotype) + configs[['quantile']]*(prediction <= phenotype)
  } else if (configs[['metric']] == 'log_loss'){
    p.hat <- exp(prediction)/(1+exp(prediction))
    residuals <- (phenotype - p.hat)
  } else if (configs[['metric']] == 'C'){
    
    residuals <- numeric(length(prediction))
    sum_approx <- function(i){ifelse(i<length(prediction),wweights[T_order][i]*sum(approxGrad(prediction[T_order][(i+1):length(prediction)]-prediction[T_order][i])*as.numeric((phenotype[,1][T_order][(i+1):length(prediction)]-phenotype[,1][T_order][i])!=0)+
                                                     0.5*approxGrad(prediction[T_order][(i+1):length(prediction)]-prediction[T_order][i])*as.numeric((phenotype[,1][T_order][(i+1):length(prediction)]-phenotype[,1][T_order][i])==0)),0)-
                              ifelse(i>1,sum(wweights[T_order][1:(i-1)]*(approxGrad(prediction[T_order][i]-prediction[T_order][1:(i-1)])*as.numeric((phenotype[,1][T_order][i]-phenotype[,1][T_order][1:(i-1)])!=0)+
                                                     0.5*approxGrad(prediction[T_order][i]-prediction[T_order][1:(i-1)])*as.numeric((phenotype[,1][T_order][i]-phenotype[,1][T_order][1:(i-1)])==0))),0) 
    }
    
    sum_approx2 <- function(i){ifelse(i<length(prediction),wweights[T_order][i]*sum(approxGrad(prediction[T_order][(i+1):length(prediction)]-prediction[T_order][i])[(phenotype[,1][T_order][(i+1):length(prediction)]-phenotype[,1][T_order][i])!=0]),0)-
        ifelse(i>1,sum(wweights[T_order][1:(i-1)]*(approxGrad(prediction[T_order][i]-prediction[T_order][1:(i-1)])*as.numeric((phenotype[,1][T_order][i]-phenotype[,1][T_order][1:(i-1)])!=0))),0) 
    }

    residuals[T_order]=parallel::mcmapply(sum_approx2,1:length(prediction),mc.cores = configs[['nCores']])
    
  } else if(configs[['metric']] == 'weightedL2'){
    time <- phenotype[,1]
    ### log(0) will produce NaN
    time[time==0] <- sort(unique(time))[2]/10
    
    residuals <- 2*IPCw*(log(time)-prediction)
  } else if(configs[['metric']] == 'Cox'){
    help_fct_0 <- function(i){sum((phenotype[,1]>=phenotype[i,1])*exp(prediction))}
    help_vct <- mcmapply(help_fct_0,1:length(prediction),mc.cores = configs[['nCores']])
    help_fct <- function(i){sum((phenotype[i,1]>=phenotype[,1])/help_vct)}
    residuals<- phenotype[,2]-exp(prediction)*mcmapply(help_fct,1:length(prediction),mc.cores = configs[['nCores']])
  } else if(configs[['metric']] == 'AFT-Weibull'){
    
    #log(0) will produce Inf
    time <- phenotype[,1]
    time[time==0] <- sort(unique(time))[2]/10
    
    lnt <- log(time)
    
    eta <- (lnt-prediction)/sigma
    
    residuals <- (exp(eta)-phenotype[,2])/sigma
    
  } else if(configs[['metric']] == 'AFT-logistic'){
    
    time <- phenotype[,1]
    time[time==0] <- sort(unique(time))[2]/10
    
    lnt <- log(time)
    
    eta <- (lnt-prediction)/sigma
    
    residuals <- (exp(eta)-phenotype[,2])/(sigma*(1+exp(eta)))
    
  } else if(configs[['metric']] == 'AFT-normal'){
    
    time <- phenotype[,1]
    time[time==0] <- sort(unique(time))[2]/10
    
    derf <- function(x){
      -x*dnorm(x)
    }
    
    lnt <- log(time)
    
    eta <- (lnt-prediction)/sigma
    
    v1 <- derf(eta)/(sigma*dnorm(eta))
    v1[phenotype[,2]==0]=0
    
    v2 <- dnorm(eta)/(sigma*(1-pnorm(eta)))
    v2[phenotype[,2]==1]=0
    
    residuals <- -v1+v2
  } else if (configs[['metric']] == 'Poisson'){
    residuals <- phenotype - exp(prediction)
  } else if (configs[['metric']] == 'negative_binomial'){
    residuals <- phenotype - (phenotype + sigma)*exp(prediction)/(exp(prediction)+sigma)
  }
  residuals
}

optimizeSigma <- function(phenotype, prediction, configs){
  if(configs[['metric']] == 'AFT-Weibull'){
    riskS <- function(sigma, phenotype, pred){
      time <- phenotype[,1]
      #### log(0) will produce NaN
      time[time==0] <- sort(unique(time))[2]/10
      
      fw <- function(x){
        exp(x - exp(x))
      }
      Sw <- function(x){
        exp(-exp(x))
      }
      
      lnt <- log(time)
      Sevent <- phenotype[,2]
      eta <- (lnt-pred) / sigma
      
      plloss <- - Sevent * log(fw(eta)/sigma) - (1 - Sevent) * log(Sw(eta))
      
      sum(plloss)
    }
  } else if(configs[['metric']] == 'AFT-logistic'){
    riskS <- function(sigma, phenotype, pred){
      time <- phenotype[,1]
      #### log(0) will produce NaN
      time[time==0] <- sort(unique(time))[2]/10
      
      fw <- function(x){
        exp(x)/(1+exp(x))^2
      }
      Sw <- function(x){
        1/(1+exp(x))
      }
      
      lnt <- log(time)
      Sevent <- phenotype[,2]
      eta <- (lnt-pred) / sigma
      
      plloss <- - Sevent * log(fw(eta)/sigma) - (1 - Sevent) * log(Sw(eta))
      
      sum(plloss)
    }
  }else if(configs[['metric']] == 'AFT-normal'){
    riskS <- function(sigma, phenotype, pred){
      time <- phenotype[,1]
      #### log(0) will produce NaN
      time[time==0] <- sort(unique(time))[2]/10
      
      fw <- function(x){
        dnorm(x)
      }
      Sw <- function(x){
        1-pnorm(x)
      }
      
      lnt <- log(time)
      Sevent <- phenotype[,2]
      eta <- (lnt-pred) / sigma
      
      plloss <- - Sevent * log(fw(eta)/sigma) - (1 - Sevent) * log(Sw(eta))
      
      sum(plloss)
      }
  } else if (configs[['metric']] == 'negative_binomial'){
    riskS <- function(sigma, phenotype, pred){
      -sum(lgamma(phenotype + sigma) - lgamma(phenotype+1) - lgamma(sigma)+sigma*log(sigma) -
             (sigma + phenotype)*log(exp(pred)+sigma) + phenotype*pred)
    }
  }
  sigma <- optimize(riskS, interval = c(0,100), phenotype = phenotype, pred= prediction)$minimum
  sigma
}


### compute evaluation metrics
computeMetric <- function(residual, phenotype, prediction, configs, phenotype2 = NULL, weights = NULL, sigma= NULL, T_order = NULL) {
  if (configs[['metric']] == 'MSEP') {
    metric <- mean(residual^2)
  } else if(configs[['metric']] == 'quantilereg'){
    metric <- sum((1-configs[['quantile']])*(prediction-phenotype)*(prediction > phenotype)+configs[['quantile']]*(phenotype - prediction)*(phenotype>=prediction))
  }else if (configs[['metric']] == 'log_loss'){
    p.hat <- exp(prediction)/(1+exp(prediction))
    metric <- - mean(phenotype * log(p.hat) + (1-phenotype) * log( 1-p.hat))
  } else if (configs[['metric']] == 'C'){
    
    sum_cindex <- function(i){ifelse(i<length(prediction),weights[T_order][i]*sum(approxLoss(prediction[T_order][(i+1):length(prediction)]-prediction[T_order][i])[(phenotype[,1][T_order][(i+1):length(prediction)]-phenotype[,1][T_order][i])!=0]),0)}
    
    metric=-sum(parallel::mcmapply(sum_cindex,1:length(prediction),mc.cores = configs[['nCores']]))
    
  } else if(configs[['metric']] == 'weightedL2'){
    time <- phenotype[,1]
    ### log(0) will produce NaN
    time[time==0] <- sort(unique(time))[2]/10
    
    metric <- mean(weights*(log(time)-prediction)^2)
  } else if(configs[['metric']] == 'Cox'){
    help_fct <- function(i){phenotype[i,2]*(prediction[i]-log(sum((phenotype[,1]>=phenotype[i,1])*exp(prediction))))}
    metric <- -sum(mcmapply(help_fct,1:length(prediction),mc.cores = configs[['nCores']]))
  } else if(configs[['metric']] == 'AFT-Weibull'){
    time <- phenotype[,1]
    #### log(0) will produce NaN
    time[time==0] <- sort(unique(time))[2]/10
    
    lfw <- function(pred){
      ordinal::dgumbel(pred,max=F,log=T)
    }
    Sw <- function(pred){
      exp(-exp(pred))
    }
    
    lnt <- log(time)
    Sevent <- phenotype[,2]
    eta <- (lnt-prediction) / sigma
    
    plloss <- - Sevent * (log(1/sigma)+lfw(eta)) - (1 - Sevent) * (-exp(eta))
    
    metric <- sum(plloss)
    
  } else if(configs[['metric']] == 'AFT-logistic'){
    time <- phenotype[,1]
    #### log(0) will produce NaN
    time[time==0] <- sort(unique(time))[2]/10
    
    fw <- function(x){
      exp(x)/(1+exp(x))^2
    }
    Sw <- function(x){
      1/(1+exp(x))
    }
    
    lnt <- log(time)
    Sevent <- phenotype[,2]
    eta <- (lnt-prediction) / sigma
    
    plloss <- - Sevent * log(fw(eta)/sigma) - (1 - Sevent) * log(Sw(eta))
    
    metric <- sum(plloss)
  } else if(configs[['metric']] == 'AFT-normal'){
    time <- phenotype[,1]
    #### log(0) will produce NaN
    time[time==0] <- sort(unique(time))[2]/10
    
    fw <- function(pred){
      dnorm(pred)
    }
    Sw <- function(pred){
      1-pnorm(pred)
    }
    
    lnt <- log(time)
    Sevent <- phenotype[,2]
    eta <- (lnt-prediction) / sigma
    
    plloss <- - Sevent * log(fw(eta)/sigma) - (1 - Sevent) * log(Sw(eta))
    
    metric <- sum(plloss)
  } else if(configs[['metric']] == 'Poisson'){
    metric <- -sum(dpois(phenotype,exp(prediction),log=T))
  } else if(configs[['metric']] == 'negative_binomial'){
    metric <- -sum(lgamma(phenotype + sigma) - lgamma(phenotype+1) - lgamma(sigma)+sigma*log(sigma) -
                     (sigma + phenotype)*log(exp(prediction)+sigma) + phenotype*prediction)
  }
  metric
}

### check outer stopping criterion, leant on checkEarlyStopping from snpnet
checkEarlyStoppingBatches <- function(metric.val, iter, b_stop, configs, m_batch_list){
  if(iter <= b_stop){
    earlyStop <- FALSE
  }else{
    max.valid.idx.lag <- 1+sum(m_batch_list[1:(iter-b_stop)])
      min.val.1 <- min(metric.val[1:max.valid.idx.lag],na.rm = T)
      min.val.2 <- min(metric.val[(max.valid.idx.lag+1):sum(!is.na(metric.val))],na.rm = T)
      snpboostLogger(sprintf('batch=%g, stopping lag=%g, min.val.1=%g min.val.2=%g', iter, max.valid.idx.lag, min.val.1, min.val.2))
      if (
        (configs[['early.stopping']]>0) &&
        (min.val.1 <= min.val.2)
      ) {
        snpboostLogger(sprintf(
          "Early stopped at iteration %d (step =%d ) with validation metric: %.14f.",
          iter, which.min(metric.val), min(metric.val, na.rm = T)
        ))
        earlyStop <- TRUE
      } else {
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
    excludeSNP = NULL,
    quantile = 0.5
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

updateConfigsWithMetric <- function(configs, metric){
  out <- configs
  if (!is.null(metric)) out[['metric']] <- metric
  out
}

### set default metric if none is specified
setDefaultMetric <- function(family){
  if (family == "gaussian") {
    metric <- 'MSEP'
  } else if (family == "binomial") {
    metric <- 'log_loss'
  } else if (family == "cox") {
    metric <- 'AFT-Weibull'
  } else if (family == "count") {
    metric <- 'Poisson'
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

## CIndex functions

### inverse probability of censoring weights
### see van der Laan & Robins (2003)
### taken from mboost::IPCweights
IPCweights <- function(x, x2= NULL, maxweight = 5) {
  
  if (!extends(class(x), "Surv"))
    stop(sQuote("x"), " is not a Surv object")
  
  if (is.null(x2)){x2 = x}
  
  event <- x[,2]
  x[,2] <- 1 - event
  km <- survival::survfit(x ~ 1)
  Ghat <- getsurv(km, times = x2[,1]) ## see github issue #54
  Ghat[x2[,2] == 0] <- 1
  w <- x2[,2] / Ghat
  w[w > maxweight] <- maxweight
  w
}

### extract survival probabilities
### taken from ipred:::getsurv
### DO NOT TOUCH HERE
getsurv <- function(obj, times)
{
  # get the survival probability for times from KM curve j'
  
  if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
  # <FIXME: methods may have problems with that>
  class(obj) <- NULL
  # </FIXME>
  lt <- length(times)
  nsurv <- times
  
  # if the times are the same, return the km-curve
  
  if(length(times) == length(obj$time)) {
    if (all(times == obj$time)) return(obj$surv)
  }
  
  # otherwise get the km-value for every element of times separatly
  
  inside <- times %in% obj$time
  for (i in (1:lt)) {
    if (inside[i])
      nsurv[i] <- obj$surv[obj$time == times[i]]
    else  {
      less <- obj$time[obj$time < times[i]]
      if (length(less) == 0)
        nsurv[i] <- 1
      else
        nsurv[i] <- obj$surv[obj$time == max(less)]
    }
  }
  nsurv
}

approxGrad <- function(x, sigma = 0.1) {    ## sigmoid function for gradient
  exp(x/sigma) / (sigma * (1 + exp(x/sigma))^2) 
}

approxLoss <- function(x, sigma = 0.1) { ## sigmoid function for loss
  1 / (1 + exp(x / sigma))
}
