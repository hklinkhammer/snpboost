#### based on snpnet function, glm part replaced by boosting

snpboost <- function(genotype.pfile, phenotype.file, phenotype, family = NULL, covariates = NULL, 
                      split.col = NULL, status.col = NULL, mem = NULL, configs = NULL,
                      p_batch = 1000, m_batch = 1000, b_max = 20000, b_stop = 2,
                      sl = 0.1, coeff_path_save = TRUE, give_residuals = FALSE) {

  ID <- ALT <- NULL

  validation <- (!is.null(split.col))
  time.start <- Sys.time()
  snpboostLogger('Start SNPboost', log.time = time.start)

  snpboostLogger('Preprocessing start..')

  ### --- Read genotype IDs --- ###
  ids <- list(); phe <- list()
  ids[['psam']] <- readIDsFromPsam(paste0(genotype.pfile, '.psam'))

  ### --- combine the specified configs with the default values --- ###
  configs <- setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, covariates, split.col, status.col, mem)

  ### --- Read phenotype file --- ###
  phe[['master']] <- readPheMaster(phenotype.file, ids[['psam']], family, covariates, phenotype, status.col, split.col, configs)

  ### --- infer family and update the configs --- ###
  if (is.null(family)) family <- inferFamily(phe[['master']], phenotype, status.col)
  configs <- updateConfigsWithFamily(configs, family)
  if (configs[['verbose']]) print(configs)
  
  ### --- so far only gaussian family is supported --- ###
  if (!family %in% c("gaussian","binomial")){
    stop("So far only gaussian and binary phenotypes are supported by snpboost.")
  }

  ### --- Define the set of individual IDs for training (and validation) set(s) --- ###
  if(is.null(split.col)){
      splits <- c('train')
      ids[['train']] <- phe[['master']]$ID
  }else{
      splits <- c('train', 'val')
      for(s in splits){
          ids[[s]] <- phe[['master']]$ID[ phe[['master']][[split.col]] == s ]
      }
  }

  ### --- Prepare the feature matrix --- ###
  features <- list()
  for(s in splits){
      phe[[s]] <- phe[['master']][match(ids[[s]], phe[['master']]$ID), ]
      rownames(phe[[s]]) <- phe[[s]]$ID
      if (length(covariates) > 0) {
          features[[s]] <- phe[[s]][, covariates, with = F]
      } else {
          features[[s]] <- NULL
      }
      if(configs[['verbose']]) snpboostLogger(sprintf("The number of individuals in %s set: %d", s, dim(phe[[s]])[1]))
  }

  ### --- Prepare the response --- ###
  response <- list() ; status <- list() ; surv <- list() ; pred <- list()
  for(s in splits){
      response[[s]] <- phe[[s]][[phenotype]]
  }

  ### --- Read genotypes --- ###
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  configs[["excludeSNP"]] <- base::intersect(configs[["excludeSNP"]], vars)
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst'))
  pgen <- list()
  for(s in splits) pgen[[s]] <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids[[s]], ids[['psam']]))
  pgenlibr::ClosePvar(pvar)

  stats <- computeStats(genotype.pfile, phe[['train']]$ID, configs = configs)

  ### --- End --- ###
  snpboostLoggerTimeDiff("Preprocessing end.", time.start, indent=1)

  snpboostLogger("Iteration 0")
    glmmod <- stats::glm(
      stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates), collapse = " + "))),
      data = phe[['train']], family = family
    )
    
    prediction <- predict(glmmod, newdata = phe[['train']])
    residual <- matrix(computeResiduals(phe[['train']][[phenotype]],prediction,family),ncol=1)
    
    prediction_val <- predict(glmmod, newdata = phe[['val']])
    residual_val <- if(validation){matrix(computeResiduals(phe[['val']][[phenotype]],prediction_val,family), ncol = 1)}else{NULL}
    
    beta <- data.frame(NA,t(glmmod$coefficients),0,row.names = NULL)
    names(beta)=c("name","intercept",covariates,"variant")

    rownames(residual) <- rownames(phe[['train']])
    colnames(residual) <- c('0')
    
    if(validation){
      rownames(residual_val) <- rownames(phe[['val']])
      colnames(residual_val) <- c('0')
    }

    if (configs[['verbose']]) snpboostLogger("  Start computing inner product for initialization ...")
    time.prod.init.start <- Sys.time()

    prod.full <- computeProduct(residual, genotype.pfile, vars, stats, configs, iter=0)
    score <- abs(prod.full[, 1])
    
    if (configs[['verbose']]) snpboostLoggerTimeDiff("  End computing inner product for initialization.", time.prod.init.start)

    nobs <- nrow(phe[['train']])
    nvars <- length(vars)-length(stats[["excludeSNP"]])-length(covariates)

    metric.train <-numeric()
    metric.val <- numeric()

    metric.train[1] <- computeMetric(residual, phe[['train']][[phenotype]], prediction, configs[['metric']])
    metric.val[1] <- computeMetric(residual_val, phe[['val']][[phenotype]], prediction_val, configs[['metric']])
    
    
    if(give_residuals){
      residual.track <- list(residual)
      residual.track.val <- list(residual_val)
    }
    earlyStopNow <- FALSE
    m_batch_list=numeric()
    
    cat("\n")
  # end of pre-processing
 
  if(! earlyStopNow){
  for (iter in 1:b_max) {
    time.iter.start <- Sys.time()
    snpboostLogger(paste0("Iteration ", iter), log.time=time.iter.start)

   ### --- Update the feature matrix --- ###
    if (configs[['verbose']]) snpboostLogger("Start updating feature matrix ...", indent=1)
    time.update.start <- Sys.time()

    sorted.score <- sort(score, decreasing = T, na.last = NA)
    
    if (length(sorted.score) > 0) {
      features.to.use <- c(covariates,names(sorted.score)[1:min(p_batch, length(sorted.score))])
      
      #features matrix set to NULL
      features.to.discard <- setdiff(colnames(features[['train']]),features.to.use)
      
      if (length(features.to.discard) > 0) {
        for(s in splits){features[[s]][,features.to.discard]= NULL}
      }
      
      features.to.add <- setdiff(features.to.use,colnames(features[['train']]))
      
      if(length(features.to.add)>0){
        for(s in splits){
          tmp.features.add <- prepareFeatures(pgen[[s]], vars, features.to.add, stats)
          if (!is.null(features[[s]])) {
            features[[s]][, colnames(tmp.features.add) := tmp.features.add]
          } else {
            features[[s]] <- tmp.features.add
          }
          rm(tmp.features.add)
        }
      }
    } else {
      break
    }

    if (configs[['verbose']]) snpboostLoggerTimeDiff("End updating feature matrix.", time.update.start, indent=2)
    
   
    ### --- Boosting --- ###
    if (configs[['verbose']]) snpboostLogger("Start Boosting ...", indent=1)
    save(list=ls(),file="test.RData")
    lm_boost_tmp <- lm_boost(y = residual, x = features[['train']],validation = validation, y_val = residual_val, x_val = features[['val']],
                             beta = beta, m_batch = m_batch, sl = sl, covariates = covariates,
                             coeff_path_save = coeff_path_save,give_residuals=give_residuals,
                             c_stop=sorted.score[min(p_batch+1, length(sorted.score))],configs=configs,
                             phe = phe, phenotype = phenotype, prediction = prediction, prediction_val = prediction_val)
    m_batch_list[iter] <- lm_boost_tmp$stop
    residual <- lm_boost_tmp$residuals
    prediction <- lm_boost_tmp$prediction
    prediction_val <- lm_boost_tmp$prediction_val
    metric.train <- c(metric.train, lm_boost_tmp$metric.train)
    residual_val <- if(validation){lm_boost_tmp$residuals_val}else{NULL}
    metric.val <- c(metric.val, lm_boost_tmp$metric.val)
    beta <- lm_boost_tmp$beta
    for(i in 1:m_batch_list[iter]){
     if(give_residuals){
        residual.track[[length(residual.track)+1]] = lm_boost_tmp$residual.track[[i]]
        if(validation){
          residual.track.val[[length(residual.track.val)+1]] = lm_boost_tmp$residual.track.val[[i]]
        }
      }
    }

    rm(lm_boost_tmp)
    
    if (configs[['verbose']]) snpboostLogger("End Boosting.", indent=1)

    ### --- Update Score --- ###
    if (configs[['verbose']]) snpboostLogger("Start updating score ...", indent=1)
    
    prod.full <- computeProduct(residual, genotype.pfile, vars, stats, configs, iter=iter)
    score <- abs(prod.full[, 1])
    rm(prod.full)

    if (configs[['verbose']]) snpboostLogger("End updating score.", indent=1)

    if (configs[['save']]) {
      save(metric.train, metric.val, beta, score, configs, residual, residual_val, m_batch_list,
           file = file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", iter, ".RData")))
    }

    time.iter.end <- Sys.time()
    snpboostLoggerTimeDiff(paste0("End iteration ", iter, '.'), time.iter.start, time.iter.end, indent=1)
    snpboostLoggerTimeDiff("The total time since start.", time.start, time.iter.end, indent=2)

    ### --- Check stopping criteria --- ####
    if(validation && checkEarlyStoppingBatches(metric.val,iter,b_stop,configs,m_batch_list))break
  }
    n_iter = iter
    n_step = if(configs[['metric']]=='MSEP' || configs[['metric']]=='log loss') which.min(metric.val)
  }
  snpboostLoggerTimeDiff("End SNP boost.", time.start)
  if(! configs[['save']]) cleanUpIntermediateFiles(configs)
  if(configs[['verbose']]) print(gc())

  out <- list(metric.train = metric.train, metric.val = metric.val, residuals = if(give_residuals){residual}else{NULL},
              residuals_val = if(give_residuals){residual_val}else{NULL}, beta = beta, configs = configs, stats = stats, residual.track = if(give_residuals){residual.track}else{NULL},
              residual.track.val = if(give_residuals){residual.track.val}else{NULL},n_iter= n_iter, n_step = n_step, m_batch_list=m_batch_list)
  out
}

