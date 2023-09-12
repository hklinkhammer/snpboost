lm_boost <- function(y,x,validation=FALSE,y_val=NULL,x_val=NULL,beta=data.frame(name=NA,intercept=0,variant=NA),m_batch=1,sl=0.1,
                     covariates=NULL,coeff_path_save=TRUE,give_residuals=FALSE,c_stop=0,configs=NULL,
                     phe=NULL,phenotype=NULL,prediction=NULL,prediction_val=NULL,
                     surv=NULL,IPCw=NULL,w_denom=NULL, T_order = NULL, wweights = NULL, sigma = NULL){

  #Initialization
  r <- y
  p <- ncol(x)
  nobs <- length(y)
  metric.train <- numeric(m_batch)
  variants <- colnames(x)[!colnames(x) %in% covariates]
  pred <- prediction
  pred_tmp <- order(pred)
  
  if(coeff_path_save | validation){
    previous.beta <- get_coefficients(coeff_path=beta,covariates=covariates)
    previous.beta[variants[!variants %in% names(previous.beta)]]=0
    previous.beta <- previous.beta[c("(Intercept)",covariates,variants)]
    
    b <- beta %>% add_row(name=rep(NA,m_batch),intercept=rep(0,m_batch),variant=rep(0,m_batch))
    
  }else{ 
    b <- beta
  }
  x_std <- as.matrix(dataPreparation::fast_scale(x%>%select(all_of(variants)),verbose=F))
  X <- as.matrix(cbind(rep(1,nrow(x)), x))
  colnames(X)=c("(Intercept)",colnames(x))
  
  if(give_residuals){
    residual.track <- list()
  }
  
  if(validation){
    r_val <- y_val
    X_val <- as.matrix(cbind(rep(1,nrow(x_val)),x_val))
    metric.val <- numeric()
    pred.val <- prediction_val
    pred.val_tmp <- order(pred.val)
    if(give_residuals){
      residual.track.val <- list()
    }
  }else{
    r_val <- NULL
    X_val <- NULL
    metric.val <- NULL
  }

  #Boosting
  
  for(m in 1:m_batch){
   #find highest correlation
    if(p==1){
      i_best <- 1
    }else{
      
      score <- abs(Rfast::Crossprod(x_std,matrix((r-mean(r))/sd(r),ncol=1)))

      i_best <- which.max(score)

      if(m>1 & score[i_best]/(nobs-1)<c_stop) break
    }
    
    #Boosting

    lm_tmp <- .lm.fit(x=X[,c("(Intercept)",covariates,variants[i_best])],y=r)

    if(coeff_path_save | (validation)){
      #update previous beta
      previous.beta[c("(Intercept)",covariates,variants[i_best])] <- previous.beta[c("(Intercept)",covariates,variants[i_best])]+sl*coef(lm_tmp)
      #save new intercept and fitted coefficient
      b[nrow(beta)+m,1] <- variants[i_best]
      b[nrow(beta)+m,-1] <- c(previous.beta[c("(Intercept)",covariates)],previous.beta[variants[i_best]])
     
    }else{
      #update previous beta
      previous.beta[c("(Intercept)",covariates,variants[i_best])] <- previous.beta[c("(Intercept)",covariates,variants[i_best])]+sl*coef(lm_tmp)
      #save new intercept and fitted coefficient
      b[nrow(beta)+m,] <- c(variants[i_best],previous.beta[c("(Intercept)",covariates)],previous.beta[variants[i_best]])
      
    }
    
    
    if(length(covariates)>0){
      pred <- pred + sl*coef(lm_tmp)[1] + sl * Rfast::mat.mult(X[,c(covariates,variants[i_best])],matrix(coef(lm_tmp)[-1]))
     
    }else{
      pred <- pred + sl*coef(lm_tmp)[1] + sl* X[,i_best+1]*coef(lm_tmp)[2]
      
    }
    
    if(configs[['family']]=="cox"){
      r <- computeResiduals(surv[['train']],pred,configs,IPCw[['train']],w_denom[['train']],T_order[['train']],wweights[['train']],sigma)
      
      if(configs[['metric']] %in% c("AFT-Weibull", "AFT-logistic", "AFT-normal")){
        sigma <- optimizeSigma(surv[['train']], pred, configs)
      }
    }else{
      r <- computeResiduals(phe[['train']][[phenotype]], pred, configs, sigma = sigma)
      if(configs[['metric']] %in% c("negative_binomial")){
        sigma <- optimizeSigma(phe[['train']][[phenotype]], pred, configs)
      }
    } 
    if(give_residuals){
      residual.track[[m]] = r
    }
    if(configs[['family']]=="cox"){
      if(configs[['metric']] == 'C'){
        if(m>1 & all(order(pred)==pred_tmp)){
          metric.train[m] = metric.train[m-1]
        }else{
          metric.train[m] = computeMetric(r, surv[['train']], pred, configs, surv[['train']], wweights[['train']],sigma, T_order[['train']])
          pred_tmp <- order(pred)
        }
      }else{
        metric.train[m] = computeMetric(r, surv[['train']], pred, configs, surv[['train']], IPCw[['train']],sigma, T_order[['train']])
        pred_tmp <- order(pred)
      }
    }else{
      metric.train[m] = computeMetric(r, phe[['train']][[phenotype]], pred, configs, sigma = sigma)
    }
    
    if(validation){
      if(length(covariates)>0){
        
        pred.val <- pred.val + sl*coef(lm_tmp)[1] + sl * Rfast::mat.mult(X_val[,c(covariates,variants[i_best])],matrix(coef(lm_tmp)[-1]))
        
      }else{
        
        pred.val <- pred.val + sl*coef(lm_tmp)[1] + sl* X_val[,i_best+1]*coef(lm_tmp)[2]
        
      }
      if(configs[['family']]=="cox"){
        
        r_val <- computeResiduals(surv[['val']],pred.val,configs,IPCw[['val']],w_denom[['val']],T_order[['val']],wweights[['val']],sigma)
        
      }else{
        
        r_val <- computeResiduals(phe[['val']][[phenotype]], pred.val, configs, sigma = sigma)
        
      }
      
      if(give_residuals){
        
        residual.track.val[[m]]=r_val
        
      }
      if(configs[['family']]=="cox"){
        if(configs[['metric']] == 'C'){
          if(m>1 & all(order(pred.val)==pred.val_tmp)){
            metric.val[m] <- metric.val[m-1]
          }else{
            metric.val[m] <- computeMetric(r_val, surv[['val']], pred.val, configs, surv[['val']], wweights[['val']], sigma, T_order[['val']])
            pred.val_tmp <- order(pred.val)
          }
        }else{
          metric.val[m] <- computeMetric(r_val, surv[['val']], pred.val, configs, surv[['val']], IPCw[['val']], sigma)
          pred.val_tmp <- order(pred.val)
        }
      }else{
        metric.val[m] <- computeMetric(r_val, phe[['val']][[phenotype]], pred.val, configs, sigma = sigma)
      }
    }
    
    rm(lm_tmp)
    
  }
  
  
  stop= if(validation){length(metric.val)}else{length(metric.train)}
  
  b <- b[1:(nrow(beta)+stop),]
  original_beta <- get_coefficients(beta,covariates=covariates)
  original_beta[variants[!variants %in% names(original_beta)]]=0
  original_beta <- original_beta[c("(Intercept)",covariates,variants)]

  new_beta <- get_coefficients(b,covariates=covariates)
  new_beta[variants[!variants %in% names(new_beta)]]=0
  new_beta <- new_beta[c("(Intercept)",covariates,variants)]
  
  pred <- prediction + X%*%(new_beta-original_beta)
 
  if(configs[['family']]=="cox"){
    r <- matrix(computeResiduals(surv[['train']],pred,configs,IPCw[['train']],w_denom[['train']],T_order[['train']],wweights[['train']], sigma),ncol=1)
    if(configs[['metric']] %in% c("AFT-Weibull", "AFT-logistic", "AFT-normal")){
      sigma <- optimizeSigma(surv[['train']], pred, configs)
    }
  }else{
    r <- matrix(computeResiduals(phe[['train']][[phenotype]], pred, configs,sigma=sigma),ncol=1)
    if(configs[['metric']] %in% c("negative_binomial")){
      sigma <- optimizeSigma(phe[['train']][[phenotype]], pred, configs)
    }
  }
  rownames(r) <- rownames(phe[['train']])
  colnames(r) <- c('0')
  
  metric.train <- if(coeff_path_save){metric.train[1:stop]}else{metric.train[stop]}

  
  if(validation){
    pred.val <- prediction_val + X_val%*%(new_beta - original_beta)
    if(configs[['family']]=="cox"){
      r_val <- matrix(computeResiduals(surv[['val']],pred.val,configs,IPCw[['val']],w_denom[['val']],T_order[['val']],wweights[['vals']],sigma),ncol=1)
    }else{
      r_val <-  matrix(computeResiduals(phe[['val']][[phenotype]],pred.val,configs, sigma = sigma),ncol=1)
    }
    rownames(r_val) <- rownames(phe[['val']])
    colnames(r_val) <- c('0')
  }else{
    r_val = NULL
    pred.val = NULL
  }
  metric.val <- if(validation){if(coeff_path_save){metric.val[1:stop]}else{metric.val[stop]}}else{NULL}
  
  b <- if(coeff_path_save){b[1:(nrow(beta)+stop),]}else{b[1+(nrow(beta)+stop),]}

  out <- list(beta=b, residuals=r, metric.train=metric.train, residuals_val=r_val, metric.val=metric.val,
              residual.track=if(give_residuals){residual.track}else{NULL}, 
              residual.track.val=if(give_residuals){residual.track.val}else{NULL},
              stop=stop, prediction = pred, prediction_val = pred.val, sigma=if(configs[['family']]=="cox" | configs[['metric']]=="negative_binomial"){sigma}else{NULL})
  suppressWarnings(rm(b,r,metric.train,r_val,metric.val,y,beta,X,y_val,x_val,stop,new_beta,original_beta,pred.val,pred,sigma))
  
  return(out)
}
