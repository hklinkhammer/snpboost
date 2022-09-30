lm_boost <- function(y,x,validation=FALSE,y_val=NULL,x_val=NULL,beta=data.frame(name=NA,intercept=0,variant=NA),m_batch=1,sl=0.1,
                     covariates=NULL,coeff_path_save=TRUE,give_residuals=FALSE,c_stop=0,configs=NULL,
                     phe=NULL,phenotype=NULL,prediction=NULL,prediction_val=NULL){

  #Initialization
  r <- y
  p <- ncol(x)
  nobs <- length(y)
  metric.train <- numeric(m_batch)
  variants <- colnames(x)[!colnames(x) %in% covariates]
  pred <- prediction
  if(coeff_path_save | validation){
    previous.beta <- get_coefficients(coeff_path=beta,covariates=covariates)
    previous.beta[variants[!variants %in% names(previous.beta)]]=0
    previous.beta <- previous.beta[c("(Intercept)",covariates,variants)]
    
    b <- beta %>% add_row(name=rep(NA,m_batch),intercept=rep(0,m_batch),variant=rep(0,m_batch))
    
  }else{ ### hat bisher keine Funktion
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
      if(score[i_best]/(nobs-1)<c_stop) break
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
    if(validation){
      if(length(covariates)>0){
        pred.val <- pred.val + sl*coef(lm_tmp)[1] + sl* X_val[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
       # if(configs[['family']]=="gaussian"){
       #   r_val <- r_val - sl*coef(lm_tmp)[1]- sl* X_val[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
       # }else if (configs[['family']]=="binomial"){
        #   p.hat.val <- exp(pred.val)/(1+exp(pred.val))
        #   r_val <- 1/length(p.hat.val)*(phe[['val']][[phenotype]]-p.hat.val)/(p.hat.val*(1-p.hat.val))
       # }
      }else{
        pred.val <- pred.val + sl*coef(lm_tmp)[1] + sl* X_val[,i_best+1]*coef(lm_tmp)[2]
       # if(configs[['family']]=="gaussian"){
       #   r_val <- r_val - sl*coef(lm_tmp)[1]- sl* X_val[,i_best+1]*coef(lm_tmp)[2]
      #  }else if (configs[['family']]=="binomial"){
       #   p.hat.val <- exp(pred.val)/(1+exp(pred.val))
       #   r_val <- 1/length(p.hat.val)*(phe[['val']][[phenotype]]-p.hat.val)/(p.hat.val*(1-p.hat.val))
       # }
      }
      r_val <- computeResiduals(phe[['val']][[phenotype]], pred.val, configs[['family']])
      
      if(give_residuals){
       residual.track.val[[m]]=r_val
      }
      metric.val[m] <- computeMetric(r_val, phe[['val']][[phenotype]], pred.val, configs[['metric']])
    }
    
    if(length(covariates)>0){
      pred <- pred + sl*coef(lm_tmp)[1] + sl* X[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
     # if(configs[['family']]=="gaussian"){
    #    r <- r - sl*coef(lm_tmp)[1]- sl* X[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
     # }else if (configs[['family']]=="binomial"){
     #   p.hat <- exp(pred)/(1+exp(pred))
    #    r <- 1/length(p.hat)*(phe[['train']][[phenotype]]-p.hat)/(p.hat*(1-p.hat))
     # }
    }else{
      pred <- pred + sl*coef(lm_tmp)[1] + sl* X[,i_best+1]*coef(lm_tmp)[2]
      #if(configs[['family']]=="gaussian"){
      #  r <- r - sl*coef(lm_tmp)[1]- sl* X[,i_best+1]*coef(lm_tmp)[2]
      #}else if (configs[['family']]=="binomial"){
       # p.hat <- exp(pred)/(1+exp(pred))
        #r <- 1/length(p.hat)*(phe[['train']][[phenotype]]-p.hat)/(p.hat*(1-p.hat))
      #}
    }
    r <- computeResiduals(phe[['train']][[phenotype]], pred, configs[['family']])
   
    if(give_residuals){
      residual.track[[m]] = r
    }
    metric.train[m] = computeMetric(r, phe[['train']][[phenotype]], pred, configs[['metric']])
    rm(lm_tmp)

  }
  
  stop= length(metric.val)
  
  b <- b[1:(nrow(beta)+stop),]
  original_beta <- get_coefficients(beta,covariates=covariates)
  original_beta[variants[!variants %in% names(original_beta)]]=0
  original_beta <- original_beta[c("(Intercept)",covariates,variants)]

  new_beta <- get_coefficients(b,covariates=covariates)
  new_beta[variants[!variants %in% names(new_beta)]]=0
  new_beta <- new_beta[c("(Intercept)",covariates,variants)]
  
  pred <- prediction + X%*%(new_beta-original_beta)
 # if(configs[['family']]=="gaussian"){
  #  r <- y - X%*%(new_beta-original_beta)
  #}else if (configs[['family']]=="binomial"){
   # p.hat <- exp(pred)/(1+exp(pred))
  #  r <- matrix(1/length(p.hat)*(phe[['train']][[phenotype]]-p.hat)/(p.hat*(1-p.hat)),ncol=1)

 # }
  r <- matrix(computeResiduals(phe[['train']][[phenotype]], pred, configs[['family']]),ncol=1)
  rownames(r) <- rownames(phe[['train']])
  colnames(r) <- c('0')
  
  metric.train <- if(coeff_path_save){metric.train[1:stop]}else{metric.train[stop]}

 
 # if(configs[['family']]=="gaussian"){
#    r_val <-  if(validation){y_val - X_val%*%(new_beta - original_beta)}else{NULL}
 # }else if (configs[['family']]=="binomial"){
 #   p.hat.val <- exp(pred.val)/(1+exp(pred.val))
 #   r_val <-  if(validation){matrix(1/length(p.hat.val)*(phe[['val']][[phenotype]]-p.hat.val)/(p.hat.val*(1-p.hat.val)),ncol=1)}else{NULL}

#  }
  
  if(validation){
    pred.val <- prediction_val + X_val%*%(new_beta - original_beta)
    r_val <-  matrix(computeResiduals(phe[['val']][[phenotype]],pred.val,configs[['family']]),ncol=1)
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
              stop=stop, prediction = pred, prediction_val = pred.val)
  suppressWarnings(rm(b,r,metric.train,r_val,metric.val,y,beta,X,y_val,x_val,stop,new_beta,original_beta,pred.val,pred))
  
  return(out)
}
