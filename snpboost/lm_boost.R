lm_boost <- function(y,x,validation=FALSE,y_val=NULL,x_val=NULL,beta=data.frame(name=NA,intercept=0,variant=NA),m_batch=1,sl=0.1,
                     covariates=NULL,coeff_path_save=TRUE,give_residuals=FALSE,c_stop=0,configs=NULL){

  #Initialization
  r <- y
  p <- ncol(x)
  nobs <- length(y)
  metric.train <- numeric(m_batch)
  variants <- colnames(x)[!colnames(x) %in% covariates]
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
      score <- abs(Rfast::Crossprod(x_std,(r-mean(r))/sd(r)))
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
        r_val <- r_val - sl*coef(lm_tmp)[1]- sl* X_val[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
      }else{
        r_val <- r_val - sl*coef(lm_tmp)[1]- sl* X_val[,i_best+1]*coef(lm_tmp)[2]
      }
      if(give_residuals){
      residual.track.val[[m]]=r_val
      }
      metric.val[m] <- computeMetric(r_val,configs[['metric']])
    }
    if(length(covariates)>0){
      r <- r - sl*coef(lm_tmp)[1]- sl* X[,c(covariates,variants[i_best])]%*%coef(lm_tmp)[-1]
    }else{
      r <- r - sl*coef(lm_tmp)[1]- sl* X[,i_best+1]*coef(lm_tmp)[2]
    }
    if(give_residuals){
      residual.track[[m]] = r
    }
    metric.train[m] = computeMetric(r,configs[['metric']])
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

  r <- y - X%*%(new_beta-original_beta)
  metric.train <- if(coeff_path_save){metric.train[1:stop]}else{metric.train[stop]}

  
  r_val <-  if(validation){y_val - X_val%*%(new_beta - original_beta)}else{NULL}
  metric.val <- if(validation){if(coeff_path_save){metric.val[1:stop]}else{metric.val[stop]}}else{NULL}
  
  b <- if(coeff_path_save){b[1:(nrow(beta)+stop),]}else{b[1+(nrow(beta)+stop),]}

  out <- list(beta=b, residuals=r, metric.train=metric.train, residuals_val=r_val, metric.val=metric.val,
              residual.track=if(give_residuals){residual.track}else{NULL}, 
              residual.track.val=if(give_residuals){residual.track.val}else{NULL},
              stop=stop)
  suppressWarnings(rm(b,r,metric.train,r_val,metric.val,y,beta,X,y_val,x_val,stop,new_beta,original_beta))
  
  return(out)
}
