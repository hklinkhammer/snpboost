library(data.table)
library(snpnet)
library(tidyverse)

predict_snpnet_hk <- function(fit, new_genotype_file, new_phenotype_file, phenotype, subset=NULL){
  if(!is.null(fit$metric.val)){
    idx <- which.max(fit$metric.val)
  }else{
    idx <- which.max(fit$metric.train)
  }
  final_betas <- fit$beta[[idx]]
  final_betas <- final_betas[final_betas!=0]
  chosen_SNPs <- cbind(str_split(names(final_betas),"_",simplify = T),
                       final_betas)
  rownames(chosen_SNPs)=rep("",nrow(chosen_SNPs))
  fwrite(chosen_SNPs,paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),sep="\t")
  
  plink2_cmd <- paste(fit$configs['plink2.path'],"--pfile",new_genotype_file,"vzs","--score",
                      paste0(fit$configs['results.dir'],"/chosen_SNPs_predict.txt"),1,2,3,"header cols=maybefid,nallele,denom,dosagesum,scoreavgs,scoresums","--out",paste0(fit$configs['results.dir'],"/PRS"))
  system(plink2_cmd, intern=F, wait=T)
  
  phe <- fread(new_phenotype_file)
  if(!is.null(subset)){
    phe <- phe[phe[[fit$configs[['split.col']]]]==subset,]
  }
  PRS <- fread(paste0(fit$configs['results.dir'],"/PRS.sscore")) %>% rename(PRS = SCORE1_SUM)
  
  phe_PRS <- left_join(phe,PRS,by="IID")
  
  pred <- matrix(fit_snpnet$a0[[idx]] + phe_PRS$PRS + if(!is.null(fit$configs[['covariates']])){as.matrix(phe_PRS%>%select(fit$configs[['covariates']])) %*% final_betas[fit$configs[['covariates']]]}else{0},ncol=1)
  
  rownames(pred)=phe_PRS$IID
  
  residuals <- phe_PRS[[phenotype]]-pred
  
  rownames(residuals)=phe_PRS$IID
  
  metric <- mean((phe_PRS[[phenotype]]-pred)^2)
  
  out <- list(prediction = pred, residuals=residuals, metric = metric)
}

#enter file paths to genotype file (prefix for pgen-files) and phenotype file (.phe)
genotype.pfile <- "ukb_cal_allChrs_XY"
phenotype.file <- "phenotypes.phe"
#enter name of phenotype (should be a column header in phe-file)
phenotype <- "height"
#enter covariates to be included in the fitting process
covariates <- NULL

configs <- list(plink2.path = "plink2", #plink2 path
                zstdcat.path = "zstdcat", #zstdcat path
                results.dir = paste0("non_imputed_data_snpnet_",phenotype),
                save = TRUE,
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.001,
                num.snps.batch = 1000,
                nlams.init = 10,
                nlams.delta = 5,
                early.stopping = TRUE,
                stopping.lag = 2,
                verbose = TRUE,
                KKT.verbose = TRUE,
                rank = TRUE,
                use.glmnetPlus = TRUE,
                mem=128000,
                niter=100,
                nCores=16,
                standardize_variant=T
                )


# fit fit_snpnet

time_snpnet_start <- Sys.time()

fit_snpnet <- snpnet(genotype.pfile = genotype.pfile,
                     phenotype.file = phenotype.file,
                     phenotype = phenotype,
                     covariates = covariates,
                     configs = configs,
                     split.col = "train_test"
                     ) 

time_snpnet_end <- Sys.time()

### predict on test set/calculate prediction for all samples (for glm)

pred_test_snpnet <- predict_snpnet_hk(fit=fit_snpnet,new_genotype_file=genotype.pfile,new_phenotype_file=phenotype.file,
                                        phenotype=phenotype,
                                        subset = "test")

pred_all_snpnet <- predict_snpnet_hk(fit=fit_snpnet,new_genotype_file=genotype.pfile,new_phenotype_file=phenotype.file,
                                       phenotype=phenotype)
### extract optimal coefficients

idx <- which.max(fit_snpnet$metric.val)

nmb_snps <- sum(fit_snpnet$beta[[idx]]!=0)

## fit glm with PRS and covariates

data <- fread(phenotype.file) %>% mutate(IID = as.character(IID))
pred <- data.table(IID=rownames(pred_all_snpnet$prediction), PRS =pred_all_snpnet$prediction)

data <- full_join(data,pred)%>%rename(PRS=PRS.V1) %>% filter(complete.cases(.))

formula_model <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep=""),"+ PRS"))

full_model <- glm(formula_model,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_full_model_test <- predict.glm(full_model, newdata = data %>% filter(train_test=="test"))

### extract MSEP and R^2 based on full model

MSEP_test_full_model=mean((data[data$train_test=="test",]%>%pull(all_of(phenotype))-pred_full_model_test)^2)

cor_squared_test_full_model= cor(data[data$train_test=="test",]%>%pull(all_of(phenotype)),pred_full_model_test)^2

### extract MSEP and R^2 based only on PRS

MSEP_only_PRS <- pred_test_snpnet$metric

cor_squared_only_PRS <- cor(data[data$train_test=="test",]%>%pull(all_of(phenotype)),pred_test_snpnet$prediction)^2

# model only covariates

formula_model_cov <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep="")))

cov_model <- glm(formula_model_cov,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_cov_model_test <- predict.glm(cov_model, newdata = data %>% filter(train_test=="test"))

### extract MSEP and R^2 based only on covariates

MSEP_test_cov_model=mean((data[data$train_test=="test",]%>%pull(all_of(phenotype))-pred_cov_model_test)^2)

cor_squared_test_cov_model= cor(data[data$train_test=="test",]%>%pull(all_of(phenotype)),pred_cov_model_test)^2

save.image(file=paste0("results/",phenotype,"_snpnet_std.RData"))
