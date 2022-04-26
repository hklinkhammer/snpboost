library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(readr)

source("lm_boost.R")
source("SNP_boost.R")
source("functions_snpboost.R")
source("functions_snpnet.R")

#enter file paths to genotype file (prefix for pgen-files) and phenotype file (.phe)
genotype.pfile <- "ukb_cal_allChrs_XY"
phenotype.file <- "phenotypes.phe"
#enter name of phenotype (should be a column header in phe-file)
phenotype <- "height"
#enter covariates to be included in the fitting process
covariates <- NULL


configs <- list(plink2.path = "plink2", #plink2 path
                zstdcat.path = "zstdcat", #zstdcat path
                results.dir = paste0("non_imputed_data_",phenotype),
                save = FALSE,
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.001,
                early.stopping = TRUE,
                verbose = TRUE,
                mem=16000,
                nCores=16,
                standardize.variant=TRUE
)


###fit fit_snpboost
time_snpboost_start <- Sys.time()

fit_snpboost <- snp_boost(genotype.pfile = genotype.pfile,
                     phenotype.file = phenotype.file,
                     phenotype = phenotype,
                     covariates = covariates,
                     configs = configs,
                     split.col = "train_test",
                     p_batch = 1000, # batch size
                     m_batch = 1000, # maximum number of boosting iterations per batch
                     b_max = 20000, # maximum number of batches
                     b_stop = 2, # outer stopping
                     sl= 0.1, # learning rate
                     coeff_path_save = TRUE,
                     give_residuals = FALSE
                     ) 


time_snpboost_end <- Sys.time()

### predict on test set/calculate prediction for all samples (for glm)

pred_test_snpboost <- predict_SNP_boost(fit_snpboost,genotype.pfile,phenotype.file,phenotype,subset="test")

pred_all_snpboost <- predict_SNP_boost(fit_snpboost,genotype.pfile,phenotype.file,phenotype)

### extract optimal coefficients

idx <- which.min(fit_snpboost$metric.val)
betas <- get_coefficients(fit_snpboost$beta,idx,covariates=fit_snpboost$configs[['covariates']])
nmb_snps <- length(betas)-1-length(fit_snpboost$configs[['covariates']])

## fit glm with PRS and covariates

data <- fread(phenotype.file) %>% mutate(IID = as.character(IID))
pred <- data.table(IID=rownames(pred_all_snpboost$prediction), PRS =pred_all_snpboost$prediction)

data <- full_join(data,pred)%>%rename(PRS=PRS.V1) %>% filter(complete.cases(.))

formula_model <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep=""),"+ PRS"))

full_model <- glm(formula_model,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_full_model_test <- predict.glm(full_model, newdata = data %>% filter(train_test %in% c("test")))

### extract MSEP and R^2 based on full model

MSEP_test_full_model=mean((data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype))-pred_full_model_test)^2)

cor_squared_test_full_model=cor(data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype)),pred_full_model_test)^2

### extract MSEP and R^2 based only on PRS

MSEP_only_PRS <- pred_test_snpboost$metric

cor_squared_only_PRS <- cor(data[data$train_test=="test",]%>%pull(all_of(phenotype)),pred_test_snpboost$prediction)^2

# model only covariates

formula_model_cov <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep="")))

cov_model <- glm(formula_model_cov,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_cov_model_test <- predict.glm(cov_model, newdata = data %>% filter(train_test=="test"))

### extract MSEP and R^2 based only on covariates

MSEP_cov_model=mean((data[data$train_test=="test",]%>%pull(all_of(phenotype))-pred_cov_model_test)^2)

cor_squared_test_cov_model=cor(data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype)),pred_cov_model_test)^2


save.image(file=paste0("results/",phenotype,"_snpboost.RData"))