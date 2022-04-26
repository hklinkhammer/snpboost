library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(batchtools)

setwd("100k_sp_01_h_50")

# help function for batchtools

reg = makeRegistry(file.dir = "100k_sp_01_h_50/Registry", seed = 1)
#reg = loadRegistry("/daten/epi/PRS/snpnet/SNP_boost/Scripts/100k_sp_1_h_10/Registry", writeable=T)


run_snpboost <- function(index,isimu){
  library(data.table)
  library(tidyverse)
  library(parallel)
  library(bigsnpr)
  library(ggplot2)
  library(cowplot)
  source("lm_boost.R")
  source("SNP_boost.R")
  source("functions_snpboost.R")
  source("functions_snpnet.R")
  
  #read arguments
  
  scene="heritability_50_sparsity_01_lr_01"
  folder=paste0("Simulations_",scene)
  n=100000
  p=100000
  n_true=100
  h_squared=0.50
  genotypes_folder="genotypes"
  genotypes=paste0(genotypes_folder,"/genotype")
  phenotypes=paste0(folder,"/pheno_simulated")
  phenotype="pheno"

  # lade Designmatrix
  load(paste0(folder,"/designmatrix.RData"))
  design <- design %>% mutate(Rsquare_testset = NA,
                              MSEP_valset = NA,
                              MSEP_testset = NA,
                              n_SNPs = NA,
                              Runtime = NA,
                              iter = NA,
                              steps = NA
                              )
  
  # Kovariaten festlegen
  covariates <- NULL
  
  
    genotype.pfile= paste0(genotypes,"_",isimu)
    phenotype.file = paste0(phenotypes,"_",isimu,".phe")
    
    configs <- list(plink2.path = "plink2",
                    zstdcat.path = "zstdcat",
                    zcat.path = "zstdcat",
                    results.dir = paste0(paste(paste0(folder,"/Simulation"),
                                        design$p_batch[index],design$m_batch[index],design$learning_rate[index],
                                        design$b_stop[index],sep="_"),"/pheno_",isimu),
                    save = FALSE,
                    prevIter = 0,
                    missing.rate = 0.1,
                    MAF.thresh = 0.001,
                    early.stopping = TRUE,
                    verbose = TRUE,
                    rank = TRUE,
                    mem=10000,
                    nCores=10,
                    standardize.variant=TRUE
    )
    
    start.time <- Sys.time()
    fit_snpboost <- try(snp_boost(genotype.pfile = genotype.pfile,
                            phenotype.file = phenotype.file,
                            phenotype = phenotype,
                            covariates = covariates,
                            configs = configs,
                            split.col = "split",
                            p_batch = design$p_batch[index], # batch size
                            m_batch = design$m_batch[index], # maximum number of boosting iterations per batch
                            b_max = 20000, # maximum number of batches
                            b_stop = design$b_stop[index], # outer stopping
                            sl= design$learning_rate[index], # learning rate
                            coeff_path_save = TRUE,
                            give_residuals = FALSE
                      )) 
    end.time <- Sys.time()
    
    results = design[index,]
    
    results$Runtime = paste(as.numeric(end.time-start.time),units(end.time-start.time))
    
    idx  <- which.min(fit_snpboost$metric.val)
    
    final_betas <- get_coefficients(fit_snpboost$beta,idx)
    chosen_SNPs <- cbind(str_split(names(final_betas)[-1],"_",simplify = T),
                         final_betas[-1])
    rownames(chosen_SNPs)=rep("",nrow(chosen_SNPs))
    fwrite(data.table(chosen_SNPs),paste0(fit_snpboost$configs[['results.dir']],"/chosen_SNPs.txt"),sep="\t")
    fwrite(list(fit_snpboost$m_stop_list),paste0(fit_snpboost$configs[['results.dir']],"/m_stops.txt"),sep="\t")
    
    results$n_SNPs=length(final_betas)-1
    results$MSEP_valset=fit_snpboost$metric.val[idx]
    
    results$MSEP_testset <- try(predict_SNP_boost(fit_snpboost,genotype.pfile,phenotype.file,phenotype,subset="test")$metric)
    
    results$Rsquare_testset <- try(predict_SNP_boost(fit_snpboost,genotype.pfile,phenotype.file,phenotype,subset="test")$Rsquare)
    
    results$iter = fit_snpboost$n_iter
    results$steps = fit_snpboost$n_step
    
    fwrite(results,paste0(fit_snpboost$configs[['results.dir']],"/results.txt"),sep="\t")
    
    return(results)
}

n_simu=100
n_design=8

arguments <- crossing(index=1:n_design,isimu=1:n_simu)
arguments <- arguments[order(arguments$isimu),]
# run batchtools

ids_done <- findDone()

clearRegistry()

ids <- batchMap(fun = run_snpboost, args = arguments, reg=reg)

ids <- ids %>% filter(!job.id %in% ids_done$job.id)

done <- submitJobs(ids=ids, reg=reg, res=list(ncpus=10,memory=2000,max.concurrent.jobs=8,walltime=6*60*60))

rm(list=ls())

