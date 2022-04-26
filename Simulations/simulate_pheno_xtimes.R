library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)

set.seed(2021)

setwd("genotypes")

# Input
n=100000
p=100000

scene = "heritability_50_sparsity_01_lr_01"
genotypes_folder = "genotypes"
n_true = 100
h_squared = 0.5
folder = paste0("Simulations_",scene)
n_simu = 100


for(i in 1:n_simu){

  fam <- paste0(genotypes_folder,"/n_",n,"_p_",p,"/fam_",i,".fam")
  pheno_names <- fread(fam)
  
  bim <- fread(paste0(genotypes_folder,"/n_",n,"_p_",p,"/bim_",i,".bim"))
  
  pheno_simulated <- list(FID= pheno_names$V1, IID = pheno_names$V2, 
                          split=c(rep("train",0.5*nrow(pheno_names)),rep("val",0.2*nrow(pheno_names)), rep("test",0.3*nrow(pheno_names))))
  
  rds <- snp_readBed2(paste0(genotypes_folder,"/genotype_general_final_filtered.bed"),ind.row=pheno_names$index,ind.col=bim$index,ncores=25)
  geno <- snp_attach(paste0(genotypes_folder,"/genotype_general_final_filtered.rds"))
  
  # simulate and save phenotype
  pheno <- snp_simuPheno(G = geno$genotypes, 
                         h2= h_squared, # heritability
                         M = n_true, # causal variants
                         ncores = 25 # Cores
  )
  pheno_simulated[["pheno"]]=pheno$pheno
  
  vars <- bim %>% rename(rsID=V2, Chr=V1, Pos = V4) %>%
    mutate(Indices = 1:nrow(.))  
  fwrite(list(Indices=bim$index[pheno$set],rsID=vars$rsID[pheno$set],betas=pheno$effects),paste0(folder,"/true_betas/true_betas_simulation_",i,".txt"),sep="\t")
  fwrite(pheno_simulated,paste0(folder,"/pheno_simulated_",i,".phe"), sep="\t")  
  
  cmd_rm <- paste("rm", paste0(genotypes_folder,"/genotype_general_final_filtered.bk"))
  system(cmd_rm)
  cmd_rm <- paste("rm", paste0(genotypes_folder,"/genotype_general_final_filtered.rds"))
  system(cmd_rm)
  rm(pheno,pheno_simulated,geno,cmd_unzip,cmd_rm)
}

rm(list=ls())


