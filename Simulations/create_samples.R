library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)

set.seed(2021)

n <- 100000
p <- 100000

nsimu <- 100

setwd("genotypes")

fam <- fread("genotype_general_final_filtered.fam")

bim <- fread("genotype_general_final_filtered.bim")

for(i in 1:nsimu){
  n_sample <- sort(sample(1:nrow(fam),n))
  p_sample <- sort(sample(1:nrow(bim),p))
  
  fam_sample <- fam[n_sample,1:2] %>% mutate(index=n_sample)
  bim_sample <- bim[p_sample,] %>% mutate(index=p_sample)
  
  fwrite(fam_sample,file=paste0("n_",n,"_p_",p,"/fam_",i,".fam"))
  fwrite(bim_sample,file=paste0("n_",n,"_p_",p,"/bim_",i,".bim"))
}
