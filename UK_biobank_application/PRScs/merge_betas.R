library(data.table)
library(tidyverse)

setwd("Output")

args <- commandArgs(T)

phenotype <- args[1]

path <- paste0(phenotype,"_pst_eff_a1_b0.5_phiauto_chr")

betas <- fread(paste0(path,1,".txt"))

for(i in 2:22){
  betas_tmp <- fread(paste0(path,i,".txt"))
  betas <- betas %>% add_row(betas_tmp)
}

rm(betas_tmp)

fwrite(betas,paste0(phenotype,"_betas.txt"),sep="\t",col.names = F)

cmd=paste0("rm ",path,"*")

system(cmd)
