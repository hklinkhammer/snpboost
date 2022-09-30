library(data.table)
library(tidyverse)

### script to prepare sumstats in PRScs format
sum_stats <- fread("plink2.height.glm.linear")

sum_stats <- sum_stats %>% select(ID,A1,AX,BETA,P) %>% rename(SNP=ID,A2=AX)

fwrite(sum_stats,"height_sumstats.txt",sep="\t")
