library(data.table)
library(tidyverse)

### script to prepare sumstats in LDAK format
sum_stats <- fread("plink2.height.glm.linear")

sum_stats <- sum_stats %>% 
  rename(Predictor=ID,A2=AX,Direction = BETA) %>% 
  filter(nchar(A1)==1 & nchar(A2)==1)

fwrite(sum_stats,"height.glm.linear",sep="\t")
