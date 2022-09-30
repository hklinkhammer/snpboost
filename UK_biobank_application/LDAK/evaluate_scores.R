library(data.table)
library(tidyverse)

setwd("/daten/epi/PRS/snpnet/LDAK")

args <- commandArgs(T)

phenotype <- args[1]

method <- args[2]

scores <- fread(paste0(phenotype,"/",method,"/",phenotype,".sscore")) %>% rename(IID='#IID')

phenotypes <- fread("../non_imputed_data/phenotypes/phenotypes_unrelated.phe")

df <- left_join(phenotypes,scores)

if(phenotype=="LDL"){phenotype="LDL_adj"}
if(phenotype=="lipo"){phenotype="lipoprotein_A"}

#print(paste("correlation", phenotype,"and SCORE1_SUM:", round(cor(df[[phenotype]][df$train_test=="test"],df$SCORE1_SUM[df$train_test=="test"])^2,4)))

### glm

df <- df %>% rename(PRS=SCORE1_AVG)

lm <- glm(paste0(phenotype, " ~ 1 + sex + age_at_2015 + ",paste0("PC",1:10,collapse="+")," + PRS"),data=df %>% filter(train_test != "test"))

summary(lm)

pred <- predict(lm,df %>% filter(train_test=="test"))

print(paste("correlation", phenotype,"and PRS:", round(cor(df[[phenotype]][df$train_test=="test"],df$PRS[df$train_test=="test"])^2,4)))
print(paste("correlation", phenotype, "and prediction with PRS and covariates:", round(cor(df[[phenotype]][df$train_test=="test"],pred)^2,4)))

