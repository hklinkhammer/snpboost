library(bigsnpr)
commandArgs_vec = commandArgs(T)

out=commandArgs_vec[1]

rds <- snp_readBed("genotype_general.bed")

genotypes <- snp_attach(rds)

## missing values are set to 0
genotypes$genotypes$code256 <- c(0,1,2,rep(0,253))

#genotypes <- snp_save(genotypes)

snp_writeBed(genotypes,"genotype_general_final_filtered.bed")
