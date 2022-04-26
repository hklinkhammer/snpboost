#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10GB

plink2="plink2"
folder="genotypes"

n=100000
p=100000
MAF=0.01
n_simu=100

##extract just white british people

${plink2} --bfile ukb_cal_allChrs_XY --keep ids_british.txt --make-bed --out ${folder}/genotype_general

genotypes=${folder}/genotype_general

## filter variants for MAF

${plink2} --bfile ${genotypes} --maf ${MAF} --make-bed --out ${genotypes}

## run R script for setting missing values to 0

Rscript prepare_genotype_simulation.R ${genotypes}

## convert into pgen-file format

${plink2} --bfile ${genotypes}_final_filtered --make-pgen vzs --out ${genotypes}_final_filtered

for (( i=1; i<=${n_simu}; i++ ))
do

${plink2} --pfile ${folder}/genotype_general_final_filtered vzs --extract ${folder}/n_${n}_p_${p}/bim_${i}.txt --keep ${folder}/n_${n}_p_${p}/fam_${i}.txt  --make-pgen vzs --out ${folder}/n_${n}_p_${p}/genotypes/genotype_${i}

done
