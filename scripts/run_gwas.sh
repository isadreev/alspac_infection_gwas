#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=10
#PBS -l walltime=24:00:00

# Run GWAS analysis using GCTA
module load apps/gcta_1.91.1beta

source config

# Continuous phenotypes
gcta64 --mlma \
	--bfile $geno_dir/data/derived/filtered/bestguess/maf0.01_info0.8/combined/data \
	--grm $geno_dir/data/derived/kinships/data \
	--pheno $files_dir/Pheno.phen \
	--covar $files_dir/Covariates_discrete.covar \
	--qcovar $files_dir/Covariates_quant.qcovar \
	--out /output_dir/out \
	--thread-num 10
	

	