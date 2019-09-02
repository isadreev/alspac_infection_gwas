#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=10
#PBS -l walltime=24:00:00

# Run GWAS analysis using GCTA
module load apps/gcta_1.91.1beta

source config

plink \
	--bfile $geno_dir/data/derived/filtered/bestguess/maf0.01_info0.8/combined/data \
	--clump-field p \
	--clump-p1 5e-8 \
	--clump-r2 0.01 \
	--clump-kb 1000 \
	--out output_dir/out \
	--clump out.mlma

