#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 16
#SBATCH --account InRoot
#SBATCH -c 1
#SBATCH -t 12:00:00

vcftools --gzvcf data/20220117_lotus_snps.vcf.gz \
--012 \
--out data/20220117_lotus_snps
