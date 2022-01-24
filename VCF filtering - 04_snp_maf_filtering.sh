#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 100g
#SBATCH -c 2
#SBATCH --account InRoot
#SBATCH -t 12:00:00

vcftools --vcf data/third_filtering_snps.recode.vcf \
--maf 0.05 \
--recode \
--recode-INFO-all \
--out data/maf_filtered_snps
