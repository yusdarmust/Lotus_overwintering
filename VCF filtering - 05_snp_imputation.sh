#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 16g
#SBATCH --account InRoot
#SBATCH -c 8
#SBATCH -t 12:00:00

java -jar share/beagle-5.1_24Aug19.3e8-1/beagle.jar \
gt=data/maf_filtered_snps.recode.vcf \
out=data/imputed_snps
