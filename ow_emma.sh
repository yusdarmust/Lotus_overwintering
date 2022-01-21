#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 20g
#SBATCH --account InRoot
#SBATCH -c 4
#SBATCH -t 12:00:00

python2.7  ~/atgwas/src/gwa.py \
	-a emmax \
	-m  8 \
	-r  data/phenotype_.csv \
	-f  data/20220117_lotus_snps.csv \
	--data_format='diploid_int'
