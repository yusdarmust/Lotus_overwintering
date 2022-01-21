# Analyses for the SNPs after filtering performed on 17.01.2022.

library(qqman)
library(magrittr)
library(ggplot2)

source('scripts/MP.R')

# Manhattan plot & important SNPs ---------------
ow14 = read.delim('data/gwas_ow14_3.csv',
                  sep = ',',
                  header = TRUE)

dfForMP = function(df) {
  df = df[!df$chromosomes %in% c(0, 7, 8), ]
  df = df[df$macs > 8, ]
  df = df[, c(1, 2, 3)]
  df$SNP = paste(df$chromosomes, df$positions, sep = '_')
  colnames(df)[c(1, 2, 3)] = c('CHR', 'BP', 'P')
  df = df[, c(4, 1, 2, 3)]
  return(df)
}

ow14 = dfForMP(ow14)
qq(ow14$P)

# Important genes
ow14$Padj = p.adjust(p = ow14$P, method = 'BH')
gene1 = ow14[as.character(170966:170976), ]
gene6 = ow14[as.character(878153:878168), ]
soi = c(gene1$SNP, gene6$SNP)
write.csv(gene1, 'results/20220117/gene1.csv')
write.csv(gene6, 'results/20220117/gene6.csv')

png('results/20220117/ow14_mp_.png', height = 4500, width = 8000, res = 1000)
MP(ow14[, 1:4],
   suggestiveline = -log10(0.05 / (nrow(ow14))),
   highlight = soi,
   col = c('#332288', '#0066CC'),
   genomewideline = FALSE)
dev.off()

# Significant SNPs ------------------------------
sig = ow14[ow14$P <= 1e-3, ]

# Make bed file for BH
bed = sig[, c(2, 3, 4)]
bed$other_pos = bed$BP - 1
bed = bed[, c(1, 4, 2, 3)]
for(i in 1:length(bed$CHR)) {
  if(bed$CHR[i] == 1) {bed$CHR[i] = 'LjG1.1_chr1'}
  if(bed$CHR[i] == 2) {bed$CHR[i] = 'LjG1.1_chr2'}
  if(bed$CHR[i] == 3) {bed$CHR[i] = 'LjG1.1_chr3'}
  if(bed$CHR[i] == 4) {bed$CHR[i] = 'LjG1.1_chr4'}
  if(bed$CHR[i] == 5) {bed$CHR[i] = 'LjG1.1_chr5'}
  if(bed$CHR[i] == 6) {bed$CHR[i] = 'LjG1.1_chr6'}
}
write.table(bed,
            sep = '\t',
            file = 'results/20220117/sig_snps_unadj.bed',
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
