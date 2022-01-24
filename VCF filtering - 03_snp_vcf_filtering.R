suppressMessages(library(dplyr))

args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]

snps = read.delim(input,
                  stringsAsFactors = FALSE,
                  header = TRUE,
                  sep = '\t')
# Remove columns with all NA
not_any_na = function(x) all(!is.na(x))
snps = select_if(snps, not_any_na)
for(i in 3:8) {
  snps[, i] = as.numeric(snps[, i])
}

# Gifu is either 0/0 or 1/1
snps = snps[snps$X.9.Gifu_.GT %in% c('0/0', '1/1'), ]
# MG020 is either 0/0 or 1/1
snps = snps[snps$X.27.MG020.GT %in% c('0/0', '1/1'), ]
# MQ >= 30
snps = snps[snps$X.5.MQ >= 30, ]
# DP >= 60
snps = snps[snps$X.4.DP >= 60, ]
# QUAL >= 60
snps = snps[snps$X.3.QUAL >= 60, ]
# Missing < 50%
missing = apply(snps[, 9:167], 1, function(x) sum(x == './.'))
snps = snps[missing <= 79, ]
positions = snps[, 1:2, drop = FALSE]
# Output
write.table(positions, file = output,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE)
