args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]

snps = read.delim(input,
                  stringsAsFactors = FALSE,
                  header = TRUE,
                  sep = '\t')

inbreeding_coefficient = function(df) {
  obs_ref = apply(df, MARGIN = 1, function(x) sum(x == '0/0'))
  obs_het = apply(df, MARGIN = 1, function(x) sum(x == '0/1'))
  obs_alt = apply(df, MARGIN = 1, function(x) sum(x == '1/1'))
  p = (2 * obs_ref + obs_het) / (2 * obs_ref + 2 * obs_het + 2 * obs_alt)
  q = (2 * obs_alt + obs_het) / (2 * obs_ref + 2 * obs_het + 2 * obs_alt)
  exp_het = 2 * p * q * 100
  ibc = 1 - (obs_het / exp_het)
  return(ibc)
}

snps$IBC = inbreeding_coefficient(snps)
snps = snps[snps$IBC >= 0, ]
positions = snps[, c(1, 2)]
write.table(positions, file = output,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE)
