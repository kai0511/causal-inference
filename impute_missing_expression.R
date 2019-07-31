library(softImpute)

setwd('~/data/')
# expr_matrix <- read.table('zscore_matrix.tsv', sep = '\t', stringsAsFactors = F, header = T)
load('zscore_matrix.RData')
imputed_matrix <- softImpute(m, rank.max = 200)
