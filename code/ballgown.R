
options(stringsAsFactors = FALSE)

library(ballgown)
library(genefilter)
library(dplyr)

pheno_data <- read.delim("data/stress_expt_meta.txt")
pheno_data <- arrange(pheno_data, sample)

express_data <- ballgown(dir("ballgown", full.names = TRUE), pData = pheno_data)

# Filter to remove low-abundance genes. One common issue with RNA-seq data is that genes often have very few or zero counts. A common step is to filter out some of these. Another approach that has been used for gene expression analysis is to apply a variance filter. Here we remove all transcripts with a variance across samples less than one:

#express_data <- subset(express_data, "rowVars(texpr(express_data)) > 1", genomesubset = TRUE)

X = gexpr(express_data)
row.max.fpkm = apply(X, 1, max)
X = X[row.max.fpkm > 5, ]
colnames(X) = gsub("FPKM.", "", colnames(X))
X = t(X)

