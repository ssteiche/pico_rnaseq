
options(stringsAsFactors = FALSE)

library(dplyr)

df = read.delim("~/Projects/llaurens/se107/annotation/SE107_augustus_gene_predictions.gtf", comment.char = "#", header = FALSE)

colnames(df) = c("chr", "source", "class", "pos1", "pos2", "score", "strand", "foo", "gene.info")


introns = df %>% 
  filter(class == "intron") %>%
  mutate(length = pos2 - pos1)

summary(introns[["length"]])


