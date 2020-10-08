
options(stringsAsFactors = FALSE)

library(ballgown)
library(genefilter)
library(dplyr)


tmap = read.delim("data/compare.stringtie_merge.gtf.tmap")
tmap = tmap %>%
       filter(ref_gene_id != "-") %>%
       group_by(qry_gene_id) %>% 
       summarise(ref_gene_id = paste(sort(unique(ref_gene_id)), collapse = ";")) %>%
       ungroup()

annot = read.delim("~/Projects/llaurens/se107/annotation/blast2go/nr_20170508/se107_nr_generic.txt")
annot = select(annot,
	Sequence.Name,
	Sequence.Description,
	Sequence.Length,
	Blast.Top.Hit.Accession,
	Mapping.Gene.Name,
	Mapping.GO.ID,
	Mapping.GO.Term,
	Mapping.GO.Category,
	Annotation.GO.ID,
	Annotation.GO.Term,
	Annotation.GO.Category,
	Enzyme.Name,
	InterPro.GO.ID,
	InterPro.GO.Term,
	InterPro.GO.Category
)

# How to merge the gene data with a matrix
# X = as.data.frame(X)
# X = mutate(X, sample = rownames(X))
# X = gather(X, gene, expression, -sample)
# X = left_join(X, tmap, by = c("gene" = "qry_gene_id"))


