
# heatmap code


options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(ggplot2)


experiment = "abiotic"
p.value.threshold = 0.01
foldchange.threshold = 5
hc = TRUE
hr = FALSE

pheno_data <- read.delim("data/stress_expt_meta.txt")
pheno_data <- arrange(pheno_data, sample)


df = read.delim("results/full_model/top_stress.txt")
df = df %>% 
		select(-starts_with("treatment")) %>%
		dplyr::filter(adj.P.Val <= p.value.threshold) %>%
		gather(sample, expression, starts_with("X")) %>%
		mutate(sample = gsub("^X", "", sample))

df = left_join(df, pheno_data, by = "sample")
df = dplyr::filter(df, type.of.stress == experiment)
df = rename(df, time = sampling.time.in.hrs)

# normalize data to the 0hr time point
df = df %>%
		group_by(gene, time) %>%
		mutate(expression = ifelse(expression == 0, 1e-16, expression)) %>%
		mutate(expression = log2(expression / mean(expression[treatment %in% c("ambient", "control")]))) %>%
		ungroup()
df = dplyr::filter(df, 
		time > 0,
		! treatment %in% c("ambient", "control")
)

#################################

kegg = read.delim("data/kegg_terms.txt")
kegg.annotations = read.delim("data/kegg_kaas_sbh_augustus.ko.txt", header = FALSE, na.strings = "")
colnames(kegg.annotations) = c("gene", "ortholog")

tmap = read.delim("data/compare.stringtie_merge.gtf.tmap")
tmap = tmap %>%
	filter(ref_gene_id != "-") %>%
	group_by(qry_gene_id) %>%
	summarise(ref_gene_id = paste(sort(unique(ref_gene_id)), collapse = ";")) %>%
	ungroup()

kegg = left_join(kegg, kegg.annotations, by = "ortholog")
kegg = left_join(kegg, tmap, by = c("gene" = "ref_gene_id"))
kegg = mutate(kegg, qry_gene_id = ifelse(is.na(qry_gene_id), gene, qry_gene_id))
kegg = select(kegg, -gene)
kegg = rename(kegg, gene = qry_gene_id)

glycolysis.genes = dplyr::filter(kegg, grepl("Glycolysis", description))
starch.genes = dplyr::filter(kegg, grepl("Starch and sucrose metabolism", description))
carbon.fix.genes = dplyr::filter(kegg, grepl("Carbon fixation in photosynthetic organisms", description))
tca.cycle.genes = dplyr::filter(kegg, grepl("TCA cycle", description))

annot = read.delim("data/se107_nr_generic.txt")
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
	       Enzyme.Name
)

extracellular = dplyr::filter(annot, 
	grepl("extracell", Annotation.GO.Term) | 
	grepl("extracell", Mapping.GO.Term) | 
	grepl("extracell", InterPro.GO.Term),
	! is.na(stringtie.gene.id)
)

#df = dplyr::filter(df, gene %in% tca.cycle.genes[["gene"]])
df = dplyr::filter(df, gene %in% extracellular[["stringtie.gene.id"]])

#################################

X <- df %>%
		dplyr::filter(treatment == "cold") %>%
		select(sample, treatment, time, gene, expression) %>%
		spread(gene, expression) %>%
		arrange(treatment, time) %>%
		select(-treatment, -time) %>%
		as.data.frame()
rownames(X) = X[["sample"]]
X = X[, -1]
X <- as.matrix(X)
# X[X > foldchange.threshold] = foldchange.threshold
# X[X < -foldchange.threshold] = -foldchange.threshold
if (hc)
{
	dd.col <- as.dendrogram(hclust(dist(t(X))))
	col.ord <- order.dendrogram(dd.col)
} else
{
	col.ord <- 1:ncol(X)
}

if (hr)
{
	dd.row <- as.dendrogram(hclust(dist(X)))
	row.ord <- order.dendrogram(dd.row)
} else 
{
	row.ord <- 1:nrow(X)
}

gene.order = colnames(X)[col.ord]

df = mutate(df, gene = factor(gene, levels = gene.order))
df = mutate(df, expression = ifelse(expression > foldchange.threshold, foldchange.threshold, expression))
df = mutate(df, expression = ifelse(expression < -foldchange.threshold, -foldchange.threshold, expression))

theme_no_gene_labels <- theme(axis.text.x = element_blank())
# scale_y_continuous(breaks=NULL)


theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour = NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

### Create plot components ###    
# Heatmap
p1 <- ggplot(df, aes(x = factor(time), y = gene)) + 
  geom_tile(aes(fill = expression)) + 
  scale_fill_gradient2(low = rgb(61, 16, 123, maxColorValue = 255), mid = "white", high = rgb(223,0,35, maxColorValue = 255)) + 
  ylab("") +
  xlab("Sampling Time (in hrs)") +
  facet_grid(. ~ treatment) +
  theme(
	  axis.title.y = element_blank(),
#	  axis.text.y  = element_blank(),
	  axis.text.y  = element_text(size = rel(0.75)),
	  axis.line    = element_blank()
  )

# if (ncol(X) > 50)
# {
# 	p1 <- p1 + theme_no_gene_labels
# } else {
# 	#make row names vertical
# }

print(p1)
