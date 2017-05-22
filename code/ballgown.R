options(stringsAsFactors = FALSE)

library(ballgown)
library(genefilter)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/new_heatmap_function.R")

pheno_data <- read.delim("data/stress_expt_meta.txt")
pheno_data <- arrange(pheno_data, sample)

express_data <- ballgown(dir("ballgown", full.names = TRUE), pData = pheno_data)

# Filter to remove low-abundance genes. One common issue with RNA-seq data is that genes often have very few or zero counts. A common step is to filter out some of these. Another approach that has been used for gene expression analysis is to apply a variance filter. Here we remove all transcripts with a variance across samples less than one:

#express_data <- subset(express_data, "rowVars(texpr(express_data)) > 1", genomesubset = TRUE)

X = gexpr(express_data)
row.max.fpkm = apply(X, 1, max)
X = X[row.max.fpkm > 5, ]
colnames(X) = gsub("FPKM.", "", colnames(X))

#################################
# PCA
#################################

X = t(X)


# PCA of the abiotic stress

X.sub = X[grep("_P\\d", rownames(X)), ]
fit = svd(scale(X.sub, scale = FALSE))

abiotic.fit <- data.frame(
					sample = rownames(X.sub), 
					PC1 = fit$u[,1], 
					PC2 = fit$u[,2], 
					PC3 = fit$u[,3]
				)
abiotic.fit <- left_join(abiotic.fit, pheno_data, by = "sample")

#c("ambient", "cold", "hot", "control", "chytrid", "rotifer")
abiotic.fit <- mutate(abiotic.fit, treatment = factor(treatment, levels = c("hot", "ambient", "cold")))

p = ggplot(abiotic.fit, aes(PC1, PC2, color = treatment)) +
    geom_point(aes(size = sampling.time.in.hrs, shape = type.of.stress)) 


# PCA of the biotic stress

X.sub = X[grep("_S\\d", rownames(X)), ]
fit = svd(scale(X.sub, scale = F))

biotic.fit <- data.frame(sample = rownames(X.sub), PC1 = fit$u[,1], PC2 = fit$u[,2], PC3 = fit$u[,3])
biotic.fit <- mutate(biotic.fit, sample = gsub("FPKM.", "", sample))
biotic.fit <- left_join(biotic.fit, pheno_data, by = "sample")

#c("ambient", "cold", "hot", "control", "chytrid", "rotifer")
biotic.fit <- mutate(biotic.fit, treatment = factor(treatment, levels = c("chytrid", "control", "rotifer")))

p = ggplot(biotic.fit, aes(PC1, PC2, color = treatment)) +
    geom_point(aes(size = sampling.time.in.hrs, shape = type.of.stress))


p = ggplot(biotic.fit, aes(PC1, PC3, color = treatment)) +
    geom_point(aes(size = sampling.time.in.hrs, shape = type.of.stress))


# data.for.3d <- select(biotic.fit, starts_with("PC"))
# data.for.3d <- as.matrix(data.for.3d)
# library(rgl)
# plot3d (
# 	data.for.3d,
# 	col = as.numeric(biotic.fit[["treatment"]]),
# 	size = ((biotic.fit[["sampling.time.in.hrs"]] / 72) * 3) + 1,
# 	type = "s",
# 	xlab = "PC1",
# 	ylab = "PC2",
# 	zlab = "PC3"
# )



# To modify the color scheme, add this line:
#    scale_color_manual(values = c("brown", "blue", "red", "grey", "lightgreen", "darkgreen"))


# Exploration:
# - find a few stress genes
# - plot the expresstion as a function of time
# - are the effects linear, plateau, early, late, etc??
# - think about optimal model


#################################
# Run Limma
#################################

x <- filter(pheno_data, type.of.stress == "abiotic", sampling.time.in.hrs <= 3) %>% arrange(treatment, sampling.time.in.hrs)
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)

treatment <- paste(x[["treatment"]], x[["sampling.time.in.hrs"]], sep = ".")
treatment <- factor(treatment)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- levels(treatment)
fit <- lmFit(X.sub, design)

cont.dif <- makeContrasts(
    Dif1hr = "(hot.1 - ambient.0) - (ambient.1 - ambient.0)",
    Dif3hr = "(hot.3 - hot.1) - (ambient.3 - ambient.1)",
	levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
hot.specific = topTableF(fit2, adjust = "BH", number = ncol(X.sub))


cont.dif <- makeContrasts(
    Dif1hr = "(cold.1 - ambient.0) - (ambient.1 - ambient.0)",
    Dif3hr = "(cold.3 - cold.1) - (ambient.3 - ambient.1)",
	levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
cold.specific = topTableF(fit2, adjust = "BH", number = ncol(X.sub))


hot.specific = mutate(hot.specific, gene = rownames(hot.specific)) %>%
				filter(adj.P.Val <= 0.05)

cold.specific = mutate(cold.specific, gene = rownames(cold.specific)) %>%
				filter(adj.P.Val <= 0.05)

de.genes <- unique(c(hot.specific[["gene"]], cold.specific[["gene"]]))


x <- filter(pheno_data, type.of.stress == "abiotic")
X.sub = X[x[["sample"]], de.genes]
X.sub = as.data.frame(X.sub)
X.sub = mutate(X.sub, sample = rownames(X.sub))
X.sub = gather(X.sub, gene, expression, starts_with("M"))
X.sub = left_join(X.sub, pheno_data, by = "sample")
X.sub = rename(X.sub, time = sampling.time.in.hrs)
X.sub = mutate(X.sub, treatment = factor(treatment, levels = c("hot", "ambient", "cold")))

pdf(file = "~/Desktop/temp_curves.pdf", height = 40, width = 4)
p = ggplot(X.sub, aes(time, expression, color = treatment)) + 
    geom_point() +
    geom_smooth(span = 0.3) +
    facet_grid(gene ~ ., scales = "free")
print(p)
dev.off()




X.sub = spread(X.sub, gene, expression)

# normalize to the 0hr

baseline = 

# For each expt:
#    - consider all time points as a linear effect
#    - consider custom time points
#    - consider first few timepoints
#    - limit to just a few timepoints
#    - try with temp as a numerical variable
#    - how to set the control? (hot vs control, and cold versus control)


# For the abiotic expt, compare:
#    - replicate, time, temp


# For the abiotic expt, compare:
#    - replicate, time, treatment


