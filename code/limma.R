
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)

source("code/ballgown.R")

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


