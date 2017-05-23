
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

source("code/ballgown.R")
source("code/load_annotation.R")

###################################################################
# Create a simple model looking for differences between treatments
###################################################################
x <- dplyr::filter(pheno_data, type.of.stress == "abiotic")
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)
X.sub = log2(X.sub + 1)

treatment <- factor(x[["treatment"]])
time      <- x[["sampling.time.in.hrs"]]

library(splines)
time = ns(time, df = 1)

design <- model.matrix(~ treatment * time)
colnames(design) <- c("intercept", "cold", "hot", "time", "cold:time", "hot:time")
fit <- lmFit(X.sub, design)
fit <- eBayes(fit)

hot.specific  = topTable(fit, coef = "hot", adjust = "BH", number = nrow(X.sub))
cold.specific = topTable(fit, coef = "cold", adjust = "BH", number = nrow(X.sub))
top.stress    = topTable(fit, coef = 5:6, adjust = "BH", number = nrow(X.sub))


#####################################################
# Follow the instructions using a natural spline
#####################################################
x <- dplyr::filter(pheno_data, type.of.stress == "abiotic")
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)
X.sub = log2(X.sub + 1)

treatment <- factor(x[["treatment"]])
time      <- x[["sampling.time.in.hrs"]]

library(splines)
time = ns(time, df = 4)

design <- model.matrix(~ treatment * time)
fit <- lmFit(X.sub, design)
fit <- eBayes(fit)

hot.specific  = topTable(fit, coef = c(8,10,12,14), adjust = "BH", number = nrow(X.sub))
cold.specific = topTable(fit, coef = c(9,11,13,15), adjust = "BH", number = nrow(X.sub))
top.stress    = topTable(fit, coef = 8:15, adjust = "BH", number = nrow(X.sub))






################################################################
# Limit analysis to the first three hours of the time course
#    - These responses are most likely to be linear
#    - likely identify the primary responses to the treatment
################################################################

x <- dplyr::filter(pheno_data, type.of.stress == "abiotic", sampling.time.in.hrs <= 9) %>% arrange(treatment, sampling.time.in.hrs)
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)
X.sub = log2(X.sub + 1)

treatment <- paste(x[["treatment"]], x[["sampling.time.in.hrs"]], sep = ".")
treatment <- factor(treatment)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- levels(treatment)
fit <- lmFit(X.sub, design)

cont.dif <- makeContrasts(
    Dif1hr = "(hot.1 - ambient.0) - (ambient.1 - ambient.0)",
    Dif3hr = "(hot.3 - hot.1) - (ambient.3 - ambient.1)",
    Dif9hr = "(hot.9 - hot.3) - (ambient.9 - ambient.3)",
	levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
hot.specific = topTable(fit2, adjust = "BH", number = ncol(X.sub))


cont.dif <- makeContrasts(
    Dif1hr = "(cold.1 - ambient.0) - (ambient.1 - ambient.0)",
    Dif3hr = "(cold.3 - cold.1) - (ambient.3 - ambient.1)",
    Dif9hr = "(cold.9 - cold.3) - (ambient.9 - ambient.3)",
	levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
cold.specific = topTableF(fit2, adjust = "BH", number = ncol(X.sub))


hot.specific = mutate(hot.specific, gene = rownames(hot.specific)) %>%
				dplyr::filter(adj.P.Val <= 0.05)

cold.specific = mutate(cold.specific, gene = rownames(cold.specific)) %>%
				dplyr::filter(adj.P.Val <= 0.05)

de.genes <- unique(c(hot.specific[["gene"]], cold.specific[["gene"]]))


###########################################################################

x <- dplyr::filter(pheno_data, type.of.stress == "abiotic")
#foo = X[x[["sample"]], de.genes]
foo = X[x[["sample"]], ]
foo = as.data.frame(foo)
foo = mutate(foo, sample = rownames(foo))
foo = gather(foo, gene, expression, matches("^[Mg]"))
foo = left_join(foo, pheno_data, by = "sample")
foo = rename(foo, time = sampling.time.in.hrs)
foo = mutate(foo, treatment = factor(treatment, levels = c("hot", "ambient", "cold")))

plot_gene <- function(df, gene.name)
{
    df = filter(df, gene == gene.name)
    p = ggplot(df, aes(time, expression, color = treatment)) + 
        geom_point() +
        geom_smooth(span = 0.3) +
        labs(title = gene.name)
    return(p)
}


pdf(file = "~/Desktop/temp_curves.pdf", height = 40, width = 4)
p = ggplot(bar, aes(time, expression, color = treatment)) + 
    geom_point() +
    geom_smooth(span = 0.3)
#    facet_grid(gene ~ ., scales = "free")
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


