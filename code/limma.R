
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

source("code/ballgown.R")
source("code/load_annotation.R")

# TODO: come up with a better name (this is the original expression data as a matrix) we want to convert to a data frame so it can be merged with the results files
df = t(X)
df = as.data.frame(df)
df = mutate(df, gene = rownames(df))



################################################################
# Limit analysis to the first few hours of the time course
#    - because only a few points, treat time as a categorical variable
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
hot.specific = topTable(fit2, adjust = "BH", number = nrow(X.sub))


cont.dif <- makeContrasts(
    Dif1hr = "(cold.1 - ambient.0) - (ambient.1 - ambient.0)",
    Dif3hr = "(cold.3 - cold.1) - (ambient.3 - ambient.1)",
    Dif9hr = "(cold.9 - cold.3) - (ambient.9 - ambient.3)",
    levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
cold.specific = topTableF(fit2, adjust = "BH", number = nrow(X.sub))


hot.specific  = mutate(hot.specific, gene = rownames(hot.specific))
cold.specific = mutate(cold.specific, gene = rownames(cold.specific))

hot.specific  = left_join(hot.specific, df, by = "gene")
cold.specific = left_join(cold.specific, df, by = "gene")

output.dir = "results/initial_timepoints"
if (! dir.exists(output.dir))
{
    dir.create(output.dir, recursive = TRUE)
}
write.table(hot.specific, file = file.path(output.dir, "hot_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cold.specific, file = file.path(output.dir, "cold_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



###################################################################
# Create a simple linear model looking for differences between treatments
###################################################################
x <- dplyr::filter(pheno_data, type.of.stress == "abiotic")
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)
X.sub = log2(X.sub + 1)

treatment <- factor(x[["treatment"]])
time      <- x[["sampling.time.in.hrs"]]

library(splines)
time = ns(time, df = 1)

# TODO: why not model.matrix(~ time + treatment) ??

design <- model.matrix(~ treatment * time)
colnames(design) <- c("intercept", "cold", "hot", "time", "cold:time", "hot:time")
fit <- lmFit(X.sub, design)
fit <- eBayes(fit)

hot.specific  = topTable(fit, coef = "hot", adjust = "BH", number = nrow(X.sub))
cold.specific = topTable(fit, coef = "cold", adjust = "BH", number = nrow(X.sub))
top.stress    = topTable(fit, coef = 5:6, adjust = "BH", number = nrow(X.sub))

hot.specific  = mutate(hot.specific, gene = rownames(hot.specific))
cold.specific = mutate(cold.specific, gene = rownames(cold.specific))
top.stress    = mutate(top.stress, gene = rownames(top.stress))


hot.specific  = left_join(hot.specific, df, by = "gene")
cold.specific = left_join(cold.specific, df, by = "gene")
top.stress    = left_join(top.stress, df, by = "gene")

output.dir = "results/simple_model"
if (! dir.exists(output.dir))
{
    dir.create(output.dir, recursive = TRUE)
}
write.table(hot.specific, file = file.path(output.dir, "hot_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cold.specific, file = file.path(output.dir, "cold_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top.stress, file = file.path(output.dir, "top_stress.txt"), sep = "\t", quote = FALSE, row.names = FALSE)






#####################################################
# Follow the instructions using a natural spline
#####################################################
x <- dplyr::filter(pheno_data, type.of.stress == "abiotic")
X.sub = X[x[["sample"]], ]
X.sub = t(X.sub)
X.sub = log2(X.sub + 1)

treatment <- factor(x[["treatment"]])
time = x[["sampling.time.in.hrs"]]
time = cbind(
    ns(time, df = 1),
    dnorm(time, mean =  9,  sd = 5),
    dnorm(time, mean = 24, sd = 12),
    pnorm(time, mean =  5,  sd = 3)
)
colnames(time) = 1:ncol(time)
y.max = apply(time, 2, max)
time = t(apply(time, 1, function(x) { (x / y.max) * 0.8 }))

design <- model.matrix(~ treatment * time)
fit <- lmFit(X.sub, design)
fit <- eBayes(fit)

hot.specific  = topTable(fit, coef = c(8,10,12,14), adjust = "BH", number = nrow(X.sub))
cold.specific = topTable(fit, coef = c(9,11,13,15), adjust = "BH", number = nrow(X.sub))
top.stress    = topTable(fit, coef = 8:15, adjust = "BH", number = nrow(X.sub))

hot.specific  = mutate(hot.specific, gene = rownames(hot.specific))
cold.specific = mutate(cold.specific, gene = rownames(cold.specific))
top.stress    = mutate(top.stress, gene = rownames(top.stress))


hot.specific  = left_join(hot.specific, df, by = "gene")
cold.specific = left_join(cold.specific, df, by = "gene")
top.stress    = left_join(top.stress, df, by = "gene")

output.dir = "results/full_model"
if (! dir.exists(output.dir))
{
    dir.create(output.dir, recursive = TRUE)
}
write.table(hot.specific, file = file.path(output.dir, "hot_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cold.specific, file = file.path(output.dir, "cold_specific.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top.stress, file = file.path(output.dir, "top_stress.txt"), sep = "\t", quote = FALSE, row.names = FALSE)




#----------------------------------------------

# functions here for downstream plotting




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

plot_gene <- function(df, gene.name, plot.title = gene.name)
{
    df = filter(df, gene == gene.name)
    p = ggplot(df, aes(time, expression, color = treatment)) + 
        geom_point() +
        geom_smooth(span = 0.3) +
        labs(title = plot.title)
    return(p)
}


