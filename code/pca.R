options(stringsAsFactors = FALSE)

library(ballgown)
library(genefilter)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/ballgown.R")

###############################
# PCA of the abiotic stress
###############################

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


p.abiotic = ggplot(abiotic.fit, aes(PC1, PC2, color = treatment)) +
    geom_point(aes(size = sampling.time.in.hrs, shape = type.of.stress)) 

# To modify the color scheme, add this line:
#    scale_color_manual(values = c("brown", "blue", "red", "grey", "lightgreen", "darkgreen"))

###############################
# PCA of the biotic stress
###############################

X.sub = X[grep("_S\\d", rownames(X)), ]
fit = svd(scale(X.sub, scale = F))

biotic.fit <- data.frame(sample = rownames(X.sub), PC1 = fit$u[,1], PC2 = fit$u[,2], PC3 = fit$u[,3])
biotic.fit <- mutate(biotic.fit, sample = gsub("FPKM.", "", sample))
biotic.fit <- left_join(biotic.fit, pheno_data, by = "sample")

#c("ambient", "cold", "hot", "control", "chytrid", "rotifer")
biotic.fit <- mutate(biotic.fit, treatment = factor(treatment, levels = c("chytrid", "control", "rotifer")))

p.biotic.PC12 = ggplot(biotic.fit, aes(PC1, PC2, color = treatment)) +
    geom_point(aes(size = sampling.time.in.hrs, shape = type.of.stress))


p.biotic.PC13 = ggplot(biotic.fit, aes(PC1, PC3, color = treatment)) +
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


########################################
# Save plots to disk
########################################

pdf(file = "PCA_abiotic.pdf")
print(p.abiotic)
dev.off()

pdf(file = "PCA_biotic_PC12.pdf")
print(p.abiotic.PC12)
dev.off()

pdf(file = "PCA_biotic_PC13.pdf")
print(p.abiotic.PC13)
dev.off()




