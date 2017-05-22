####################################################################
## Quality Control
####################################################################
PCA.plot <- function(X, sample.names, sample.index, components = 1:2, ...)
{
	x.svd <- svd(scale(X, scale = FALSE))

	# Sometimes long sample labels go beyond the borders of the plot. 
	# let's expand the borders of the plot by 10%
	x.min <- min(x.svd$u[, components[1]])
	x.max <- max(x.svd$u[, components[1]])
	x <- x.max - x.min
	x.min <- x.min - (x * 0.1)
	x.max <- x.max + (x * 0.1)
	
	pc1 <- percent_var_explained(x.svd, components[1])
	pc2 <- percent_var_explained(x.svd, components[2])
	pc1 <- round(pc1, 1)
	pc2 <- round(pc2, 1)
	
	if (is.null(sample.names))
	{
		type = "p"
	} else {
		type = "n"
	}
	
	plot(
		x.svd$u[, components], 
		type = type, 
		pch = 20,
		col = sample.index,
		xlab = paste0("PC", components[1], " (", pc1, "%)"), 
		ylab = paste0("PC", components[2], " (", pc2, "%)"), 
		cex.lab = 1.2, 
		xlim = c(x.min, x.max)
	)
	abline(h = 0, col = grey(0.85))
	abline(v = 0, col = grey(0.85))
	if (! is.null(sample.names))
	{
		text(x.svd$u[, 1:2], labels = sample.names, col = sample.index, ...)
	}
}

Norm <- function(x) { return(sqrt(sum(x^2))) }

percent_var_explained <- function(fit, component)
{
	d <- fit$d / sqrt(max(1, ncol(fit$u) - 1))
	percent.explained <- d^2 / sum(d^2)
	percent.explained <- percent.explained * 100
	return(percent.explained[component])
}

Dendrogram.plot <- function(X.temp, sample.names, sample.index, hang = 0.05, ...)
{
	d   <- dist(X.temp, method = "euclidean")
	fit <- hclust(d, method = "ward.D2")
	myplclust(fit, lab = sample.names, lab.col = sample.index, ylab = "Euclidean Distance", xlab = "", sub = "", cex = 1, hang = hang, ...)
}

Plot.PCA.3D <- function(X.temp, sample.index)
{
	library(rgl)
	x.svd <- svd(scale(X.temp, scale = FALSE))
	plot3d (
		x.svd$u[, 1:3],
		col = sample.index,
		size = 2,
		type = "s",
		xlab = "PC1",
		ylab = "PC2",
		zlab = "PC3"
		)
}



####################################################################
# Limma
####################################################################
Run.Limma <- function(X.temp, design, coefficient)
{
	X.temp[X.temp == 0] <- 1e-16
	X.temp <- log2(X.temp)	
	fit <- lmFit(t(X.temp), design)
	fit <- eBayes(fit)
	top.genes <- topTable(fit, coef = coefficient, lfc = log2(1.5), adjust.method = "BH", number = ncol(X.temp))
	return(top.genes)
}



####################################################################
# GSEA
####################################################################
Run.gsea <- function(X.temp, sample.index, control.title, treated.title, gsea.folder = "gsea")
{	
	#gene.set    <- "go_all"
	#gmt.file    <- "~/gsea/c5.all.v4.0.symbols.gmt"
	gene.set    <- "kegg"
	gmt.file    <- "~/gsea/c2.cp.kegg.v4.0.symbols.gmt"
	
	if (! paste0("./", gsea.folder) %in% list.dirs()) dir.create(gsea.folder)

	writeGCT(X.temp, file = paste0(gsea.folder, "/data.gct"))
	writeCLS(sample.index, control.title, treated.title, file = paste0(gsea.folder, "/data.cls"))
	runGSEA(gct.file = paste0(gsea.folder, "/data.gct"), cls.file = paste0(gsea.folder, "/data.cls"), gene.set = gmt.file, output.folder = paste(getwd(), gsea.folder, sep = "/"), folder.prefix = gene.set, verbose = TRUE) 

	sample.pathways <- getCoreGenesAsList(control.title, folder.prefix = gene.set, p.value = 0.05, gsea.folder = gsea.folder)
	sample.pathways <- c(sample.pathways, getCoreGenesAsList(treated.title, folder.prefix = gene.set, p.value = 0.05, gsea.folder = gsea.folder))
	core.genes <- unlist(lapply(sample.pathways, function(x) x[,1]))

	pathways.per.core.gene <- sapply(unique(core.genes), function(x) paste(names(core.genes)[core.genes == x], collapse = " // ")) 
	pathways.per.core.gene <- data.frame(gene = names(pathways.per.core.gene), gsea.pathways = pathways.per.core.gene)
	return(pathways.per.core.gene)
}


####################################################################
# Compile Results
####################################################################
Compile.Results <- function()
{
  results <- data
  if (exists("top.genes"))
  {
    results <- merge(results, top.genes, by.x = "symbol", by.y = "row.names", all.x = TRUE)
  }
  if (exists("pathways.per.core.gene"))
  {
    results <- merge(results, pathways.per.core.gene, by.x = "symbol", by.y = "gene", all.x = TRUE)
  }
  write.table(results, "results/FullData.txt", sep = "\t", quote = FALSE, row.names = FALSE)
}
