plot_heatmap <- function(X, threshold = 1, hc = TRUE, hr = TRUE)
{
	X <- as.matrix(X)
	X[X > threshold] <- threshold
	X[X < -threshold] <- -threshold

	if (hc)
	{
		dd.col <- as.dendrogram(hclust(dist(X)))
		col.ord <- order.dendrogram(dd.col)
	} else
	{
		col.ord <- 1:ncol(X)
	}
	
	if (hr)
	{
		dd.row <- as.dendrogram(hclust(dist(t(X))))
		row.ord <- order.dendrogram(dd.row)
	} else 
	{
		row.ord <- 1:nrow(X)
	}

	xx <- X
	#xx <- scale(xx, scale = FALSE)
	xx <- xx[col.ord, row.ord]
	xx_names <- attr(xx, "dimnames")
	df <- as.data.frame(xx)
	colnames(df) <- xx_names[[2]]
	df$sample <- xx_names[[1]]
	df$sample <- with(df, factor(sample, levels = sample, ordered = TRUE))

	mdf <- melt(df, id.vars = "sample")

	theme_no_gene_labels <- theme(axis.text.x = element_blank())
	# scale_y_continuous(breaks=NULL)


	theme_none <- theme(
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  panel.background = element_blank(),
	  axis.title.x = element_text(colour=NA),
	  axis.title.y = element_blank(),
	  axis.text.x = element_blank(),
	  axis.text.y = element_blank(),
	  axis.line = element_blank()
	  #axis.ticks.length = element_blank()
	)

	### Create plot components ###    
	# Heatmap
	p1 <- ggplot(mdf, aes(x = variable, y = sample)) + 
	  geom_tile(aes(fill = value)) + scale_fill_gradient2() + ylab("")

	if (ncol(X) > 50)
	{
		p1 <- p1 + theme_no_gene_labels
	} else {
		#make row names vertical
	}

	print(p1)
	
}

