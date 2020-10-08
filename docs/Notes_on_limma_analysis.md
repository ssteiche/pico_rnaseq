
# Notes on the statistical analysis of the Desmodesmus RNAseq data

The primary interest is the comparison with temperature stress versus control, and stress of pests versus the control. The researchers used two different temperature extremes and two different pests to cover their bases. They were also not sure when to collect their samples, so they opted to do a time course instead. 

The goals of the project are to 1) identify biomarkers that might alert the researchers when the algae are stressed and 2) understand the biology of the stress response.

1. The main deliverable is a list of genes and pathways that are differentially expressed when these algae are stressed. 
	- list of genes/pathways that are differentially expressed in all conditions
	- temp specific genes/pathways
	- pest specific genes/pathways
	- hot/cold specific genes/pathways
	- Rot/FD specific genes/pathways

2. It would be nice to separate genes based on when they are activated (early versus late). The primary responses are likely to be signal transduction/transcriptional pathways that sense the stress and turn on the necessary genes to deal with the problem. The secondary and tertiary responses are those genes that help the algae overcome the stress. The early genes are likely to be good early warning signs that the algae are stressed, allowing the researchers to respond before there is a loss in productivity. The late responding genes are also likely to be useful for diagnostics/prediction, but would offer more biological insight. Identification of the genes that are turned on in response to the stress helps understand how the algae plan to deal with the stress.

## How to do the stat test?

### Find DE genes using the first few timepoints, as these are likely to show a linear response

1. Separate data into abiotic and biotic responses
2. Select only the first few timepoints
3. Treat time as a categorical variable
4. Find DE genes significant for each treatment


This strategy gives 4842 genes for hot and 4802 for cold. 2556 genes overlap. At p-value <= 0.01

And gives **??** genes for Rot and **??** for FD. **??** genes overlap


### Look for a linear response using the full data

It's possible that the expression of some genes is linear with time. If so it would be dead simple to look for and identify these genes.

This strategy gives 413 genes for hot, 5566 for cold, and 6409 for combined response. At p-value <= 0.01

And gives **??** genes for Rot and **??** for FD. **??** genes overlap



### Fit a linear model to a range of different response curves

This should help us identify early, mid, and late responses as well as linear responses.

This strategy gives 1969 genes for hot, 121 for cold, and 1633 for combined response. At p-value <= 0.01



