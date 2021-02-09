# Replicate statistical analyses in limma.R for comparison

library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(stringr)
library(readr)
library(tximport)
library(DESeq2)

################################################################
# Consider turning the following into a stand alone module/script
################################################################
# Create named vector of salmon quant file locations
dir <- "salmon_quants"
samples <- str_extract(list.files(dir), ".*(?=_quant)")
files <- file.path(dir, paste0(samples,"_quant"), "quant.sf")
names(files) <- samples

# Make table of transcripts and corresponding genes
gff <- read_table2("mikado.loci.gff3", col_names = FALSE, comment = "#")

gffl <- gff %>% subset(gff$X3 == "mRNA") %>% 
  separate(X9, c("id","parent"), sep = ";") 
gffl <- gffl %>% mutate(TXNAME = gsub("ID=", "", gffl$id), 
                        GENEID = gsub("Parent=", "", gffl$parent))

tx2gene <- select(gffl, c(TXNAME, GENEID))

# Read in count data, using 'countsfromabundance' option to make compatible with limma
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")

################################################################
# DEseq2 implementation
################################################################

pheno_data   <- read.delim("data/stress_expt_meta.txt")
t <- paste(pheno_data[["treatment"]], pheno_data[["sampling.time.in.hrs"]], sep = ".")
pheno_data   <- arrange(pheno_data, sample) %>% 
  mutate(treatment = t)

dds <- DESeqDataSetFromTximport(txi,
                                   colData = pheno_data,
                                   design = ~ 0 + treatment)

dds <- DESeq(dds)
res <- results(dds)
as.data.frame(colData(dds))
hotDiff3hr <- results(dds, 
                      contrast=list(c("treatmenthot.3", "treatmenthot.1"), 
                                    c("treatmentambient.3", "treatmentambient.1")), 
                      listValues=c(1, -1),
                      alpha = 0.05)
write.csv(as.data.frame(hotDiff3hr), file = "salmon_results/hotDiff3hr_deseq.csv")

summary(hotDiff3hr)
plotMA(hotDiff3hr)
plotCounts(dds, gene=which.min(hotDiff3hr$padj), intgroup="treatment")
plotDispEsts(dds)
plot(metadata(hotDiff3hr)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(hotDiff3hr)$lo.fit, col="red")
abline(v=metadata(hotDiff3hr)$filterTheta)

resultsNames(hotDiff3hr)

coldDiff1hr <- results(dds, 
                       contrast=list(c("treatmentcold.1", "treatmentambient.0"), 
                                     c("treatmentambient.1", "treatmentambient.0")), 
                       listValues=c(1, -1),
                       alpha = 0.05)
write.csv(as.data.frame(coldDiff1hr), file = "salmon_results/coldDiff1hr_deseq.csv")

coldDiff3hr <- results(dds, 
                      contrast=list(c("treatmentcold.3", "treatmentcold.1"), 
                                    c("treatmentambient.3", "treatmentambient.1")), 
                      listValues=c(1, -1),
                      alpha = 0.05)
write.csv(as.data.frame(coldDiff3hr), file = "salmon_results/coldDiff3hr_deseq.csv")

coldDiff9hr <- results(dds, 
                       contrast=list(c("treatmentcold.9", "treatmentcold.3"), 
                                     c("treatmentambient.9", "treatmentambient.3")), 
                       listValues=c(1, -1),
                       alpha = 0.05)
write.csv(as.data.frame(coldDiff9hr), file = "salmon_results/coldDiff9hr_deseq.csv")

# Compare to cold_specific genes identified by limma approach
both <- left_join(as.data.frame(coldDiff3hr) %>% mutate(gene = rownames(as.data.frame(coldDiff3hr))), 
                  as.data.frame(coldDiff9hr) %>% mutate(gene = rownames(as.data.frame(coldDiff9hr))),
                  by = 'gene')
bfilt <- both %>% filter(padj.x <= 0.05 | padj.y <= 0.05) %>% arrange(gene)

limc <- read_tsv("salmon_results/initial_timepoints/cold_specific.txt")
limc <- limc %>% filter(adj.P.Val <= 0.05) %>% arrange(gene)

comp <- inner_join(bfilt, limc, by = 'gene')
lven <- anti_join(limc, bfilt)
rven <- anti_join(bfilt, limc)
################################################################################
# simple model
################################################################################
treatment <- factor(x[["treatment"]])
time      <- x[["sampling.time.in.hrs"]]

ddsm <- DESeqDataSetFromTximport(txi,
                                 colData = pheno_data,
                                 design = ~ sampling.time.in.hrs + treatment)

################################################################################
results(dds, contrast=c("treatment","hot.1","ambient.0"))

ddsMF <- dds
levels(ddsMF$treatment)
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
# filt <- which(pheno_data$type.of.stress == "abiotic" & pheno_data$sampling.time.in.hrs <= 9)
# ddsf <- dds[,filt]
# ddsf$treatment <- relevel(ddsf$treatment, ref = "ambient")
# 
# x <- pheno_data[filt,]
# treatment <- paste(x[["treatment"]], x[["sampling.time.in.hrs"]], sep = ".")
# treatment <- factor(treatment)
# 
# design <- model.matrix(~ 0 + treatment)
# colnames(design) <- levels(treatment)
# 
# ddsf <- DESeq(ddsf, full = "design", betaPrior = FALSE, modelMatrixType = "standard")

# This removed no rows
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

