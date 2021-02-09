library(stringr)
library(readr)
library(tximport)
library(dplyr)
library(tidyr)

# Create named vector of salmon qunt file locations
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

# Match input format for limma.R script
Z <- txi$counts
# This did not seem to remove 
fpkm.threshold = 5
row.max.fpkm = apply(X, 1, max)
Z = Z[row.max.fpkm > fpkm.threshold, ]
Z <- t(Z)
