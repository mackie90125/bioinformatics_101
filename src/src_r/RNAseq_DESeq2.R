###############################################################################################
###
###
###
###############################################################################################

library(tidyverse)

# install libraries from Bioconductor
# BiocManager::install("DESeq2")
# BiocManager::install("airway")
library(DESeq2)
library(airway)


###############################################################################################
### get data from airway package and re-format

data(airway)
head(airway)

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "/Users/ryan/test_analysis/bioinformatics_101/data/sample_info.csv",
            sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "/Users/ryan/test_analysis/bioinformatics_101/data/counts_data.csv",
            sep = ',', col.names = T, row.names = T, quote = F)

###############################################################################################
### DESeq2

# ensure row names in sample_info match the columns in countsData
all(colnames(countsData) %in% rownames(sample_info))
all(colnames(countsData) == rownames(sample_info))

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countsData, colData = sample_info, design = ~ dexamethasone)
dds

# recommended pre-filtering step
# remove rows with less than 10 reads across all samples
keep_rows <- rowSums(counts(dds)) >= 10
dds <- dds[keep_rows, ]
dds

# set reference level for design variable
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# collapse all technical replicates before running DESeq function
# use collapseReplicates() from DESeq2 library

# run DESeq
dds <- DESeq(dds)

# save results
res <- results(dds)
head(res)

summary(res)

# change significance cut-off level
res01 <- results(dds, alpha = 0.01)
summary(res01)

# contrasts
# this is used to define pairs to compare if our design has more than two groups
resultsNames(dds)

# volcano plot
# blue indicates genes that are significanlty differentially expressed
plotMA(res01)