###############################################################################################
### manipulate and graph gene expression data
###
### data in csv file has already been processed from raw reads to fpkm
###############################################################################################

library(tidyverse)
library(GEOquery)

###############################################################################################
### load and reformat data ###

# load data
d.raw <- read.csv("/Users/ryan/test_analysis/bioinformatics_101/data/GSE183947_fpkm.csv")

# get metadata using GEOquery package
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
d.meta <- pData(phenoData(gse[[1]])) %>%

  select(c(1,10,11,17)) %>%

  rename(tissue = characteristics_ch1, metastasis = characteristics_ch1.1, samples = description) %>%

  mutate(metastasis = str_remove(metastasis, "metastasis: "),
         tissue = if_else(grepl("tumor", tissue), "tumor", "normal"))

# convert to long format and merge with meta-data
d.all <- rename(d.raw, gene = X) %>%

  pivot_longer(!gene, names_to = "samples", values_to = "FPKM") %>%

  arrange(samples, gene) %>%

  left_join(., d.meta, by = "samples")


###############################################################################################
### explore data ###

# BRCA gene expression summary statistics
filter(d.all, gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM))

# bar plot
filter(d.all, gene == "BRCA1") %>%
  ggplot(aes(samples, FPKM, fill = tissue)) +
  geom_col()

# density plot
filter(d.all, gene == "BRCA1") %>%
  ggplot(aes(FPKM, fill = tissue)) +
  geom_density(alpha = 0.4)

# box plot
filter(d.all, gene == "BRCA1") %>%
  ggplot(aes(metastasis, FPKM)) +
  geom_boxplot(fill = "lightblue")
  # geom_jitter(width = 0.1)

# box plot
filter(d.all, gene == "BRCA1" | gene == "BRCA2") %>%
  pivot_wider(names_from = gene, values_from = FPKM) %>%
  ggplot(aes(BRCA1, BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, size = 0.5)

# heatmap
# ordered by TP53 expression
genes_interest <- c("BRCA1","BRCA2","TP53","ALK","MYCN")
filter(d.all, gene %in% genes_interest) %>%
  pivot_wider(names_from = gene, values_from = FPKM) %>%
  mutate(samples = factor(samples)) %>%
  mutate(samples = fct_reorder(samples, TP53)) %>%
  pivot_longer(genes_interest, names_to = "gene", values_to = "FPKM") %>%
  ggplot(aes(samples, gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue")