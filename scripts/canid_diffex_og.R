setwd("~/github/dom")

library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(clusterProfiler)

# import counts -----------------------------------------------------------

txi_og <- read_tsv("outputs/tx2gene/eggnog_tx2gene_og.tsv")
canids <- read_csv("inputs/tximport_samples/canid_samples.csv")
canid_counts_og <- tximport(files = canids$filepath, type = "salmon", tx2gene = txi_og)

# differential abundance -------------------------------------------------

ddsTxi_og <- DESeqDataSetFromTximport(canid_counts_og,
                                      colData = canids,
                                      design = ~ species + diet)
dds_og <- DESeq(ddsTxi_og)
resultsNames(dds_og) # lists the coefficients
res <- results(dds_og, name="diet_Wild_vs_Chow")
res <- res[order(res$pvalue), ]

res_sig <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("@", ortholog)) %>%
  filter(padj < .05) 

res_og <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("@", ortholog))

# enrichment --------------------------------------------------------------

## separate to the last go ortholog if multiple are observed
res_sig$ortholog <- gsub(".*,","", res_sig$ortholog)

## separate to the last go ortholog if multiple are observed
res_og$ortholog <- gsub(".*,","", res_go$ortholog)
res_og$ko <- gsub("K", "", res_og$ortholog)

res_ind <- res_sig %>% 
  filter(log2FoldChange > 0) %>%
  mutate(KO = gsub("K", "", ortholog))
enriched_ind <- enrichgo(gene = res_ind$KO, species = "ko")
dotplot(enriched_ind)
View(enriched_ind@result)

res_rep <- res_sig %>% 
  filter(log2FoldChange < 0) %>%
  mutate(KO = gsub("K", "", ortholog))
enriched_rep <- enrichGO(gene = res_rep$Ko)
dotplot(enriched_rep)
View(enriched_rep@result)
