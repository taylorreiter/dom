setwd("~/github/dom")

library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(clusterProfiler)

# import counts -----------------------------------------------------------

txi <- read_tsv("outputs/tx2gene/eggnog_tx2gene_kegg.tsv")
canids <- read_csv("inputs/tximport_samples/canid_samples.csv")
canid_counts <- tximport(files = canids$filepath, type = "salmon", tx2gene = txi)

# differential abundance -------------------------------------------------

ddsTxi <- DESeqDataSetFromTximport(canid_counts,
                                   colData = canids,
                                   design = ~ species + diet)
dds <- DESeq(ddsTxi)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="diet_Wild_vs_Chow")
res <- res[order(res$pvalue), ]

res_sig <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("K", ortholog)) %>%
  filter(padj < .05) 

res_kegg <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("K", ortholog))

# enrichment --------------------------------------------------------------

## separate to the last kegg ortholog if multiple are observed
res_sig$ortholog <- gsub(".*,","", res_sig$ortholog)

## separate to the last kegg ortholog if multiple are observed
res_kegg$ortholog <- gsub(".*,","", res_kegg$ortholog)
res_kegg$ko <- gsub("K", "", res_kegg$ortholog)

res_ind <- res_sig %>% 
  filter(log2FoldChange > 0) %>%
  mutate(KO = gsub("K", "", ortholog))
enriched_ind <- enrichKEGG(gene = res_ind$KO, species = "ko")
dotplot(enriched_ind)
View(enriched_ind@result)

res_rep <- res_sig %>% 
  filter(log2FoldChange < 0) %>%
  mutate(KO = gsub("K", "", ortholog))
enriched_rep <- enrichKEGG(gene = res_rep$KO, 
                           organism = "eco",
                           keyType = "kegg")
dotplot(enriched_rep)
View(enriched_rep@result)
