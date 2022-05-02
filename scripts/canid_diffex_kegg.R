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

# kegg pathway enrichment -------------------------------------------------

## contrast is: log2 fold change (MLE): diet Wild vs Chow 
## this means that for these results, "repressed" means the value is 
## lower with the "wild" diet.

plotCounts(dds, gene = "K03293", intgroup = "diet", 
           normalized = TRUE, transform = TRUE)

## download ortholog to pathway table
url <- "http://rest.kegg.jp/link/pathway/ko"
download.file(url = url, destfile = "inputs/kegg_pathways.tsv")
pathways <- read_tsv("inputs/kegg_pathways.tsv", col_names = c("KO", "path"))
pathways <- pathways %>%
  mutate(KO = gsub("ko:", "", KO)) %>%
  filter(grepl(pattern = "ko", path)) %>%
  mutate(path = gsub("path:", "", path)) %>%
  select(path, KO)

## format results
## separate to the last kegg ortholog if multiple are observed
res_sig$ortholog <- gsub(".*,","", res_sig$ortholog)

## separate to the last kegg ortholog if multiple are observed
res_kegg$ortholog <- gsub(".*,","", res_kegg$ortholog)
res_kegg$ko <- gsub("K", "", res_kegg$ortholog)

res_ind <- res_sig %>% 
  filter(log2FoldChange > 0)

enriched_ind <- enricher(gene = res_ind$ortholog, TERM2GENE = pathways)
dotplot(enriched_ind) +
  ggtitle("'Induced' in ")

res_rep <- res_sig %>% 
  filter(log2FoldChange < 0)
enriched_rep <- enricher(gene = res_rep$ortholog, TERM2GENE = pathways)
dotplot(enriched_rep) +
  ggtitle("'Repressed'")
