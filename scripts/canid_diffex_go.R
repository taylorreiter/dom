setwd("~/github/dom")

library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(topGO)

# import counts -----------------------------------------------------------

txi_go <- read_tsv("outputs/tx2gene/eggnog_tx2gene_go.tsv")
canids <- read_csv("inputs/tximport_samples/canid_samples.csv")
canid_counts_go <- tximport(files = canids$filepath, type = "salmon", tx2gene = txi_go)

# differential abundance -------------------------------------------------

ddsTxi_go <- DESeqDataSetFromTximport(canid_counts_go,
                                      colData = canids,
                                      design = ~ species + diet)
dds_go <- DESeq(ddsTxi_go)
resultsNames(dds_go) # lists the coefficients
res <- results(dds_go, name="diet_Wild_vs_Chow")
res <- res[order(res$pvalue), ]

res_sig <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("GO", ortholog)) %>%
  filter(padj < .05) %>%
  mutate(id = seq(1:nrow(.)))

res_go <- res %>%
  as.data.frame %>%
  mutate(ortholog = rownames(.)) %>%
  filter(grepl("GO", ortholog))

# enrichment example --------------------------------------------------------------
# how to construct topGO object for enrichment:
# 1. Make a tab-separated file with GeneID and GO terms. GO terms need to be i
#    one column but should be separated by a comma. Export this file to a tsv
# 2. Read in the exported file using the readMappings() function. Name this obj
#    geneID2GO
# 3. Using the logFC table, make a named vector of p values
# 4. Define a topDiffGenes function selector (p < .05)
# 5. run new() function


head(res_go)
geneID2GO <- res_go %>%
  mutate(row_num = row.names(.)) %>%
  dplyr::select(row_num, ortholog)

write_tsv(geneID2GO, "outputs/diff_abund/geneID2GO_canid.tsv")
geneID2GO <- readMappings("outputs/diff_abund/geneID2GO_canid.tsv")
# define a function to select top genes
topDiffGenes <- function(allScore) {return(allScore < 0.05)}

# make a vector of pvalues named by row num
res_go_genes <- res_go$pvalue
names(res_go_genes) <- rownames(res_go)

## CC (cellular compartment) ontology
cc_go <- new("topGOdata", description= "Canid Diet", ontology = "CC", 
              allGenes = res_go_genes,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO, 
              geneSel = topDiffGenes)

cc_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
go_cc_res <- GenTable(GOdata, classicFisher = cc_fisher,
                  ranksOf = "classicFisher", topNodes = 50)
head(cc_go_res)
kable(cc_go_res)

## MF (molecular function) ontology
mf_go <- new("topGOdata", description= "Canid Diet", ontology = "MF", 
             allGenes = res_go_genes,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO, 
             geneSel = topDiffGenes)

mf_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
go_mf_res <- GenTable(GOdata, classicFisher = mf_fisher,
                      ranksOf = "classicFisher", topNodes = 50)
head(mf_go_res)
kable(mf_go_res)

## BP (biological process) ontology
bp_go <- new("topGOdata", description= "Canid Diet", ontology = "BP", 
             allGenes = res_go_genes,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO, 
             geneSel = topDiffGenes)

bp_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
go_bp_res <- GenTable(GOdata, classicFisher = bp_fisher,
                      ranksOf = "classicFisher", topNodes = 50)
head(bp_go_res)
kable(bp_go_res)
