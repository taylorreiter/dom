setwd("~/github/dom")

library(readr)
library(dplyr)

eggnog <- read_tsv("outputs/atlas_out/Genecatalog/annotations/eggNog.tsv.gz")

## TX2GENE GO
#tx2gene_go <- eggnog %>%
#  select(query_name, GO_terms) %>%
#  rename(TXNAME = query_name, GENEID = GO_terms)
#tx2gene_go$GENEID <- ifelse(is.na(tx2gene_go$GENEID), tx2gene_go$TXNAME, tx2gene_go$GENEID)
#write_tsv(tx2gene_go, "outputs/tx2gene/eggnog_tx2gene_go.tsv")

## TX2GENE KEGG
tx2gene_kegg <- eggnog %>%
  select(Query, KO) %>%
  rename(TXNAME = Query, GENEID = KO)
tx2gene_kegg$GENEID <- ifelse(is.na(tx2gene_kegg$GENEID), tx2gene_kegg$TXNAME, tx2gene_kegg$GENEID)
write_tsv(tx2gene_kegg, "outputs/tx2gene/eggnog_tx2gene_kegg.tsv")

## TX2GENE OG
tx2gene_og <- eggnog %>%
  select(Query, EggNog) %>%
  rename(TXNAME = Query, GENEID = EggNog)
tx2gene_og$GENEID <- ifelse(is.na(tx2gene_og$GENEID), tx2gene_og$TXNAME, tx2gene_og$GENEID)
write_tsv(tx2gene_og, "outputs/tx2gene/eggnog_tx2gene_og.tsv")
