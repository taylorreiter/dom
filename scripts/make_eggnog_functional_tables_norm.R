# This script generates tables of counts for each of the eggnog ortholog group
# categories. To do this, I first use tximport to import all of the gene-level
# ortholog counts and summarize them to per-sample ortholog-level counts.
# Then, we create a DESeq object as if we were doing differential abundance
# analysis, and extract counts from this object. We do this because counts 
# need to be normalized for both library size and for gene length, and this is 
# performed automatically (and correctly :) ) by DESeq. 

library(readr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(tximport)

# read in file that contains eggnog gene ID to category mapping
eggnog <- read_tsv("outputs/atlas_out/Genecatalog/annotations/eggNog.tsv.gz")

# for eggnogs with multiple categories, separate into different columns
eggnog_category <- separate(eggnog, col = FunctionalCategory, sep = c(1,2,3,4,5), remove = F,
                            into = c("category_1", "category_2", "category_3", 
                                     "category_4", "category_5", "category_6"))


# make a 1-to-1 mapping of eggnog category to Ortholog Group Identifier
eggnog_category <- eggnog_category %>% 
  select(Query, EggNog, category_1, category_2, category_3, category_4, category_5, category_6) %>%
  distinct() %>%
  pivot_longer(cols = c(-Query, -EggNog), names_to = "tmp",values_to = "category", values_drop_na = T) %>%
  filter(category != "") %>%
  select(EggNog, category) %>%
  distinct()

# read table that maps gene ID to ortholog for ortholog groups
# txi_og <- read_tsv("outputs/tx2gene/eggnog_tx2gene_og.tsv")
# only keep genes that were annotated as orthologs
# txi_og <- filter(txi_og, grepl("@", GENEID))
# read table that maps sample name to file name that contains gene counts. 
# obs <- read_csv("inputs/tximport_samples/observational_samples.csv")
# summarize gene counts to ortholog counts
# obs_counts_og <- tximport(files = obs$filename, type = "salmon", tx2gene = txi_og)
# ddsTxi_og <- DESeqDataSetFromTximport(obs_counts_og,
#                                      colData = obs,
#                                      design = ~ status + diet_type + physiology)
# dds_og <- DESeq(ddsTxi_og)
dds_og <- readRDS("outputs/diff_abund/observational_species/dds_og.RDS")
dds_og <- estimateSizeFactors(dds_og)
dds_og_counts <- counts(dds_og, normalized=TRUE)
colnames(dds_og_counts) <- colData(dds_og)$sample
dds_og_counts <- as.data.frame(dds_og_counts)
dds_og_counts$EggNog <- rownames(dds_og_counts)

# subset eggnog by category, and use this subset to separate counts per category
# and write this to a file

for(category_i in unique(eggnog_category$category)){
  eggnog_sub_category <- filter(eggnog_category, category == category_i)
  category_counts <- filter(dds_og_counts, EggNog %in% eggnog_sub_category$EggNog)
  write_tsv(category_counts, paste0("outputs/functional_tables/eggnog/observational_eggnog_category_",
                                    category_i, "_norm_counts.tsv"))
}
