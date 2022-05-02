# This script generates tables of counts for each CAZy group for each sample.
# I first used tximport to import all of the gene-level
# ortholog counts and summarize them to per-sample CAZy-level counts.
# Because CAZys make up a small proportion of total orthologs, I did not
# perform normalization as I did with e.g. EggNog categories.

library(readr)
library(dplyr)
library(tidyr)
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
txi_og <- read_tsv("outputs/tx2gene/eggnog_tx2gene_og.tsv")
# only keep genes that were annotated as orthologs
txi_og <- filter(txi_og, grepl("@", GENEID))
# read table that maps sample name to file name that contains gene counts. 
obs <- read_csv("inputs/tximport_samples/observational_samples.csv")
# summarize gene counts to ortholog counts
#obs_counts_og <- tximport(files = obs$filename, type = "salmon", tx2gene = txi_og)
obs_counts_og <- readRDS("outputs/diff_abund/observational_species/tximport_og.RDS")
obs_counts_og <- obs_counts_og$counts
colnames(obs_counts_og) <- obs$sample
obs_counts_og <- as.data.frame(obs_counts_og)
obs_counts_og$EggNog <- rownames(obs_counts_og) 

# subset eggnog by category, and use this subset to separate counts per category
# and write this to a file
for(category_i in unique(eggnog_category$category)){
  eggnog_sub_category <- filter(eggnog_category, category == category_i)
  category_counts <- filter(obs_counts_og, EggNog %in% eggnog_sub_category$EggNog)
  write_tsv(category_counts, paste0("outputs/functional_tables/eggnog/observational_eggnog_category_",
                                    category_i, "_raw_counts.tsv"))
}

