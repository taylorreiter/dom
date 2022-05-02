# This script generates tables of counts for each of the eggnog ortholog group
# categories. To do this, I first use tximport to import all of the gene-level
# ortholog counts and summarize them to per-sample ortholog-level counts.

library(readr)
library(dplyr)
library(tidyr)
library(tximport)

# read in file that contains eggnog gene ID to category mapping
eggnog <- read_tsv("outputs/atlas_out/Genecatalog/annotations/eggNog.tsv.gz")

cazy_tx2gene <- eggnog %>%
  select(Query, CAZy) %>%
  filter(!is.na(CAZy)) %>%
  dplyr::rename(TXNAME = Query, GENEID = CAZy)


# observational study -----------------------------------------------------

obs <- read_csv("inputs/tximport_samples/observational_samples.csv")
# summarize gene counts to ortholog counts
cazy_count <- tximport(files = obs$filename, type = "salmon", tx2gene = cazy_tx2gene)
cazy_count <- cazy_count$counts
colnames(cazy_count) <- obs$sample
cazy_count <- as.data.frame(cazy_count)
cazy_count$CAZy <- rownames(cazy_count)
write_tsv(cazy_count, "outputs/functional_tables/cazy/cazy_raw_counts_observational.tsv")


# canid study -------------------------------------------------------------

canid <- read_csv("inputs/tximport_samples/canid_samples.csv")
# summarize gene counts to ortholog counts
cazy_count <- tximport(files = canid$filepath, type = "salmon", tx2gene = cazy_tx2gene)
cazy_count <- cazy_count$counts
colnames(cazy_count) <- canid$sample
cazy_count <- as.data.frame(cazy_count)
cazy_count$CAZy <- rownames(cazy_count)
write_tsv(cazy_count, "outputs/functional_tables/cazy/cazy_raw_counts_canid.tsv")

# mouse study -------------------------------------------------------------
mouse<- read_csv("inputs/tximport_samples/mouse_samples.csv")
# summarize gene counts to ortholog counts
cazy_count <- tximport(files = mouse$filepath, type = "salmon", tx2gene = cazy_tx2gene)
cazy_count <- cazy_count$counts
colnames(cazy_count) <- mouse$SampleID
cazy_count <- as.data.frame(cazy_count)
cazy_count$CAZy <- rownames(cazy_count)
write_tsv(cazy_count, "outputs/functional_tables/cazy/cazy_raw_counts_mouse.tsv")


