# This script generates tables of counts for BRITE hierarchies or KEGG pathways
# annotated from KEGG
# For example: BRITE hierarchy 01504 encodes antimicrobial resistance genes. 
# I first used tximport to import all of the gene-level KEGG ortholog counts and 
# summarize them to per-sample counts. 
# Then, I create a DESeq object as if we were doing differential abundance
# analysis, and extract counts from this object. We do this because counts 
# need to be normalized for both library size and for gene length, and this is 
# performed automatically by DESeq. Then, these are filtered to orthologs
# that are classified as antimicrobial resistance in the BRITE hierarchy.

library(readr)
library(dplyr)
library(tidyr)
library(tximport)
library(DESeq2)

# import raw counts -------------------------------------------------------

txi_kegg <- read_tsv("outputs/tx2gene/eggnog_tx2gene_kegg.tsv")
txi_kegg <- filter(txi_kegg, grepl("K", GENEID))
info <- read_csv("inputs/tximport_samples/all_samples.csv")
all_counts <- tximport(files = info$filename, type = "salmon", tx2gene = txi_kegg)
ddsTxi_kegg <- DESeqDataSetFromTximport(all_counts,
                                        colData = info,
                                        design = ~ species)
dds_kegg <- DESeq(ddsTxi_kegg)
dds_kegg <- estimateSizeFactors(dds_kegg)
dds_kegg_counts <- counts(dds_kegg, normalized=TRUE)
colnames(dds_kegg_counts) <- colData(dds_kegg)$sample
dds_kegg_counts <- as.data.frame(dds_kegg_counts)
dds_kegg_counts$KO <- rownames(dds_kegg_counts)

# import annotation file --------------------------------------------------
eggnog <- read_tsv("outputs/atlas_out/Genecatalog/annotations/eggNog.tsv.gz")


# KO to KEGG Pathway map --------------------------------------------------
kegg_pathways <- eggnog %>%
  select("KO", "KEGG_Pathway") %>%
  separate(col = KEGG_Pathway, sep = ",", 
           into = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9",
                    "p10", "p11", "p12", "p13", "p14", "p15", "p16", "p17",
                    "p18", "p19", "p20", "p21", "p22", "p23", "p24", "p25",
                    "p26", "p27", "p28", "p29", "p30", "p31", "p32", "p33",
                    "p34", "p35", "p36", "p37", "p38", "p39", "p40", "p41",
                    "p42", "p43", "p44", "p45", "p46", "p47", "p48", "p49",
                    "p50", "p51", "p52", "p53", "p54", "p55", "p56", "p57",
                    "p58", "p59", "p60", "p61", "p62", "p63", "p64", "p65",
                    "p66", "p67", "p68", "p69", "p70", "p71", "p72", "p73",
                    "p74", "p75", "p76", "p77", "p78", "p79", "p80", "p81",
                    "p82", "p83", "p84", "p85", "p86", "p87", "p88", "p89",
                    "p90", "p91", 
                    "p92", "p93", "p94", "p95", "p96", "p97", "p98", "p99",
                    "p100", "p101", "p102", "p103", "p104", "p105", "p106", 
                    "p107", "p108", "p109", 'p110', 'p111', 'p112', 'p113', 
                    'p114', 'p115', 'p116', 'p117', 'p118', 'p119', 'p120', 
                    'p121', 'p122', 'p123', 'p124', 'p125', 'p126', 'p127', 
                    'p128', 'p129', 'p130', 'p131', 'p132', 'p133', 'p134', 
                    'p135', 'p136', 'p137', 'p138', 'p139', 'p140', 'p141', 
                    'p142', 'p143', 'p144', 'p145', 'p146', 'p147', 'p148', 
                    'p149', 'p150', 'p151', 'p152', 'p153', 'p154', 'p155', 
                    'p156', 'p157', 'p158', 'p159', 'p160', 'p161', 'p162', 
                    'p163', 'p164', 'p165', 'p166', 'p167', 'p168', 'p169', 
                    'p170', 'p171', 'p172', 'p173', 'p174', 'p175', 'p176', 
                    'p177', 'p178', 'p179', 'p180', 'p181', 'p182', 'p183', 
                    'p184', 'p185', 'p186', 'p187', 'p188', 'p189', 'p190', 
                    'p191', 'p192', 'p193', 'p194', 'p195', 'p196', 'p197', 
                    'p198', 'p199', 'p200', "p201", "p202", "p203", "p204",
                    'p205', "p206", "p207", "p208", "p209", "p210", "p211",
                    "p212", "p213")) %>%
  pivot_longer(cols = -KO, names_to = "tmp", values_to = "KEGG_Pathway") %>%
  select(KO, KEGG_Pathway)

kegg_pathways <- kegg_pathways %>%
  filter(!is.na(KEGG_Pathway)) %>%
  distinct()

write_tsv(kegg_pathways, "outputs/functional_tables/ko_to_pathway_map.tsv")

# filter to specific sets of counts ---------------------------------------

## ANTIMICROBIAL RESISTANCE
# eggnog_amr <- eggnog  %>%
#   select(KO, BRITE) %>%
#   filter(grepl("01504", BRITE)) %>%
#   distinct()
# 
# norm_amr <- dds_kegg_counts %>%
#   filter(KO %in% eggnog_amr$KO)
# write_tsv(norm_amr, "outputs/functional_tables/antibiotic_resistance/all_amr_norm_counts.tsv")

## VITAMINS
# 09108 Metabolism of cofactors and vitamins
# 00730 Thiamine metabolism [PATH:ko00730]
# 00740 Riboflavin metabolism [PATH:ko00740]
# 00750 Vitamin B6 metabolism [PATH:ko00750]
# 00760 Nicotinate and nicotinamide metabolism [PATH:ko00760]
# 00770 Pantothenate and CoA biosynthesis [PATH:ko00770]
# 00780 Biotin metabolism [PATH:ko00780]
# 00785 Lipoic acid metabolism [PATH:ko00785]
# 00790 Folate biosynthesis [PATH:ko00790]
# 00670 One carbon pool by folate [PATH:ko00670]
# 00830 Retinol metabolism [PATH:ko00830]
# 00860 Porphyrin and chlorophyll metabolism [PATH:ko00860]
# 00130 Ubiquinone and other terpenoid-quinone biosynthesis [PATH:ko00130]

# tmp <- eggnog[(grepl("ko00730", eggnog$KEGG_Pathway)), ]
# 
# vitamin_kegg_pathway <- c("ko00730", "ko00740", "ko00750", "ko00760", "ko00770", 
#                           "ko00780", "ko00785", "ko00790", "ko00670", "ko00830", 
#                           "ko00860", "ko00130")
# 
# eggnog_vitamin <- kegg_pathways %>% 
#   filter(KEGG_Pathway %in% vitamin_kegg_pathway)
# 
# norm_vitamin <- dds_kegg_counts %>%
#   filter(KO %in% eggnog_vitamin$KO)
# write_tsv(norm_vitamin, "outputs/functional_tables/vitamins/all_vitamins_norm_counts.tsv")
# 
# ## METABOLISM: koO11OO
# tmp <- kegg_pathways[grepl("01100", kegg_pathways$KEGG_Pathway), ]
# eggnog_metabolism<- kegg_pathways %>% 
#   filter(KEGG_Pathway == "ko01100")
# 
# norm_metabolism <- dds_kegg_counts %>%
#   filter(KO %in% eggnog_metabolism$KO)
# write_tsv(norm_metabolism, "outputs/functional_tables/metabolism/all_ko01100_norm_counts.tsv")

## STARCH AND SUCROSE METABOLISM
eggnog_starch <- eggnog  %>%
  select(KO, KEGG_Pathway) %>%
  filter(grepl("map00500", KEGG_Pathway)) %>%
  distinct()

norm_starch <- dds_kegg_counts %>%
  filter(KO %in% eggnog_starch$KO)
write_tsv(norm_starch, "outputs/functional_tables/starch_and_sucrose_metabolism/all_starch_and_sucrose_norm_counts.tsv")
