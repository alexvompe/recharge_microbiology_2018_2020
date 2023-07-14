library(reshape2)

#Transpose ASV table (create for seq_table_1, 2, and 3)
asvtab = readRDS("seq_table_3.rds")
asvtabcorrected = t(asvtab)
write.csv(asvtabcorrected, "corrected_seq_table_3.csv")

#Save taxonomy table as csv for downstream
taxtab = readRDS("seq3_taxa.rds")
write.csv(taxtab, "seq3_taxa.csv")