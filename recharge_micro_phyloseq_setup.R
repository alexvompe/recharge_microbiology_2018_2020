library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid) 
library(reshape2)
library(phyloseq)
library(readxl)
library(microbiome)
library(decontam)

#import data
asv_mat1 = read_excel("Seq1_metadata.xlsx", sheet = "asv matrix")
tax_table1 = read_excel("Seq1_metadata.xlsx", sheet = "taxonomy table")
samples_df1 = read_excel("Seq1_metadata.xlsx", sheet = "samples")

#make phyloseq object
asv_mat1 = asv_mat1 %>%
  tibble::column_to_rownames("asv")
tax_table1 = tax_table1 %>% 
  tibble::column_to_rownames("asv")
samples_df1 = samples_df1 %>% 
  tibble::column_to_rownames("sample")
asv_mat1 = as.matrix(asv_mat1)
tax_table1 = as.matrix(tax_table1)
OTU1 = otu_table(asv_mat1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax_table1)
samples1 = sample_data(samples_df1)
ps1 = phyloseq(OTU1, TAX1, samples1)
ps1

#import data
asv_mat2 = read_excel("Seq2_metadata.xlsx", sheet = "asv matrix")
tax_table2 = read_excel("Seq2_metadata.xlsx", sheet = "taxonomy table")
samples_df2 = read_excel("Seq2_metadata.xlsx", sheet = "samples")

#make phyloseq object
asv_mat2 = asv_mat2 %>%
  tibble::column_to_rownames("asv")
tax_table2 = tax_table2 %>% 
  tibble::column_to_rownames("asv")
samples_df2 = samples_df2 %>% 
  tibble::column_to_rownames("sample")
asv_mat2 = as.matrix(asv_mat2)
tax_table2 = as.matrix(tax_table2)
OTU2 = otu_table(asv_mat2, taxa_are_rows = TRUE)
TAX2 = tax_table(tax_table2)
samples2 = sample_data(samples_df2)
ps2 = phyloseq(OTU2, TAX2, samples2)
ps2

#import data
asv_mat3 = read_excel("Seq3_metadata.xlsx", sheet = "asv matrix")
tax_table3 = read_excel("Seq3_metadata.xlsx", sheet = "taxonomy table")
samples_df3 = read_excel("Seq3_metadata.xlsx", sheet = "samples")

#make phyloseq object
asv_mat3 = asv_mat3 %>%
  tibble::column_to_rownames("asv")
tax_table3 = tax_table3 %>% 
  tibble::column_to_rownames("asv")
samples_df3 = samples_df3 %>% 
  tibble::column_to_rownames("sample")
asv_mat3 = as.matrix(asv_mat3)
tax_table3 = as.matrix(tax_table3)
OTU3 = otu_table(asv_mat3, taxa_are_rows = TRUE)
TAX3 = tax_table(tax_table3)
samples3 = sample_data(samples_df3)
ps3 = phyloseq(OTU3, TAX3, samples3)
ps3

#Merge Phyloseq objects for downstream pre-processing and analyses
microbes = merge_phyloseq(ps1, ps2, ps3)
microbes

#Change ranks
colnames(tax_table(microbes)) = c("Kingdom", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species")
microbes #n=759
rank_names(microbes)

#Remove mitochondria, chloroplasts, any other non-bacterial/archaeal seqs
microbes = microbes %>% subset_taxa(Family!= "Mitochondria" | is.na(Family))
microbes = microbes %>% subset_taxa(Order!= "Chloroplast" | is.na(Order))
microbes = microbes %>% subset_taxa(Kingdom!= "Eukaryota" | is.na(Kingdom))
microbes = microbes %>% subset_taxa(Kingdom!= "NA" | is.na(Kingdom))
microbes

#Remove NA problem samples
ps1 = subset_samples(microbes, Nutrients!="NA")
NoSampleNA = subset_samples(ps1, Herbivory!="NA")
NoSampleNA = subset_samples(NoSampleNA, Run!="NA")
NoProblem = subset_samples(NoSampleNA, Nutrients!="multi")
OnlyCorals = subset_samples(NoProblem, Coral!="Hal")

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(OnlyCorals)$is.neg = sample_data(OnlyCorals)$Coral == "control"
contamdf.prev = isContaminant(OnlyCorals, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant) #34 contam ASVs
head(which(contamdf.prev$contaminant))

#Remove contaminants
physeq.noncont = prune_taxa(!contamdf.prev$contaminant, OnlyCorals)

#Remove controls as we have used them for decontam
physeq.noncont = subset_samples(physeq.noncont, Coral!="control")

#Remove singletons
NoSingle = prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)

#Let's get the final numbers and remove low-read samples
NoSingle #748 samples, 19702 taxa
NoSingleover1000 = prune_samples(sample_sums(NoSingle)
                                 >=1000, NoSingle)

#Extract the sequences of the filtered object
write.csv(tax_table(NoSingleover1000), "allrunsseqs.csv")

#Change ASV names for more simple work
taxa_names(microbes) = paste0("ASV", seq(ntaxa(microbes)))
tax_table(microbes) = cbind(tax_table(microbes),
                            rownames(tax_table(microbes)))
head(taxa_names(microbes))

#Export filtered sequence and taxonomy tables for full phyloseq object
saveRDS(otu_table(microbes), "filtered_seq_table.rds")
saveRDS(sample_data(microbes), "filtered_sample_data.rds")
saveRDS(tax_table(microbes), "filtered_tax_table.rds")