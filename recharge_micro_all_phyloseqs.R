library(phyloseq)
library(vegan)
library(readxl)
library(plyr)
library(dplyr)

##Filtered phyloseq
sample_filt = readRDS("filtered_sample_data.rds")
seqtab_filt = readRDS("filtered_seq_table.rds")
taxtab_filt = readRDS("filtered_tax_table.rds")
tree_filt = read_tree("seq.nwk")

microbes.filt = phyloseq(sample_filt, seqtab_filt, taxtab_filt, tree_filt)

#Make 'Date' a factor and order the months sequentially
desired_order = list("Jul18","Nov18","Mar19","Aug19","Nov19","Mar20","Aug20")
sample_data(microbes.filt)$Date = factor(sample_data(microbes.filt)$Date, 
                                         levels = desired_order)

#Save as RDS for downstream access and check file
saveRDS(microbes.filt, "filtered_ps.rds")
microbes.filt = readRDS("filtered_ps.rds")

#Save the filtered sample data as a csv file for mortality and bleaching data integration
write.csv(sample_data(microbes.filt), 
          "filtered sample data_mortality_bleaching integration.csv")

##Rarefied phyloseq
#Rarefying to 1000 reads post filtering
microbes.rarefied = rarefy_even_depth(microbes.filt, 
                                      rngseed=1, 
                                      sample.size=1000, 
                                      replace=F)

#Save as RDS for downstream access and check file
saveRDS(microbes.rarefied, "rarefied_ps.RDS")
microbes.rarefied = readRDS("rarefied_ps.RDS")

##Filtered phyloseq at Family level
families = tax_glom(microbes.filt, taxrank = "Family", NArm = FALSE)
sample_data(families)$Date = factor(sample_data(families)$Date, levels=desired_order)

#Save as RDS for downstream access and check file
saveRDS(families, "families_ps.rds")
families = readRDS("families_ps.rds")

##Rarefied phyloseq at Family level
families_rare = tax_glom(microbes.rarefied, taxrank = "Family", NArm = FALSE)

#Save as RDS for downstream access and check file
saveRDS(families_rare, "rarefied_fams_ps.rds")
families_rare = readRDS("rarefied_fams_ps.rds")

##Rarefied phyloseq with top 20 families + "other" taxa (for stacked bar)
Top19Fams = names(sort(taxa_sums(families_rare), TRUE)[1:19])
Top19Fams

#Replace the taxonomy of non-top families with "Other" manually to get 20 taxa w Other
write.csv(tax_table(families_rare), "taxother_rarefied.csv")

#Assemble new phyloseq object with "other" in taxtable
tax_table_other = read_excel("taxother_rarefied.xlsx")
tax_table_other = tax_table_other %>% 
  tibble::column_to_rownames("asv")
tax_table_other = as.matrix(tax_table_other)
TAX = tax_table(tax_table_other)
Rarefied_w_other = phyloseq(TAX, sample_data(families_rare), 
                         otu_table(families_rare))

#Save as RDS for downstream access and check file
saveRDS(Rarefied_w_other, "ps_for_stacked_bar.rds")
Rarefied_w_other = readRDS("ps_for_stacked_bar.rds")

##Phyloseq with mortality and bleaching data
survival_sample = read_excel("survival sample data.xlsx")
survival_sample = survival_sample %>% 
  tibble::column_to_rownames("sample")
survival_sample = sample_data(survival_sample)
survival = phyloseq(otu_table(families), tax_table(families), 
                         survival_sample, phy_tree(families))

#Filtered ps object with no NAs
survival_filt = subset_samples(survival, percent_dead!="NA")

#Make percent_dead and percent bleached continuous variables
sample_data(survival_filt)$percent_dead = as.numeric(sample_data(survival_filt)$percent_dead)
sample_data(survival_filt)$percent_bleached = as.numeric(sample_data(survival_filt)$percent_bleached)
sample_data(survival_filt)$Date = factor(sample_data(survival_filt)$Date,
                                              levels=desired_order)

#Make bins for mortality and bleaching data
sample_df_mortality = as.data.frame(sample_data(survival_filt))

#Mortality bins: acr and plob: 0%, 0%<i<25%, and >25%; pver: none and some.
sample_df_mortality$mort_bin = NA
for (i in 1:nrow(sample_df_mortality)) {
  if(sample_df_mortality$percent_dead[i]==0){
    sample_df_mortality$mort_bin[i] = "No Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>0 & sample_df_mortality$percent_dead[i]<25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>0 & sample_df_mortality$percent_dead[i]<25 
          & sample_df_mortality$Coral[i]=="Plob"){
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>=25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$mort_bin[i] = "High Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>=25 
          & sample_df_mortality$Coral[i]=="Plob"){
    sample_df_mortality$mort_bin[i] = "High Mortality"
  }
  else {
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
}

length(sample_df_mortality$Coral[sample_df_mortality$Coral == "Pver" 
                                 & sample_df_mortality$percent_bleached != 0])
#only 5 Pver corals bleached, no plob corals bleached, so bleaching analyses only on Acr.

sample_df_mortality$bleach_bin = NA
for (i in 1:nrow(sample_df_mortality)) {
  if(sample_df_mortality$percent_bleached[i]==0){
    sample_df_mortality$bleach_bin[i] = "No Bleaching"
  }
  else if(sample_df_mortality$percent_bleached[i]>0 & sample_df_mortality$percent_bleached[i]<25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$bleach_bin[i] = "Some Bleaching"
  }
  else if(sample_df_mortality$percent_bleached[i]>=25 & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$bleach_bin[i] = "High Bleaching"
  }
  else {
    sample_df_mortality$bleach_bin[i] = "NA"
  }
}

sample_df_mortality = sample_data(sample_df_mortality)
sample_data(survival_filt) = sample_df_mortality

#Save the phyloseq object for ease of access
saveRDS(survival_filt, "survival_fams_bins_ps.rds")
survival_filt = readRDS("survival_fams_bins_ps.rds")

##Rarefied phyloseq with mortality and bleaching data
survival_sample = read_excel("survival sample data.xlsx")
survival_sample = survival_sample %>% 
  tibble::column_to_rownames("sample")
survival_sample = sample_data(survival_sample)
survival = phyloseq(otu_table(families_rare), tax_table(families_rare), 
                    survival_sample, phy_tree(families_rare))

#Filtered ps object with no NAs
survival_rare_filt = subset_samples(survival, percent_dead!="NA")

#Make percent_dead and percent bleached continuous variables
sample_data(survival_rare_filt)$percent_dead = as.numeric(sample_data(survival_rare_filt)$percent_dead)
sample_data(survival_rare_filt)$percent_bleached = as.numeric(sample_data(survival_rare_filt)$percent_bleached)
sample_data(survival_rare_filt)$Date = factor(sample_data(survival_rare_filt)$Date,
                                         levels=desired_order)

#Make bins for mortality and bleaching data
sample_df_mortality = as.data.frame(sample_data(survival_rare_filt))

#Mortality bins: acr and plob: 0%, 0%<i<25%, and >25%; pver: none and some.
sample_df_mortality$mort_bin = NA
for (i in 1:nrow(sample_df_mortality)) {
  if(sample_df_mortality$percent_dead[i]==0){
    sample_df_mortality$mort_bin[i] = "No Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>0 & sample_df_mortality$percent_dead[i]<25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>0 & sample_df_mortality$percent_dead[i]<25 
          & sample_df_mortality$Coral[i]=="Plob"){
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>=25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$mort_bin[i] = "High Mortality"
  }
  else if(sample_df_mortality$percent_dead[i]>=25 
          & sample_df_mortality$Coral[i]=="Plob"){
    sample_df_mortality$mort_bin[i] = "High Mortality"
  }
  else {
    sample_df_mortality$mort_bin[i] = "Some Mortality"
  }
}

length(sample_df_mortality$Coral[sample_df_mortality$Coral == "Pver" 
                                 & sample_df_mortality$percent_bleached != 0])
#only 5 Pver corals bleached, no plob corals bleached, so bleaching analyses only on Acr.

sample_df_mortality$bleach_bin = NA
for (i in 1:nrow(sample_df_mortality)) {
  if(sample_df_mortality$percent_bleached[i]==0){
    sample_df_mortality$bleach_bin[i] = "No Bleaching"
  }
  else if(sample_df_mortality$percent_bleached[i]>0 & sample_df_mortality$percent_bleached[i]<25 
          & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$bleach_bin[i] = "Some Bleaching"
  }
  else if(sample_df_mortality$percent_bleached[i]>=25 & sample_df_mortality$Coral[i]=="Acr"){
    sample_df_mortality$bleach_bin[i] = "High Bleaching"
  }
  else {
    sample_df_mortality$bleach_bin[i] = "NA"
  }
}

sample_df_mortality = sample_data(sample_df_mortality)
sample_data(survival_rare_filt) = sample_df_mortality

#Save the rarefied families phyloseq object for ease of access
saveRDS(survival_rare_filt, "survival_rare_fams_bins_ps.rds")
survival_rare_filt = readRDS("survival_rare_fams_bins_ps.rds")