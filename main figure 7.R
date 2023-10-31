library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ANCOMBC)
library(plyr)
library(dplyr)
library(Rmisc)

##Read in data
families = readRDS("families_ps.rds")
acr = subset_samples(families, Coral == "Acr")
plob = subset_samples(families, Coral == "Plob")
pver = subset_samples(families, Coral == "Pver")

##ANCOM-BC by Date for Aropora
out = ancombc(phyloseq = acr, formula = "Date", # test for DA by covariate of interest
              p_adj_method = "BH", zero_cut = 0.90,
              group = "Date", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE) # conserve = FALSE if you're expecting large # of DA taxa
res=out$res
res_global = out$res_global

out_df = data.frame(
  Species = row.names(res$beta),
  beta = unlist(res$beta),
  se = unlist(res$se),
  W = unlist(res$W),
  p_val = unlist(res$p_val),
  q_val = unlist(res$q_val),
  diff_abn = unlist(res$diff_abn))
fdr_out_df = out_df %>%
  dplyr::filter(q_val < 0.05) # filter significant ones
dim(fdr_out_df) # how many taxa are DA
fdr_out_df$Species
unique(fdr_out_df$Species)# gives you the DA taxa names

#Visualization
theme_set(theme_bw())
p = ggplot(data=fdr_out_df, aes(x=reorder(Species, W, sum), y=W))+
  geom_point(size=6, color="darkgreen", alpha=0.5) + xlab("Taxa")+
  # geom_hline(yintercept = 0, linetype = "dotted", size = 1)+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
write.csv(tax_table(family_data),"acropora families for DA.csv")

##Heatmap
#transform acr to relative abundance by date 
da_acr = subset_taxa(acr, Family %in% c("Endozoicomonadaceae","Oxalobacteraceae",
                                        "Burkholderiaceae","Weeksellaceae","Saccharospirillaceae","Xanthobacteraceae",
                                        "Xanthomonadaceae","Vibrionaceae","Cellvibrionaceae","Hyphomonadaceae","Rubritaleaceae",
                                        "Sphingomonadaceae","Pirellulaceae","Flavobacteriaceae","Microtrichaceae","Halieaceae",
                                        "Rhodobacteraceae"))

da_acr = subset_samples(da_acr, sample_sums(da_acr)>0)
taxa_sums(da_acr)

relative1 = transform_sample_counts(da_acr, 
                                    function(x) {x/sum(x)})
taxa_sums(relative1)

##Main figure manually computed from rel abund avgs by category by date and
##compiled into 'diff abund average figure.csv'
merge = merge_samples(da_acr, "Date")
rel = transform_sample_counts(merge, 
                              function(x) {x/sum(x)})
jul18 = subset_samples(rel, Date == 1)
df = as.data.frame(taxa_sums(jul18))
df = tibble::rownames_to_column(df, "taxa")

df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(`[taxa=="ASV1"]
other = (`taxa_sums`[taxa=="ASV13"] + `taxa_sums(jul18)`[taxa=="ASV473"]+
  `taxa_sums(jul18)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(jul18)`[taxa=="ASV10455"]+`taxa_sums(jul18)`[taxa=="ASV584"]+
  `taxa_sums(jul18)`[taxa=="ASV406"]+`taxa_sums(jul18)`[taxa=="ASV593"]+
  `taxa_sums(jul18)`[taxa=="ASV13767"]+`taxa_sums(jul18)`[taxa=="ASV165"]+
  `taxa_sums(jul18)`[taxa=="ASV328"])/7
bad = (`taxa_sums(jul18)`[taxa=="ASV30"]+`taxa_sums(jul18)`[taxa=="ASV255"]+
         `taxa_sums(jul18)`[taxa=="ASV65"]+`taxa_sums(jul18)`[taxa=="ASV25"]+
         `taxa_sums(jul18)`[taxa=="ASV139"]+`taxa_sums(jul18)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

mar19 = subset_samples(rel, Date == 2)
df = as.data.frame(taxa_sums(mar19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar19)`[taxa=="ASV1"]
other = (`taxa_sums(mar19)`[taxa=="ASV13"] + `taxa_sums(mar19)`[taxa=="ASV473"]+
           `taxa_sums(mar19)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(mar19)`[taxa=="ASV10455"]+`taxa_sums(mar19)`[taxa=="ASV584"]+
             `taxa_sums(mar19)`[taxa=="ASV406"]+`taxa_sums(mar19)`[taxa=="ASV593"]+
             `taxa_sums(mar19)`[taxa=="ASV13767"]+`taxa_sums(mar19)`[taxa=="ASV165"]+
             `taxa_sums(mar19)`[taxa=="ASV328"])/7
bad = (`taxa_sums(mar19)`[taxa=="ASV30"]+`taxa_sums(mar19)`[taxa=="ASV255"]+
         `taxa_sums(mar19)`[taxa=="ASV65"]+`taxa_sums(mar19)`[taxa=="ASV25"]+
         `taxa_sums(mar19)`[taxa=="ASV139"]+`taxa_sums(mar19)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

mar19 = subset_samples(rel, Date == 3)
df = as.data.frame(taxa_sums(mar19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar19)`[taxa=="ASV1"]
other = (`taxa_sums(mar19)`[taxa=="ASV13"] + `taxa_sums(mar19)`[taxa=="ASV473"]+
           `taxa_sums(mar19)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(mar19)`[taxa=="ASV10455"]+`taxa_sums(mar19)`[taxa=="ASV584"]+
             `taxa_sums(mar19)`[taxa=="ASV406"]+`taxa_sums(mar19)`[taxa=="ASV593"]+
             `taxa_sums(mar19)`[taxa=="ASV13767"]+`taxa_sums(mar19)`[taxa=="ASV165"]+
             `taxa_sums(mar19)`[taxa=="ASV328"])/7
bad = (`taxa_sums(mar19)`[taxa=="ASV30"]+`taxa_sums(mar19)`[taxa=="ASV255"]+
         `taxa_sums(mar19)`[taxa=="ASV65"]+`taxa_sums(mar19)`[taxa=="ASV25"]+
         `taxa_sums(mar19)`[taxa=="ASV139"]+`taxa_sums(mar19)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

aug19 = subset_samples(rel, Date == 4)
df = as.data.frame(taxa_sums(aug19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug19)`[taxa=="ASV1"]
other = (`taxa_sums(aug19)`[taxa=="ASV13"] + `taxa_sums(aug19)`[taxa=="ASV473"]+
           `taxa_sums(aug19)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(aug19)`[taxa=="ASV10455"]+`taxa_sums(aug19)`[taxa=="ASV584"]+
             `taxa_sums(aug19)`[taxa=="ASV406"]+`taxa_sums(aug19)`[taxa=="ASV593"]+
             `taxa_sums(aug19)`[taxa=="ASV13767"]+`taxa_sums(aug19)`[taxa=="ASV165"]+
             `taxa_sums(aug19)`[taxa=="ASV328"])/7
bad = (`taxa_sums(aug19)`[taxa=="ASV30"]+`taxa_sums(aug19)`[taxa=="ASV255"]+
         `taxa_sums(aug19)`[taxa=="ASV65"]+`taxa_sums(aug19)`[taxa=="ASV25"]+
         `taxa_sums(aug19)`[taxa=="ASV139"]+`taxa_sums(aug19)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

nov19 = subset_samples(rel, Date == 5)
df = as.data.frame(taxa_sums(nov19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(nov19)`[taxa=="ASV1"]
other = (`taxa_sums(nov19)`[taxa=="ASV13"] + `taxa_sums(nov19)`[taxa=="ASV473"]+
           `taxa_sums(nov19)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(nov19)`[taxa=="ASV10455"]+`taxa_sums(nov19)`[taxa=="ASV584"]+
             `taxa_sums(nov19)`[taxa=="ASV406"]+`taxa_sums(nov19)`[taxa=="ASV593"]+
             `taxa_sums(nov19)`[taxa=="ASV13767"]+`taxa_sums(nov19)`[taxa=="ASV165"]+
             `taxa_sums(nov19)`[taxa=="ASV328"])/7
bad = (`taxa_sums(nov19)`[taxa=="ASV30"]+`taxa_sums(nov19)`[taxa=="ASV255"]+
         `taxa_sums(nov19)`[taxa=="ASV65"]+`taxa_sums(nov19)`[taxa=="ASV25"]+
         `taxa_sums(nov19)`[taxa=="ASV139"]+`taxa_sums(nov19)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

mar20 = subset_samples(rel, Date == 6)
df = as.data.frame(taxa_sums(mar20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar20)`[taxa=="ASV1"]
other = (`taxa_sums(mar20)`[taxa=="ASV13"] + `taxa_sums(mar20)`[taxa=="ASV473"]+
           `taxa_sums(mar20)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(mar20)`[taxa=="ASV10455"]+`taxa_sums(mar20)`[taxa=="ASV584"]+
             `taxa_sums(mar20)`[taxa=="ASV406"]+`taxa_sums(mar20)`[taxa=="ASV593"]+
             `taxa_sums(mar20)`[taxa=="ASV13767"]+`taxa_sums(mar20)`[taxa=="ASV165"]+
             `taxa_sums(mar20)`[taxa=="ASV328"])/7
bad = (`taxa_sums(mar20)`[taxa=="ASV30"]+`taxa_sums(mar20)`[taxa=="ASV255"]+
         `taxa_sums(mar20)`[taxa=="ASV65"]+`taxa_sums(mar20)`[taxa=="ASV25"]+
         `taxa_sums(mar20)`[taxa=="ASV139"]+`taxa_sums(mar20)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

aug20 = subset_samples(rel, Date == 7)
df = as.data.frame(taxa_sums(aug20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                                     "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                                     "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug20)`[taxa=="ASV1"]
other = (`taxa_sums(aug20)`[taxa=="ASV13"] + `taxa_sums(aug20)`[taxa=="ASV473"]+
           `taxa_sums(aug20)`[taxa=="ASV453"])/3
neutral = (`taxa_sums(aug20)`[taxa=="ASV10455"]+`taxa_sums(aug20)`[taxa=="ASV584"]+
             `taxa_sums(aug20)`[taxa=="ASV406"]+`taxa_sums(aug20)`[taxa=="ASV593"]+
             `taxa_sums(aug20)`[taxa=="ASV13767"]+`taxa_sums(aug20)`[taxa=="ASV165"]+
             `taxa_sums(aug20)`[taxa=="ASV328"])/7
bad = (`taxa_sums(aug20)`[taxa=="ASV30"]+`taxa_sums(aug20)`[taxa=="ASV255"]+
         `taxa_sums(aug20)`[taxa=="ASV65"]+`taxa_sums(aug20)`[taxa=="ASV25"]+
         `taxa_sums(aug20)`[taxa=="ASV139"]+`taxa_sums(aug20)`[taxa=="ASV14"])/6
endo
other
neutral
bad
detach(df)

##ANCOM-BC by Date for Porites
out = ancombc(phyloseq = plob, formula = "Date", # test for DA by covariate of interest
              p_adj_method = "BH", zero_cut = 0.90,
              group = "Date", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE) # conserve = FALSE if you're expecting large # of DA taxa
res=out$res
res_global = out$res_global

out_df = data.frame(
  Species = row.names(res$beta),
  beta = unlist(res$beta),
  se = unlist(res$se),
  W = unlist(res$W),
  p_val = unlist(res$p_val),
  q_val = unlist(res$q_val),
  diff_abn = unlist(res$diff_abn))
fdr_out_df = out_df %>%
  dplyr::filter(q_val < 0.05) # filter significant ones
dim(fdr_out_df) # how many taxa are DA
fdr_out_df$Species
unique(fdr_out_df$Species)# gives you the DA taxa names

#Visualization
theme_set(theme_bw())
p = ggplot(data=fdr_out_df, aes(x=reorder(Species, W, sum), y=W))+
  geom_point(size=6, color="darkgreen", alpha=0.5) + xlab("Taxa")+
  # geom_hline(yintercept = 0, linetype = "dotted", size = 1)+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
write.csv(tax_table(family_data),"porites families for DA.csv")

##Heatmap
#transform plob to relative abundance by date 
da_plob = subset_taxa(plob, Family %in% c("Endozoicomonadaceae","Saccharospirillaceae",
                                          "Oxalobacteraceae","Vibrionaceae","Puniceicoccaceae",
                                          "Crocinitomicaceae","Comamonadaceae","Microtrichaceae",
                                          "Moraxellaceae","Ilumatobacteraceae"))

da_plob = subset_samples(da_plob, sample_sums(da_plob)>0)

relative2 = transform_sample_counts(da_plob, 
                                    function(x) {x/sum(x)})

#Main text summary figure:
merge = merge_samples(da_plob, "Date")
rel = transform_sample_counts(merge, 
                              function(x) {x/sum(x)})
jul18 = subset_samples(rel, Date == 1)
df = as.data.frame(taxa_sums(jul18))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(jul18)`[taxa=="ASV1"]
other = (`taxa_sums(jul18)`[taxa=="ASV13"] + `taxa_sums(jul18)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(jul18)`[taxa=="ASV327"]+`taxa_sums(jul18)`[taxa=="ASV593"]+
             `taxa_sums(jul18)`[taxa=="ASV230"]+`taxa_sums(jul18)`[taxa=="ASV328"]+
             `taxa_sums(jul18)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(jul18)`[taxa=="ASV65"]+`taxa_sums(jul18)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

mar19 = subset_samples(rel, Date == 2)
df = as.data.frame(taxa_sums(mar19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar19)`[taxa=="ASV1"]
other = (`taxa_sums(mar19)`[taxa=="ASV13"] + `taxa_sums(mar19)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(mar19)`[taxa=="ASV327"]+`taxa_sums(mar19)`[taxa=="ASV593"]+
             `taxa_sums(mar19)`[taxa=="ASV230"]+`taxa_sums(mar19)`[taxa=="ASV328"]+
             `taxa_sums(mar19)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(mar19)`[taxa=="ASV65"]+`taxa_sums(mar19)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

mar19 = subset_samples(rel, Date == 3)
df = as.data.frame(taxa_sums(mar19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar19)`[taxa=="ASV1"]
other = (`taxa_sums(mar19)`[taxa=="ASV13"] + `taxa_sums(mar19)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(mar19)`[taxa=="ASV327"]+`taxa_sums(mar19)`[taxa=="ASV593"]+
             `taxa_sums(mar19)`[taxa=="ASV230"]+`taxa_sums(mar19)`[taxa=="ASV328"]+
             `taxa_sums(mar19)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(mar19)`[taxa=="ASV65"]+`taxa_sums(mar19)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

aug19 = subset_samples(rel, Date == 4)
df = as.data.frame(taxa_sums(aug19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug19)`[taxa=="ASV1"]
other = (`taxa_sums(aug19)`[taxa=="ASV13"] + `taxa_sums(aug19)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(aug19)`[taxa=="ASV327"]+`taxa_sums(aug19)`[taxa=="ASV593"]+
             `taxa_sums(aug19)`[taxa=="ASV230"]+`taxa_sums(aug19)`[taxa=="ASV328"]+
             `taxa_sums(aug19)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(aug19)`[taxa=="ASV65"]+`taxa_sums(aug19)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

nov19 = subset_samples(rel, Date == 5)
df = as.data.frame(taxa_sums(nov19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(nov19)`[taxa=="ASV1"]
other = (`taxa_sums(nov19)`[taxa=="ASV13"] + `taxa_sums(nov19)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(nov19)`[taxa=="ASV327"]+`taxa_sums(nov19)`[taxa=="ASV593"]+
             `taxa_sums(nov19)`[taxa=="ASV230"]+`taxa_sums(nov19)`[taxa=="ASV328"]+
             `taxa_sums(nov19)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(nov19)`[taxa=="ASV65"]+`taxa_sums(nov19)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

mar20 = subset_samples(rel, Date == 6)
df = as.data.frame(taxa_sums(mar20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar20)`[taxa=="ASV1"]
other = (`taxa_sums(mar20)`[taxa=="ASV13"] + `taxa_sums(mar20)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(mar20)`[taxa=="ASV327"]+`taxa_sums(mar20)`[taxa=="ASV593"]+
             `taxa_sums(mar20)`[taxa=="ASV230"]+`taxa_sums(mar20)`[taxa=="ASV328"]+
             `taxa_sums(mar20)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(mar20)`[taxa=="ASV65"]+`taxa_sums(mar20)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

aug20 = subset_samples(rel, Date == 7)
df = as.data.frame(taxa_sums(aug20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV13", "ASV41", "ASV327", 
                                     "ASV593", "ASV230", "ASV328", "ASV4476", 
                                     "ASV65", "ASV313"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug20)`[taxa=="ASV1"]
other = (`taxa_sums(aug20)`[taxa=="ASV13"] + `taxa_sums(aug20)`[taxa=="ASV41"])/2
neutral = (`taxa_sums(aug20)`[taxa=="ASV327"]+`taxa_sums(aug20)`[taxa=="ASV593"]+
             `taxa_sums(aug20)`[taxa=="ASV230"]+`taxa_sums(aug20)`[taxa=="ASV328"]+
             `taxa_sums(aug20)`[taxa=="ASV4476"])/5
bad = (`taxa_sums(aug20)`[taxa=="ASV65"]+`taxa_sums(aug20)`[taxa=="ASV13"])/2
endo
other
neutral
bad
detach(df)

##ANCOM-BC by Date for Pocillopora
out = ancombc(phyloseq = pver, formula = "Date", # test for DA by covariate of interest
              p_adj_method = "BH", zero_cut = 0.90,
              group = "Date", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE) # conserve = FALSE if you're expecting large # of DA taxa
res=out$res
res_global = out$res_global

out_df = data.frame(
  Species = row.names(res$beta),
  beta = unlist(res$beta),
  se = unlist(res$se),
  W = unlist(res$W),
  p_val = unlist(res$p_val),
  q_val = unlist(res$q_val),
  diff_abn = unlist(res$diff_abn))
fdr_out_df = out_df %>%
  dplyr::filter(q_val < 0.05) # filter significant ones
dim(fdr_out_df) # how many taxa are DA
fdr_out_df$Species
unique(fdr_out_df$Species)# gives you the DA taxa names

#Visualization
theme_set(theme_bw())
p = ggplot(data=fdr_out_df, aes(x=reorder(Species, W, sum), y=W))+
  geom_point(size=6, color="darkgreen", alpha=0.5) + xlab("Taxa")+
  # geom_hline(yintercept = 0, linetype = "dotted", size = 1)+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
write.csv(tax_table(family_data),"pocillopora families for DA.csv")

##Heatmap
#transform pver to relative abundance by date 
da_pver = subset_taxa(pver, Family %in% c("NS9 marine group","Endozoicomonadaceae",
                                          "Cryomorphaceae","Cellvibrionaceae",
                                          "AEGEAN-169 marine group","Woeseiaceae","Halieaceae",
                                          "Saprospiraceae","Pseudoalteromonadaceae",
                                          "Alteromonadaceae","Nostocaceae","SAR116 clade",
                                          "SAR86 clade","Phormidiaceae","Vibrionaceae",
                                          "Moraxellaceae",
                                          "Cyanobiaceae","Amoebophilaceae","Flavobacteriaceae",
                                          "Rhodobacteraceae"))

da_pver = subset_samples(da_pver, sample_sums(da_pver)>0)

#Main text average relative abundance figure:
merge = merge_samples(da_pver, "Date")
rel = transform_sample_counts(merge, 
                              function(x) {x/sum(x)})
jul18 = subset_samples(rel, Date == 1)
df = as.data.frame(taxa_sums(jul18))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(jul18)`[taxa=="ASV1"]
other = (`taxa_sums(jul18)`[taxa=="ASV5"] + `taxa_sums(jul18)`[taxa=="ASV473"]+
           `taxa_sums(jul18)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(jul18)`[taxa=="ASV239"]+`taxa_sums(jul18)`[taxa=="ASV368"]+
             `taxa_sums(jul18)`[taxa=="ASV327"]+`taxa_sums(jul18)`[taxa=="ASV584"]+
             `taxa_sums(jul18)`[taxa=="ASV324"]+`taxa_sums(jul18)`[taxa=="ASV86"]+
             `taxa_sums(jul18)`[taxa=="ASV289"]+`taxa_sums(jul18)`[taxa=="ASV235"])/8
bad = (`taxa_sums(jul18)`[taxa=="ASV255"]+`taxa_sums(jul18)`[taxa=="ASV61"]+
         `taxa_sums(jul18)`[taxa=="ASV65"]+`taxa_sums(jul18)`[taxa=="ASV1051"]+
         `taxa_sums(jul18)`[taxa=="ASV150"]+`taxa_sums(jul18)`[taxa=="ASV379"]+
         `taxa_sums(jul18)`[taxa=="ASV14"]+`taxa_sums(jul18)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

merge = merge_samples(da_pver, "Date")
rel = transform_sample_counts(merge, 
                              function(x) {x/sum(x)})
nov18 = subset_samples(rel, Date == 2)
df = as.data.frame(taxa_sums(nov18))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(nov18)`[taxa=="ASV1"]
other = (`taxa_sums(nov18)`[taxa=="ASV5"] + `taxa_sums(nov18)`[taxa=="ASV473"]+
           `taxa_sums(nov18)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(nov18)`[taxa=="ASV239"]+`taxa_sums(nov18)`[taxa=="ASV368"]+
             `taxa_sums(nov18)`[taxa=="ASV327"]+`taxa_sums(nov18)`[taxa=="ASV584"]+
             `taxa_sums(nov18)`[taxa=="ASV324"]+`taxa_sums(nov18)`[taxa=="ASV86"]+
             `taxa_sums(nov18)`[taxa=="ASV289"]+`taxa_sums(nov18)`[taxa=="ASV235"])/8
bad = (`taxa_sums(nov18)`[taxa=="ASV255"]+`taxa_sums(nov18)`[taxa=="ASV61"]+
         `taxa_sums(nov18)`[taxa=="ASV65"]+`taxa_sums(nov18)`[taxa=="ASV1051"]+
         `taxa_sums(nov18)`[taxa=="ASV150"]+`taxa_sums(nov18)`[taxa=="ASV379"]+
         `taxa_sums(nov18)`[taxa=="ASV14"]+`taxa_sums(nov18)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

mar19 = subset_samples(rel, Date == 3)
df = as.data.frame(taxa_sums(mar19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar19)`[taxa=="ASV1"]
other = (`taxa_sums(mar19)`[taxa=="ASV5"] + `taxa_sums(mar19)`[taxa=="ASV473"]+
           `taxa_sums(mar19)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(mar19)`[taxa=="ASV239"]+`taxa_sums(mar19)`[taxa=="ASV368"]+
             `taxa_sums(mar19)`[taxa=="ASV327"]+`taxa_sums(mar19)`[taxa=="ASV584"]+
             `taxa_sums(mar19)`[taxa=="ASV324"]+`taxa_sums(mar19)`[taxa=="ASV86"]+
             `taxa_sums(mar19)`[taxa=="ASV289"]+`taxa_sums(mar19)`[taxa=="ASV235"])/8
bad = (`taxa_sums(mar19)`[taxa=="ASV255"]+`taxa_sums(mar19)`[taxa=="ASV61"]+
         `taxa_sums(mar19)`[taxa=="ASV65"]+`taxa_sums(mar19)`[taxa=="ASV1051"]+
         `taxa_sums(mar19)`[taxa=="ASV150"]+`taxa_sums(mar19)`[taxa=="ASV379"]+
         `taxa_sums(mar19)`[taxa=="ASV14"]+`taxa_sums(mar19)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

aug19 = subset_samples(rel, Date == 4)
df = as.data.frame(taxa_sums(aug19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug19)`[taxa=="ASV1"]
other = (`taxa_sums(aug19)`[taxa=="ASV5"] + `taxa_sums(aug19)`[taxa=="ASV473"]+
           `taxa_sums(aug19)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(aug19)`[taxa=="ASV239"]+`taxa_sums(aug19)`[taxa=="ASV368"]+
             `taxa_sums(aug19)`[taxa=="ASV327"]+`taxa_sums(aug19)`[taxa=="ASV584"]+
             `taxa_sums(aug19)`[taxa=="ASV324"]+`taxa_sums(aug19)`[taxa=="ASV86"]+
             `taxa_sums(aug19)`[taxa=="ASV289"]+`taxa_sums(aug19)`[taxa=="ASV235"])/8
bad = (`taxa_sums(aug19)`[taxa=="ASV255"]+`taxa_sums(aug19)`[taxa=="ASV61"]+
         `taxa_sums(aug19)`[taxa=="ASV65"]+`taxa_sums(aug19)`[taxa=="ASV1051"]+
         `taxa_sums(aug19)`[taxa=="ASV150"]+`taxa_sums(aug19)`[taxa=="ASV379"]+
         `taxa_sums(aug19)`[taxa=="ASV14"]+`taxa_sums(aug19)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

nov19 = subset_samples(rel, Date == 5)
df = as.data.frame(taxa_sums(nov19))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(nov19)`[taxa=="ASV1"]
other = (`taxa_sums(nov19)`[taxa=="ASV5"] + `taxa_sums(nov19)`[taxa=="ASV473"]+
           `taxa_sums(nov19)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(nov19)`[taxa=="ASV239"]+`taxa_sums(nov19)`[taxa=="ASV368"]+
             `taxa_sums(nov19)`[taxa=="ASV327"]+`taxa_sums(nov19)`[taxa=="ASV584"]+
             `taxa_sums(nov19)`[taxa=="ASV324"]+`taxa_sums(nov19)`[taxa=="ASV86"]+
             `taxa_sums(nov19)`[taxa=="ASV289"]+`taxa_sums(nov19)`[taxa=="ASV235"])/8
bad = (`taxa_sums(nov19)`[taxa=="ASV255"]+`taxa_sums(nov19)`[taxa=="ASV61"]+
         `taxa_sums(nov19)`[taxa=="ASV65"]+`taxa_sums(nov19)`[taxa=="ASV1051"]+
         `taxa_sums(nov19)`[taxa=="ASV150"]+`taxa_sums(nov19)`[taxa=="ASV379"]+
         `taxa_sums(nov19)`[taxa=="ASV14"]+`taxa_sums(nov19)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

mar20 = subset_samples(rel, Date == 6)
df = as.data.frame(taxa_sums(mar20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(mar20)`[taxa=="ASV1"]
other = (`taxa_sums(mar20)`[taxa=="ASV5"] + `taxa_sums(mar20)`[taxa=="ASV473"]+
           `taxa_sums(mar20)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(mar20)`[taxa=="ASV239"]+`taxa_sums(mar20)`[taxa=="ASV368"]+
             `taxa_sums(mar20)`[taxa=="ASV327"]+`taxa_sums(mar20)`[taxa=="ASV584"]+
             `taxa_sums(mar20)`[taxa=="ASV324"]+`taxa_sums(mar20)`[taxa=="ASV86"]+
             `taxa_sums(mar20)`[taxa=="ASV289"]+`taxa_sums(mar20)`[taxa=="ASV235"])/8
bad = (`taxa_sums(mar20)`[taxa=="ASV255"]+`taxa_sums(mar20)`[taxa=="ASV61"]+
         `taxa_sums(mar20)`[taxa=="ASV65"]+`taxa_sums(mar20)`[taxa=="ASV1051"]+
         `taxa_sums(mar20)`[taxa=="ASV150"]+`taxa_sums(mar20)`[taxa=="ASV379"]+
         `taxa_sums(mar20)`[taxa=="ASV14"]+`taxa_sums(mar20)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

aug20 = subset_samples(rel, Date == 7)
df = as.data.frame(taxa_sums(aug20))
df = tibble::rownames_to_column(df, "taxa")
df$taxa = factor(df$taxa, levels = c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                     "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                     "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))
df = arrange(df, taxa)
attach(df)
endo = `taxa_sums(aug20)`[taxa=="ASV1"]
other = (`taxa_sums(aug20)`[taxa=="ASV5"] + `taxa_sums(aug20)`[taxa=="ASV473"]+
           `taxa_sums(aug20)`[taxa=="ASV93"])/3
neutral = (`taxa_sums(aug20)`[taxa=="ASV239"]+`taxa_sums(aug20)`[taxa=="ASV368"]+
             `taxa_sums(aug20)`[taxa=="ASV327"]+`taxa_sums(aug20)`[taxa=="ASV584"]+
             `taxa_sums(aug20)`[taxa=="ASV324"]+`taxa_sums(aug20)`[taxa=="ASV86"]+
             `taxa_sums(aug20)`[taxa=="ASV289"]+`taxa_sums(aug20)`[taxa=="ASV235"])/8
bad = (`taxa_sums(aug20)`[taxa=="ASV255"]+`taxa_sums(aug20)`[taxa=="ASV61"]+
         `taxa_sums(aug20)`[taxa=="ASV65"]+`taxa_sums(aug20)`[taxa=="ASV1051"]+
         `taxa_sums(aug20)`[taxa=="ASV150"]+`taxa_sums(aug20)`[taxa=="ASV379"]+
         `taxa_sums(aug20)`[taxa=="ASV14"]+`taxa_sums(aug20)`[taxa=="ASV127"])/8
endo
other
neutral
bad
detach(df)

##Main Text averages figure 7:
df = read.csv("diff abund average figure.csv", header = TRUE)
df = subset(df, df$Classification!="Neutral/Unknown")
head(df)

df$Date = trimws(df$Date, "left")
df$Date = factor(df$Date, levels = c("Jul18","Nov18","Mar19","Aug19","Nov19",
                                     "Mar20","Aug20"))
df$Classification = factor(df$Classification, levels = c("Detrimental",
                                                         "Other Beneficial",
                                                         "Endozoicomonadaceae"))
acr = subset(df, df$Coral=="Acr")
plob = subset(df, df$Coral=="Plob")
pver = subset(df, df$Coral=="Pver")

p1 = ggplot(acr, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p1_alt = ggplot(acr, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "black")+
  theme(legend.position = "none")

p2 = ggplot(plob, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p2_alt = ggplot(plob, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "#E69F00")+
  theme(legend.position = "none")

p3 = ggplot(pver, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p3_alt = ggplot(pver, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "#56B4E9")+
  theme(legend.position = "none")

##ratios over time with mortality
survival = readRDS("survival_fams_bins_ps.rds")
df_survival = data.frame(sample_data(survival))
acr_mort = subset(df_survival, df_survival$Coral=="Acr")
acr_mort = summarySEwithin(acr_mort, measurevar="percent_dead", withinvars="Date",
                           na.rm=FALSE, conf.interval=.95)

p4 = ggplot(acr, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="black", shape = 17)+
  geom_line(color = "black")+
  theme_bw()+
  geom_point(data=acr_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=acr_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

plob_mort = subset(df_survival, df_survival$Coral=="Plob")
plob_mort = summarySEwithin(plob_mort, measurevar="percent_dead", withinvars="Date",
                           na.rm=FALSE, conf.interval=.95)

p5 = ggplot(plob, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="#E69F00", shape = 17)+
  geom_line(color = "#E69F00")+
  geom_point(data=plob_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=plob_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

pver_mort = subset(df_survival, df_survival$Coral=="Pver")
pver_mort = summarySEwithin(pver_mort, measurevar="percent_dead", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

p6 = ggplot(pver, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="#56B4E9", shape = 17)+
  geom_line(color = "#56B4E9")+
  geom_point(data=pver_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=pver_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

p = ggarrange(p4, p1, p5, p2, p6, p3, ncol=1, nrow=6, align = "v",
              heights = c(1,0.85,1,0.85,1,0.85))
ggsave(plot = p, "workshop figure 7 legend_experiment.tiff", units="mm",
       scale = 0.9, height = 200, width = 180, dpi = 1000)

##Stats
#Acr
unique(log(acr$Ratio))
ratios = c(7.907402,7.906351,3.869733,3.618751,4.779950,5.003201)
cor.test(ratios, acr_mort$percent_dead, method = "pearson")

#Plob
unique(log(plob$Ratio))
ratios = c(4.7945142,-0.6619048,-0.6186997,0.5672443,1.6807686,1.2601039)
cor.test(ratios, plob_mort$percent_dead, method = "pearson")

#Pver
unique(log(pver$Ratio))
ratios = c(2.170288,5.065308,7.352671)
pver_mort_mid = subset(pver_mort, pver_mort$Date==c("Aug19","Nov19","Aug20"))
cor.test(ratios, pver_mort_mid$percent_dead, method = "pearson")