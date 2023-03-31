library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ANCOMBC)
library(plyr)
library(dplyr)

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

relative1 = transform_sample_counts(da_acr, 
                                    function(x) {x/sum(x)})

family_order_acr = as.vector(c("ASV1", "ASV13", "ASV473", "ASV453", "ASV10455",
                               "ASV584", "ASV406", "ASV593", "ASV13767", "ASV165", "ASV328",
                               "ASV30", "ASV255", "ASV65", "ASV25", "ASV139","ASV14"))

#plot
theme_set(theme_bw())
p1=plot_heatmap(relative1, taxa.label = "Family", low="white", 
                high="black", na.value = "white", taxa.order = rev(family_order_acr))+
  facet_grid(.~Date, scales="free_x") +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  guides(fill=guide_legend(title="Abundance Proportion"))+
  theme(strip.text.x = element_text(size = 20))
p1$data$Date = factor(p1$data$Date, levels = desired_order)

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

family_order_plob = as.vector(c("ASV1", "ASV13", "ASV41", "ASV327", "ASV593", "ASV230", "ASV328", "ASV4476", "ASV65", "ASV313"))

#plot
theme_set(theme_bw())
p2=plot_heatmap(relative2, taxa.label = "Family", low="white", 
                high="#E69F00", na.value = "white", taxa.order = rev(family_order_plob))+
  facet_grid(.~Date, scales="free_x") +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size = 10))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  guides(fill=guide_legend(title="Abundance Proportion"))+
  theme(strip.text.x = element_text(size = 20))
p2$data$Date = factor(p2$data$Date, levels = desired_order)

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

relative3 = transform_sample_counts(da_pver, 
                                    function(x) {x/sum(x)})

family_order_pver = as.vector(c("ASV1", "ASV5", "ASV473", "ASV93", "ASV239", "ASV368", "ASV327", 
                                "ASV584", "ASV324", "ASV86", "ASV289", "ASV235", "ASV255", "ASV61", "ASV65",
                                "ASV1051", "ASV150", "ASV379","ASV14","ASV127"))

#Plot
theme_set(theme_bw())
p3=plot_heatmap(relative3, taxa.label = "Family", low="white", 
                high="#56B4E9", na.value = "white", taxa.order = rev(family_order_pver))+
  facet_grid(.~Date, scales="free_x") +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size = 8.49))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  guides(fill=guide_legend(title="Abundance Proportion"))+
  theme(strip.text.x = element_text(size = 20))
p3$data$Date = factor(p3$data$Date, levels = desired_order)

##Make the figure
p=ggarrange(p1,p2,p3, ncol=1)
ggsave(plot=p, filename="DA heatmap.tiff", scale=2,
       width=180,height=200,units="mm")