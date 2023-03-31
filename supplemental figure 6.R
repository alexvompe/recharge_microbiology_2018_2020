library(tidyverse)
library(phyloseq)
library(ggpubr)

##Colorblind-friendly palette
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

##Load phyloseq and make dfs for each alpha metric
families_rare = readRDS("rarefied_fams_ps.rds")

df_shannon = data.frame(sample_data(families_rare), 
                        estimate_richness(families_rare, measures = "Shannon"))

df_observed = data.frame(sample_data(families_rare), 
                         estimate_richness(families_rare, measures = "Observed"))

source("estimate_richness_wPD.R")
df_faith = data.frame(sample_data(families_rare), 
                      estimate_richness_wPD(families_rare, measures = "FaithPD"))

##Make Panel A
p1 = ggplot(data = df_shannon, aes(x=Date, y=Shannon, color=Coral))+
  theme_bw()+
  facet_grid(.~Coral)+
  geom_point(alpha=0.6)+
  geom_boxplot(alpha=0.6, linewidth=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Shannon Diversity")+
  theme(legend.position="none")

##Make Panel B
p2 = ggplot(data = df_observed, aes(x=Date, y=Observed, color=Coral))+
  theme_bw()+
  facet_grid(.~Coral)+
  geom_point(alpha=0.6)+
  geom_boxplot(alpha=0.6, linewidth=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Observed Richness")+
  theme(legend.position="none")

##Make Panel C
p3 = ggplot(data = df_faith, aes(x=Date, y=FaithPD, color=Coral))+
  theme_bw()+
  facet_grid(.~Coral)+
  geom_point(alpha=0.6)+
  geom_boxplot(alpha=0.6, linewidth=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Faith's Phylogenetic Diversity")+
  theme(legend.position="none")

##Make the figure
p = ggarrange(p1,p2,p3, ncol=1)
ggsave(plot=p, "family alpha div by other metrics.tiff", units = "mm",
       width = 180, height = 200, scale = 2)