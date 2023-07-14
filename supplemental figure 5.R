library(phyloseq)
library(tidyverse)

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

microbes.rarefied = rarefy_even_depth(microbes.filt, 
                                      rngseed=1, 
                                      sample.size=1000, 
                                      replace=F)
##Colorblind-friendly palette
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

p = plot_richness(microbes.rarefied, x="Date", color="Coral", measures = "Shannon")+
  geom_boxplot(linewidth=1.2, alpha=0.6)+
  facet_grid(.~Coral)+
  theme_bw()+
  scale_color_manual(values=cbPalette)+
  theme(text = element_text(size = 30))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(plot=p, "supplemental figure 5.tiff", units = "mm", height = 100,
       width = 180, scale = 2)