library(ggplot2)
library(ggpubr)
library(phyloseq)
library(plyr)
library(dplyr)

##NMDS ordination of Pocillopora samples at Jul18
families = readRDS("families_ps.rds")
pver = subset_samples(families, Coral=="Pver" & Date=="Jul18")
ord_pver = ordinate(pver, "NMDS", "bray")#stress = 0.19

ordination_df = plot_ordination(pver, ord_pver, type="samples", 
                                color="Coral", justDF = TRUE)

p = ggplot(ordination_df, aes(x=NMDS1, y=NMDS2))+
  geom_point(alpha = 0.6, color = "#56B4E9", size = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(plot=p, "figure s10 intermediate.tiff", units="mm", 
       height = 100, width = 100, scale=1)