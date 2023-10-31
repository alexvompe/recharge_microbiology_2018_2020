library(ggplot2)
library(ggpubr)
library(phyloseq)
library(plyr)
library(dplyr)

families = readRDS("families_ps.rds")

#NMDS ordination of Porites samples at Jul18 (Panel A)
plob = subset_samples(families, Coral=="Plob" & Date=="Jul18")
ord_plob = ordinate(plob, "NMDS", "bray")#stress = 0.15

ordination_df = plot_ordination(plob, ord_plob, type="samples", 
                                color="Coral", justDF = TRUE)

p1 = ggplot(ordination_df, aes(x=NMDS1, y=NMDS2))+
  geom_point(alpha = 0.6, color = "#E69F00", size = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#NMDS ordination of Pocillopora samples at Jul18 (Panel B)
pver = subset_samples(families, Coral=="Pver" & Date=="Jul18")
ord_pver = ordinate(pver, "NMDS", "bray")#stress = 0.15

ordination_df = plot_ordination(pver, ord_pver, type="samples", 
                                color="Coral", justDF = TRUE)

p2 = ggplot(ordination_df, aes(x=NMDS1, y=NMDS2))+
  geom_point(alpha = 0.6, color = "#56B4E9", size = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p = ggarrange(p1,p2, ncol = 2, align = "hv")

ggsave(plot=p, "figure s1 intermediate.tiff", units="mm", 
       height = 100, width = 180, scale=1, dpi = 700)