library(phyloseq)
library(tidyverse)
library(ggpubr)

families = readRDS("families_ps.rds")
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

#Split analyses by species
acr = subset_samples(families, Coral == "Acr")
plob = subset_samples(families, Coral == "Plob")
pver = subset_samples(families, Coral == "Pver")

ord_acr = ordinate(acr, "PCoA", "unifrac", weighted=TRUE)
ord_plob = ordinate(plob, "PCoA", "unifrac", weighted=TRUE)
ord_pver = ordinate(pver, "PCoA", "unifrac", weighted=TRUE)

##Source ordinations for PCoA1 regressions by date
theme_set(theme_bw())
ordination_acr = plot_ordination(acr, ord_acr, type="samples", 
                                 color="Coral")+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), linewidth=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), linewidth=10, level=0.001)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=cbPalette) + 
  scale_fill_manual(values=cbPalette)+
  theme(legend.position="none")

ordination_plob = plot_ordination(plob, ord_plob, type="samples", 
                                  color="Coral")+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), linewidth=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), linewidth=10, level=0.001)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#E69F00", "#E69F00")) + 
  scale_fill_manual(values=c("#E69F00", "#E69F00"))+
  theme(legend.position="none")

ordination_pver = plot_ordination(pver, ord_pver, type="samples", 
                                  color="Coral")+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), linewidth=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), linewidth=10, level=0.001)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#56B4E9", "#56B4E9")) + 
  scale_fill_manual(values=c("#56B4E9", "#56B4E9"))+
  theme(legend.position="none")

p = ggarrange(ordination_acr, ordination_plob, ordination_pver, ncol=3)
ggsave(plot=p, filename="supplemental figure 11.tiff", scale=4,
       width=180,height=100,units="mm")