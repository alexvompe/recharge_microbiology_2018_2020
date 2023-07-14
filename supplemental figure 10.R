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

ordination_acr = plot_ordination(acr, ord_acr, type="samples", 
                                 color="Coral", justDF = TRUE)
ordination_plob = plot_ordination(plob, ord_plob, type="samples", 
                                  color="Coral", justDF = TRUE)
ordination_pver = plot_ordination(pver, ord_pver, type="samples", 
                                  color="Coral", justDF = TRUE)
theme_set(theme_bw())

p1=ggplot(ordination_acr, aes(Axis.1, Axis.2, color="black", fill="black")) + 
  geom_point(size=3, alpha=0.4)+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), lwd=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), lwd=10, level=0)+
  facet_grid(.~Date)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=cbPalette) + 
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [55.6%]")+
  ylab("PCoA Axis 2 [8.6%]")+
  xlim(-1.05,1.2)+
  ylim(-0.65,0.45)+
  geom_hline(yintercept=0, linetype=2, color="red", size=1, alpha=0.5)+
  geom_vline(xintercept=0, linetype=2, color="red", size=1, alpha=0.5)

p2=ggplot(ordination_plob, aes(Axis.1, Axis.2, color="#E69F00", fill="#E69F00")) + 
  geom_point(size=3, alpha=0.4)+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), lwd=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), lwd=10, level=0)+
  facet_grid(.~Date)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#E69F00", "#E69F00")) + 
  scale_fill_manual(values=c("#E69F00", "#E69F00"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [26.5%]")+
  ylab("PCoA Axis 2 [13.7%]")+
  xlim(-1.05,1.2)+
  ylim(-0.65,0.45)+
  geom_hline(yintercept=0, linetype=2, color="red", size=1, alpha=0.5)+
  geom_vline(xintercept=0, linetype=2, color="red", size=1, alpha=0.5)

p3=ggplot(ordination_pver, aes(Axis.1, Axis.2, color="#56B4E9", fill="#56B4E9")) + 
  geom_point(size=3, alpha=0.4)+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), lwd=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), lwd=10, level=0)+
  facet_grid(.~Date)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#56B4E9", "#56B4E9")) + 
  scale_fill_manual(values=c("#56B4E9", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [66.7%]")+
  ylab("PCoA Axis 2 [7.8%]")+
  xlim(-1.05,1.2)+
  ylim(-0.65,0.45)+
  geom_hline(yintercept=0, linetype=2, color="red", size=1, alpha=0.5)+
  geom_vline(xintercept=0, linetype=2, color="red", size=1, alpha=0.5)

p = ggarrange(p1,p2,p3,ncol=1)
ggsave(plot=p, filename="supplemental figure 10.tiff", scale=3,
       width=180,height=200,units="mm")