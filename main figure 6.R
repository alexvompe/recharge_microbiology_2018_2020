library(ggplot2)
library(ggpubr)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)

##Figure 6 panel A
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")
families = readRDS("families_ps.rds")
acr = subset_samples(families, Coral=="Acr")
plob = subset_samples(families, Coral=="Plob")
pver = subset_samples(families, Coral=="Pver")

ord_families = ordinate(families, "PCoA", "unifrac", weighted=TRUE)

ordination_df = plot_ordination(families, ord_families, type="samples", 
                                color="Coral", justDF = TRUE)

p1=ggplot(ordination_df, aes(Axis.1, Axis.2, color=Coral, fill=Coral)) +
  theme_bw()+
  geom_point(size=3, alpha=0.4)+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Coral), linewidth=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Coral), linewidth=10, level=0.001)+
  facet_grid(.~Date)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=cbPalette) + 
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [60.6%]")+
  ylab("PCoA Axis 2 [5.8%]")

##Figure 6 panel B
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
disp_date_acr = betadisper(dist_uni_acr, sampledf$Date, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_date_acr$distances,Date_disp=disp_date_acr$group, 
                    sampledf)

theme_set(theme_bw())
p2=ggplot(data = df_acr, aes(x=Date, y=Distance_to_centroid))+
  geom_boxplot(color="black", lwd=1.2)+
  geom_point(alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Weighted UniFrac Distance to Centroid")

dist_uni_plob = phyloseq::distance(plob, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(plob))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
disp_date_plob = betadisper(dist_uni_plob, sampledf$Date, type = "centroid")
df_plob = data.frame(Distance_to_centroid=disp_date_plob$distances,Date_disp=disp_date_plob$group,
                     sampledf)

theme_set(theme_bw())
p3=ggplot(data = df_plob, aes(x=Date, y=Distance_to_centroid))+
  geom_boxplot(color="#E69F00", lwd=1.2)+
  geom_point(color="#E69F00", alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Weighted UniFrac Distance to Centroid")

dist_uni_pver = phyloseq::distance(pver, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(pver))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
disp_date_pver = betadisper(dist_uni_pver, sampledf$Date, type = "centroid")
df_pver = data.frame(Distance_to_centroid=disp_date_pver$distances,Date_disp=disp_date_pver$group,
                     sampledf)

theme_set(theme_bw())
p4=ggplot(data = df_pver, aes(x=Date, y=Distance_to_centroid))+
  geom_boxplot(color="#56B4E9", lwd=1.2)+
  geom_point(color="#56B4E9", alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Weighted UniFrac Distance to Centroid")

p_1 = ggarrange(p2,p3,p4, ncol=3)

##Figure 6 panel C
ord_acr = ordinate(acr, "PCoA", "unifrac", weighted=TRUE)
ord_plob = ordinate(plob, "PCoA", "unifrac", weighted=TRUE)
ord_pver = ordinate(pver, "PCoA", "unifrac", weighted=TRUE)

theme_set(theme_bw())
ordination_acr = plot_ordination(acr, ord_acr, type="samples", 
                                 color="Coral", justDF = TRUE)
ordination_plob = plot_ordination(plob, ord_plob, type="samples", 
                                  color="Coral", justDF = TRUE)
ordination_pver = plot_ordination(pver, ord_pver, type="samples", 
                                  color="Coral", justDF = TRUE)

ordination_acr$Date = factor(ordination_acr$Date, ordered = TRUE, 
                             c("Jul18", "Nov18", "Mar19", "Aug19", "Nov19", "Mar20", "Aug20"))
#Acr stats
model = lm(Axis.1 ~ Date, data = ordination_acr)
summary(model) #adjusted R^2 = 0.1863, p = 2.404e-07

p5=ggplot(ordination_acr, aes(Date, Axis.1, color="black", fill="black", group=1)) + 
  geom_point(size=3, alpha=0.4)+
  geom_smooth(aes(as.numeric(Date), Axis.1), method = "loess")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=cbPalette) + 
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  ylab("PCoA Axis 1 [55.6%]")+
  xlab("Date")

#plob stats
model = lm(Axis.1 ~ Date, data = ordination_plob)
summary(model) #adjusted R^2 = 0.2432, p = 1.719e-08

ordination_plob$Date = factor(ordination_plob$Date, ordered = TRUE, 
                              c("Jul18", "Nov18", "Mar19", "Aug19", "Nov19", "Mar20", "Aug20"))

p6=ggplot(ordination_plob, aes(Date, Axis.1, color="#E69F00", fill="#E69F00", group=1)) + 
  geom_point(size=3, alpha=0.4)+
  geom_smooth(aes(as.numeric(Date), Axis.1), method = "loess")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#E69F00", "#E69F00")) + 
  scale_fill_manual(values=c("#E69F00", "#E69F00"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  ylab("PCoA Axis 1 [26.5%]")+
  xlab("Date")

#pver stats
model = lm(Axis.1 ~ Date, data = ordination_pver)
summary(model) #adjusted R^2 = 0.2272, p = 4.344e-12

ordination_pver$Date = factor(ordination_pver$Date, ordered = TRUE, 
                              c("Jul18", "Nov18", "Mar19", "Aug19", "Nov19", "Mar20", "Aug20"))

p7=ggplot(ordination_pver, aes(Date, Axis.1, color="#56B4E9", fill="#56B4E9", group=1)) + 
  geom_point(size=3, alpha=0.4)+
  geom_smooth(aes(as.numeric(Date), Axis.1), method = "loess")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=c("#56B4E9", "#56B4E9")) + 
  scale_fill_manual(values=c("#56B4E9", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(legend.position="none")+
  ylab("PCoA Axis 1 [66.7%]")+
  xlab("Date")

p_2 = ggarrange(p5,p6,p7, ncol=3)

##Assembling Figure 6. Additional annotations added in powerpoint.
p = ggarrange(p1,p_1,p_2, ncol=1)
ggsave(plot=p, filename="figure 6.tiff", scale=3,
       width=180,height=200,units="mm")