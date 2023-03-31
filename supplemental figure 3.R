library(tidyverse)
library(phyloseq)
library(ggpubr)
library(vegan)

##Colorblind-friendly palettes
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

cbPalette_herb=c("#196F3D",
                 "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
                 "#009E73", "#0072B2", "#D55E00", 
                 "#CC79A7", "#999999", "#FF468F", "#89472F", 
                 "#F0E442", "#FF4040", "#66CCCC", "#808080", 
                 "#B4CEFF")

##Make Panel A
families_rare = readRDS("rarefied_fams_ps.rds")

p1=plot_richness(families_rare, x="Herbivory", 
                measures="Shannon", color = "Coral") + 
  theme_bw()+
  facet_grid(.~Coral)+
  geom_boxplot(alpha=0.6, lwd=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Shannon Diversity Index")+
  theme(legend.position="none")

#Stats
shannon=estimate_richness(families_rare, measures = "Shannon")
sample=sample_data(families_rare)
df=data.frame(shannon,sample)
df=na.omit(df)

shapiro.test(df$Shannon)#non-normal, what about log transform?
shapiro.test(log(df$Shannon))#non-normal, proceed with pairwise WRST

pairwise.wilcox.test(df$Shannon[df$Coral=="Acr"], df$Herbivory[df$Coral=="Acr"])#no effect
pairwise.wilcox.test(df$Shannon[df$Coral=="Plob"], df$Herbivory[df$Coral=="Plob"])#no effect
pairwise.wilcox.test(df$Shannon[df$Coral=="Pver"], df$Herbivory[df$Coral=="Pver"])#no effect

##Make Panel B
families = readRDS("families_ps.rds")

acr = subset_samples(families, Coral == "Acr")
plob = subset_samples(families, Coral == "Plob")
pver = subset_samples(families, Coral == "Pver")

dist_uni_families = phyloseq::distance(families, method="unifrac", 
                                       weighted=TRUE)
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
dist_uni_plob = phyloseq::distance(plob, method="unifrac", weighted=TRUE)
dist_uni_pver = phyloseq::distance(pver, method="unifrac", weighted=TRUE)

ord_families = ordinate(families, "PCoA", "unifrac", weighted=TRUE)
ordination_families = plot_ordination(families, ord_families, type="samples", justDF = TRUE)

theme_set(theme_bw())
p2=ggplot(ordination_families, aes(Axis.1, Axis.2, fill = Herbivory, color = Herbivory))+ 
  geom_point(size=3, alpha=0.4)+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, aes(fill=Herbivory), linewidth=1.2)+
  stat_ellipse(geom = "polygon", type="euclid", 
               aes(fill=Herbivory), linewidth=10, level=0.001)+
  facet_grid(.~Date)+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_manual(values=cbPalette_herb)+ 
  scale_fill_manual(values=cbPalette_herb)+
  xlab("PCoA Axis 1 [60.6%]")+
  ylab("PCoA Axis 2 [5.8%]")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  geom_hline(yintercept=0, linetype=2, color="red", linewidth=1, alpha=0.5)+
  geom_vline(xintercept=0, linetype=2, color="red", linewidth=1, alpha=0.5)

#Stats
stat.data.r = as(sample_data(families), "data.frame")
adonis2(dist_uni_families ~ Herbivory, data = stat.data.r)#p = 0.046, R2 = 0.010
adonis2(dist_uni_families ~ Herbivory*Date, data = stat.data.r)#ns

stat.data.r = as(sample_data(acr), "data.frame")
adonis2(dist_uni_acr ~ Herbivory, data = stat.data.r)#ns

stat.data.r = as(sample_data(plob), "data.frame")
adonis2(dist_uni_plob ~ Herbivory, data = stat.data.r)#p = 0.036, R2 = 0.031

stat.data.r = as(sample_data(pver), "data.frame")
adonis2(dist_uni_pver ~ Herbivory, data = stat.data.r)#p = 0.049, R2 = 0.026

##Make Panel C
sampledf = data.frame(sample_data(acr))
sampledf$Herbivory = factor(sampledf$Herbivory, levels=c("1x1", "2x2", "3x3", "open"))
disp_herb_acr = betadisper(dist_uni_acr, sampledf$Herbivory, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_herb_acr$distances,Herbivory=disp_herb_acr$group)

sampledf = data.frame(sample_data(plob))
sampledf$Herbivory = factor(sampledf$Herbivory, levels=c("1x1", "2x2", "3x3", "open"))
disp_herb_plob = betadisper(dist_uni_plob, sampledf$Herbivory, type = "centroid")
df_plob = data.frame(Distance_to_centroid=disp_herb_plob$distances,Herbivory=disp_herb_plob$group)

sampledf = data.frame(sample_data(pver))
sampledf$Herbivory = factor(sampledf$Herbivory, levels=c("1x1", "2x2", "3x3", "open"))
disp_herb_pver = betadisper(dist_uni_pver, sampledf$Herbivory, type = "centroid")
df_pver = data.frame(Distance_to_centroid=disp_herb_pver$distances,Herbivory=disp_herb_pver$group)

p_1=ggplot(data = df_acr, aes(x=Herbivory, y=Distance_to_centroid))+
  geom_boxplot(color="black", lwd=1.2)+
  geom_point(alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylim(0,0.8)+
  ylab("Weighted UniFrac Distance to Centroid")

p_2=ggplot(data = df_plob, aes(x=Herbivory, y=Distance_to_centroid))+
  geom_boxplot(color="#E69F00", lwd=1.2)+
  geom_point(color="#E69F00", alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylim(0,0.8)+
  ylab("Weighted UniFrac Distance to Centroid")

p_3=ggplot(data = df_pver, aes(x=Herbivory, y=Distance_to_centroid))+
  geom_boxplot(color="#56B4E9", lwd=1.2)+
  geom_point(color="#56B4E9", alpha=0.6)+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylim(0,0.8)+
  ylab("Weighted UniFrac Distance to Centroid")

p3 = ggarrange(p_1, p_2, p_3, ncol = 3)

#Stats
shapiro.test(df_acr$Distance_to_centroid)#non-normal, pairwise Wilcoxon rank sum
pairwise.wilcox.test(df_acr$Distance_to_centroid, df_acr$Herbivory)

shapiro.test(df_plob$Distance_to_centroid)#non-normal, pairwise Wilcoxon rank sum
pairwise.wilcox.test(df_plob$Distance_to_centroid, df_plob$Herbivory)

shapiro.test(df_pver$Distance_to_centroid)#non-normal, pairwise Wilcoxon rank sum
pairwise.wilcox.test(df_pver$Distance_to_centroid, df_pver$Herbivory)

##Make the herbivory figure
p = ggarrange(p1,p2,p3, ncol=1)
ggsave(plot=p, "figure S3 intermediate.tiff", units = "mm", 
       height = 200, width = 180, scale = 3)