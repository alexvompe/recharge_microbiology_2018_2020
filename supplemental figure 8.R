library(tidyverse)
library(phyloseq)
library(ggpubr)
library(rstatix)

##Colorblind-friendly palette
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

##Make Panel A
families_rare = readRDS("rarefied_fams_ps.rds")

#Remove Jul18, as this is when enrichment started
families_nutrients = subset_samples(families_rare, Date != "Jul18")

p1=plot_richness(families_nutrients, x="Nutrients", 
                 measures="Shannon", color = "Coral") + 
  theme_bw()+
  facet_grid(.~Coral)+
  geom_boxplot(alpha=0.6, lwd=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Shannon Diversity Index")+
  theme(legend.position="none")

#Stats
shannon=estimate_richness(families_nutrients, measures = "Shannon")
sample=sample_data(families_nutrients)
df=data.frame(shannon,sample)
df=na.omit(df)

shapiro.test(df$Shannon)#non-normal, what about log transform?
shapiro.test(log(df$Shannon))#non-normal, proceed with pairwise WRST

pairwise.wilcox.test(df$Shannon[df$Coral=="Acr"], df$Nutrients[df$Coral=="Acr"])#no effect
pairwise.wilcox.test(df$Shannon[df$Coral=="Plob"], df$Nutrients[df$Coral=="Plob"])#no effect
pairwise.wilcox.test(df$Shannon[df$Coral=="Pver"], df$Nutrients[df$Coral=="Pver"])#no effect

##Make Panel B
acr = subset_samples(families_nutrients, Coral == "Acr")
plob = subset_samples(families_nutrients, Coral == "Plob")
pver = subset_samples(families_nutrients, Coral == "Pver")

mixed_df1 = data.frame(estimate_richness(acr, measures = "Shannon"), sample_data(acr))
mixed_df2 = data.frame(estimate_richness(plob, measures = "Shannon"), sample_data(plob))
mixed_df3 = data.frame(estimate_richness(pver, measures = "Shannon"), sample_data(pver))

p_1 = ggline(mixed_df1, x = "Date", y = "Shannon", linetype = "Nutrients", shape = "Nutrients", 
             color = "black", add = "mean_ci", linewidth = 2, size = 1.2,
             ylab = "Shannon Diversity Index") + font("xy.text", size=30) +
  font("xlab", size = 30) + font("ylab", size = 30) + font ("legend.text", size = 20) + 
  font("legend.title", size = 20) + rotate_x_text()

p_2 = ggline(mixed_df2, x = "Date", y = "Shannon", linetype = "Nutrients", shape = "Nutrients",
             color = "#E69F00", add = "mean_ci", linewidth = 2, size = 1.2,
             ylab = "Shannon Diversity Index") + font("xy.text", size=30) +
  font("xlab", size = 30) + font("ylab", size = 30) + font ("legend.text", size = 20) + 
  font("legend.title", size = 20) + rotate_x_text()

p_3 = ggline(mixed_df3, x = "Date", y = "Shannon", linetype = "Nutrients", shape = "Nutrients",
             color = "#56B4E9", add = "mean_ci", linewidth = 2, size = 1.2,
             ylab = "Shannon Diversity Index") + font("xy.text", size=30) +
  font("xlab", size = 30) + font("ylab", size = 30) + font ("legend.text", size = 20) + 
  font("legend.title", size = 20) + rotate_x_text()

p2 = ggarrange(ggpar(p_1, ylim=c(0,4), legend = "none"),ggpar(p_2, ylim=c(0,4), legend = "none"),
               ggpar(p_3, ylim=c(0,4), legend = "none"), ncol = 3)

#Stats
y = estimate_richness(pver, measures = "Shannon")
x = sample_data(pver)

mixed_df = data.frame(x,y)

mixed_df_res = anova_test(data=mixed_df, Shannon ~ Herbivory + Nutrients*Date)
get_anova_table(mixed_df_res)

res = aov(Shannon ~ Nutrients*Date + Herbivory, data = mixed_df)
summary(res)
tukey = TukeyHSD(res)

tukey$`Nutrients:Date`%>%as_tibble(rownames = "comparators")%>%dplyr::filter(`p adj` < 0.05)

a=t.test(mixed_df$Shannon[mixed_df$Date=="Mar19"] ~ mixed_df$Nutrients[mixed_df$Date=="Mar19"],
         alternative = "two.sided")$p.value

b=t.test(mixed_df$Shannon[mixed_df$Date=="Mar20"] ~ mixed_df$Nutrients[mixed_df$Date=="Mar20"],
         alternative = "two.sided")$p.value

p.vals=c(a,b)
p.adjust(p.vals, method="holm")

##Make Panel C
families = readRDS("families_ps.rds")
families_nutrients = subset_samples(families, Date != "Jul18")
acr = subset_samples(families_nutrients, Coral == "Acr")
plob = subset_samples(families_nutrients, Coral == "Plob")
pver = subset_samples(families_nutrients, Coral == "Pver")

acr_df=plot_ordination(acr, ordinate(acr, method="PCoA", 
                                     distance="unifrac", 
                                     weighted=TRUE), type="samples", 
                       color="Nutrients", justDF=TRUE)

a=mean(acr_df$Axis.1[acr_df$Nutrients=="Ambient"])
b=mean(acr_df$Axis.2[acr_df$Nutrients=="Ambient"])
c=mean(acr_df$Axis.1[acr_df$Nutrients=="Nutrient"])
d=mean(acr_df$Axis.2[acr_df$Nutrients=="Nutrient"])

p_4=ggplot(acr_df, aes(x=Axis.1, y=Axis.2, fill=Nutrients, color=Nutrients))+
  geom_point(aes(color=Nutrients, shape=Nutrients), size = 2)+
  theme_bw()+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, linewidth=1.2, aes(linetype=Nutrients))+
  geom_point(aes(x=a, y=b), size=6, color="black")+
  geom_point(aes(x=c, y=d), size=6, color="black", shape=17)+
  scale_shape_manual(values = c(19, 17))+
  scale_color_manual(values = c("black", "black"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme(text = element_text(size = 30), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [54%]")+
  ylab("PCoA Axis 2 [8.9%]")

plob_df=plot_ordination(plob, ordinate(plob, method="PCoA", 
                                       distance="unifrac", 
                                       weighted=TRUE), type="samples", 
                        color="Nutrients", justDF = TRUE)

a=mean(plob_df$Axis.1[plob_df$Nutrients=="Ambient"])
b=mean(plob_df$Axis.2[plob_df$Nutrients=="Ambient"])
c=mean(plob_df$Axis.1[plob_df$Nutrients=="Nutrient"])
d=mean(plob_df$Axis.2[plob_df$Nutrients=="Nutrient"])

p_5=ggplot(plob_df, aes(x=Axis.1, y=Axis.2, fill=Nutrients, color=Nutrients))+
  geom_point(aes(color=Nutrients, shape=Nutrients), size = 2)+
  theme_bw()+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, linewidth=1.2, aes(linetype=Nutrients))+
  geom_point(aes(x=a, y=b), size=6, color="#E69F00")+
  geom_point(aes(x=c, y=d), size=6, color="#E69F00", shape=17)+
  scale_shape_manual(values = c(19, 17))+
  scale_color_manual(values = c("#E69F00", "#E69F00"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme(text = element_text(size = 30), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [21.9%]")+
  ylab("PCoA Axis 2 [15.5%]")

pver_df=plot_ordination(pver, ordinate(pver, method="PCoA", 
                                       distance="unifrac", 
                                       weighted=TRUE), type="samples", 
                        color="Nutrients", justDF = TRUE)

a=mean(pver_df$Axis.1[pver_df$Nutrients=="Ambient"])
b=mean(pver_df$Axis.2[pver_df$Nutrients=="Ambient"])
c=mean(pver_df$Axis.1[pver_df$Nutrients=="Nutrient"])
d=mean(pver_df$Axis.2[pver_df$Nutrients=="Nutrient"])

p_6=ggplot(pver_df, aes(x=Axis.1, y=Axis.2, fill=Nutrients, color=Nutrients))+
  geom_point(aes(color=Nutrients, shape=Nutrients), size = 2)+
  theme_bw()+
  stat_ellipse(geom = "polygon", type="norm", 
               alpha=0, linewidth=1.2, aes(linetype=Nutrients))+
  geom_point(aes(x=a, y=b), size=6, color="#56B4E9")+
  geom_point(aes(x=c, y=d), size=6, color="#56B4E9", shape=17)+
  scale_shape_manual(values = c(19, 17))+
  scale_color_manual(values = c("#56B4E9", "#56B4E9"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme(text = element_text(size = 30), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(legend.position="none")+
  xlab("PCoA Axis 1 [66.5%]")+
  ylab("PCoA Axis 2 [7.8%]")

p3 = ggarrange(p_4, p_5, p_6, ncol=3)

#Stats
dist_uni_acr_nut = phyloseq::distance(acr, method = "unifrac", weighted = TRUE)
sampledf = data.frame(sample_data(acr))
disp_nut_acr = betadisper(dist_uni_acr_nut, sampledf$Nutrients, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_nut_acr$distances,Nutrients=disp_nut_acr$group)

shapiro.test(df_acr$Distance_to_centroid)#non-normal, Wilcoxon rank sum 
pairwise.wilcox.test(df_acr$Distance_to_centroid, df_acr$Nutrients) #ns

dist_uni_plob_nut = phyloseq::distance(plob, method = "unifrac", weighted = TRUE)
sampledf = data.frame(sample_data(plob))
disp_nut_plob = betadisper(dist_uni_plob_nut, sampledf$Nutrients, type = "centroid")
df_plob = data.frame(Distance_to_centroid=disp_nut_plob$distances,Nutrients=disp_nut_plob$group)

shapiro.test(df_plob$Distance_to_centroid)#non-normal, Wilcoxon rank sum
pairwise.wilcox.test(df_plob$Distance_to_centroid, df_plob$Nutrients) #ns

sampledf = data.frame(sample_data(pver))
dist_uni_pver_nut = phyloseq::distance(pver, method = "unifrac", weighted = TRUE)
disp_nut_pver = betadisper(dist_uni_pver_nut, sampledf$Nutrients, type = "centroid")
df_pver = data.frame(Distance_to_centroid=disp_nut_pver$distances,Nutrients=disp_nut_pver$group)

shapiro.test(df_pver$Distance_to_centroid)#non-normal, Wilcoxon rank sum
pairwise.wilcox.test(df_pver$Distance_to_centroid, df_pver$Nutrients) #significant

stat.data.r = as(sample_data(acr), "data.frame")
adonis2(dist_uni_acr_nut ~ Nutrients, data = stat.data.r)#ns

stat.data.r = as(sample_data(plob), "data.frame")
adonis2(dist_uni_plob_nut ~ Nutrients, data = stat.data.r)#ns

stat.data.r = as(sample_data(pver), "data.frame")
adonis2(dist_uni_pver_nut ~ Nutrients, data = stat.data.r)#ns

##Complete nutrient results figure for supplemental
p = ggarrange(p1, p2, p3, ncol=1)
ggsave(plot=p, filename="figure S1 intermediate.tiff", scale=2.5,
       width=180,height=200,units="mm")