library(tidyverse)
library(vegan)
library(ggpubr)
library(phyloseq)

##Color gradient of mortality and bleaching for PCoA1 by Date regression
survival_filt = readRDS("survival_fams_bins_ps.rds")
acr = subset_samples(survival_filt, Coral == "Acr")
plob = subset_samples(survival_filt, Coral == "Plob")
pver = subset_samples(survival_filt, Coral == "Pver")

ord_acr = ordinate(acr, method = "PCoA", distance = "unifrac", 
                   weighted = TRUE)
ord_plob = ordinate(plob, method = "PCoA", distance = "unifrac", 
                    weighted = TRUE)
ord_acr_df = plot_ordination(acr, ord_acr, type="samples", 
                             color="Date", justDF = TRUE)
#59.2 and 10.1 % variance explained by 1st 2 axes

ord_plob_df = plot_ordination(plob, ord_plob, type="samples", 
                             color="Date", justDF = TRUE)
#26.5 and 14.5 % variance explained by 1st 2 axes

#Dispersion differences by binary mortality and bleaching####MOVE TO SUPP
p2=ggplot(ord_acr_df, aes(x=Axis.1, y=Axis.2, color=percent_dead>0,
                          linetype=percent_dead>0)) + 
  theme_bw()+
  geom_point(size=3, alpha = 0.6)+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_acr_df$percent_dead==0]), 
                 y=mean(ord_acr_df$Axis.2[ord_acr_df$percent_dead==0])), size=8,
             color="black")+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_acr_df$percent_dead>0]), 
                 y=mean(ord_acr_df$Axis.2[ord_acr_df$percent_dead>0])), size=8,
             color="darkred")+
  stat_ellipse(type="norm", linewidth=2)+
  scale_color_manual(values=c("black","darkred"), "Mortality")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  xlab("PCoA Axis 1 [59.2%]")+
  ylab("PCoA Axis 2 [10.1%]")+
  theme(legend.position="none")

p1=ggplot(ord_acr_df, aes(x=Axis.1, y=Axis.2, color=percent_bleached>0,
                          linetype=percent_bleached>0)) + 
  theme_bw()+
  geom_point(size=3, alpha = 0.6)+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_acr_df$percent_bleached==0]), 
                 y=mean(ord_acr_df$Axis.2[ord_acr_df$percent_bleached==0])), size=8,
             color="black")+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_acr_df$percent_bleached>0]), 
                 y=mean(ord_acr_df$Axis.2[ord_acr_df$percent_bleached>0])), size=8,
             color="darkred")+
  stat_ellipse(type="norm", linewidth=2)+
  scale_color_manual(values=c("black","darkred"), "Bleaching")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  xlab("PCoA Axis 1 [59.2%]")+
  ylab("PCoA Axis 2 [10.1%]")+
  theme(legend.position="none")

##Keep this one in main text Figure 8
p5=ggplot(ord_plob_df, aes(x=Axis.1, y=Axis.2, color=percent_dead>0,
                          linetype=percent_dead>0)) + 
  theme_bw()+
  geom_point(alpha = 0.6)+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_plob_df$percent_dead==0]), 
                 y=mean(ord_acr_df$Axis.2[ord_plob_df$percent_dead==0])), size=5,
             color="black")+
  geom_point(aes(x=mean(ord_acr_df$Axis.1[ord_plob_df$percent_dead>0]), 
                 y=mean(ord_acr_df$Axis.2[ord_plob_df$percent_dead>0])), size=5,
             color="darkred")+
  stat_ellipse(type="norm")+
  scale_color_manual(values=c("black","darkred"), "Mortality")+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("PCoA Axis 1 [26.5%]")+
  ylab("PCoA Axis 2 [14.5%]")+
  theme(legend.position="none")

##Change to dual axis plots
##new
survival_rare_filt = readRDS("survival_rare_fams_bins_ps.rds")
acr_shannon = subset_samples(survival_rare_filt, Coral=="Acr")
df_acr_shannon = data.frame(sample_data(acr_shannon),
                            estimate_richness(acr_shannon,
                                              measures = "Shannon"))

p1=ggplot(df_acr_shannon, aes(x=Date, y=Shannon))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="black")+
  stat_summary(fun=mean, geom="line", color="black", aes(group = 1))+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", 
                     sec.axis = sec_axis(trans = ~.*8, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_acr = betadisper(dist_uni_acr, sampledf$Date, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_date_acr$distances,Date_disp=disp_date_acr$group, 
                    sampledf)
p2=ggplot(df_acr, aes(x=Date, y=Distance_to_centroid))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="black", shape = 17)+
  stat_summary(fun=mean, geom="line", linetype = "dotdash",
               color="black", aes(group = 1))+
  stat_summary(data=df_acr, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Dispersion", 
                     sec.axis = sec_axis(trans = ~.*50, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p = ggarrange(p1,p2, ncol = 1, align = "hv")

#Pver dispersion by mortality prevalence panel
dist_uni_pver = phyloseq::distance(pver, "unifrac", weighted=TRUE)
stat.data.r = as(sample_data(pver), "data.frame")

##check mortality prevalence
stat.data.r$binary_mortality = 0
for (i in 1:nrow(stat.data.r)) {
  if(stat.data.r$percent_dead[i]>0){
    stat.data.r$binary_mortality[i] = 1
  }
  else {
    stat.data.r$binary_mortality[i] = 0
  }
}

stat.data.r$binary_mortality = factor(stat.data.r$binary_mortality)

#Dispersion
disp_mort_pver = betadisper(dist_uni_pver, stat.data.r$binary_mortality, type = "centroid")
stat.data.r$mort_disp = disp_mort_pver$distances

mar=subset(stat.data.r, Date=="Mar19")
aug=subset(stat.data.r, Date=="Aug20")

p6 = ggplot(aug, aes(x=binary_mortality, y=log(mort_disp)))+
  theme_bw()+
  geom_point(alpha = 0.8, color = "#56B4E9")+
  geom_boxplot(color = "#56B4E9", alpha = 0)+
  labs(x="August 2020 Mortality Prevalence", y="log(Dispersion)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
length(mar$binary_mortality[mar$binary_mortality==1]) #only one coral with any mortality in Mar19! All trends are for Aug20
#stats
pairwise.wilcox.test(log(aug$mort_disp), aug$binary_mortality)

p = ggarrange(p1, p2, ncol=1, align = "hv")
ggsave(plot=p, "main figure 8 experimental.tiff",height = 100,
       width = 180, units = "mm", scale=1, dpi = 700)

##Make above dual axis plots for each coral species
survival_rare_filt = readRDS("survival_rare_fams_bins_ps.rds")
acr_shannon = subset_samples(survival_rare_filt, Coral=="Acr")
plob_shannon = subset_samples(survival_rare_filt, Coral=="Plob")
pver_shannon = subset_samples(survival_rare_filt, Coral=="Pver")

df_acr_shannon = data.frame(sample_data(acr_shannon),
                            estimate_richness(acr_shannon,
                                              measures = "Shannon"))
df_plob_shannon = data.frame(sample_data(plob_shannon),
                            estimate_richness(plob_shannon,
                                              measures = "Shannon"))
df_pver_shannon = data.frame(sample_data(pver_shannon),
                             estimate_richness(pver_shannon,
                                               measures = "Shannon"))
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_acr = betadisper(dist_uni_acr, sampledf$Date, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_date_acr$distances,Date_disp=disp_date_acr$group, 
                    sampledf)

dist_uni_plob = phyloseq::distance(plob, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(plob))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_plob = betadisper(dist_uni_plob, sampledf$Date, type = "centroid")
df_plob = data.frame(Distance_to_centroid=disp_date_plob$distances,Date_disp=disp_date_plob$group, 
                    sampledf)

dist_uni_pver = phyloseq::distance(pver, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(pver))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_pver = betadisper(dist_uni_pver, sampledf$Date, type = "centroid")
df_pver = data.frame(Distance_to_centroid=disp_date_pver$distances,Date_disp=disp_date_pver$group, 
                    sampledf)

p1=ggplot(df_acr_shannon, aes(x=Date, y=Shannon))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="black")+
  stat_summary(fun=mean, geom="line", color="black", aes(group = 1))+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_acr_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", 
                     sec.axis = sec_axis(trans = ~.*8, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p2=ggplot(df_acr, aes(x=Date, y=Distance_to_centroid))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="black", shape = 17)+
  stat_summary(fun=mean, geom="line", linetype = "dotdash",
               color="black", aes(group = 1))+
  stat_summary(data=df_acr, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_acr, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Dispersion", 
                     sec.axis = sec_axis(trans = ~.*50, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p3=ggplot(df_plob_shannon, aes(x=Date, y=Shannon))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="#E69F00")+
  stat_summary(fun=mean, geom="line", color="#E69F00", aes(group = 1))+
  stat_summary(data=df_plob_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_plob_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_plob_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_plob_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", 
                     sec.axis = sec_axis(trans = ~.*8, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p4=ggplot(df_plob, aes(x=Date, y=Distance_to_centroid))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="#E69F00", shape = 17)+
  stat_summary(fun=mean, geom="line", linetype = "dotdash",
               color="#E69F00", aes(group = 1))+
  stat_summary(data=df_plob, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_plob, aes(x=Date, y=percent_bleached/50, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_plob, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_plob, aes(x=Date, y=percent_dead/50, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Dispersion", 
                     sec.axis = sec_axis(trans = ~.*50, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p5=ggplot(df_pver_shannon, aes(x=Date, y=Shannon))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="#56B4E9")+
  stat_summary(fun=mean, geom="line", color="#56B4E9", aes(group = 1))+
  stat_summary(data=df_pver_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_pver_shannon, aes(x=Date, y=percent_bleached/8, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_pver_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_pver_shannon, aes(x=Date, y=percent_dead/8, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", 
                     sec.axis = sec_axis(trans = ~.*8, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p6=ggplot(df_pver, aes(x=Date, y=Distance_to_centroid))+
  theme_bw()+
  stat_summary(fun=mean, geom="point", color="#56B4E9", shape = 17)+
  stat_summary(fun=mean, geom="line", linetype = "dotdash",
               color="#56B4E9", aes(group = 1))+
  stat_summary(data=df_pver, aes(x=Date, y=percent_bleached/25, group=1),fun=mean, 
               geom="point", color="darkred")+
  stat_summary(data=df_pver, aes(x=Date, y=percent_bleached/25, group=1),fun=mean, 
               geom="line", color="darkred")+
  stat_summary(data=df_pver, aes(x=Date, y=percent_dead/25, group=1),fun=mean, 
               geom="point", shape = 4, color="darkred")+
  stat_summary(data=df_pver, aes(x=Date, y=percent_dead/25, group=1),fun=mean, 
               geom="line", linetype = 2, color="darkred")+
  scale_y_continuous(name = "Dispersion", 
                     sec.axis = sec_axis(trans = ~.*25, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p = ggarrange(p1,p2,p3,p4,p5,p6, ncol = 1, align = "hv")
ggsave(plot = p, "recent experimental figure 8.tiff", units = "mm", width = 180,
       height = 300, dpi = 1000, scale = 1)

##Try to add the 3rd y axis with code:
library(cowplot)
library(patchwork)
library(Rmisc)

survival_rare_filt = readRDS("survival_rare_fams_bins_ps.rds")
acr_shannon = subset_samples(survival_rare_filt, Coral=="Acr")
plob_shannon = subset_samples(survival_rare_filt, Coral=="Plob")
pver_shannon = subset_samples(survival_rare_filt, Coral=="Pver")

survival_filt = readRDS("survival_fams_bins_ps.rds")
acr = subset_samples(survival_filt, Coral == "Acr")
plob = subset_samples(survival_filt, Coral == "Plob")
pver = subset_samples(survival_filt, Coral == "Pver")

df_acr_shannon = data.frame(sample_data(acr_shannon),
                            estimate_richness(acr_shannon,
                                              measures = "Shannon"))
df_plob_shannon = data.frame(sample_data(plob_shannon),
                             estimate_richness(plob_shannon,
                                               measures = "Shannon"))
df_pver_shannon = data.frame(sample_data(pver_shannon),
                             estimate_richness(pver_shannon,
                                               measures = "Shannon"))
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_acr = betadisper(dist_uni_acr, sampledf$Date, type = "centroid")
df_acr_dispersion = data.frame(Distance_to_centroid=disp_date_acr$distances,Date_disp=disp_date_acr$group, 
                    sampledf)

dist_uni_plob = phyloseq::distance(plob, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(plob))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_plob = betadisper(dist_uni_plob, sampledf$Date, type = "centroid")
df_plob_dispersion = data.frame(Distance_to_centroid=disp_date_plob$distances,Date_disp=disp_date_plob$group, 
                     sampledf)

dist_uni_pver = phyloseq::distance(pver, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(pver))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Aug20"))
disp_date_pver = betadisper(dist_uni_pver, sampledf$Date, type = "centroid")
df_pver_dispersion = data.frame(Distance_to_centroid=disp_date_pver$distances,Date_disp=disp_date_pver$group, 
                     sampledf)

df_acr = data.frame(df_acr_shannon, df_acr_dispersion)
df_plob = data.frame(df_plob_shannon, df_plob_dispersion)
df_pver = data.frame(df_pver_shannon, df_pver_dispersion)
df_acr_bleaching = summarySEwithin(df_acr, measurevar="percent_bleached", withinvars="Date",
                         na.rm=FALSE, conf.interval=.95)
df_acr_bleaching$N = NULL
df_acr_bleaching$sd = NULL
df_acr_bleaching$se = NULL
df_acr_bleaching$ci = NULL

df_acr_mortality = summarySEwithin(df_acr, measurevar="percent_dead", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_acr_mortality$N = NULL
df_acr_mortality$sd = NULL
df_acr_mortality$se = NULL
df_acr_mortality$ci = NULL

df_acr_shannon = summarySEwithin(df_acr, measurevar="Shannon", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_acr_shannon$N = NULL
df_acr_shannon$sd = NULL
df_acr_shannon$se = NULL
df_acr_shannon$ci = NULL

df_acr_dispersion = summarySEwithin(df_acr, measurevar="Distance_to_centroid", withinvars="Date",
                                 na.rm=FALSE, conf.interval=.95)
df_acr_dispersion$N = NULL
df_acr_dispersion$sd = NULL
df_acr_dispersion$se = NULL
df_acr_dispersion$ci = NULL

df_acr_master = data.frame(df_acr_bleaching,df_acr_mortality,
                           df_acr_shannon,df_acr_dispersion)
df_acr_master$Date.1 = NULL
df_acr_master$Date.2 = NULL
df_acr_master$Date.3 = NULL

df_plob_bleaching = summarySEwithin(df_plob, measurevar="percent_bleached", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_plob_bleaching$N = NULL
df_plob_bleaching$sd = NULL
df_plob_bleaching$se = NULL
df_plob_bleaching$ci = NULL

df_plob_mortality = summarySEwithin(df_plob, measurevar="percent_dead", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_plob_mortality$N = NULL
df_plob_mortality$sd = NULL
df_plob_mortality$se = NULL
df_plob_mortality$ci = NULL

df_plob_shannon = summarySEwithin(df_plob, measurevar="Shannon", withinvars="Date",
                                 na.rm=FALSE, conf.interval=.95)
df_plob_shannon$N = NULL
df_plob_shannon$sd = NULL
df_plob_shannon$se = NULL
df_plob_shannon$ci = NULL

df_plob_dispersion = summarySEwithin(df_plob, measurevar="Distance_to_centroid", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_plob_dispersion$N = NULL
df_plob_dispersion$sd = NULL
df_plob_dispersion$se = NULL
df_plob_dispersion$ci = NULL

df_plob_master = data.frame(df_plob_bleaching,df_plob_mortality,
                           df_plob_shannon,df_plob_dispersion)
df_plob_master$Date.1 = NULL
df_plob_master$Date.2 = NULL
df_plob_master$Date.3 = NULL

df_pver_bleaching = summarySEwithin(df_pver, measurevar="percent_bleached", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_pver_bleaching$N = NULL
df_pver_bleaching$sd = NULL
df_pver_bleaching$se = NULL
df_pver_bleaching$ci = NULL

df_pver_mortality = summarySEwithin(df_pver, measurevar="percent_dead", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_pver_mortality$N = NULL
df_pver_mortality$sd = NULL
df_pver_mortality$se = NULL
df_pver_mortality$ci = NULL

df_pver_shannon = summarySEwithin(df_pver, measurevar="Shannon", withinvars="Date",
                                  na.rm=FALSE, conf.interval=.95)
df_pver_shannon$N = NULL
df_pver_shannon$sd = NULL
df_pver_shannon$se = NULL
df_pver_shannon$ci = NULL

df_pver_dispersion = summarySEwithin(df_pver, measurevar="Distance_to_centroid", withinvars="Date",
                                     na.rm=FALSE, conf.interval=.95)
df_pver_dispersion$N = NULL
df_pver_dispersion$sd = NULL
df_pver_dispersion$se = NULL
df_pver_dispersion$ci = NULL

df_pver_master = data.frame(df_pver_bleaching,df_pver_mortality,
                            df_pver_shannon,df_pver_dispersion)
df_pver_master$Date.1 = NULL
df_pver_master$Date.2 = NULL
df_pver_master$Date.3 = NULL

p1=ggplot(df_acr_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="black")+
  geom_line(color="black", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="black", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="black")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,4),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p2=ggplot(df_plob_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="#E69F00")+
  geom_line(color="#E69F00", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="#E69F00", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="#E69F00")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,4),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p3=ggplot(df_pver_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="#56B4E9")+
  geom_line(color="#56B4E9", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="#56B4E9", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="#56B4E9")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,4),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p = ggarrange(p1,p2,p3, ncol=1, align="hv")

##add dispersion axis in illustrator, it is 1/2 of left y-axis
ggsave(plot = p, "recent experimental figure 8.tiff", dpi = 1500, units = "mm",
       width = 180, height = 200, scale = 0.8)

##Supplemental Figure S13: same thing but only with corals present since Jul18
acr_since_start = df_acr$ID[df_acr$Date=="Jul18"]
df_new_acr = df_acr[df_acr$ID %in% acr_since_start,]

plob_since_start = df_plob$ID[df_plob$Date=="Jul18"]
df_new_plob = df_plob[df_plob$ID %in% plob_since_start,]

pver_since_start = df_pver$ID[df_pver$Date=="Jul18"]
df_new_pver = df_pver[df_pver$ID %in% pver_since_start,]

df_acr_bleaching = summarySEwithin(df_new_acr, measurevar="percent_bleached", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_acr_bleaching$N = NULL
df_acr_bleaching$sd = NULL
df_acr_bleaching$se = NULL
df_acr_bleaching$ci = NULL

df_acr_mortality = summarySEwithin(df_new_acr, measurevar="percent_dead", withinvars="Date",
                                   na.rm=FALSE, conf.interval=.95)
df_acr_mortality$N = NULL
df_acr_mortality$sd = NULL
df_acr_mortality$se = NULL
df_acr_mortality$ci = NULL

df_acr_shannon = summarySEwithin(df_new_acr, measurevar="Shannon", withinvars="Date",
                                 na.rm=FALSE, conf.interval=.95)
df_acr_shannon$N = NULL
df_acr_shannon$sd = NULL
df_acr_shannon$se = NULL
df_acr_shannon$ci = NULL

df_acr_dispersion = summarySEwithin(df_new_acr, measurevar="Distance_to_centroid", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_acr_dispersion$N = NULL
df_acr_dispersion$sd = NULL
df_acr_dispersion$se = NULL
df_acr_dispersion$ci = NULL

df_acr_master = data.frame(df_acr_bleaching,df_acr_mortality,
                           df_acr_shannon,df_acr_dispersion)
df_acr_master$Date.1 = NULL
df_acr_master$Date.2 = NULL
df_acr_master$Date.3 = NULL

df_plob_bleaching = summarySEwithin(df_new_plob, measurevar="percent_bleached", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_plob_bleaching$N = NULL
df_plob_bleaching$sd = NULL
df_plob_bleaching$se = NULL
df_plob_bleaching$ci = NULL

df_plob_mortality = summarySEwithin(df_new_plob, measurevar="percent_dead", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_plob_mortality$N = NULL
df_plob_mortality$sd = NULL
df_plob_mortality$se = NULL
df_plob_mortality$ci = NULL

df_plob_shannon = summarySEwithin(df_new_plob, measurevar="Shannon", withinvars="Date",
                                  na.rm=FALSE, conf.interval=.95)
df_plob_shannon$N = NULL
df_plob_shannon$sd = NULL
df_plob_shannon$se = NULL
df_plob_shannon$ci = NULL

df_plob_dispersion = summarySEwithin(df_new_plob, measurevar="Distance_to_centroid", withinvars="Date",
                                     na.rm=FALSE, conf.interval=.95)
df_plob_dispersion$N = NULL
df_plob_dispersion$sd = NULL
df_plob_dispersion$se = NULL
df_plob_dispersion$ci = NULL

df_plob_master = data.frame(df_plob_bleaching,df_plob_mortality,
                            df_plob_shannon,df_plob_dispersion)
df_plob_master$Date.1 = NULL
df_plob_master$Date.2 = NULL
df_plob_master$Date.3 = NULL

df_pver_bleaching = summarySEwithin(df_new_pver, measurevar="percent_bleached", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_pver_bleaching$N = NULL
df_pver_bleaching$sd = NULL
df_pver_bleaching$se = NULL
df_pver_bleaching$ci = NULL

df_pver_mortality = summarySEwithin(df_new_pver, measurevar="percent_dead", withinvars="Date",
                                    na.rm=FALSE, conf.interval=.95)
df_pver_mortality$N = NULL
df_pver_mortality$sd = NULL
df_pver_mortality$se = NULL
df_pver_mortality$ci = NULL

df_pver_shannon = summarySEwithin(df_new_pver, measurevar="Shannon", withinvars="Date",
                                  na.rm=FALSE, conf.interval=.95)
df_pver_shannon$N = NULL
df_pver_shannon$sd = NULL
df_pver_shannon$se = NULL
df_pver_shannon$ci = NULL

df_pver_dispersion = summarySEwithin(df_new_pver, measurevar="Distance_to_centroid", withinvars="Date",
                                     na.rm=FALSE, conf.interval=.95)
df_pver_dispersion$N = NULL
df_pver_dispersion$sd = NULL
df_pver_dispersion$se = NULL
df_pver_dispersion$ci = NULL

df_pver_master = data.frame(df_pver_bleaching,df_pver_mortality,
                            df_pver_shannon,df_pver_dispersion)
df_pver_master$Date.1 = NULL
df_pver_master$Date.2 = NULL
df_pver_master$Date.3 = NULL

p1=ggplot(df_acr_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="black")+
  geom_line(color="black", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="black", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="black")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,9),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p2=ggplot(df_plob_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="#E69F00")+
  geom_line(color="#E69F00", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="#E69F00", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="#E69F00")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,9),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p3=ggplot(df_pver_master, aes(x=Date, y=Shannon))+
  theme_bw()+
  geom_point(color="#56B4E9")+
  geom_line(color="#56B4E9", aes(group=1))+
  geom_point(aes(x=Date, y=Distance_to_centroid*2), color="#56B4E9", shape = 17)+
  geom_line(aes(x=Date, y=Distance_to_centroid*2, group = 1), 
            linetype = "dotdash", color="#56B4E9")+
  geom_point(aes(x=Date, y=percent_bleached/10),color="darkred")+
  geom_line(aes(x=Date, y=percent_bleached/10, group=1),color="darkred")+
  geom_point(aes(x=Date, y=percent_dead/10),shape = 4, color="darkred")+
  geom_line(aes(x=Date, y=percent_dead/10, group=1),linetype = 2, 
            color="darkred")+
  scale_y_continuous(name = "Shannon Diversity", limits = c(0,9),
                     sec.axis = sec_axis(trans = ~.*10, 
                                         name = "Avg. % Colony Affected"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y.right = element_text(colour="darkred"))

p = ggarrange(p1,p2,p3, ncol=1, align="hv")
ggsave(plot = p, "Figure S13.tiff", dpi = 1500, units = "mm",
       width = 180, height = 200, scale = 0.8)