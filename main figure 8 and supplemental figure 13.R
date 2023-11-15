library(tidyverse)
library(vegan)
library(ggpubr)
library(phyloseq)
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