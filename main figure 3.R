library(ggplot2)
library(Rmisc)
library(ggpubr)

#Read in mortality and bleaching data. These data are for all corals ever
#sampled for this project, tracked since the start of the time series.
mortality = read.csv("porites integration_v3.csv", header = TRUE)
mortality$Date = factor(mortality$Date, levels = c("Jul18", "Nov18",
                                                   "Mar19", "Aug19",
                                                   "Nov19", "Mar20", "Aug20"))
mortality_acr = subset(mortality, Coral=="Acr")
mortality_pver = subset(mortality, Coral=="Pver")
mortality_plob = subset(mortality, Coral=="Plob")

#Mortality stats
model = lm(percent_dead ~ Coral, data = mortality)
TukeyHSD(aov(model))
mortality_acr_plob = subset(mortality, Coral=="Acr" | Coral=="Plob")
t.test(percent_dead ~ Coral, data = mortality_acr_plob, alternative = "two.sided")
mortality_aug20 = subset(mortality, Date=="Aug20")
model = lm(percent_dead ~ Coral, data = mortality_aug20)
TukeyHSD(aov(model))
mortality_Mar19 = subset(mortality, Date=="Mar19")
model = lm(percent_dead ~ Coral, data = mortality_Mar19)
TukeyHSD(aov(model))

#Mortality
stats_acr = summarySEwithin(mortality_acr, measurevar="percent_dead", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

stats_pver = summarySEwithin(mortality_pver, measurevar="percent_dead", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

stats_plob = summarySEwithin(mortality_plob, measurevar="percent_dead", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

theme_set(theme_bw())
p1 = ggplot()+
  geom_point(data=mortality_acr, aes(x=Date, y=percent_dead), position="jitter", 
             size=4, alpha=0.2, color="black")+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="black", size = 8)+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line", linetype=2, color="black", linewidth = 2)+
  geom_errorbar(data=stats_acr, width=0.3, size=3, aes(x=Date, y=percent_dead,
                                                       ymin=percent_dead-ci,
                                                       ymax=percent_dead+ci), color="black")+
  geom_point(data=mortality_plob, aes(x=Date, y=percent_dead), position="jitter", 
             size=4, alpha=0.2, color="#E69F00")+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="#E69F00", size = 8)+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line", linetype=2, color="#E69F00", linewidth = 2)+
  geom_errorbar(data=stats_plob, width=0.3, size=3, aes(x=Date, y=percent_dead,
                                                        ymin=percent_dead-ci,
                                                        ymax=percent_dead+ci), color="#E69F00")+
  geom_point(data=mortality_pver, aes(x=Date, y=percent_dead), position="jitter", 
             size=4, alpha=0.2, color="#56B4E9")+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="#56B4E9", size = 8)+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line", linetype=2, color="#56B4E9", linewidth = 2)+
  geom_errorbar(data=stats_pver, width=0.3, size=3, aes(x=Date, y=percent_dead,
                                                        ymin=percent_dead-ci,
                                                        ymax=percent_dead+ci), color="#56B4E9")+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("% of Colony Dead")+
  theme(legend.position="none")

#Bleaching
#Make datasets without any fully dead colonies, so that mortality doesn't confound bleaching levels
bleaching_acr = subset(mortality_acr, percent_dead != 100)
bleaching_plob = subset(mortality_plob, percent_dead != 100)
bleaching_pver = subset(mortality_pver, percent_dead != 100)
bleaching_acr = na.omit(bleaching_acr)
bleaching_plob = na.omit(bleaching_plob)
bleaching_pver = na.omit(bleaching_pver)

#What proportion of corals had signs of bleaching in 2019 vs 2020?
a = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="Mar19"])
b = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="Mar19" & 
                                            bleaching_acr$percent_bleached > 0])
b/a #2019: 53%

c = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="Aug20"])
d = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="Aug20" & 
                                            bleaching_acr$percent_bleached > 0])
d/c #2020: 2%

#Bleaching stats
model = lm(percent_bleached ~ Date, data = bleaching_acr)
TukeyHSD(aov(model))

#Make plot with 95% confidence intervals
stats_acr = summarySEwithin(bleaching_acr, measurevar="percent_bleached", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

stats_plob = summarySEwithin(bleaching_plob, measurevar="percent_bleached", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

stats_pver = summarySEwithin(bleaching_pver, measurevar="percent_bleached", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

theme_set(theme_bw())
p2 = ggplot()+
  geom_point(data=mortality_acr, aes(x=Date, y=percent_bleached), position="jitter", 
             size=4, alpha=0.2, color="black")+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="black", size = 8)+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", linetype=2, color="black", size = 2)+
  geom_errorbar(data=stats_acr, width=0.3, size=3, aes(x=Date, y=percent_bleached,
                                                       ymin=percent_bleached-ci,
                                                       ymax=percent_bleached+ci), color="black")+
  geom_point(data=mortality_plob, aes(x=Date, y=percent_bleached), position="jitter", 
             size=4, alpha=0.2, color="#E69F00")+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="#E69F00", size = 8)+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", linetype=2, color="#E69F00", size = 2)+
  geom_errorbar(data=stats_plob, width=0.3, size=3, aes(x=Date, y=percent_bleached,
                                                        ymin=percent_bleached-ci,
                                                        ymax=percent_bleached+ci), color="#E69F00")+
  geom_point(data=mortality_pver, aes(x=Date, y=percent_bleached), position="jitter", 
             size=4, alpha=0.2, color="#56B4E9")+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="#56B4E9", size = 8)+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", linetype=2, color="#56B4E9", size = 2)+
  geom_errorbar(data=stats_pver, width=0.3, size=3, aes(x=Date, y=percent_bleached,
                                                        ymin=percent_bleached-ci,
                                                        ymax=percent_bleached+ci), color="#56B4E9")+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("% of Colony Bleached")+
  theme(legend.position="none")

#Overall plot
p = ggarrange(p2,p1,ncol=2)
ggsave(plot=p, filename="Figure 3.tiff", scale=2,
       width=180,height=100,units="mm")