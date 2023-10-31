library(ggplot2)
library(Rmisc)
library(ggpubr)

#Read in mortality and bleaching data. These data are for all corals ever
#sampled for this project, tracked since the start of the time series.
mortality = read.csv("porites integration_v3.csv", header = TRUE)
mortality$Date = trimws(mortality$Date, "left")
mortality$Date = factor(mortality$Date, levels = c("Jul18", "Nov18",
                                                   "Mar19", "May19","Aug19",
                                                   "Nov19", "Mar20", "Aug20"))
mortality = na.omit(mortality)
mortality_acr = subset(mortality, Coral=="Acr")
mortality_pver = subset(mortality, Coral=="Pver")
mortality_plob = subset(mortality, Coral=="Plob")

#Make datasets without any fully dead colonies, so that mortality doesn't confound bleaching levels
bleaching_acr = subset(mortality_acr, percent_dead != 100)
bleaching_plob = subset(mortality_plob, percent_dead != 100)
bleaching_pver = subset(mortality_pver, percent_dead != 100)
bleaching = subset(mortality, percent_dead != 100)

#Stats
all_test_bleaching = aov(percent_bleached ~ Coral, data = bleaching)
summary(all_test_bleaching)
TukeyHSD(all_test_bleaching)
all_test_mortality = aov(percent_dead ~ Coral, data = mortality)
summary(all_test_mortality)
TukeyHSD(all_test_mortality)

bleaching_may19 = subset(bleaching, Date=="May19")
mhw1_bleaching = aov(percent_bleached ~ Coral, data = bleaching_may19)
TukeyHSD(mhw1_bleaching)
mortality_aug19 = subset(mortality, Date=="Aug19")
mhw1_mortality = aov(percent_dead ~ Coral, data = mortality_aug19)
TukeyHSD(mhw1_mortality)

mortality_aug20 = subset(mortality, Date=="Aug20")
mhw2_mortality = aov(percent_dead ~ Coral, data = mortality_aug20)
TukeyHSD(mhw2_mortality)

model = aov(percent_dead ~ Date, data = mortality_acr)
TukeyHSD(model)

model = aov(percent_dead ~ Date, data = mortality_plob)
TukeyHSD(model)

model = aov(percent_dead ~ Date, data = mortality_pver)
TukeyHSD(model)

#What proportion of corals had signs of bleaching in May 2019?
a = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="May19"])
b = length(bleaching_acr$percent_bleached[bleaching_acr$Date=="May19" & 
                                            bleaching_acr$percent_bleached > 0])
b/a #2019: 58%

#Mortality plot
stats_acr = summarySEwithin(mortality_acr, measurevar="percent_dead", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

stats_pver = summarySEwithin(mortality_pver, measurevar="percent_dead", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

stats_plob = summarySEwithin(mortality_plob, measurevar="percent_dead", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

theme_set(theme_bw())
p1 = ggplot()+
  geom_point(data=mortality_acr, aes(x=Date, y=percent_dead), position="jitter", 
             alpha=0.2, color="black")+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="black")+
  stat_summary(data=mortality_acr, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line", color="black")+
  geom_errorbar(data=stats_acr, width=0.1,aes(x=Date, y=percent_dead,
                                                       ymin=percent_dead-ci,
                                                       ymax=percent_dead+ci), color="black")+
  geom_point(data=mortality_plob, aes(x=Date, y=percent_dead), position="jitter", 
             alpha=0.2, color="#E69F00")+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="#E69F00")+
  stat_summary(data=mortality_plob, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line", color="#E69F00")+
  geom_errorbar(data=stats_plob, width=0.1,aes(x=Date, y=percent_dead,
                                                        ymin=percent_dead-ci,
                                                        ymax=percent_dead+ci), color="#E69F00")+
  geom_point(data=mortality_pver, aes(x=Date, y=percent_dead), position="jitter", 
             alpha=0.2, color="#56B4E9")+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="point", color="#56B4E9")+
  stat_summary(data=mortality_pver, aes(x=Date, y=percent_dead, group=1),fun=mean, 
               geom="line",color="#56B4E9")+
  geom_errorbar(data=stats_pver, width=0.1, aes(x=Date, y=percent_dead,
                                                        ymin=percent_dead-ci,
                                                        ymax=percent_dead+ci), color="#56B4E9")+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("% of Colony Dead")+
  theme(legend.position="none")

#Bleaching plot
stats_acr = summarySEwithin(bleaching_acr, measurevar="percent_bleached", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

stats_plob = summarySEwithin(bleaching_plob, measurevar="percent_bleached", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

stats_pver = summarySEwithin(bleaching_pver, measurevar="percent_bleached", withinvars="Date",
                             na.rm=FALSE, conf.interval=.95)

theme_set(theme_bw())
p2 = ggplot()+
  geom_point(data=bleaching_acr, aes(x=Date, y=percent_bleached), position="jitter", 
             alpha=0.2, color="black")+
  stat_summary(data=bleaching_acr, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="black")+
  stat_summary(data=bleaching_acr, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", color="black")+
  geom_errorbar(data=stats_acr, width=0.1,aes(x=Date, y=percent_bleached,
                                                       ymin=percent_bleached-ci,
                                                       ymax=percent_bleached+ci), color="black")+
  geom_point(data=bleaching_plob, aes(x=Date, y=percent_bleached), position="jitter", 
             alpha=0.2, color="#E69F00")+
  stat_summary(data=bleaching_plob, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="#E69F00")+
  stat_summary(data=bleaching_plob, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", color="#E69F00")+
  geom_errorbar(data=stats_plob, width=0.1, aes(x=Date, y=percent_bleached,
                                                        ymin=percent_bleached-ci,
                                                        ymax=percent_bleached+ci), color="#E69F00")+
  geom_point(data=bleaching_pver, aes(x=Date, y=percent_bleached), position="jitter", 
             alpha=0.2, color="#56B4E9")+
  stat_summary(data=bleaching_pver, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="point", color="#56B4E9")+
  stat_summary(data=bleaching_pver, aes(x=Date, y=percent_bleached, group=1),fun=mean, 
               geom="line", color="#56B4E9")+
  geom_errorbar(data=stats_pver, width=0.1,aes(x=Date, y=percent_bleached,
                                                        ymin=percent_bleached-ci,
                                                        ymax=percent_bleached+ci), color="#56B4E9")+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("% of Colony Bleached")+
  theme(legend.position="none")

#Overall plot
p = ggarrange(p2,p1,ncol=1, align = "hv")
ggsave(plot=p, filename="Figure 3 workshop.tiff", scale=0.6,
       width=150,height=200,units="mm", dpi = 1000)