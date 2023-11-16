library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(Rmisc)

##Main Text averages figure 7:
df = read.csv("diff abund average figure.csv", header = TRUE)
df = subset(df, df$Classification!="Neutral/Unknown")
head(df)

df$Date = trimws(df$Date, "left")
df$Date = factor(df$Date, levels = c("Jul18","Nov18","Mar19","Aug19","Nov19",
                                     "Mar20","Aug20"))
df$Classification = factor(df$Classification, levels = c("Detrimental",
                                                         "Other Beneficial",
                                                         "Endozoicomonadaceae"))
acr = subset(df, df$Coral=="Acr")
plob = subset(df, df$Coral=="Plob")
pver = subset(df, df$Coral=="Pver")

p1 = ggplot(acr, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p1_alt = ggplot(acr, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "black")+
  theme(legend.position = "none")

p2 = ggplot(plob, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p2_alt = ggplot(plob, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "#E69F00")+
  theme(legend.position = "none")

p3 = ggplot(pver, aes(Date, Classification, color = Classification,
                     size = avg_rel_abund))+
  geom_point(alpha = 0.7)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_color_manual(values = c("#5D3FD3","darkgreen","#173518"))+
  theme(legend.position = "none")

p3_alt = ggplot(pver, aes(Date, Classification, fill = log(avg_rel_abund)))+
  geom_tile()+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient(low = "white", high = "#56B4E9")+
  theme(legend.position = "none")

##ratios over time with mortality
survival = readRDS("survival_fams_bins_ps.rds")
df_survival = data.frame(sample_data(survival))
acr_mort = subset(df_survival, df_survival$Coral=="Acr")
acr_mort = summarySEwithin(acr_mort, measurevar="percent_dead", withinvars="Date",
                           na.rm=FALSE, conf.interval=.95)

p4 = ggplot(acr, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="black", shape = 17)+
  geom_line(color = "black")+
  theme_bw()+
  geom_point(data=acr_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=acr_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

plob_mort = subset(df_survival, df_survival$Coral=="Plob")
plob_mort = summarySEwithin(plob_mort, measurevar="percent_dead", withinvars="Date",
                           na.rm=FALSE, conf.interval=.95)

p5 = ggplot(plob, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="#E69F00", shape = 17)+
  geom_line(color = "#E69F00")+
  geom_point(data=plob_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=plob_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

pver_mort = subset(df_survival, df_survival$Coral=="Pver")
pver_mort = summarySEwithin(pver_mort, measurevar="percent_dead", withinvars="Date",
                            na.rm=FALSE, conf.interval=.95)

p6 = ggplot(pver, aes(x=Date, y=log(Ratio), group=1))+
  geom_point(size=2, color="#56B4E9", shape = 17)+
  geom_line(color = "#56B4E9")+
  geom_point(data=pver_mort, aes(x=Date, y=percent_dead/4), color="darkred")+
  geom_line(data=pver_mort, aes(x=Date, y=percent_dead/4, group=1),
            color="darkred")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4), limits = c(-1,10))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_hline(yintercept=0, linetype = "dashed", alpha = 0.5)

p = ggarrange(p4, p1, p5, p2, p6, p3, ncol=1, nrow=6, align = "v",
              heights = c(1,0.85,1,0.85,1,0.85))
ggsave(plot = p, "figure 7.tiff", units="mm",
       scale = 0.9, height = 200, width = 180, dpi = 1000)

##Stats
#Acr
unique(log(acr$Ratio))
ratios = c(7.907402,7.906351,3.869733,3.618751,4.779950,5.003201)
cor.test(ratios, acr_mort$percent_dead, method = "pearson")

#Plob
unique(log(plob$Ratio))
ratios = c(4.7945142,-0.6619048,-0.6186997,0.5672443,1.6807686,1.2601039)
cor.test(ratios, plob_mort$percent_dead, method = "pearson")

#Pver
unique(log(pver$Ratio))
ratios = c(2.170288,5.065308,7.352671)
pver_mort_mid = subset(pver_mort, pver_mort$Date==c("Aug19","Nov19","Aug20"))
cor.test(ratios, pver_mort_mid$percent_dead, method = "pearson")