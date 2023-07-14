library(scales)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(Rmisc)
library(phyloseq)

##Figure 4 Panels A & B alpha points + 95% CIs
families_rare = readRDS("rarefied_fams_ps.rds")

acr = subset_samples(families_rare, Coral == "Acr")
plob = subset_samples(families_rare, Coral == "Plob")
pver = subset_samples(families_rare, Coral == "Pver")

shannon_acr=estimate_richness(acr, measures = "Shannon")
sample_acr=sample_data(acr)
df_acr=data.frame(shannon_acr,sample_acr)
ci_acr = summarySEwithin(df_acr, measurevar="Shannon", withinvars="Date",
                         na.rm=FALSE, conf.interval=.95)

ci_acr

shannon_plob=estimate_richness(plob, measures = "Shannon")
sample_plob=sample_data(plob)
df_plob=data.frame(shannon_plob,sample_plob)
ci_plob = summarySEwithin(df_plob, measurevar="Shannon", withinvars="Date",
                          na.rm=FALSE, conf.interval=.95)

ci_plob

shannon_pver=estimate_richness(pver, measures = "Shannon")
sample_pver=sample_data(pver)
df_pver=data.frame(shannon_pver,sample_pver)
ci_pver = summarySEwithin(df_pver, measurevar="Shannon", withinvars="Date",
                          na.rm=FALSE, conf.interval=.95)

ci_pver

p = ggplot()+
  theme_bw()+
  geom_errorbar(data=ci_acr, width=0.1, linewidth=1.2, aes(x=Date, y=Shannon,
                                                      ymin=Shannon-ci,
                                                      ymax=Shannon+ci), color="black")+
  geom_point(data=ci_acr, size=4, color="black", aes(x=Date, y=Shannon, group=1))+
  geom_errorbar(data=ci_plob, width=0.1, linewidth=1.2, aes(x=Date, y=Shannon,
                                                       ymin=Shannon-ci,
                                                       ymax=Shannon+ci), color="#E69F00")+
  geom_point(data=ci_plob, size=4, color="#E69F00", aes(x=Date, y=Shannon, group=1))+
  geom_errorbar(data=ci_pver, width=0.1, linewidth=1.2, aes(x=Date, y=Shannon,
                                                       ymin=Shannon-ci,
                                                       ymax=Shannon+ci), color="#56B4E9")+
  geom_point(data=ci_pver, size=4, color="#56B4E9", aes(x=Date, y=Shannon, group=1))+
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Shannon Diversity Index")

ggsave(plot=p, filename="richness by date_points.tiff", scale=2,
       width=180,height=100,units="mm") #overlayed temp and alpha div in
#Panels A & B overlayed with Figure 2 Panels B and C in Adobe Illustrator and
#Powerpoint.

##Figure 4 panel C
#Colorblind-friendly palette
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
                       "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
                       "#009E73", "#0072B2", "#D55E00", 
                       "#CC79A7", "#999999", "#FF468F", "#89472F", 
                       "#F0E442", "#FF4040", "#66CCCC", "#808080", 
                       "#B4CEFF")
                       
p=plot_richness(families_rare, x="Date", 
                measures="Shannon", color = "Coral") + 
  facet_grid(.~Coral)+
  theme_bw()+
  geom_boxplot(alpha=0.6, linewidth=1.2) + 
  scale_color_manual(values=cbPalette) +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Shannon Diversity Index")+
  theme(legend.position="none")

ggsave(plot=p, filename="shannon by date.tiff", scale=2,
       width=180,height=100,units="mm")