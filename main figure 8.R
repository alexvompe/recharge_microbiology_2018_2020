library(tidyverse)
library(ggpubr)
library(phyloseq)

##Color gradient of mortality and bleaching for PCoA1 by Date regression
survival_filt = readRDS("survival_fams_bins_ps.rds")
acr = subset_samples(survival_filt, Coral == "Acr")

ord_acr = ordinate(acr, method = "PCoA", distance = "unifrac", 
                   weighted = TRUE)
ord_acr_df = plot_ordination(acr, ord_acr, type="samples", 
                             color="Date", justDF = TRUE) 
#59.2 and 10.1 variance explained by 1st 2 axes%

#Dispersion differences by binary mortality and bleaching->
p1=ggplot(ord_acr_df, aes(x=Axis.1, y=Axis.2, color=percent_dead>0,
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

p2=ggplot(ord_acr_df, aes(x=Axis.1, y=Axis.2, color=percent_bleached>0,
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

p3=ggplot(ord_acr_df, aes(x=Date, y=Axis.1, color=percent_dead, group=1)) + 
  theme_bw()+
  geom_point(size=3, position="jitter")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_smooth(aes(as.numeric(Date), Axis.1, group=1), method = "loess",
              color="black", linewidth=1.2)+
  scale_color_gradientn(colours = c("#C0E5FE","#60B2FE", 
                                             "blue", "darkblue","black", 
                                             "#400000","darkred", "red"), 
                                             "Percent Bleached")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  ylab("PCoA Axis 1 [59.2%]")+
  xlab("Date")+
  theme(legend.position="none")

p4=ggplot(ord_acr_df, aes(x=Date, y=Axis.1, color=percent_bleached, group=1)) + 
  theme_bw()+
  geom_point(size=3, position="jitter")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_smooth(aes(as.numeric(Date), Axis.1, group=1), method = "loess",
              color="black", linewidth=1.2)+
  scale_color_gradientn(colours = c("#C0E5FE","#60B2FE", 
                                             "blue", "darkblue","black", 
                                             "#400000","darkred", "red"), 
                                             "Percent Bleached")+
  theme(text = element_text(size = 26), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  ylab("PCoA Axis 1 [59.2%]")+
  xlab("Date")+
  theme(legend.position="none")

p = ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
ggsave(plot=p, "main figure 8 intermediate.tiff",height = 200,
       width = 180, units = "mm", scale=2)