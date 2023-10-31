library(ggplot2)
library(vegan)
library(scales)
library(grid) 
library(reshape2)
library(phyloseq)
library(readxl)
library(microbiome)
library(nlme)
library(lme4)
library(lmerTest)
library(devtools)
library(stats)
library(ANCOMBC)
library(ggpubr)
library(Rmisc)
library(dplyr)
library(rstatix)

##Temperature: THERM 1 and 2 (LTER 0 and LTER 2) correlation:
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(lubridate)

## Code adapted with permission from Dr. Kelly Speare

# Data --------------------------------------------------------------------
# data from Moorea Coral Reef LTER core time series, bottom mounted termistors:
# http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.1035
# Leichter, J, K. Seydel and C. Gotschalk of Moorea Coral Reef LTER. 2018. MCR LTER: Coral Reef: Benthic Water Temperature, ongoing since 2005. knb-lter-mcr.1035.11

LTER0<-read.csv("MCR_LTER00_BottomMountThermistors_20220304.csv", header=TRUE)
LTER1<-read.csv("MCR_LTER01_BottomMountThermistors_20220304.csv", header=TRUE)
LTER2<-read.csv("MCR_LTER02_BottomMountThermistors_20220304.csv", header=TRUE)
LTER3<-read.csv("MCR_LTER03_BottomMountThermistors_20220304.csv", header=TRUE)
LTER4<-read.csv("MCR_LTER04_BottomMountThermistors_20220304.csv", header=TRUE)
LTER5<-read.csv("MCR_LTER05_BottomMountThermistors_20220304.csv", header=TRUE)
LTER6<-read.csv("MCR_LTER06_BottomMountThermistors_20220304.csv", header=TRUE)

temperature <- rbind(LTER0, LTER2, LTER3, LTER4, LTER5, LTER6)
# format time and date
temperature$time_use = ymd_hms(temperature$time_local)
temperature$day = format(temperature$time_use, '%Y-%m-%d')

# check metadata
temperature$sensor_depth_m=factor(as.factor(temperature$sensor_depth_m))
temperature = temperature %>% subset(reef_type_code=='FOR') %>% subset(sensor_depth_m=="10")

temperature=subset(temperature, site=='LTER00' | site=='LTER02')

# median temperature by day - mean of all 6 LTER sites
lter.day.site.1.2 <- temperature %>% group_by(day, site) %>% summarise(temp_c = median(temperature_c)) %>% ungroup()
lter.day.site.1.2$day <- ymd(lter.day.site.1.2$day)

lter.day.site.1.2.w<-spread(lter.day.site.1.2, site, temp_c)
cor.test(lter.day.site.1.2.w$LTER00, lter.day.site.1.2.w$LTER02, 
         method = c("pearson", "kendall", "spearman"))

##Microbiome Alpha Diversity (more nutrient enrichment and herbivore reduction
##statistics can be found in 'supplemental figure 8.R' and
##'supplemental figure 9.R'):
families_rare = readRDS("rarefied_fams_ps.rds")

acr = subset_samples(families_rare, Coral == "Acr")
acr
plob = subset_samples(families_rare, Coral == "Plob")
plob
pver = subset_samples(families_rare, Coral == "Pver")
pver

##PRWST for Figure 4 Significance codes
shannon=estimate_richness(families_rare, measures = "Shannon")
sample=sample_data(families_rare)
df=data.frame(shannon,sample)
df=na.omit(df)

shapiro.test(df$Shannon)#non-normal, what about log transform?
shapiro.test(log(df$Shannon))#non-normal, proceed with pairwise WRST

pairwise.wilcox.test(df$Shannon[df$Coral=="Acr"], df$Date[df$Coral=="Acr"])
pairwise.wilcox.test(df$Shannon[df$Coral=="Plob"], df$Date[df$Coral=="Plob"])
pairwise.wilcox.test(df$Shannon[df$Coral=="Pver"], df$Date[df$Coral=="Pver"])

#LMEM and ANCOVA general form (plug in acr, plob, or pver)
mixed_df = data.frame(sample_data(pver), 
                      estimate_richness(pver, measures = "Shannon"))

#Make Tag nominal categorical variable
mixed_df$Tag = factor(mixed_df$Tag)

#Make the full model
full_model = lmer(Shannon~Date + Herbivory + Nutrients + Herbivory*Nutrients + (1|Tag), 
                  data = mixed_df,
                  na.action=na.omit)
summary(full_model)
anova(full_model)

##Beta Diversity (more nutrient enrichment and herbivore reduction
##statistics can be found in 'supplemental figure 8.R' and
##'supplemental figure 9.R'):
###Dispersion
families = readRDS("families_ps.rds")
acr = subset_samples(families, Coral=="Acr")
plob = subset_samples(families, Coral=="Plob")
pver = subset_samples(families, Coral=="Pver")

##LMEM and ANCOVA general form (plug in acr, plob, or pver)
dist_uni_pver = phyloseq::distance(pver, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(pver))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
disp_date_pver = betadisper(dist_uni_pver, sampledf$Date, type = "centroid")
df_pver = data.frame(Distance_to_centroid=disp_date_pver$distances,Date_disp=disp_date_pver$group, 
                     sampledf)
df_pver$Tag = factor(df_pver$Tag)
model = lmer(data = df_pver, Distance_to_centroid ~ Date + Herbivory +
               Nutrients + Nutrients*Herbivory + (1|Tag), na.action=na.omit)
summary(model)
anova(model)

mixed_df_res = anova_test(data=df_acr, Distance_to_centroid ~ Date + Herbivory*Nutrients)
get_anova_table(mixed_df_res)

##PWRST for Figure 6 Panel B. Use dataframes from above:
shapiro.test(df_acr$Distance_to_centroid)#non-normal, pairwise WRST
pairwise.wilcox.test(df_acr$Distance_to_centroid, df_acr$Date)

shapiro.test(df_plob$Distance_to_centroid)#non-normal, pairwise WRST
pairwise.wilcox.test(df_plob$Distance_to_centroid, df_plob$Date)#ns

shapiro.test(df_pver$Distance_to_centroid)#non-normal, pairwise WRST
pairwise.wilcox.test(df_pver$Distance_to_centroid, df_pver$Date)
#statistics for Figure 6 Panel C in 'main figure 6.R'

###PERMANOVA
#plug in acr, plob, and pver for results by each coral
dist_uni_acr = phyloseq::distance(acr, "unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
stat.data.r = as(sample_data(acr), "data.frame")
adonis(dist_uni_acr ~ Date, data = stat.data.r)
adonis(dist_uni_acr ~ Herbivory, data = stat.data.r)

#remove Jul18 for nutrient analyses
acr_nut = subset_samples(acr, Date != "Jul18")
dist_uni_acr_nut = phyloseq::distance(acr_nut, "unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr_nut))
sampledf$Date = factor(sampledf$Date, levels=c("Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
stat.data.r = as(sample_data(acr_nut), "data.frame")
adonis(dist_uni_acr_nut ~ Nutrients, data = stat.data.r)#ns

###Connecting microbiome to phenotype (note: statistically comparable bleaching
###only observed in A. retusa):
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid) 
library(reshape2)
library(phyloseq)
library(readxl)
library(ggsignif)
library(plyr)
library(grid)
library(ggpubr)
library(microbiome)
library(ggthemes)
library(lubridate)
library(tidyr)
library(magrittr)
library(randomForest)
library(knitr)
library(kableExtra)
library(cplm)
library(pscl)
library(grid) 
library(nlme)
library(lme4)
library(lmerTest)
library(devtools)
library(stats)
library(Rmisc)

##Alpha diversity:
survival_rare_filt = readRDS("survival_rare_fams_bins_ps.rds")
acr = subset_samples(survival_rare_filt, Coral == "Acr")
plob = subset_samples(survival_rare_filt, Coral == "Plob")
pver = subset_samples(survival_rare_filt, Coral == "Pver")

#LMMs (substitute acr, plob, pver for results by species)
mixed_df = data.frame(sample_data(acr), 
                      estimate_richness(acr, measures = "Shannon"))

model = lmer(percent_bleached ~ Shannon + Nutrients +
               Herbivory+ Nutrients*Herbivory + (1|ID), 
             data = mixed_df, na.action = na.omit)
summary(model)
anova(model)

model = lmer(percent_dead ~ Shannon + Nutrients +
               Herbivory+ Nutrients*Herbivory + (1|ID), 
             data = mixed_df, na.action = na.omit)
summary(model)
anova(model)

##Beta dispersion:
survival_filt = readRDS("survival_fams_bins_ps.rds")
acr = subset_samples(survival_filt, Coral == "Acr")
plob = subset_samples(survival_filt, Coral == "Plob")
pver = subset_samples(survival_filt, Coral == "Pver")

#LMMs (substitute acr, plob, pver for results by species)
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
sampledf = data.frame(sample_data(acr))
sampledf$Date = factor(sampledf$Date, levels=c("Jul18", "Nov18", "Mar19", "Aug19",
                                               "Nov19", "Mar20", "Aug20"))
disp_date_acr = betadisper(dist_uni_acr, sampledf$Date, type = "centroid")
df_acr = data.frame(Distance_to_centroid=disp_date_acr$distances,Date_disp=disp_date_acr$group, 
                     sampledf)

model = lmer(percent_bleached ~ Distance_to_centroid + Nutrients + Herbivory +
               Nutrients*Herbivory + (1|ID), 
             data = df_acr, na.action = na.omit)
summary(model)
anova(model)

model = lmer(percent_dead ~ Distance_to_centroid + Nutrients + Herbivory +
               Nutrients*Herbivory + (1|ID), 
             data = df_acr, na.action = na.omit)
summary(model)
anova(model)

##Beta dissimilarity (PERMANOVA) (sub in acr, plob, pver for results by species):
dist_uni_acr = phyloseq::distance(acr, method="unifrac", weighted=TRUE)
stat.data.r = as(sample_data(acr), "data.frame")

#Bleaching and mortality Severity
adonis2(dist_uni_acr~percent_bleached, data = stat.data.r)
adonis2(dist_uni_acr~percent_dead, data = stat.data.r)

#Bleaching and mortality Prevalence
stat.data.r$binary_mortality = 0
for (i in 1:nrow(stat.data.r)) {
  if(stat.data.r$percent_dead[i]>0){
    stat.data.r$binary_mortality[i] = 1
  }
  else {
    stat.data.r$binary_mortality[i] = 0
  }
}

stat.data.r$binary_bleaching = 0
for (i in 1:nrow(stat.data.r)) {
  if(stat.data.r$percent_bleached[i]>0){
    stat.data.r$binary_bleaching[i] = 1
  }
  else {
    stat.data.r$binary_bleaching[i] = 0
  }
}
stat.data.r$binary_mortality = factor(stat.data.r$binary_mortality)
stat.data.r$binary_bleaching = factor(stat.data.r$binary_bleaching)

adonis2(dist_uni_acr ~ binary_mortality, data = stat.data.r)
adonis2(dist_uni_acr ~ binary_bleaching, data = stat.data.r)

#PWRST for FIgure 8F
dist_uni_pver = phyloseq::distance(pver, "unifrac", weighted=TRUE)
stat.data.r = as(sample_data(pver), "data.frame")

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
disp_mort_pver = betadisper(dist_uni_pver, stat.data.r$binary_mortality, type = "centroid")
stat.data.r$mort_disp = disp_mort_pver$distances

aug=subset(stat.data.r, Date=="Aug20")

pairwise.wilcox.test(aug$mort_disp, aug$binary_mortality)