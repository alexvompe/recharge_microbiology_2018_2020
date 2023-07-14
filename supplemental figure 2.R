library(lubridate)
library(scales)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

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

# median temperature by day - mean of all 6 LTER sites
lter.day.site.1.2 <- temperature %>% group_by(day, site) %>% summarise(temp_c = median(temperature_c)) %>% ungroup()
lter.day.site.1.2$day <- ymd(lter.day.site.1.2$day)

lter.day.site.1.2.w<-spread(lter.day.site.1.2, site, temp_c)

#Test correlation of LTER 0 and LTER 2 temperature data
lter00_lter02_daily_correlation<-ggplot(lter.day.site.1.2.w, aes(x=LTER00, y=LTER02))+
  geom_point()+
  geom_abline(intercept=0, slope=1, color="red")+
  stat_smooth()
ggsave("lter00_lter02_daily_correlation.pdf", width=4, height=4, units="in")

cor.test(lter.day.site.1.2.w$LTER00, lter.day.site.1.2.w$LTER02, method = c("pearson", "kendall", "spearman"))