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

LTER0=read.csv("MCR_LTER00_BottomMountThermistors_20220304.csv", header=TRUE)
LTER1=read.csv("MCR_LTER01_BottomMountThermistors_20220304.csv", header=TRUE)
LTER2=read.csv("MCR_LTER02_BottomMountThermistors_20220304.csv", header=TRUE)
LTER3=read.csv("MCR_LTER03_BottomMountThermistors_20220304.csv", header=TRUE)
LTER4=read.csv("MCR_LTER04_BottomMountThermistors_20220304.csv", header=TRUE)
LTER5=read.csv("MCR_LTER05_BottomMountThermistors_20220304.csv", header=TRUE)
LTER6=read.csv("MCR_LTER06_BottomMountThermistors_20220304.csv", header=TRUE)

temperature = rbind(LTER0, LTER2, LTER3, LTER4, LTER5, LTER6)
#saveRDS(temperature,file='temperature_combined.Rdata')
# temperature = readRDS(file='temperature_combined.Rdata') # does not work from saved file for some reason

# Figure 1 - Temperature trends -------------------------------------------
# format time and date
temperature$time_use = ymd_hms(temperature$time_local)
temperature$day = format(temperature$time_use, '%Y-%m-%d')

# check metadata
temperature$sensor_depth_m=factor(as.factor(temperature$sensor_depth_m))
temperature = temperature %>% subset(reef_type_code=='FOR') %>% subset(sensor_depth_m=="10")

temperature=subset(temperature, site=='LTER00' | site=='LTER02')
temperature.0=subset(temperature, site=="LTER00")

temperature.0$day=as.Date(temperature.0$day, '%Y-%m-%d')

temperature.0$day=ymd(temperature.0$day)

# median temperature by day - mean of all 6 LTER sites
lter.day = temperature.0 %>% group_by(day) %>% summarise(temp_c = mean(temperature_c)) %>% ungroup()

lter.day$day=as.Date(lter.day$day, '%Y-%m-%d')

lter.day$day = ymd(lter.day$day)

lter.day$day = as.Date(lter.day$day, '%Y-%m-%d')

lter.time.seq = data.frame(day=unique(lter.day$day),ind=seq(1:length(unique(lter.day$day))))
lter.time.seq$year = year(lter.time.seq$day)
lter.time.seq$week.num = week(lter.time.seq$day)
lter.time.seq.week = lter.time.seq %>% group_by(week.num,year) %>% summarise(day=min(day)) %>% ungroup()

lter.day.temp=left_join(lter.day, lter.time.seq, by='day')
lter.week.temp= lter.day.temp %>% group_by(year,week.num) %>% summarise(temp_c = mean(temp_c)) %>% ungroup()
lter.week.temp = lter.week.temp[with(lter.week.temp, order(year,week.num)),]

mma_ref = 29

# calculate accumulated heat stress
lter.week.temp$hotspot = lter.week.temp$temp_c - mma_ref
lter.week.temp$hotspot[lter.week.temp$hotspot < 0] = 0

lter.week.temp$cumstress = NA

# 12 week running sum 
lter.out = lter.week.temp 
for(i in 13:nrow(lter.out)){
  lter.out$cumstress[i] = sum(lter.out$hotspot[(i-12):i],na.rm=T)
}

lter.out$cumstress[is.na(lter.out$cumstress)] = 0
lter.out = left_join(lter.out, lter.time.seq.week, by=c('year','week.num'))

lter.out.2018.2020=subset(lter.out, day >= '2018-07-01' & day <= '2020-08-31')

##Figure 1 Plot Panel 2
p2=ggplot(lter.out.2018.2020, aes(x=day, y=cumstress))+
  geom_line(color="black", linewidth=1.2)+
  theme_bw()+
  scale_x_date(breaks = date_breaks("months"),labels = date_format("%b%y"))+
  labs(x="Date", y="Accumulated Thermal Stress")+
  theme(axis.text.x = element_text(colour="black", angle=45, vjust=1, hjust=1), 
        axis.text.y = element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(text = element_text(size = 30))

##Figure 1 plot history of disturbance Panel 1
lter.out_hist_subset = subset(lter.out, day >= '2009-12-01' & day <= '2021-01-01')
p1=ggplot(lter.out_hist_subset, aes(x=day, y=cumstress))+
  geom_line(color="black", size=1.2)+
  theme_bw()+
  scale_x_date(breaks = date_breaks("years"),labels = date_format("%Y"))+
  labs(x="Year", y="Accumulated Thermal Stress")+
  theme(axis.text.x = element_text(colour="black", angle=45, vjust=1, hjust=1), 
        axis.text.y = element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(text = element_text(size = 30))

############################################################################################# 
#only plotting data from LTER 0

thermTemp.0=temperature.0

# format time and date
#left the time_local and time_utc columns alone, unaltered. created new date column to reformat
thermTemp.0$date_use = ymd_hms(thermTemp.0$time_local)
thermTemp.0$day_mo_yr = format(thermTemp.0$date_use, '%Y-%m-%d')
thermTemp.0$day = format(thermTemp.0$date_use, '%d') #making new factor for each day (number but as a factor)
thermTemp.0$month = format(thermTemp.0$date_use, '%m') #making new factor for each month (number but as a factor)
thermTemp.0$year = format(thermTemp.0$date_use, '%y') 

thermTemp.0$day=as.factor(thermTemp.0$day)
thermTemp.0$month=as.factor(thermTemp.0$month)
thermTemp.0$year=as.factor(thermTemp.0$year)
thermTemp.0$date_use=as.Date(thermTemp.0$date_use, '%Y-%m-%d')

# subsetting data from August 2018 to July 2020. this will become the line for the high thermal stress year
thermTemp_2018_2020=subset(thermTemp.0, date_use>="2018-07-01" & date_use<="2020-08-31")
thermTemp_2018_2020_mean=ddply(thermTemp_2018_2020, .(day, month, year), summarize,
                                mean_daily_temp=mean(temperature_c),
                                sd=NA)
thermTemp_2018_2020_mean$timeframe=as.factor("2018to2020")

thermTemp_2018_2020_mean$date= as.Date(with(thermTemp_2018_2020_mean, paste(month, day, year,sep="-")), format="%m-%d-%Y")
years_2018_2020_mean=thermTemp_2018_2020_mean[c(1:3,7,6,4,5)]
#making the upper and lower (mean +- sd) columns
years_2018_2020_mean$temp_upper=years_2018_2020_mean$mean_daily_temp + years_2018_2020_mean$sd
years_2018_2020_mean$temp_lower=years_2018_2020_mean$mean_daily_temp - years_2018_2020_mean$sd
colnames(years_2018_2020_mean)[6]="temp"

# creating a new df with only data up to Dec 31 2017. this will become the mean line and SD
# is >700 rows long because each temp value is listed twice so that it plots over 2 years
thermTemp_toDec2017=subset(thermTemp.0, date_use<"2017-12-31")
thermTemp_toDec2017_mean=ddply(thermTemp_toDec2017, .(day, month), summarize,
                                mean_daily_temp=mean(temperature_c),
                                sd=sd(temperature_c))

#creating a factor column for this year
thermTemp_toDec2017_mean$timeframe=as.factor("mean")

thermTemp_toDec2017_mean=merge(thermTemp_2018_2020_mean[c(1:3,7)], thermTemp_toDec2017_mean, by=c("day", "month"))
years_toDec2017_mean=thermTemp_toDec2017_mean[c(1:4,7,5,6)] #reordering columns
#making the upper and lower (mean +- sd) columns
years_toDec2017_mean$temp_upper=years_toDec2017_mean$mean_daily_temp + years_toDec2017_mean$sd
years_toDec2017_mean$temp_lower=years_toDec2017_mean$mean_daily_temp - years_toDec2017_mean$sd
colnames(years_toDec2017_mean)[6]="temp"

data=rbind(years_toDec2017_mean, years_2018_2020_mean)

data=merge(data, thermTemp.0[c(7,9:13)], by=c("day", "month", "year"))

##Figure 1 Panel 3
p3=ggplot(data)+
  geom_ribbon(aes(x=date_use, ymin=temp_lower, ymax=temp_upper, fill=timeframe), alpha=0.2)+
  geom_hline(yintercept=29, linetype=2, color="black", size = 1.2)+
  geom_line(aes(x=date_use, y=temp, color=timeframe))+
  scale_color_manual(values = c("purple", "black")) +
  scale_fill_manual(values = c("purple", "pink"))+
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%b%y"), 
               limits=as.Date(c('2018-07-01', '2020-08-31')))+
  labs(x="Date", y=expression("Temperature " ( degree*C)))+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black", angle=45, vjust=1, hjust=1), 
        axis.text.y = element_text(colour="black"))+
  theme(text = element_text(size = 30))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

##Make figure 1. Additional annotations added in powerpoint.
p = ggarrange(p1,p2,p3, ncol=1)
ggsave(plot=p, filename="updated main figure 1.tiff", scale=3,
       width=180,height=200,units="mm")