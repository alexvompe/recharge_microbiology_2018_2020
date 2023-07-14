library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)

#Read in data
Rarefied_w_Other = readRDS("ps_for_stacked_bar.rds")

#First, merge samples by covariates of interest
merge1 = merge_samples(Rarefied_w_Other, "Coral")
merge2 = merge_samples(Rarefied_w_Other, "Date")
variable1 = as.character(get_variable(Rarefied_w_Other, "Coral"))
variable2 = as.character(get_variable(Rarefied_w_Other, "Date"))
sample_data(Rarefied_w_Other)$Coralbydate <- mapply(paste0, variable1, variable2, 
                                                    collapse = "_")
merge3 = merge_samples(Rarefied_w_Other, "Coralbydate")

#Transform sample counts with the merged object
relative3 = transform_sample_counts(merge3, 
                                    function(x) {x/sum(x)} )

#For plots that need faceting, make manual modifications
sample_data(relative3)
x = data.frame(sample_data(relative3))
x$Date = NULL
x$Coral = NULL
x$Date = c("Aug19","Aug20","Jul18","Mar19","Mar20","Nov18","Nov19",
           "Aug19","Aug20","Jul18","Mar19","Mar20","Nov18","Nov19",
           "Aug19","Aug20","Jul18","Mar19","Mar20","Nov18","Nov19")#This ordering
#was done when Date was an unordered factor object. The order might be different
#for both coral and date with the supplied phyloseq. Check x$Date and x$Coral
#before assigning the list.

x$Coral = c("Acr","Acr","Acr","Acr","Acr","Acr","Acr",
            "Plob","Plob","Plob","Plob","Plob","Plob","Plob",
            "Pver","Pver","Pver","Pver","Pver","Pver","Pver")
sample_data(relative3) = x

#Set taxa order for plot
taxa_order = list("Amoebophilaceae", "Cyanobiaceae", "Cyclobacteriaceae", "Endozoicomonadaceae", "Entomoplasmatales Incertae Sedis", "Flavobacteriaceae", "Moraxellaceae", "NA_Alphaproteobacteria", "NA_Bacilli", "NA_Campylobacterales", "NA_Gammaproteobacteria", "NA_Proteobacteria", "Oxalobacteraceae", "Pirellulaceae", "Rhodobacteraceae", "Simkaniaceae", "Sphingomonadaceae", "Vibrionaceae", "Xenococcaceae", "Other")
desired_order = list("Jul18","Nov18","Mar19","Aug19","Nov19","Mar20","Aug20")
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#196F3D",
              "#922B21", "#0055CC", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "#D55E00", 
              "#CC79A7", "#999999", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "#808080", 
              "#B4CEFF")

#Relative Abundance Plot
theme_set(theme_bw())
p=plot_bar(relative3, fill = "Family", 
           x="Date",facet_grid=~Coral)+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=cbPalette) + 
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Date") + 
  ylab("Relative Abundance")
p$data$Date = factor(p$data$Date, levels = desired_order)
p$data$Family = factor(p$data$Family, levels = taxa_order)

ggsave(plot=p, filename="community composition by date.tiff", scale=2,
       width=180,height=100,units="mm")