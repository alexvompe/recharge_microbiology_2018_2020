library(tidyverse)
library(ggmap)

##Figure 1: annotated map of Moorea
#Terrain with labels
map = get_stamenmap(bbox = c(left = -149.93, 
                             bottom = -17.6, 
                             right = -149.75, 
                             top = -17.465),
                    maptype = "terrain", zoom=13)

p = ggmap(map)+theme_bw()+
  labs(x="Longitude", y="Latitude")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_point(aes(x=-149.817650,y=-17.473100), size=6)+
  geom_point(aes(x=-149.839050,y=-17.475333), size=6)+
  geom_point(aes(x=-149.806624,y=-17.471851), size=6)

ggsave(plot=p, "terrain moorea.tiff", units="mm", height=200,
       width=360, scale=0.5)

#Terrain no labels
map = get_stamenmap(bbox = c(left = -149.93, 
                             bottom = -17.6, 
                             right = -149.75, 
                             top = -17.465),
                    maptype = "terrain", zoom=13)

p = ggmap(map)+theme_bw()+
  labs(x="Longitude", y="Latitude")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot=p, "terrain moorea_no labels.tiff", units="mm", height=200,
       width=360, scale=0.5)

#Watercolor
map = get_stamenmap(bbox = c(left = -149.93, 
                             bottom = -17.6, 
                             right = -149.75, 
                             top = -17.465),
                    maptype = "watercolor", zoom=13)

p = ggmap(map)+theme_bw()+
  labs(x="Longitude", y="Latitude")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_point(aes(x=-149.817650,y=-17.473100), size=6)+
  geom_point(aes(x=-149.839050,y=-17.475333), size=6)+
  geom_point(aes(x=-149.806624,y=-17.471851), size=6)

ggsave(plot=p, "Watercolor moorea.tiff", units="mm", height=200,
       width=360, scale=0.5)

#Watercolor no labels
map = get_stamenmap(bbox = c(left = -149.93, 
                             bottom = -17.6, 
                             right = -149.75, 
                             top = -17.465),
                    maptype = "watercolor", zoom=13)

p = ggmap(map)+theme_bw()+
  labs(x="Longitude", y="Latitude")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot=p, "Watercolor moorea_no labels.tiff", units="mm", height=200,
       width=360, scale=0.5)