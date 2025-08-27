rm(list = ls())

## SET WORKING DIRECTORY
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

save.dir <- "outputs"
dir.create(save.dir, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(xtable)
library(colorspace)
library(eurostat)
library(sf)
library(tmap)

myggtheme <- theme(panel.grid.major.x = element_blank() ,
                   panel.grid.major.y = element_line( size=.01, color="lightgrey"), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "lightgrey"),
                   legend.position = "bottom")

# PRELIMINARY ACTIONS for MAPS
if (file.exists("auxiliary/geodata_NUTS2.rda")){
  load("auxiliary/geodata_NUTS2.rda")
} else {
  geodata_NUTS2 <- get_eurostat_geospatial(nuts_level = 2, year = 2021, cache = TRUE, update_cache = TRUE)
  save(geodata_NUTS2, file="auxiliary/geodata_NUTS2.rda")
}

## 1. MAP PERFROMANCE ------------------------------------------------------------
load("auxiliary/best_CRPS_avg.rda")

df_map <- collect_sub_avg %>%
  filter(horizon == 0) %>%
  ungroup() %>%
  select(nuts, value_rel)

map_data <- left_join(geodata_NUTS2, df_map, by = c("geo"="nuts"))

ggplot(map_data) +
  geom_sf(aes(fill=value_rel), color="grey", size=0.01) +
  coord_sf(xlim=c(-11,30), ylim=c(35,67)) +
  colorspace::scale_fill_continuous_diverging(palette="Blue-Red", mid=1, na.value = "white", name="") +
  scale_colour_manual(values=NA) + 
  theme_void() +
  guides(colour = "none") +
  theme(legend.position = "bottom") 

ggsave(paste0(save.dir, "/map_CRPS_MFREGvsHM.pdf"), height=8, width=10)

## 2. MAP COVID ----------------------------------------------------------------
## data
load("auxiliary/best_preds_lat_fullsample.rda")

Y <- Zlat_sub %>%
  mutate(date = as.Date(as.character(date))) %>%
  select(date, nuts, p50) %>%
  filter(date %in% as.Date(c("2019-10-01", "2020-01-01", "2020-04-01", "2020-07-01"))) %>%
  mutate(date = zoo::as.yearqtr(date))

eurMapDf2020 <- merge(geodata_NUTS2, unique(Y$date)) %>%
  rename(date=y)
eurMapDataDf <- left_join(eurMapDf2020, Y, by = c("geo" = "nuts", "date"))

eurMapDataDf %>%
  mutate(date = as.character(date)) %>%
  ggplot()+
  geom_sf(aes(fill=p50, colour="lightgrey"), size=0.01) +
  facet_wrap(~date, ncol=4) +
  coord_sf(xlim=c(-11,30), ylim=c(35,67)) +
  colorspace::scale_fill_continuous_diverging(palette="Blue-Red 3", mid=0, na.value ="white", name="GVA growth\n", rev=T, p1=0.8, p2=0.8) +
  scale_colour_manual(values=NA) + 
  theme_void() +
  guides(colour = "none") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(paste0(save.dir, "/map_GVA_20192020_wide.pdf"), height=3, width=8)



