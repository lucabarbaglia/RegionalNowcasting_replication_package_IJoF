rm(list = ls())

## SET WORKING DIRECTORY
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

save.dir <- "outputs"
dir.create(save.dir, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(zoo)
library(lubridate)

load("01_estimation/mf_data_2025-05.rda") # original data
load("01_estimation/insample_df.rda") # nowcasted quarterly data
Zann_df$qf <- factor(Zann_df$qf, levels = c("single","multi"))
countries <- sort(unique(Zann_df$cn))

sl_dates <- as.character(seq(as.Date("2004-01-01"), as.Date("2100-12-01"), by = "month"))

# ----------------------------------------------------------------------------------------------
# annualized measures
# ----------------------------------------------------------------------------------------------
Zann_all <- list()
for(cc in 1:length(countries)){
  Zann_all[[cc]] <- Zann_df %>%
    subset(date %in% sl_dates) %>%
    # subset(qf %in% "multi") %>%
    subset(cn %in% countries[cc]) %>%
    
    ggplot(aes(x = as.Date(date), group = nuts, y = p50)) +
    geom_line(alpha = 0.35) +
    geom_point(aes(y = value), size = 0.3, shape = 1) +
    
    facet_nested(.~cn + qf) + 
    
    ylab("") + xlab("") + 
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_minimal() + 
    theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(0.2, "cm"),
          plot.margin = margin(0, 0, -0.25, 0, "cm"),
          strip.background = element_rect(fill = "grey95", linewidth = 0, color = "grey95"),
          axis.text.x = element_text(angle = 0))  
}

# outputs Figure 4 (insample_all.pdf; main paper)
pdf(file = paste0(save.dir,"/insample_all.pdf"), width = 10, height = 8.5)
print(plot_grid(plotlist = Zann_all, ncol = 3))
dev.off()

tmp_all <- Zann_df %>% subset(!(cn %in% "SK")) %>% subset(qf %in% "multi") %>% select(!qf)
tmp_SK <- Zann_df %>% subset(cn %in% "SK") %>% select(!qf)

Zann_sub <- bind_rows(tmp_all,tmp_SK)
save(Zann_sub, file = "auxiliary/best_preds_fullsample.rda")

# ----------------------------------------------------------------------------------------------
# quarterly latent states
# ----------------------------------------------------------------------------------------------
Zlat_df$qf <- factor(Zlat_df$qf, levels = c("single","multi"))

Zlat_all <- list()
for(cc in 1:length(countries)){
  Zlat_all[[cc]] <- Zlat_df %>%
    subset(date %in% sl_dates) %>%
    subset(cn %in% countries[cc]) %>%
    
    ggplot(aes(x = as.Date(date), group = nuts, y = p50)) +
    geom_line(alpha = 0.35) +

    facet_nested(.~cn + qf) + 
    
    ylab("") + xlab("") + 
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_minimal() + 
    theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(0.2, "cm"),
          plot.margin = margin(0, 0, -0.25, 0, "cm"),
          strip.background = element_rect(fill = "grey95", linewidth = 0, color = "grey95"),
          axis.text.x = element_text(angle = 0))  
}

pdf(file = paste0(save.dir,"/insample_all_lat.pdf"), width = 10, height = 8.5)
print(plot_grid(plotlist = Zlat_all, ncol = 3))
dev.off()

tmp_all <- Zlat_df %>% subset(!(cn %in% "SK")) %>% subset(qf %in% "multi") %>% select(!qf)
tmp_SK <- Zlat_df %>% subset(cn %in% "SK") %>% select(!qf)

Zlat_sub <- bind_rows(tmp_all,tmp_SK)
save(Zlat_sub, file = "auxiliary/best_preds_lat_fullsample.rda")


