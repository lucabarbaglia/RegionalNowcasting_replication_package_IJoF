rm(list = ls())

## SET WORKING DIRECTORY
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

save.dir <- "outputs"
dir.create(save.dir, showWarnings = FALSE)

library(forecast)
library(zoo)
library(dplyr)
library(tidyr)
library(reshape2)
library(R.utils)
library(lubridate)
library(xtable)
library(readxl)
library(ggplot2)
library(ggh4x)
library(lemon)
library(scales)
library(grid)
library(gridExtra) 
library(cowplot)

outsig <- function(x,digs,lab){
  paste0(format(round(x,digs),nsmall=digs),lab)
}

col_blue <- "navy"
col_red <- "firebrick"

mode <- "weights"

set.digits <- 2
tile_alpha <- 0.3

## DATA LOAD -----
load("01_estimation/collect_df.rda")
load("01_estimation/mf_data_2025-05.rda")

wghs_ls <- list()
for(cc in 1:length(gva_ls)){
  tmp <- gva_ls[[cc]]$regraw
  tmp <- tmp[,!grepl("ZZ",colnames(tmp))]
  if(NCOL(tmp) == 1) next
  
  labs <- colnames(tmp)
  tmp_ls <- data.frame("nuts" = labs, "weights" = apply(tmp / apply(tmp, 1, sum),2,mean))
  if(mode == "raw"){
    tmp_ls$weights <- 1/length(tmp_ls$weights)
  }
  
  wghs_ls[[cc]] <- tmp_ls
}
wghs_df <- do.call("bind_rows", wghs_ls)

collect_df <- collect_df %>%
  subset(metric %in% c("crps","error"))
collect_avg <- collect_df %>%
  group_by(horizon,country,nuts,model,states,qf,csr,me,metric) %>%
  left_join(wghs_df, by = "nuts") %>%
  mutate(value_raw = value) %>%
  mutate(value = value * weights)

collect_avg$csr <- factor(collect_avg$csr,levels=c("none","loose","tight"))
levels(collect_avg$csr) <- paste0("CSR:", c("no","loose","tight"))
collect_avg$me <- factor(collect_avg$me, levels = c("loose","tight"))
levels(collect_avg$me) <- paste0("ME:", c("loose","tight"))
collect_avg$model <- factor(gsub("mfreg-","",collect_avg$model))
levels(collect_avg$model) <- c("AR","HM","MF-DFM-hom","RW","MF-DFM-SV")

collect_avg$uid <- paste0(collect_avg$country,"_h",collect_avg$horizon,"_",collect_avg$model,"_",collect_avg$qf,"_",collect_avg$csr,"_",collect_avg$me)

# --------------------------------------------------------------------------------------------------
# main results tables, averaged for countries
# --------------------------------------------------------------------------------------------------
bench_df <- collect_avg %>%
  subset(model %in% "HM") %>%
  rename("bench" = "value", "bench_raw" = "value_raw") %>% ungroup() %>%
  select(c("horizon","country","nuts","metric", "bench_raw","bench"))

collect_bench <- left_join(collect_avg, bench_df, 
                           by = c("horizon","country","nuts","metric")) %>% ungroup() %>%
  group_by(horizon, country, model, states, qf, csr, me, metric) %>% 
  summarize(value = sum(value), bench = sum(bench), 
            value_raw = mean(value_raw), bench_raw = mean(bench_raw)) %>%
  subset(states %in% "ffbs") %>%
  mutate(value_rel = value / bench, value_raw_rel = value_raw / bench_raw)

collect_bench <- collect_bench %>% ungroup() %>%
  group_by(horizon,country,metric) %>%
  
  mutate(best = ifelse(value == min(value),outsig(value_rel,digs = set.digits,""),"")) %>%
  mutate(print = ifelse(best == "", outsig(value_rel,digs = set.digits,""),"")) %>%
  mutate(print = ifelse(best != "" & print != "","",print)) %>%
  
  mutate(fill = value_rel)

collect_bench$qf <- factor(collect_bench$qf, levels = c(12:1))

plot_ls <- list()

# output Table 1 (CRPS_h0_csrtighttable.pdf; main paper)
# output Table A.1--A.3 (CRPS_h0_csrloosetable.pdf, CRPS_h0_csrnotable.pdf, CRPS_h-1_csrtighttable.pdf; appendix), plus additional tables that were not included in the paper 
for(hh in c(-1,0)){
  for(icsr in c("CSR:no","CSR:tight","CSR:loose")){
    plot_ls <- collect_bench %>%
      subset(csr %in% icsr) %>%
      subset(horizon == hh) %>%
      subset(metric %in% "crps") %>%
      
      ggplot() +
      
      geom_tile(aes(x=country,y=qf,fill=fill),alpha=tile_alpha) +
      geom_text(aes(x=country,y=qf,label=best),size=3,fontface="bold") +
      geom_text(aes(x=country,y=qf,label=print),size=3) +
      
      xlab("") + ylab("Factors (Q)") +
      facet_nested(rows = vars(model,csr, me),
                   scales = "free_y",
                   switch = "y",
                   space = "free_y") +  # allows varying width
      
      scale_fill_gradient2(midpoint=1,low=col_blue,high=col_red, mid = "grey98",na.value = "grey80", limits = c(0,1.25), oob=squish) +
      scale_color_gradient2(midpoint=1,low=col_blue,high=col_red, mid = "grey98",na.value = "grey80", limits = c(0,1.25), oob=squish) +
      scale_x_discrete(position = "top") +
      
      coord_cartesian(expand = FALSE) +
      theme_minimal() +
      theme(strip.placement = "outside",
            strip.background.y = element_rect(fill = "grey90", linewidth = 0, color = "grey90"),
            legend.position="none",legend.key.width = unit(1.5,"cm"),legend.key.height = unit(0.1,"cm"),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(hjust=1),
            axis.text.x = element_text(angle=0, vjust = 0),
            axis.line=element_blank(),
            plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
    
    # save outputs mentioned above
    pdf(file = paste0(save.dir,"/","CRPS_h",hh,"_",gsub(":","",icsr),"table.pdf"), width = 8, height = 6)
    print(plot_ls)
    dev.off()
  }
}

sl_best <- collect_bench %>% 
  subset(metric %in% "crps") %>%
  subset(best != "")
sl_best$uid <- paste0(sl_best$country,"_h",sl_best$horizon,"_",sl_best$model,"_",sl_best$qf,"_",sl_best$csr,"_",sl_best$me)
sl_best$uid2 <- paste0(sl_best$country,"_h",sl_best$horizon,"_",sl_best$model,"_",sl_best$csr,"_",sl_best$me)
sl_best <- sl_best %>% subset(horizon == 0)

subset_regs <- sl_best$uid
subset2_regs <- sl_best$uid2
save(sl_best, file = "01_estimation/bestmodels.rda")

# --------------------------------------------------------------------------------------------------
# boxplots: over the cross-section and time, best performing models
# --------------------------------------------------------------------------------------------------
collect_sub <- collect_avg %>%
  subset(metric %in% "crps") %>%
  subset(uid %in% subset_regs)
bench_df <- collect_avg %>%
  subset(model %in% "HM") %>%
  rename("bench" = "value", "bench_raw" = "value_raw") %>% ungroup() %>%
  select(c("horizon","country","nuts","vintage","metric","bench_raw","bench"))

collect_sub <- left_join(collect_sub, bench_df, 
          by = c("horizon","vintage","country","nuts","metric")) %>% ungroup() %>%
  mutate(value_rel = value / bench, value_raw_rel = value_raw / bench_raw)

pp_boxplot <- collect_sub %>%
  ggplot() +
  geom_boxplot(aes(x = as.Date(vintage), y = value_raw_rel, group = vintage), 
               outlier.size = 1, outlier.shape = 1, linewidth = 0.3, fill = "grey60") +
  geom_hline(yintercept = c(1), linewidth = 0.75) +
  geom_hline(yintercept = c(0), linewidth = 0.3) +
  facet_wrap(.~country, ncol = 7, scales = "free") +
  xlab("") + ylab("CRPS (ratio to benchmark)") +
  coord_cartesian(expand = TRUE, ylim = c(0,2.5), clip = "off") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        panel.spacing.y = unit(0.5, "cm"), panel.spacing.x = unit(0.25,"cm"))

# outputs Figure 3a (CRPS_boxplot_years.pdf; main paper)
pdf(file = paste0(save.dir,"/","CRPS_boxplot_years.pdf"), width = 10, height = 3.5)
print(pp_boxplot)
dev.off()

save(collect_sub, file = "auxiliary/best_CRPS_years.rda")

pp_boxplot <- collect_sub %>%
  ggplot() +
  
  geom_point(aes(x = as.Date(vintage), y = value_raw, group = vintage, size = weights), 
              shape = 1, alpha = 0.4) +
  
  geom_hline(yintercept = c(0), linewidth = 0.3) +
  facet_wrap(.~country, ncol = 7, scales = "free") +
  xlab("") + ylab("CRPS (level)") +

  scale_size_continuous(range = c(0.1, 3)) + 

  coord_cartesian(expand = TRUE, clip = "off") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        panel.spacing.y = unit(0.5, "cm"), panel.spacing.x = unit(0.25,"cm"))

# outputs Figure 3b (CRPSraw_boxplot_years.pdf; main paper)
pdf(file = paste0(save.dir,"/","CRPSraw_boxplot_years.pdf"), width = 10, height = 3.5)
print(pp_boxplot)
dev.off()

# --------------------------------------------------------------------------------------------------
# boxplots: across number of factors, averaged over the holdout
# --------------------------------------------------------------------------------------------------
collect_avg$uid2 <- paste0(collect_avg$country,"_h",collect_avg$horizon,"_",
                           collect_avg$model,"_",collect_avg$csr,"_",collect_avg$me)

collect_sub <- collect_avg %>%
  subset(metric %in% "crps") %>%
  subset(uid2 %in% subset2_regs)
bench_df <- collect_avg %>%
  subset(model %in% "HM") %>%
  rename("bench" = "value", "bench_raw" = "value_raw") %>% ungroup() %>%
  select(c("horizon","country","nuts","vintage","metric","bench_raw","bench"))

collect_sub <- left_join(collect_sub, bench_df, 
                         by = c("horizon","vintage","country","nuts","metric")) %>% ungroup() %>%
  group_by(horizon,uid,uid2,country,nuts,model,states,qf,csr,me,metric) %>%
  summarize(value = mean(value), bench = mean(bench), 
            value_raw = mean(value_raw), bench_raw = mean(bench_raw),
            weights = mean(weights)) %>%
  mutate(value_rel = value / bench, value_raw_rel = value_raw / bench_raw)
collect_sub$best <- 0
collect_sub$best[collect_sub$uid %in% subset_regs] <- 1
collect_sub$best <- factor(collect_sub$best)
collect_sub$qf <- factor(collect_sub$qf)

cn_num <- 5
sl_cn <- list("g1" = unique(collect_sub$country)[1:cn_num],
              "g2" = unique(collect_sub$country)[(cn_num+1):length(unique(collect_sub$country))])

pp_box_ls <- list()
pp_raw_ls <- list()

collect_sub_avg <- collect_sub %>%
  subset(best == 1)
save(collect_sub_avg, file = "auxiliary/best_CRPS_avg.rda")

for(i in 1:length(sl_cn)){
  pp_box_ls[[i]] <- collect_sub %>%
    subset(country %in% sl_cn[[i]]) %>%
    ggplot() +
    geom_boxplot(aes(x = qf, y = value_raw_rel, group = qf, fill = best),
                 outlier.size = 1, outlier.shape = 1, linewidth = 0.3) +
    
    geom_hline(yintercept = c(1), linewidth = 0.75) +
    geom_hline(yintercept = c(0), linewidth = 0.3) +
    facet_grid(.~country, scales = "free", space = "free") +
    
    scale_fill_manual(values = c("white","grey60")) +
    
    xlab("Factors (Q)") + ylab("") +
    coord_cartesian(expand = TRUE, clip = "off", ylim = c(0,1.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          strip.placement = "outside",
          panel.spacing.y = unit(0.25, "cm"), panel.spacing.x = unit(0.5,"cm"))
  
  pp_raw_ls[[i]] <- collect_sub %>%
    subset(country %in% sl_cn[[i]]) %>%
    ggplot() +
    geom_point(aes(x = qf, y = value_raw, group = qf, 
                   color = best, size = weights), 
               shape = 1, alpha = 0.4) +
    facet_grid(.~country, scales = "free", space = "free") +
    geom_hline(yintercept = c(0), linewidth = 0.3) +
    
    scale_color_manual(values = c("grey60","black")) +
    scale_size_continuous(range = c(0.1, 3)) + 

    xlab("Factors (Q)") + ylab("") +
    coord_cartesian(expand = TRUE, clip = "off", ylim = c(0,0.065)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          strip.placement = "outside",
          panel.spacing.y = unit(0.25, "cm"), panel.spacing.x = unit(0.5,"cm"))
}

pp_box <- plot_grid(plotlist = pp_box_ls,
                           ncol = 1,           # or ncol = 2 for horizontal
                           align = "v",
                           axis = "l")         # keep left axis line
pp_box <- ggdraw() +
  draw_label("CRPS (ratio to benchmark)", x = 0.015, y = 0.5, angle = 90, vjust = 0.5, size = 11) +
  draw_plot(pp_box, x = 0.01, y = 0, width = 0.98, height = 1)

# outputs Figure 1a (CRPS_boxplot_factors.pdf; main paper)
pdf(file = paste0(save.dir,"/","CRPS_boxplot_factors.pdf"), width = 10, height = 3.5)
print(pp_box)
dev.off()

pp_raw <- plot_grid(plotlist = pp_raw_ls,
                    ncol = 1,           # or ncol = 2 for horizontal
                    align = "v",
                    axis = "l")         # keep left axis line
pp_raw <- ggdraw() +
  draw_label("CRPS (level)", x = 0.015, y = 0.5, angle = 90, vjust = 0.5, size = 11) +
  draw_plot(pp_raw, x = 0.01, y = 0, width = 0.98, height = 1)

# outputs Figure 1b (CRPSraw_boxplot_factors.pdf; main paper)
pdf(file = paste0(save.dir,"/","CRPSraw_boxplot_factors.pdf"), width = 10, height = 3.5)
print(pp_raw)
dev.off()

# --------------------------------------------------------------------------------------------------
# plots of how data updates the latent states
# --------------------------------------------------------------------------------------------------
load("01_estimation/insample_df.rda")

# load("01_estimation/preds_df.rda")
# sl_subs <- unique(Zann_df$uid)
# sl_subs <- c(sl_subs[!grepl("Qf1_",sl_subs)],sl_subs[grepl("SK_",sl_subs)])
# 
# tmp <- preds_df %>%
#   pivot_wider(names_from = moment, values_from = value) %>%
#   mutate(vintage = substr(vintage,1,4)) %>%
#   mutate(uid = paste0(country,"_",model,"-",states,"_Qf",qf,"_",csr,"_me",me)) %>%
#   subset(uid %in% sl_subs)
# tmp2 <- Zann_df %>% 
#   select(date, nuts, uid, p16, p50, p84, value) %>%
#   rename("p16fs" = "p16", "p50fs" = "p50", "p84fs" = "p84", "true" = "value")
# preds_sub <- tmp %>% left_join(tmp2, by = c("date" = "date", "nuts" = "nuts", "uid" = "uid"))
# preds_sub$date <- as.character(as.Date(as.character(preds_sub$date)) %m+% months(2))
# 
# sl_cn <- c("DE","IT","FR","ES")
# sl_nuts_big <- sl_nuts_median <- c()
# for(cc in sl_cn){
#   tmp <- collect_sub %>% subset(country %in% cc) %>% subset(best %in% 1)
#   sl_nuts_big <- c(sl_nuts_big, tmp$nuts[which(tmp$weights == max(tmp$weights))])
#   sl_nuts_median <- c(sl_nuts_median, tmp$nuts[which.min(abs(tmp$weights - median(tmp$weights)))])
# }
# 
# sl_nuts <- c(sl_nuts_big, sl_nuts_median)
# sl_start <- c("2019-12-01")
# 
# preds_final <- preds_sub %>% subset(vintage == "2023") %>%
#   subset(nuts %in% sl_nuts) %>%
#   subset(date %in% seq(as.Date(sl_start),as.Date("2023-12-01"),by = "month")) %>%
#   mutate(vintage = "Full sample") %>%
#   mutate("nuts_type" = ifelse(nuts %in% sl_nuts_big, "Largest", "Median"))
# preds_plot <- preds_sub %>%
#   subset(nuts %in% sl_nuts) %>%
#   subset(date %in% seq(as.Date(sl_start),as.Date("2023-12-01"),by = "month")) %>%
#   mutate(vintage_date = paste0(vintage,"-12-01")) %>%
#   mutate(type = ifelse(vintage_date == date, "Nowcast", 
#                        ifelse(as.character(as.Date(vintage_date) %m-% months(12)) == date, "Backcast", ""))) %>%
#   mutate("nuts_type" = ifelse(nuts %in% sl_nuts_big, "Largest", "Median"))
# 
# preds_plot$type <- factor(preds_plot$type, levels = c("Backcast","Nowcast",""))
# save(preds_plot, preds_final, file = "01_estimation/preds_sub.rda")

# load collected data files
load("01_estimation/preds_sub.rda")

vintage_cols <- c("#e7298a", "#33a02c","#1f78b4","#d95f02")
pp_vintage <- preds_plot %>%
  subset(!(vintage %in% "2019")) %>%
  
  ggplot(aes(x = as.Date(date), y = p50, color = vintage)) +
  
  geom_ribbon(aes(ymin = p16fs, ymax = p84fs),
              fill = "black", color = "black", alpha = 0, linewidth = 0.5, data = preds_final, linetype = "dotdash") +
  geom_ribbon(aes(ymin = p16, ymax = p84, fill = vintage), alpha = 0.2, linewidth = NA) +
  geom_text(aes(label = nuts), x = as.Date("2021-01-01"), y = 0.1, stat = "identity", color = "black", size = 4) +
  
  geom_line() + 
  geom_point(aes(y = p50, color = vintage, shape = type), size = 1.5) +
  geom_point(aes(y = true), color = "black", size = 2, shape = 4) +
  
  ylab("GVA") + xlab("Date") +
  scale_shape_manual(values = c(1,19,NA), name = "") +
  scale_color_manual(values = vintage_cols, name = "Vintage") +
  scale_fill_manual(values = vintage_cols, name = "Vintage") +
  
  geom_hline(yintercept = 0) +
  facet_wrap(nuts_type ~ country, switch = "y", ncol = 4) +
  coord_cartesian(clip = "off", expand = FALSE) +
   
  theme_minimal() +
  theme(strip.placement = "outside", 
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, vjust = 0),
        panel.grid.minor = element_blank())

# outputs Figure 6 (nuts_realtime.pdf; main paper)
pdf(file = paste0(save.dir,"/nuts_realtime.pdf"), width = 8, height = 4.5)
print(pp_vintage)
dev.off()

# --------------------------------------------------------------------------------------------------
# plots VAR state equation diagnostics
# --------------------------------------------------------------------------------------------------
load("01_estimation/diag_df.rda")
sl_diag <- gsub("-SV","-hom",sl_best$uid)

diag_df <- diag_df %>% 
  mutate(model = gsub("-sv","-SV",gsub("mfreg","MF-DFM",model))) %>%
  mutate(uid = paste0(country,"_h0_",model,"_",qf,"_CSR:",csr,"_ME:",me))
  
diag_logdet <- diag_df %>% 
  subset(type %in% c("logdetSig")) %>%
  subset(uid %in% sl_diag) %>%
  pivot_wider(names_from = moment, values_from = value)

pp_logdet <- diag_logdet %>%
  ggplot(aes(x = as.Date(vintage))) +
  geom_line(aes(y = p50)) +
  xlab("") + ylab(expression(log(det(H)))) +
  geom_ribbon(aes(ymin = p16, ymax = p84), alpha = 0.15) +
  facet_wrap(~country, scales = "free", ncol = 7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

# plot the eigenvalues
n_theta <- 400
circle <- data.frame(
  x = cos(seq(0, 2*pi, length.out = n_theta)),
  y = sin(seq(0, 2*pi, length.out = n_theta))
)

diag_eigen <- diag_df %>% 
  subset(type %in% c("maxeigReA","maxeigImA")) %>%
  subset(uid %in% sl_diag) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  select(country,vintage,moment,maxeigReA,maxeigImA) %>%
  melt(id = c("country","moment","vintage")) %>%
  pivot_wider(names_from = c(variable, moment), values_from = value)

disc_cols <- c("#4eb3d3","#2b8cbe","#0868ac","#084081","#fd8d3c","#e31a1c","#bd0026","#800026")

pp_eigen <- diag_eigen %>%
  ggplot() +
  geom_path(data = circle, aes(x, y), linewidth = .5) +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  
  geom_point(aes(x = maxeigReA_p50, y = maxeigImA_p50, color = vintage), size = 1.5, alpha = 0.5) +
  
  scale_color_manual(values = disc_cols, name = "Max. eigenvalue\nVintage:") +
  scale_x_continuous(breaks = c(-1,0,1)) + scale_y_continuous(breaks = c(-1,0,1)) +
  
  facet_rep_wrap(~country, ncol = 7) + xlab("Real") + ylab("Imaginary") +
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = FALSE, clip = "off") +
  theme_minimal() + theme(legend.position = "bottom")

# outputs Figure 7a (paras_eigen.pdf, main paper)
pdf(file = paste0(save.dir,"/paras_eigen.pdf"), width = 8, height = 4)
print(pp_eigen)
dev.off()

# outputs Figure 7b (paras_logdet.pdf, main paper)
pdf(file = paste0(save.dir,"/paras_logdet.pdf"), width = 8, height = 3)
print(pp_logdet)
dev.off()


# --------------------------------------------------------------------------------------------------
# plots unconditional variance of regions
# --------------------------------------------------------------------------------------------------
tmp <- diag_df %>% 
  subset(type %in% c("var")) %>%
  subset(uid %in% sl_diag) %>%
  pivot_wider(names_from = moment, values_from = value) %>%
  mutate(vintage = substr(vintage, 1,4))

pp_ucv <- tmp %>%
  subset(vintage %in% as.character(2019:2023)) %>%
  ggplot() +
  geom_boxplot(aes(x = log(p50), y = vintage), outliers = FALSE) +
  facet_wrap(~country, scales = "free_x", ncol = 5) +
  xlab("Log (unconditional variance)") + ylab("Vintage") +
  coord_cartesian(expand = TRUE, clip = "off") +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.5, "cm"),
        plot.margin = unit(c(0,0.5,0,0.1),"cm"))

# outputs Figure A.2 (paras_uncondvar.pdf; appendix)
pdf(file = paste0(save.dir,"/paras_uncondvar.pdf"), width = 9, height = 3.5)
print(pp_ucv)
dev.off()
