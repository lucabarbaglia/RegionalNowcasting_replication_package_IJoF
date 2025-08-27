rm(list=ls())

library(Matrix)
library(zoo)
library(forecast)
library(scoringRules)

library(dplyr)
library(tidyr)
library(reshape2)

## DATA LOAD -----
load("mf_data_2025-05.rda")
source("aux_func.R")
save.dir <- "results"
error.dir <- "errors"

# algorithm setup
P <- 7
set.conf <- c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)
sl_quants <- c("16%","50%","84%")

set.seed(35256)

# set grid for the forecast exercise
countries <- names(gva_ls)
countries <- countries[countries != "PT"]

grid_h <- 8

csr_grid <- c("tight","loose","none")
me_grid <- c("tight","loose")

run.grid_ls <- list()
for (cc in countries){
  nuts2 <- colnames(gva_ls[[cc]]$reg)
  nuts2 <- nuts2[!grepl("ZZ",nuts2)]
  
  nuts2raw <- colnames(gva_ls[[cc]]$regraw)
  nuts2raw <- nuts2[!grepl("ZZ",nuts2raw)]
  
  qf_max <- min(length(nuts2raw),12)
  if(qf_max == 1){
    qf_seq <- 1
  }else{
    if(qf_max <= 12){
      qf_seq <- c(1,seq(2, qf_max, by = 1))
    }else{
      qf_seq <- c(1,seq(2, 10, by = 2), seq(12, qf_max, by=4))
    }  
  }
  
  run.grid.mfreg <- expand.grid("h"=1:grid_h,
                                "model"="mfreg",
                                "sv"=c(FALSE,TRUE),
                                "country"=cc,
                                "qf"= qf_seq,
                                "csr"=csr_grid,
                                "me"=me_grid,
                                stringsAsFactors = FALSE)
  run.grid.other <- expand.grid("h"=1:grid_h,
                                "model"=c("ar","rw","hm"),
                                "sv" = FALSE,
                                "country"=cc,
                                "qf"= NA,
                                "csr"=NA,
                                "me"=NA,
                                stringsAsFactors = FALSE)
  run.grid_ls[[cc]] <- rbind(run.grid.mfreg,
                             run.grid.other)
}
run.grid <- do.call("rbind",run.grid_ls)
rownames(run.grid) <- 1:nrow(run.grid)

states <- "ffbs"
run <- 1
collect_ls <- list()
preds_ls <- list()
diag_ls <- list()
missing_ls <- list()

for(states in c("ffbs")){
  for(run in 1:nrow(run.grid)){
    # set grid for the forecast exercise
    run.spec <- run.grid[run,]
    run.h <- run.spec$h
    run.model <- run.spec$model
    run.sv <- run.spec$sv
    run.cc <- run.spec$country
    run.qf <- run.spec$qf
    run.csr <- run.spec$csr
    run.me <- run.spec$me
    
    if(run.model != "mfreg"){
      states_run <- NA
    }else{
      states_run <- states
    }
    spec <- paste0("h",run.h,
                   "_",run.model,"-",ifelse(run.sv,"sv","hom"),"-",states_run,
                   "_",run.cc,
                   "_Qf",run.qf,
                   "_csr",run.csr,
                   "_me",run.me)
    file.dir <- paste0(save.dir,"/",spec,".rda")
    logs.dir <- paste0(error.dir,"/",run,"-",spec,".txt")
    
    if(file.exists(file.dir)){
      load(file.dir)
    }else{
      if(file.exists(logs.dir)){
        spec_error <- readLines(logs.dir)[1]
      }else{
        spec_error <- NA
      }
      missing_ls[[spec]] <- cbind(run,spec,spec_error)
      next
    }
    
    # ---
    metr_ls <- list()
    for(metr in c("crps","lps","error")){
      tmp <- melt(out$fcmetr[[metr]]) %>% rename("type"="Var1","nuts"="Var2") %>%
        separate(type, into=c("horizon","date"), sep="_") %>%
        mutate(metric = metr, horizon = ifelse(horizon=="nc",0,-1)) %>%
        mutate("holdout" = run.h,
               "model" = paste0(run.model,"-",ifelse(run.sv,"sv","hom")),
               "states" = states_run,
               "country" = run.cc,
               "qf" = run.qf,
               "csr" = run.csr,
               "me" = run.me)
      vintage <- tmp %>% subset(horizon == 0) %>% select(date) %>% unique() %>% as.character()
      metr_ls[[metr]] <- tmp %>% 
        mutate("vintage" = vintage) %>%
        select(vintage,holdout,date,horizon,country,nuts,model,states,qf,csr,me,metric,value)
    }
    metr_df <- do.call("bind_rows",metr_ls)
    collect_ls[[spec]] <- metr_df
    
    # ---
    if(!is.null(out$Za)){
      tmp <- melt(out$Za[sl_quants,,]) %>% na.exclude() %>%
        rename("moment"="Var1","date"="Var2","nuts"="Var3") %>%
        mutate("holdout" = run.h,
               "model" = paste0(run.model,"-",ifelse(run.sv,"sv","hom")),
               "states" = states_run,
               "country" = run.cc,
               "qf" = run.qf,
               "csr" = run.csr,
               "me" = run.me,
               "vintage" = vintage) %>%
        select(vintage,holdout,date,country,nuts,model,states,qf,csr,me,moment,value)
      levels(tmp$moment) <- paste0("p",gsub("%","",levels(tmp$moment)))
      preds_ls[[spec]] <- tmp
    }
    
    if(!is.null(out$sediag)){
      # state equation diagnostics
      tmp <- melt(out$sediag[sl_quants,]) %>%
        rename("moment"="Var1","type"="Var2") %>%
        mutate("holdout" = run.h,
               "model" = paste0(run.model,"-",ifelse(run.sv,"sv","hom")),
               "states" = states_run,
               "country" = run.cc,
               "nuts" = NA,
               "qf" = run.qf,
               "csr" = run.csr,
               "me" = run.me,
               "vintage" = vintage) %>%
        select(vintage,holdout,country,nuts,model,states,qf,csr,me,moment,type,value)
      
      # regional unconditional moments
      tmp2 <- melt(out$uncm) %>%
        rename("moment"="Var1","type"="Var2","nuts"="Var3") %>%
        mutate("holdout" = run.h,
               "model" = paste0(run.model,"-",ifelse(run.sv,"sv","hom")),
               "states" = states_run,
               "country" = run.cc,
               "qf" = run.qf,
               "csr" = run.csr,
               "me" = run.me,
               "vintage" = vintage) %>%
        select(vintage,holdout,country,nuts,model,states,qf,csr,me,moment,type,value)
      tmp <- bind_rows(tmp,tmp2)
      levels(tmp$moment) <- paste0("p",gsub("%","",levels(tmp$moment)))
      
      diag_ls[[spec]] <- tmp
    }
  }
}

if(length(missing_ls)!=0){
  missings <- do.call("rbind",missing_ls) %>% 
    as.data.frame() %>% 
    separate(spec, into = c("horizon","model","var","states","country","facs","csr","me"))  
}else{
  missings <- NA
}

collect_df <- do.call("bind_rows",collect_ls)
save(collect_df, file = "collect_df.rda")

preds_df <- do.call("bind_rows",preds_ls)
save(preds_df, file = "preds_df.rda")

diag_df <- do.call("bind_rows",diag_ls)
save(diag_df, file = "diag_df.rda")


