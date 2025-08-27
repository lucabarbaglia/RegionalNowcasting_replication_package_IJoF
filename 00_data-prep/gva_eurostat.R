rm(list=ls())

# ------- SET WORKING DIRECTORY ------- #
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# ------- PACKAGES ------- #
library(eurostat)
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)

# ------- DATA ------- #
gva_raw <- get_eurostat(tolower("NAMA_10R_2GVAGR"))
gdp_raw <- get_eurostat(tolower("NAMA_10R_2GDP"))
colnames(gdp_raw) <- c("freq","unit","geo","time","values")
colnames(gva_raw) <- c("freq","na_item","unit","geo","time","values")

gva_reg <- gva_raw %>%
  subset(unit %in% "I15") %>%
  subset(na_item %in% "B1GQ") %>% ## Gross domestic product at market prices [B1GQ]
  select(geo,time,values)
gdp_reg <- gdp_raw %>%
  subset(unit %in% "MIO_EUR") %>%
  select(geo,time,values)

gva_cnraw <- get_eurostat(tolower("NAMQ_10_A10"))
gva_cn <- gva_cnraw %>% 
  subset(unit %in% "CLV_I15") %>% 
  subset(nace_r2 %in% "TOTAL") %>% 
  subset(s_adj %in% "SCA") %>% rename("time" = "TIME_PERIOD") %>%
  select(geo,time,values)

# --------------------------------------------------------------------------------------------------------
# re-labeling
gva_cn <- gva_cn[nchar(gva_cn$geo)==2,]
gva_reg <- gva_reg[nchar(gva_reg$geo)==4,] # subset to nuts-2 regions only
gdp_reg <- gdp_reg[nchar(gdp_reg$geo)==4,] # subset to nuts-2 regions only

gva_reg$nuts <- gva_reg$geo
gva_reg$geo <- substr(gva_reg$nuts,1,2)
gdp_reg$nuts <- gdp_reg$geo
gdp_reg$geo <- substr(gdp_reg$nuts,1,2)

list.cn <- c("AT", "BE", "CY", "DE", "EE", "EL", "ES", "FI", "FR", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "SI", "SK") ## EA19 (excl. PT because of NUTS change)

gva_ls <- list()
for(cn in 1:length(list.cn)){
  gva_cn_sub <- gva_cn %>% 
    subset(geo %in% list.cn[cn]) %>%
    select(time,values) %>%
    arrange(time)
  gva_reg_sub <- gva_reg %>% 
    subset(geo %in% list.cn[cn]) %>% 
    select(time,nuts,values) %>%
    #spread(nuts, values, -time)
    pivot_wider(names_from = nuts, values_from = values) %>%
    arrange(time)
  gdp_reg_sub <- gdp_reg %>% 
    subset(geo %in% list.cn[cn]) %>% 
    select(time,nuts,values) %>%
    pivot_wider(names_from = nuts, values_from = values) %>%
    arrange(time)
  
  start_q <- as.numeric(unlist(strsplit(as.character(min(gva_cn_sub$time)),"-"))[1:2])
  start_y <- as.numeric(substr(min(gva_reg_sub$time),1,4))
  
  gdp_reg_sub <- data.frame("time" = gva_reg_sub$time) %>%
    left_join(gdp_reg_sub, by = "time")
  
  # exclude all countries with less than 2 NUTS regions
  if(ncol(gdp_reg_sub) < 3){
    next
  }
  
  if(any(is.na(gdp_reg_sub))){
    for(j in 2:ncol(gdp_reg_sub)) {
      tmp_gdp <- gdp_reg_sub[,j]
      if(any(is.na(tmp_gdp))){
        ind_1st <- min(which(!is.na(tmp_gdp)))
        tmp_gdp[1:((ind_1st)-1)] <- tmp_gdp[ind_1st]
        gdp_reg_sub[,j] <- tmp_gdp
      }else{
        next
      }
    }  
  }
  
  nat <- diff(log(ts(gva_cn_sub[,-1],start=start_q,frequency=4)))
  reg <- diff(log(ts(gva_reg_sub[,-1],start=start_y)))
  regraw <- ts(gdp_reg_sub[,-1],start=start_y)
  
  # NUTS change 2021-2024 (see xlsx)
  if(list.cn[cn] == "NL"){
    colnames(reg)[colnames(reg) == "NL35"] <- "NL31"
    colnames(reg)[colnames(reg) == "NL36"] <- "NL33"
    reg <- reg[,sort(colnames(reg))]
    
    colnames(regraw)[colnames(regraw) == "NL35"] <- "NL31"
    colnames(regraw)[colnames(regraw) == "NL36"] <- "NL33"
    regraw <- regraw[,sort(colnames(regraw))]
    
    # add missings from earlier vintage (they are not provided anymore by eurostat)
    # same script for mf_data_OLD.rda, run at different date (earlier "vintage")
    load("mf_data_OLD.rda") 
    reg[1:10,"NL31"] <- gva_ls_OLD[["NL"]]$reg[1:10,"NL31"]
    reg[1:10,"NL33"] <- gva_ls_OLD[["NL"]]$reg[1:10,"NL33"]
  }
  
  gva_ls[[list.cn[cn]]] <- list("nat"=nat,
                                "reg"=reg,
                                "regraw"=regraw)
}

save(gva_ls,file="mf_data_new.rda") # save current data-file
 
