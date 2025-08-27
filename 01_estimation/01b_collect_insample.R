rm(list=ls())

library(reshape2)
library(dplyr)
library(tidyr)

save.dir <- "insample"
file_list <- list.files(save.dir)

# storage
Zlat_ls <- list()
Zann_ls <- list()

run <- 1
for(run in 1:length(file_list)){
  load(paste0(save.dir,"/",file_list[run]))
  tmp <- gsub("\\.rda","",unlist(strsplit(file_list[run],"_")))
  
  mod <- tmp[2]
  cn <- tmp[3]
  qf_num <- tmp[4]
  qf <- ifelse(gsub("Qf","",tmp[4]) == 1, "single", "multi")
  csr <- gsub("csr","",tmp[5])
  me <- gsub("csr","",tmp[6])
  
  Zlat_df <- melt(out$Z_lat) %>%
    rename("moment" = "Var1","date" = "Var2", "nuts" = "Var3") %>%
    mutate(moment = paste0("p",gsub("%","",moment))) %>%
    pivot_wider(names_from = moment, values_from = value) %>%
    mutate(model = mod, 
           cn = cn, 
           qf = qf,
           qf_num = qf_num,
           csr = csr, 
           me = me)
  
  Zann_df <- melt(out$Z_ann) %>%
    rename("moment" = "Var1","date" = "Var2", "nuts" = "Var3") %>%
    mutate(moment = paste0("p",gsub("%","",moment))) %>%
    pivot_wider(names_from = moment, values_from = value)
  tmp <- melt(out$Y) %>%
    rename("date" = "Var1", "nuts" = "Var2")
  Zann_df <- left_join(Zann_df, tmp, by = c("date","nuts")) %>%
    mutate(model = mod, 
           cn = cn, 
           qf = qf, 
           qf_num = qf_num,
           csr = csr, 
           me = me)
  
  # -------------------------------
  # storage
  Zlat_ls[[run]] <- Zlat_df
  Zann_ls[[run]] <- Zann_df
}

Zlat_df <- do.call("bind_rows",Zlat_ls) %>%
  mutate(uid = paste0(cn,"_",model,"_",qf_num,"_",csr,"_",me))
Zann_df <- do.call("bind_rows",Zann_ls) %>%
  mutate(uid = paste0(cn,"_",model,"_",qf_num,"_",csr,"_",me))

save(Zlat_df, Zann_df, file = "insample_df.rda")


