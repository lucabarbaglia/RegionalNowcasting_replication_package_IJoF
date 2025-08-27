rm(list=ls())
run <- commandArgs(trailingOnly = TRUE)
# run <- 1 # select specification from grid if run locally 

## SET WORKING DIRECTORY
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(Matrix)
library(zoo)
library(forecast)
library(scoringRules)

## DATA LOAD -----
load("mf_data_2025-05.rda")
source("aux_func.R")
save.dir <- "results"
error.dir <- "errors"
dir.create(save.dir,showWarnings = FALSE)
dir.create(error.dir,showWarnings = FALSE)

# algorithm setup
nburn <- 1000                       # burn observations
nsave <- 3000                       # number of draws to be stored
nthin <- 1                          # thinning
toplot <- NULL                      # plot summary during sampling (each 'tplot' draw or NULL for no output)
check.stab <- TRUE                  # check for stationarity of the system
eigen.crit <- 1                     # criterion for stability of the system (try 1.05)
max.try <- 20                       # iterations of tries for stable draw
states <- "ffbs"                    # simulation algorithm for the states
ident_facs <- FALSE                 # identification of factors (upper triangular block of zeroes with unit diag)

# model setup
P <- 7
set.conf <- c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)

# set.seed(35256)

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

# set grid for the forecast exercise
run.spec <- run.grid[run,]
run.h <- run.spec$h
run.model <- run.spec$model
run.sv <- run.spec$sv
run.cc <- run.spec$country
run.qf <- run.spec$qf
run.csr <- run.spec$csr
run.me <- run.spec$me

if(run.model != "mfreg") states <- NA
spec <- paste0("h",run.h,
               "_",run.model,"-",ifelse(run.sv,"sv","hom"),"-",states,
               "_",run.cc,
               "_Qf",run.qf,
               "_csr",run.csr,
               "_me",run.me)
file.dir <- paste0(save.dir,"/",spec,".rda")
logs.dir <- paste0(error.dir,"/",run,"-",spec)
if(file.exists(file.dir)) stop("File already exists.")

# ----------------------------------------------------------------------
# select data
Q.lat <- run.qf              # number of latent factors
itr_scl <- 1                 # scale intertemporal restriction
reg_wghs <- TRUE             # set weights in regional average (if FALSE: equal weights)
sv <- run.sv                 # setting for SV, insample test
csr <- run.csr
me_prior <- run.me

start <- c(2002,1)
end.m <- c(2023,12)
end.q <- c(2023,4)
end.a <- c(2023,1)

# construct design matrices
data_nat <- window(ts.union(ts(NA,start=start,end=end.q,frequency=4),
                            window(gva_ls[[run.cc]]$nat,start=start))[,-1],end=end.q)
data_reg <- window(gva_ls[[as.character(run.cc)]]$reg,start=start, end = end.a)
data_regraw <- window(gva_ls[[as.character(run.cc)]]$regraw,start=start, end = end.a)

# exclude overseas french territories
if(run.cc == "FR"){
  data_reg <- data_reg[,!grepl("FRY",colnames(data_reg))]
  data_regraw <- data_regraw[,!grepl("FRY",colnames(data_regraw))]
}
# exclude "ZZ" nuts category
if (sum(substr(colnames(data_reg),3,4)=="ZZ")>0){
  data_reg <- data_reg[, -which(substr(colnames(data_reg),3,4)=="ZZ")]
}
data_regraw <- data_regraw[,colnames(data_reg)]
omega <- data_regraw/apply(data_regraw,1,sum) # obtain weights for cross-section

date.labs <- as.character(zoo::as.Date.ts(data_nat))
date.labs.p <- date.labs[(P+1):length(date.labs)]

YQ0 <- as.matrix(window(data_nat,start=start(data_reg)))
YA0 <- as.matrix(data_reg)
reg.labs <- colnames(YA0)
rownames(YA0) <- as.character(as.Date.ts(data_reg))

tdiff_hl <- (nrow(YQ0)/4)-nrow(YA0)
tdiff_hl_q <- 4*tdiff_hl
YA0 <- rbind(YA0,matrix(0,tdiff_hl,ncol(YA0)))

# some key dimensions
N <- ncol(YA0)
Q <- Q.lat + NCOL(YQ0) - 1 # (because of regional GVA being excluded)
Q.obs <- Q - Q.lat

# ----------------------------------------------------------
# create raw data objects and POOS design
if(ncol(YQ0) == 1){
  if(!is.na(Q.lat)){
    P0_sig2 <- c(rep(mean(apply(YA0[1:4,],2,sd)),Q.lat))^2 # scale prior on initial state 
  }else{
    P0_sig2 <- 1
  }
}else{
  P0_sig2 <- c(rep(mean(apply(YA0[1:4,],2,sd)),Q.lat),apply(YQ0[1:12,-sl_GVA,drop=FALSE],2,sd))^2 # scale prior on initial state  
}

if(run.cc == "FR"){
  P0_sig2 <- P0_sig2 / 4
}

Z <- cbind(kronecker(YA0,matrix(1,4,1)),YQ0)
omega_big <- kronecker(omega,matrix(1,4,1))

colnames(Z) <- c(colnames(YA0),"GDP")
T.full <- nrow(Z)
id_obs <- rep(0,T.full); id_obs[seq(4,T.full-tdiff_hl_q,by=4)] <- 1

Z.pc <- princomp(Z[,1:N])$scores
Z.og <- Z <- Z[(P+1):nrow(Z),]
id_obs <- id_obsNA <- id_obs[(P+1):length(id_obs)]

# ----------------------------------------------------------
# construct holdout stuff for out-of-sample exercise
rownames(Z.pc) <- rownames(YQ0) <- date.labs
rownames(Z) <- rownames(Z.og) <- names(id_obs) <- date.labs.p

# select training sample (only needed when computing pseudo out-of-samples)
sl.train <- 1:((T.full) - ((grid_h - run.h) * 4))
sl.train.p <- 1:((T.full - P) - ((grid_h - run.h) * 4))
Z <- Z[sl.train.p,]

# weights for the cross-sectional average
omega_big <- omega_big[sl.train.p,]
if(!reg_wghs){
  omega_big[] <- 1/N
}

Z.pc <- Z.pc[sl.train,]
YQ0 <- YQ0[sl.train,]

# identify when quarterly values are observed
id_obs <- id_obs[sl.train.p]
id_obsNA <- id_obsNA[sl.train.p]

# consider release calendar (exclude final two yearly observations)
id_obs[length(id_obs)] <- 0 # final year
id_obs[length(id_obs)-4] <- 0 # penultimate year

Z.og <- Z.og[sl.train.p,]

# selection for holdout sample (only needed when computing pseudo out-of-samples)
sl.true <- nrow(data_reg)-grid_h+run.h
sl.true <- c(sl.true-1,sl.true) # backcast and nowcast
true.gva <- data_reg[sl.true,]

dates_ncbc <- as.character(as.Date.ts(data_reg)[sl.true])
rownames(true.gva) <- paste0(c("bc","nc"),"_",dates_ncbc)

# ---------------
# estimate corresponding specification
if(run.model == "mfreg"){
  outputFile <- file(paste0(logs.dir, ".txt"))
  tryCatch({
    source("est_fcfunc.R")
  }, error = function(e) {
    writeLines(as.character(e), outputFile)
  })
  
  Z_ann_store <- Z_store
  
  fc_mu <- apply(Z_ann_store,c(2,3),median)[c(dim(Z_ann_store)[2]-4,dim(Z_ann_store)[2]),]
  fc_sd <- apply(Z_ann_store,c(2,3),sd)[c(dim(Z_ann_store)[2]-4,dim(Z_ann_store)[2]),]
  fc_mcmc <- Z_ann_store[,c(dim(Z_ann_store)[2]-4,dim(Z_ann_store)[2]),]
}else{
  fc_mcmc <- array(NA,dim=c(nsave,2,N))
  dimnames(fc_mcmc) <- list(paste0("draw",1:nsave),dates_ncbc,reg.labs)
  for(i in 1:N){
    Yraw <- YA0[1:(nrow(YA0)-grid_h+run.h-2),i]
    if(run.model == "ar") sim_bench <- bar1(Yraw,rw=FALSE,hm=FALSE,draws=nsave)
    if(run.model == "rw") sim_bench <- bar1(Yraw,rw=TRUE,hm=FALSE,draws=nsave)
    if(run.model == "hm") sim_bench <- bar1(Yraw,rw=FALSE,hm=TRUE,draws=nsave)
    
    fc_mcmc[,,i] <- sim_bench
  }
  fc_mu <- apply(fc_mcmc, c(2,3), median)
  fc_sd <- apply(fc_mcmc, c(2,3), sd)
}

# ---------------
# evaluation
lps <- dnorm(true.gva,fc_mu,fc_sd,log=TRUE)
error <- abs(true.gva-fc_mu)
crps <- true.gva*NA

for(j in 1:2){
  for(i in 1:N){
    crps[j,i] <- crps_sample(y=true.gva[j,i],dat=fc_mcmc[,j,i])
  }
}

# save additional moments
if(run.model == "mfreg"){
  sediag_post <- apply(sediag_store, c(2), quantile, probs = set.conf, na.rm=TRUE)
  uncm_post <- apply(uncm_store,c(2,3),quantile,probs = set.conf, na.rm=TRUE)  
  
  Z_ann_post <- apply(Z_ann_store, c(2,3), quantile, probs = set.conf, na.rm=TRUE)

  out <- list("Za"=Z_ann_post, # annualized GVA series on quarterly frequency
              "fa"=NULL, # annualized factors on quarterly frequency
              "reg"=reg.labs,
              
              # forecast metrics
              "fcmetr"=list("crps"=crps,
                            "lps"=lps,
                            "error"=error),
              
              # VAR/measurement diagnostics for the state equation
              "sediag" = sediag_post,
              "uncm" = uncm_post
  )
}else{
  out <- list("Za"=NULL, # annualized GVA series on quarterly frequency
              "fa"=NULL, # annualized factors on quarterly frequency
              "reg"=reg.labs,
              
              # forecast metrics
              "fcmetr"=list("crps"=crps,
                            "lps"=lps,
                            "error"=error),
              
              # VAR/measurement diagnostics for the state equation
              "sediag" = NULL,
              "uncm" = NULL
  )
}

save(out,file=file.dir)


