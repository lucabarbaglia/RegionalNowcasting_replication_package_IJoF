# ----------------------------------------------------------
# design matrices
if(Q.obs > 0){
  Yfraw <- cbind(Z.pc[,1:Q.lat,drop=FALSE] + matrix(rnorm(nrow(Z.pc) * Q.lat, mean = 0, sd = sqrt(1e-3)), nrow(Z.pc), Q.lat), 
                 YQ0[,-1])
}else{
  Yfraw <- cbind(Z.pc[,1:Q.lat,drop=FALSE] + matrix(rnorm(nrow(Z.pc) * Q.lat, mean = 0, sd = sqrt(1e-3)), nrow(Z.pc), Q.lat))
}

Yflags <- embed(Yfraw,P+1)
T <- nrow(Z)

# ----------------------------------------------------------
# state equation objects
Yf <- Yf_ann <- Yflags[,1:Q,drop=FALSE]

Xf <- cbind(Yflags[,-c(1:Q)],1)
Yf0 <- Yfraw[1:P,,drop=FALSE]
A_K <- ncol(Xf)

A_OLS <- solve(crossprod(Xf)+0.1*diag(ncol(Xf)))%*%crossprod(Xf,Yf)
Sigma <- 1e-2*diag(Q)
Sigma.big <- Matrix(0,Q*P,Q*P,sparse=TRUE)
Sigma.big[1:Q,1:Q] <- Sigma

Sigma.big_ls <- list()
for(tt in 1:T){
  Sigma.big_ls[[tt]] <- Sigma.big
}

# construct several needed variants of the loadings matrices
beta <- matrix(0,N+Q.obs,Q)
for(i in 1:(N+Q.obs)){
  if(i<=N){
    btmp <- try(solve(crossprod(Yf,Yf))%*%crossprod(Yf,Z[,i]), silent = TRUE)
    if(is(btmp,"try-error")){
      btmp <- try(solve(crossprod(Yf,Yf) + diag(ncol(Yf)))%*%crossprod(Yf,Z[,i]), silent = TRUE)
      if(is(btmp,"try-error")){
        btmp <- 0
      }
    }
    beta[i,] <- btmp
  }else{
    beta[i,Q.lat-N+i] <- 1
  }
}

if(ident_facs){
  beta[1:Q.lat,1:Q.lat] <- diag(Q.lat)
}

beta.big <- Matrix(rep(beta,P),N+Q.obs,Q*P,sparse=TRUE)
beta.big_ls <- list()
if(csr!="none"){
  for(tt in 1:T){
    beta.big.t <- rbind(beta.big,c(apply(as.matrix(omega_big[tt,]*beta[1:N,1:Q.lat]),2,sum),rep(0,Q*(P-1))))
    beta.big_ls[[tt]] <- beta.big.t
  }
}else{
  for(tt in 1:T){
    beta.big_ls[[tt]] <- beta.big
  }
}

theta_beta <- lambda_beta <- array(1,dim=c(N+Q.obs,Q))
nu_beta <- array(1,dim=c(N+Q.obs,Q))
tau_beta <- 1
zeta_beta <- 1

# measurement errors and cross-sectional restriction
sig2m <- c(rep(1e-3,N),rep(0,Q.obs))
if(csr!="none") sig2csr <- 1e-3 else sig2csr <- NULL

# MF restriction (intertemporal)
itr <- 1/itr_scl*c(1/4,1/2,3/4,1,3/4,1/2,1/4)
ITR <- matrix(0,nrow=1,Q*P)
for(i in 1:Q){
  if(i <= Q.lat) ITR[1,seq(i,Q*P,by=Q)] <- itr else ITR[1,seq(i,Q*P,by=Q)] <- c(rep(1,Q-Q.lat),rep(0,(Q-Q.lat)*(P-1)))
}

# setup for loadings matrices and identities
if(csr!="none"){
  H1 <- Matrix(rbind(kronecker(matrix(1,nrow=N),ITR),matrix(0,Q.obs+1,Q*P)),sparse=TRUE)
  H0 <- Matrix(0,nrow=1+ifelse(Q.obs>0,1,0),ncol=Q*P,sparse=TRUE)
  H1[nrow(H1),1:Q.lat] <- 1
  H0[nrow(H0),1:Q.lat] <- 1
}else{
  H1 <- Matrix(rbind(kronecker(matrix(1,nrow=N),ITR),matrix(0,Q.obs,Q*P)),sparse=TRUE)
  H0 <- Matrix(0,nrow=1,ncol=Q*P,sparse=TRUE)
}

# objects used for filtering/smoothing
s0 <- rep(0,Q*P)
P0 <- Matrix(0,Q*P,Q*P,sparse=TRUE)
P0[1:Q,1:Q] <- P0_sig2*diag(Q)

sl.lat <- NULL
for(i in 1:Q.lat){
  sl.lat <- sort(c(sl.lat,seq(i,Q*P,by=Q)))
}

st <- matrix(0,T,Q*P)
sst <- matrix(0,T,Q*P)
Pt <- list()
for(tt in 1:T){
  Pt[[tt]] <- Matrix(0,Q*P,Q*P)
}

# ----------------------------------------------------------
# VAR coefficients and covariance matrix (exogenous variables)
A_draw <- A_prior <- matrix(0,A_K,Q)
A_draw[1:Q,1:Q] <- 0.3*diag(Q)

A.big <- Matrix(rbind(t(A_draw),cbind(diag(Q*(P-1)),matrix(0,Q*(P-1),Q+1))),sparse=TRUE)
theta_A <- 1e-3*A_draw^0 # prior on VAR coefficients
id_A <- seq(1,A_K)

B0_draw <- Matrix(diag(Q),sparse=TRUE)
theta_B0 <- matrix(1,Q,Q)
id_B0 <- lower.tri(diag(Q))

# priors on VAR coefficients
lambda_A <- array(1,dim=c(A_K,Q))
nu_A <- array(1,dim=c(A_K,Q))
tau_A <- 1
zeta_A <- 1

# priors on covariance matrix
lambda_B0 <- array(1,dim=c(Q,Q))
nu_B0 <- array(1,dim=c(Q,Q))
tau_B0 <- 1
zeta_B0 <- 1

# priors on structural shocks in state equation
a_0 <- 3
b_0 <- 0.3

if(sv){
  library(stochvol)
  sv_sig2 <- 0.1
  sv_draw <- sv_latent <- sv_priors <- list()
  for (qq in seq_len(Q)){
    sv_draw[[qq]] <- list(mu = -5, phi = 0.95, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[qq]] <- rep(0,T)
    
    sv_priors[[qq]] <- specify_priors(
      mu = sv_normal(mean=-5,sd=0.1),
      phi = sv_beta(shape1 = 5, shape2 = 1.5),
      sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*sv_sig2)),
      nu = sv_infinity(),
      rho = sv_constant(0)
    )
  }
}

# priors on the measurement errors
if(me_prior == "tight"){
  a_m0 <- 1000
  b_m0 <- 0.01
}else{
  a_m0 <- 100
  b_m0 <- 0.01
}

# prior on cross-sectional restriction (tightness can be controlled here)
if(csr=="tight"){
  a_csr <- 500
  b_csr <- 0.05
  # plot(density(1/rgamma(3000,a_csr,b_csr)), main="Prior density for CSR (tight)")
}else if(csr=="loose"){
  a_csr <- 50
  b_csr <- 0.05
  # plot(density(1/rgamma(3000,a_csr,b_csr)), main="Prior density for CSR (loose)")
}else{
  a_csr <- 1
  b_csr <- 1
  # plot(density(1/rgamma(3000,a_csr,b_csr)), main="Prior density for CSR (none)")
}

H <- matrix(-3,T,Q)
eps <- matrix(0,T,Q)

# ----------------------------------------------------------
# sampler setup
ntot <- nburn + nthin*nsave
thinset <- seq(nburn+1,ntot,by=nthin)
if(is.null(toplot)){
  plotset <- 0
}else{
  plotset <- seq(toplot,ntot,by=toplot)
}

savecnt <- 0

# storage objects
A_store <- array(NA,dim=c(nsave,A_K,Q))
Sigma_store <- array(NA,dim=c(nsave,Q,Q))

beta_store <- array(0,dim=c(nsave,N+Q.obs,Q))
f_store <- array(NA,dim=c(nsave,T,Q.lat))
Z_store <- Zlat_store <- array(NA,dim=c(nsave,T,N))

uncm_store <- array(NA,dim=c(nsave,2,N))
sediag_store <- array(NA,dim=c(nsave,4))

# labeling (outputs with dates)
dimnames(f_store) <- list(paste0("draw",1:nsave),dimnames(Z.og)[[1]],paste0("fac",1:Q.lat))
dimnames(Z_store) <- dimnames(Zlat_store) <- list(paste0("draw",1:nsave),dimnames(Z.og)[[1]],reg.labs)
dimnames(beta_store) <- list(paste0("draw",1:nsave),reg.labs,paste0("fac",1:Q.lat))

dimnames(uncm_store) <- list(paste0("draw",1:nsave),c("mean","var"),reg.labs)
dimnames(sediag_store) <- list(paste0("draw",1:nsave),c("maxeigReA","maxeigImA","lneigenA","logdetSig"))

# start sampling loop
pb <- txtProgressBar(min = 0, max = ntot, style = 3)
time.st <- Sys.time()
irep <- 1

for(irep in 1:ntot){
  # Step 1a: State equation coefficients
  eigen.check <- TRUE
  check.count <- 0
  while(eigen.check){
    for(nn in 1:Q){
      normalizer <- exp(-H[,nn]/2)
      YY <- Yf[,nn]*normalizer
      if(nn!=1){
        XXo <- cbind(Xf,-Yf[,1:(nn-1)])
        XX <- XXo*normalizer
        VA_inv <- diag(A_K+(nn-1))/c(theta_A[,nn],theta_B0[nn,1:(nn-1)])
        a_prior <- c(A_prior[,nn],rep(0,nn-1))
        k_n <- ncol(XX)
      }else{
        XXo <- Xf
        XX <- XXo*normalizer
        VA_inv <- diag(A_K)/theta_A[,nn]
        a_prior <- A_prior[,nn]
        k_n <- A_K
      }
      
      VA_post <- solve(crossprod(XX) + VA_inv)
      A_post <- VA_post %*% (VA_inv %*% a_prior + crossprod(XX,YY))
      a_draw <- A_post + t(chol(VA_post))%*%rnorm(k_n)
      A_draw[,nn] <- a_draw[id_A]
      if(nn>1) B0_draw[nn,1:(nn-1)] <- a_draw[-id_A]
      eps[,nn] <- eps_nn <- Yf[,nn] - XXo%*%a_draw
    }
    
    # map back to reduced form
    B0i_draw <- solve(B0_draw)
    B_draw <- t(B0i_draw%*%t(A_draw))
    A.big[1:Q,] <- t(B_draw)
    
    # re-sampling if eigenvalue criterion is not fulfilled
    check.count <- check.count+1
    if(irep < nburn*0.1){
      eigen.check <- FALSE
      eigenvalues <- eigen(A.big[1:(A_K-1),1:(A_K-1)])$value
      ln_eigen <- sqrt(Re(eigenvalues)^2 + Im(eigenvalues)^2)
      eigen.check <- max(abs(sqrt(Re(eigenvalues)^2 + Im(eigenvalues)^2))) > eigen.crit
      
      A.eigen <- Re(eigenvalues[which.max(ln_eigen)])
      A.eigen.Im <- Im(eigenvalues[which.max(ln_eigen)])
    }else{
      eigenvalues <- eigen(A.big[1:(A_K-1),1:(A_K-1)])$value
      ln_eigen <- sqrt(Re(eigenvalues)^2 + Im(eigenvalues)^2)
      eigen.check <- max(ln_eigen) > eigen.crit
      
      A.eigen <- Re(eigenvalues[which.max(ln_eigen)])
      A.eigen.Im <- Im(eigenvalues[which.max(ln_eigen)])
    }
    if(check.count == max.try){
      eigen.check <- FALSE
      message("\nNo stable draw, A eigen Re=",round(A.eigen,digits=3),", Im=",round(A.eigen.Im,digits=3),", lngth=",round(max(ln_eigen),digits = 3))
    }
  }
  
  # "structural" variances for the factors
  if(sv){
    for(nn in 1:Q){
      svdraw_nn <- svsample_fast_cpp(eps[,nn],startpara=sv_draw[[nn]],startlatent=sv_latent[[nn]],priorspec=sv_priors[[nn]])
      sv_draw[[nn]][c("mu", "phi", "sigma", "nu")] <- as.list(svdraw_nn$para[, c("mu", "phi", "sigma", "nu")])
      H[,nn] <- sv_latent[[nn]] <- svdraw_nn$latent
    }
    
    H[H>0] <- 0
    H[H<(-8)] <- -8
    for(tt in 1:T){
      Sigma <- B0i_draw %*% (exp(H[tt,]) * diag(Q)) %*% t(B0i_draw)
      Sigma.big_ls[[tt]][1:Q,1:Q] <- Sigma
    }
  }else{
    for(nn in 1:Q){
      a_p <- a_0 + T/2
      b_p <- b_0 + crossprod(eps[,nn])/2
      H[,nn] <- -log(rgamma(1,a_p,b_p))
    }
    Sigma <- B0i_draw %*% (exp(H[1,]) * diag(Q)) %*% t(B0i_draw)
    for(tt in 1:T){
      Sigma.big_ls[[tt]][1:Q,1:Q] <- Sigma
    }
  }
  
  # Step 1b: Shrinkage prior
  hs_draw_A <- get.hs(bdraw=as.numeric(A_draw-A_prior),lambda.hs = lambda_A, nu.hs = nu_A, tau.hs = tau_A, zeta.hs = zeta_A)
  lambda_A <- matrix(hs_draw_A$lambda,A_K,Q)
  nu_A <- matrix(hs_draw_A$nu,A_K,Q)
  tau_A <- hs_draw_A$tau
  zeta_A <- hs_draw_A$zeta
  theta_A <- matrix(hs_draw_A$psi,A_K,Q)
  theta_A[theta_A<1e-8] <- 1e-8
  theta_A[theta_A>1] <- 1
  
  hs_draw_B0 <- get.hs(bdraw=as.matrix(B0_draw)[id_B0],lambda.hs = lambda_B0[id_B0], nu.hs = nu_B0[id_B0], tau.hs = tau_B0, zeta.hs = zeta_B0)
  lambda_B0[id_B0] <- hs_draw_B0$lambda
  nu_B0[id_B0] <- hs_draw_B0$nu
  tau_B0 <- hs_draw_B0$tau
  zeta_B0 <- hs_draw_B0$zeta
  theta_B0[id_B0] <- hs_draw_B0$psi
  theta_B0[theta_B0<1e-8] <- 1e-8
  theta_B0[theta_B0>1] <- 1
  
  # Step 2a: Loadings
  if(ident_facs){
    for(nn in 2:N){
      # from latent representation
      # yy <- Z[,nn]/sqrt(sig2m[nn])
      # xx <- Yf/sqrt(sig2m[nn])
      
      # from observables
      if(nn <= Q.lat){
        yy <- (Z.og[id_obs == 1,nn] - Yf_ann[id_obs == 1,nn])/sqrt(sig2m[nn])
        xx <- Yf_ann[id_obs == 1,1:(nn-1),drop=FALSE]/sqrt(sig2m[nn])
      }else{
        yy <- Z.og[id_obs == 1,nn]/sqrt(sig2m[nn])
        xx <- Yf_ann[id_obs == 1,,drop=FALSE]/sqrt(sig2m[nn])
      }
      Q_nn <- ncol(xx)
      
      bvinv <- solve(crossprod(xx) + theta_beta[nn,1:Q_nn]*diag(Q_nn))
      bm <- bvinv %*% (crossprod(xx,yy))
      beta[nn,1:Q_nn] <- bm + t(chol(bvinv)) %*% rnorm(Q_nn)
    }  
  }else{
    for(nn in 1:N){
      # from latent representation
      # yy <- Z[,nn]/sqrt(sig2m[nn])
      # xx <- Yf/sqrt(sig2m[nn])
      
      # from observables
      yy <- Z.og[id_obs == 1,nn]/sqrt(sig2m[nn])
      xx <- Yf_ann[id_obs == 1,,drop=FALSE]/sqrt(sig2m[nn])
      
      bvinv <- solve(crossprod(xx) + theta_beta[nn,]*diag(Q))
      bm <- bvinv %*% (crossprod(xx,yy))
      beta[nn,] <- bm + t(chol(bvinv)) %*% rnorm(Q)
    }  
  }
  
  # shrinkage on loadings
  hs_draw_beta <- get.hs(bdraw=as.numeric(beta[1:N,]),lambda.hs=lambda_beta[1:N,],nu.hs=nu_beta[1:N,],tau.hs=tau_beta,zeta.hs=zeta_beta)
  lambda_beta[1:N,] <- matrix(hs_draw_beta$lambda,N,Q)
  nu_beta[1:N,] <- matrix(hs_draw_beta$nu,N,Q)
  tau_beta <- hs_draw_beta$tau
  zeta_beta <- hs_draw_beta$zeta
  theta_beta[1:N,] <- matrix(hs_draw_beta$psi,N,Q)
  theta_beta[theta_beta > 1] <- 1
  theta_beta[theta_beta < 1e-8] <- 1e-8
  # theta_beta <- theta_beta^0 # in case no shrinkage
  
  # build cross-sectional and intertemporal loadings matrix
  beta.big <- Matrix(rep(beta,P),N+Q.obs,Q*P,sparse=TRUE)
  if(csr!="none"){
    for(tt in 1:T){
      if(id_obs[tt]==1){
        beta.big.t <- rbind(beta.big,c(apply(as.matrix(omega_big[tt,]*beta[1:N,1:Q.lat,drop=FALSE]),2,sum),rep(0,Q*(P-1))))
      }else{
        beta.big.t <- matrix(c(apply(as.matrix(omega_big[tt,]*beta[1:N,1:Q.lat]),2,sum),rep(0,Q*(P-1))), nrow = 1)
      }
      beta.big_ls[[tt]] <- beta.big.t
    }
  }else{
    for(tt in 1:T){
      beta.big_ls[[tt]] <- beta.big
    }
  }
  
  # ---
  # Step 2b: Measurement errors
  
  # from latent representation
  # eta <- Z[,1:N] - t(beta[1:N,] %*% t(Yf))
  # for(nn in 1:N){
  #   sig2m[nn] <- 1/rgamma(1,a_m0 + T/2, b_m0 + crossprod(eta[,nn])/2)
  # }
  
  # from observables
  eta <- Z.og[id_obs==1,1:N] - t(beta[1:N,] %*% t(Yf_ann[id_obs==1,,drop=FALSE]))
  for(nn in 1:N){
    sig2m[nn] <- 1/rgamma(1,a_m0 + sum(id_obs)/2, b_m0 + crossprod(eta[,nn])/2)
  }
  
  # Step 2d: cross-sectional restriction (in expectations)
  if(csr!="none"){
    eps_csr <- Z.og[,N+1] - apply(omega_big * t(beta %*% t(Yf)),1,sum)
    sig2csr <- 1/rgamma(1,a_csr+T/2,b_csr+crossprod(eps_csr)/2)
  }

  # Step 3a: Filtering
  for(tt in 1:T){
    Sigma.big <- Sigma.big_ls[[tt]]
    A <- A.big[1:(A_K-1),1:(A_K-1)]
    if(tt > 1){
      stt <- A.big %*% c(st[tt-1,],1)
      Ptt <- A %*% Pt[[tt-1]] %*% t(A) + Sigma.big
    }else{
      stt <- A.big %*% c(s0,1)
      Ptt <- A %*% P0 %*% t(A) + Sigma.big
    }
    
    # kalman updating steps (checks whether we observe something or not)
    if(csr!="none"){
      if(id_obs[tt]==1){
        H1beta <- beta.big_ls[[tt]]*H1
        Kg <- Ptt %*% t(H1beta) %*% solve(H1beta %*% Ptt %*% t(H1beta) + diag(N+Q.obs+1)*c(sig2m,sig2csr))
        st[tt,] <- as.numeric(stt + Kg %*% (c(Z.og[tt,1:N],Z.og[tt,N+1]/itr_scl) - H1beta %*% stt))
        Pt[[tt]] <- (diag(Q*P) - Kg %*% H1beta) %*% Ptt
      }else{
        H0beta <- beta.big_ls[[tt]]*H0
        if(Q.obs == 0){
          Kg <- Ptt %*% t(H0beta) %*% solve(H0beta %*% Ptt %*% t(H0beta) + diag(Q.obs+1)*c(sig2csr))
          st[tt,] <- as.numeric(stt + Kg %*% (c(Z.og[tt,N+1]/itr_scl) - H0beta %*% stt))
        }else{
          Kg <- Ptt %*% t(H0beta) %*% solve(H0beta %*% Ptt %*% t(H0beta) + diag(Q.obs+1)*c(sig2m[(N+1):ncol(Z.og)],sig2csr))
          st[tt,] <- as.numeric(stt + Kg %*% (c(Z.og[tt,(N+1):ncol(Z.og)],Z.og[tt,N+1]/itr_scl) - H0beta %*% stt)) # recode this in case of other observables
        }
        Pt[[tt]] <- (diag(Q*P) - Kg %*% H0beta) %*% Ptt
      }
    }else{
      if(id_obs[tt]==1){
        H1beta <- beta.big_ls[[tt]]*H1
        Kg <- Ptt %*% t(H1beta) %*% solve(H1beta %*% Ptt %*% t(H1beta) + diag(N+Q.obs)*sig2m)
        st[tt,] <- as.numeric(stt + Kg %*% (Z.og[tt,1:N] - H1beta %*% stt))
        Pt[[tt]] <- (diag(Q*P) - Kg %*% H1beta) %*% Ptt
      }else{
        st[tt,] <- as.numeric(stt) # no update of the states
        Pt[[tt]] <- Ptt
      }
    }
  }
  
  # Step 3b: Smoothing
  Pt.d <- forceSymmetric(Pt[[T]][sl.lat,sl.lat]) # stability
  sst[T,sl.lat] <- as.numeric(st[T,sl.lat]) # as.numeric(st[T,sl.lat] + t(chol(Pt.d)) %*% rnorm(Q.lat*P))
  sst[T,-sl.lat] <- st[T,-sl.lat]
  
  for(tt in (T-1):1){
    sst.d <- sst[tt+1,]
    st.d <- st[tt,]
    Pt.d <- Pt[[tt]]
    
    Sigma.big <- Sigma.big_ls[[tt]]
    f <- A %*% Pt.d %*% t(A) + Sigma.big
    inv_f <- try(Matrix::solve(f[sl.lat,sl.lat]),silent=TRUE)
    
    if(!is(inv_f,"try-error")){
      cfe <- sst.d - (A.big %*% c(st.d,1))
      sst.m <- st.d[sl.lat] + (Pt.d %*% t(A))[sl.lat,sl.lat] %*% inv_f %*% cfe[sl.lat]
      
      if(states == "ffbs"){
        sst.v <- Pt.d[sl.lat,sl.lat] - Pt.d[sl.lat,sl.lat] %*% t(A[sl.lat,sl.lat]) %*% inv_f %*% A[sl.lat,sl.lat] %*% Pt.d[sl.lat,sl.lat] 
      }
    }else{
      sst.m <- st.d[sl.lat]
      if(states == "ffbs"){
        sst.v <- Pt.d[sl.lat,sl.lat]
      }
    }
    
    sst[tt,sl.lat] <- as.numeric(sst.m)
    sst[tt,-sl.lat] <- st.d[-sl.lat]
    
    # vs. backward sampling (RTS sampler of Koop et al. 2020 above)
    if(states == "ffbs"){
      chol_sst.v <- try(chol(sst.v), silent = TRUE)
      if(is(chol_sst.v,"try-error")){
        chol_sst.v <- try(chol(sst.v + 1e-8 * diag(ncol(sst.v)), pivot=TRUE), silent = TRUE)
        if(is(chol_sst.v,"try-error")){
          chol_sst.v <- sst.v*0
        }
      }
      sst[tt,sl.lat] <- as.numeric(sst.m + t(chol_sst.v) %*% rnorm(ncol(sst.v))) 
    }
  }
  
  # Step 3c: construct state space data
  Y0 <- rbind(matrix(s0,P,Q),sst[,1:Q,drop=FALSE])
  Y0_mu <- matrix(apply(Y0,2,mean),nrow=nrow(Y0),ncol=ncol(Y0),byrow = TRUE)
  Y0_sd <- matrix(apply(Y0,2,sd),nrow=nrow(Y0),ncol=ncol(Y0),byrow = TRUE)
  Y0 <- (Y0 - Y0_mu)/Y0_sd + Y0_mu # standardize
  
  # fix sign
  for(qq in 1:Q.lat){
    beta_sign <- sign(beta[qq,qq])
    beta[,qq] <- beta[,qq] * beta_sign
    Y0[,qq] <- Y0[,qq] * beta_sign
  }
  
  # design matrices
  Yflags <- embed(as.matrix(Y0),P+1)
  Yf <- Yf_ann <- Yflags[,1:Q,drop=FALSE]
  Xf <- cbind(Yflags[,-c(1:Q)],1)
  
  for(qq in 1:Q.lat){
    for(t in P:(T+P-1)){
      Yf_ann[t-P+1,qq] <- (1/itr_scl)*(Y0[t,qq]/4 + 
                                     Y0[t-1,qq]/2 +
                                     3*Y0[t-2,qq]/4 + 
                                     Y0[t-3,qq] + 
                                     3*Y0[t-4,qq]/4 + 
                                     Y0[t-5,qq]/2 + 
                                     Y0[t-6,qq]/4)
    }
  }
  
  me_draw <- t(matrix(rep(sqrt(sig2m[1:N]),T)*rnorm(T*N),N,T))
  Z <- t(beta[1:N,,drop=FALSE] %*% t(Yf)) + me_draw
  Z_ann <- t(beta %*% t(Yf_ann)) + me_draw
  
  # Final: Storage
  if(irep %in% thinset){
    # storage
    savecnt <- savecnt + 1
    A_store[savecnt,,] <- as.matrix(B_draw)
    Sigma_store[savecnt,,] <- as.matrix(Sigma)
    
    beta_store[savecnt,,] <- beta
    f_store[savecnt,,] <- as.matrix(Yf_ann[,1:Q.lat])
    Z_store[savecnt,,] <- as.matrix(Z_ann)
    Zlat_store[savecnt,,] <- as.matrix(Z)
    
    # compute desired moments
    A_m <- A.big[,ncol(A.big),drop=FALSE]
    A_c <- A.big[,-ncol(A.big)]
    
    uncm <- beta[1:N,,drop=FALSE] %*% (solve(diag(A_K-1) - A_c) %*% A_m)[1:Q]
    # Sigma.hh <- matrix(solve(diag((A_K-1)^2) - kronecker(A_c,A_c)) %*% c(as.matrix(Sigma.big)), A_K-1, A_K-1)
    Sigma.hh <- Sigma.big
    for(hh in 1:100){
      Sigma.hh <- A_c %*% Sigma.hh %*% t(A_c) + Sigma.big
    }
    uncv <- diag(beta[1:N,,drop=FALSE] %*% Sigma.hh[1:Q,1:Q,drop=FALSE] %*% t(beta[1:N,,drop=FALSE]))
    sediag <- c(A.eigen,A.eigen.Im,max(ln_eigen),log(det(Sigma)))
    
    uncm_store[savecnt,,] <- t(cbind(uncm,uncv))
    sediag_store[savecnt,] <- sediag
  }
  
  if(irep %in% plotset){
    # if(Q.lat==1) par(mfrow=c(1,1),mar=c(4,4,1,1)) else par(mfrow=pdim(Q.lat),mar=c(4,4,1,1))
    # for(qq in 1:Q.lat) ts.plot(cbind(Yf[,qq]), ylab=paste0("irep = ",irep),col=c("red","black"))
    
    # par(mar=c(4,4,1,1), mfrow = c(2,2))
    # ts.plot(Yf)
    # ts.plot(cbind(Z.og[,N+1],apply(omega_big * t(beta %*% t(Yf)),1,sum)), lwd = c(2,1))
    # ts.plot(Z)
    # ts.plot(Z_ann)

    # par(mar=c(4,4,1,1), mfrow = c(3,3))
    # Z.og_NA <- Z.og
    # Z.og_NA[id_obs == 0,] <- NA
    # for(ii in 1:N){
    #   plot(Z_ann[,ii], type = "l")
    #   points(Z.og_NA[,ii])
    # }
    # 
    # message("\nTime per draw (last ",toplot," draw(s)): ",
    #         round(as.numeric(difftime(Sys.time(),time.st,units="secs")),digits=2)/toplot," secs. 
    #         A eigen Re=", round(A.eigen,digits=3),", Im=",round(A.eigen.Im,digits=3), ", lngth=",round(max(ln_eigen),digits = 3)," (",check.count," tries)")
    # time.st <- Sys.time()
  }
  
  setTxtProgressBar(pb, irep)
}

if(!is.null(toplot)){dev.off()}



