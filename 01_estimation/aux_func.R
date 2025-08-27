# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# bayesian univariate autoregression
bar1 <- function(Yraw,rw=TRUE,hm=FALSE,draws=3000){
  Ylag <- embed(Yraw,dimension=2)
  yy <- Ylag[,1]
  if(rw) xx <- Ylag[,2]
  if(!rw) xx <- cbind(Ylag[,2],1)
  if(hm) xx <- matrix(1, nrow = nrow(Ylag), ncol = 1)
  K <- NCOL(xx)
  
  a0 <- 0.01
  b0 <- 0.01
  phi0 <- rep(10,K)
  
  if(rw) phi <- 0 else phi <- c(0,0)
  if(hm) phi <- 0
  sig2 <- 0.1
  
  yfc_store <- sig2_store <- phi_store <- matrix(NA,draws,2)
  for(irep in 1:(draws+500)){
    if(!rw){
      v_post <- solve(crossprod(xx)/sig2 + diag(K) / phi0)
      phi_post <- v_post %*% (crossprod(xx,yy)/sig2)
      phi <- as.numeric(phi_post + t(chol(v_post)) %*% rnorm(K)) 
    }
    
    if(rw) eps <- yy - phi*xx
    if(!rw) eps <- yy - xx %*% phi
    sig2 <- 1/rgamma(1,a0+length(yy)/2,b0+crossprod(eps)/2)
    
    if(irep > 500){
      phi_store[irep-500,] <- phi
      sig2_store[irep-500,] <- sig2
      
      if(hm){
        yfc_store[irep-500,] <- phi + rnorm(2,0,sqrt(sig2))
      }else{
        # first quarter
        if(rw) yfc_store[irep-500,1] <- yfc <- phi*yy[length(yy)] + rnorm(1,0,sqrt(sig2))
        if(!rw) yfc_store[irep-500,1] <- yfc <- sum(phi*c(yy[length(yy)],1)) + rnorm(1,0,sqrt(sig2))
        
        # second quarter
        if(rw) yfc_store[irep-500,2] <- phi*yfc + rnorm(1,0,sqrt(sig2))
        if(!rw) yfc_store[irep-500,2] <- sum(phi*c(yfc,1)) + rnorm(1,0,sqrt(sig2)) 
      }
    }
  }
  return(yfc_store)
}

pdim <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  set <- cbind(factors,x/factors)
  return(c(set[floor(nrow(set)/2),]))
}

# --------------------------------------------------------------------------------------
# factor imputation
factorize <- function(data,n_fac,id=1){
  
  # Author: Stephane Surprenant
  # Creation: 12/10/2017
  #
  # This function estimates factors and loadings using principal component
  # analysis. The model is X = F^0 Lambda^0' + e, with X is T by N data matrix,
  # F is the the T by r factor matrix, Lambda is the N by r loadings matrix and
  # e is the T by r errors matrix. (Supercript k indicates k factors used.)
  #
  # Under F^k'F^k/T = I, maximizing the trace of F^k'(XX')F^k, \hat{F}^k = 
  # sqrt(T) times the eigenvectors of the k largest eigenvalues of XX'. Under 
  # Lambda'Lambda/N = I, we have Lambda_k_hat = sqrt(N) times the eignevectors
  # of X'X. It also implies \hat{F}^k = X \hat{Lambda}^k/N. (Bai & Ng, 2002).
  # The second choice is implemented here.
  #
  # Bai, J. and Serena Ng. 2002. Determining the number of factors in
  #     approximate factor models. Econometrica, 70(1), p. 191-221.
  #
  # ========================================================================= #
  #                               ESTIMATION
  # ========================================================================= #
  
  X <- data  # Reassign
  r <- n_fac # Reassign
  
  bign <- dim(X)[2]
  bigt <- dim(X)[1]
  
  svd <- svd(crossprod(X))         # svd$u = svd$v (symmetry); svd$d = eigenvalues
  lambda <- svd$u[,1:r]*sqrt(bign) # r th column times r th biggest eigenvalue
  f_hat <- X%*%lambda/bign         # factors
  e_hat <- X - f_hat%*%t(lambda)   # errors
  mse <- sum(e_hat^2)/(bign*bigt)
  
  if(id!=3){lambda = lambda %*% diag(r)} #does nothing
  
  if(id==3){
    F.F = pracma::sqrtm(crossprod(f_hat))
    
    lambda =  lambda %*% F.F$B
    f_hat = f_hat%*%F.F$Binv
    
    #rotation to PC3 identification
    lambda=lambda %*% solve(lambda[1:r,])
    f_hat = tcrossprod(f_hat,lambda[1:r,])
  }
  
  results <- list(factors = f_hat, lambda = lambda, mse = mse)
  
  return(results)
}

EM_sw <- function(data, n, it_max=50,std.opt=TRUE){
  #
  # Author: Stephane Surprenant
  # Creation: 16/11/2017
  #
  # Description:
  # This function inputs missing values in 'data' by using an
  # Expectation-Maximization algorithm. The data must is standerdized and
  # centered at zero. The algorithm begins by replacing the missing values
  # by 0 (the unconditional mean of the available observations). Then, it
  # estimates a factor model through principal component. Missing values
  # are then replaced by the predicted common component of the series.
  #
  # INPUTS: data: a matrix of data with time series as column vectors;
  #         n: number of factors considered
  # OUTPUT: X: a balanced panel
  #
  # Note: Requires package pracma (for repmat) and factor.R for PCA.
  #
  # ========================================================================= #
  #                               ESTIMATION
  # ========================================================================= #
  
  X <- as.matrix(data)
  
  # 1. Standerdize data
  size <- dim(X)
  if(std.opt==TRUE){
    a <- standard(X)
    X0 <- a$Y
    std <- a$std
    mu <- a$mean
  }else{
    X0 <- X
    std <- rep(1,ncol(X))
    mu <- rep(0,ncol(X)) 
  }
  # 2. First replacement and PCA
  X0[is.na(X)] <- 0 # Unconditional mean of standard data is zero
  X0[is.na(X0)] <- 0 # Unconditional mean of standard data is zero
  
  a <- factorize(X0, n_fac=n) # PCA
  f_1 <- a$factors
  L_1 <- a$lambda
  x_hat <- X0
  cc <- f_1%*%t(L_1) # Common component 
  x_hat[is.na(X)] <- cc[is.na(X)]
  x_hat <- x_hat*repmat(std,size[1],1) + repmat(mu,size[1],1)
  
  # 3. Initialize algorithm
  it <- 0      # Initialize iterations
  err <- 999   # Initialize error criterion to arbitrarily high value
  # it_max <- 50 # Selection maximum iterations
  con <- 1e-6  # Set convergence criterion
  old_f <- array(0, dim=c(size[1],1))
  
  while(it < it_max && err > con){
    # 1. Standerdize dataset
    a <- standard(x_hat) 
    std_new <- a$std
    mu_new <- a$mean
    x_new <- a$Y
    
    old_f <- f_1 # Save factors from previous iteration
    
    # 2. Extract factor through PCA
    a <- factorize(x_new, n_fac=n)
    f_1 <- a$factors
    L_1 <- a$lambda
    mse <- a$mse
    cc <- f_1%*%t(L_1) # Common component 
    x_hat <- x_new
    x_hat[is.na(X)] <- cc[is.na(X)]
    x_hat <- x_hat*repmat(std_new,size[1],1) + repmat(mu_new,size[1],1)
    
    mean_old_f2 <- apply(old_f^2, c(2), mean)
    mean_new_f2 <- apply(f_1^2, c(2), mean)
    
    err <- abs(mean(mean_old_f2 - mean_new_f2))
    
    it <- it + 1;
  }
  x_hat[!is.na(X)] <- X[!is.na(X)]
  results <- list(data=x_hat, factors=f_1, lambda=L_1, iterations=it, mse=mse)
  return(results)
}

# Subfunction to standerdize dataset ======================================== #
standard <- function(Y){
  require(pracma)
  Y <- as.matrix(Y)
  size <- dim(Y)
  
  mean_y <-apply(Y, c(2), mean, na.rm=TRUE)
  sd_y <- apply(Y, c(2), sd, na.rm=TRUE)
  # print(Y)
  # print(mean_y)
  # print(size[1])
  # print(repmat(mean_y, size[1],1))
  
  Y0 <- (Y - pracma::repmat(mean_y, size[1],1))/pracma::repmat(sd_y, size[1],1)
  return(list(Y=Y0, mean=mean_y, std=sd_y))
}