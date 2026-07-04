#
#
#     Internal fitting functions for the pnam R package to fit panel network 
#     autocorrelation model. Internally, the models follow the work of 
#     Millo and Piras (2012). splm: Spatial Panel data models in R. Journal of Statistical Software
#     where first, the rho, lambda, and theta values are found by numerical optimization
#     with optim(), then the values for beta and sigma2 are estimated via GLS. 
#
#     last updated on 07-02-2026 (pnam version: 0.0.1)
#
#

#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters (supplied by pnam())
#' @param optim.method The user-provided optim search method (supplied by pnam())
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.random <- function(X,  #the model formula
                            Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                            W,        #this should be a list of k NT x NT block-diagonal network matrices
                            n,       #a vector that indexes the unit (respondents)
                            t,  #a vector that indexes the time (respondents)
                            optim.control, #the list of optim control options
                            optim.method){    
  
  nt <- n*t #the total number of panel data points NT
  k <- length(W) + 1 #the number of searching parameters (length of the rho vector,theta) 
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  #an r side wrapper here to source to c++
  pnam_ll_rule_searchR <- function(par,Y,X,W,n,t){
    rho <- par[1:(length(par)-1)] #the vector of possible rho values (in case of an mstar model)
    theta <- par[length(par)] #the last element will be the theta parameter value
    return(pnam_ll_random_find(Y,X,W,n,t,rho,theta))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values (for the starting theta value!)
  mle.values <- optim(par = starting,
                      fn = pnam_ll_rule_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      n = n,
                      t = t,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = c(rep(length(W),-0.99999),0),
                      upper = c(rep(length(W),0.99999),Inf),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1:(length(mle.values$par)-1)]
  theta <- mle.values$par[length(mle.values$par)]
  net.film <- netfilter(W,rho = rho)
  omega.est <- omega(theta=theta,n,t)
  inv.omega.est <- solve(omega.est)
  beta.est <- gls_compute_beta(net.film,omegainv=inv.omega.est,X=X,Y=Y,random=TRUE)
  vt.est <- gls_compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2(vt.est,omegainv=inv.omega.est,random=TRUE)
  ll.est <- pnam_ll_random(netA = net.film,sigma2 = sigma2.est,
                         v=vt.est,omega=omega.est,n=n,t=t)
  noneta.vcov <- solve(pnamhess)
  rho.vcov <- as.matrix(noneta.vcov[1:length(W),1:length(W)])
  colnames(rho.vcov) <- rownames(rho.vcov) <- names(W)
  names(rho) <- names(W)
  theta.vcov <- noneta.vcov[length(W)+1,length(W)+1]
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%inv.omega.est%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  k <- length(rho) + length(beta.est) + 2 #plus 2 for theta and sigma2
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  results <- list(rho = rho,
                  vcov.rho = rho.vcov,
                  theta = theta,
                  vcov.theta = theta.vcov,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = vt.est,
                  fitted.values = fitted,
                  error.dis = list(theta = theta,
                                   idiosyncratic = sigma2.est,
                                   unit = (theta*sigma2.est)),
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}


#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters
#' @param optim.method The user-provided optim search method
#' @param dyanmic.lag The type of dynamic lag structure (supplied by pnam())
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.dynam <- function(X,  #the model formula
                           Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                           W,   #this should be a list of T network N x N matrices
                           n,       #a vector that indexes the unit (respondents)
                           t, #a vector that indexes the time (respondents)
                           dyanmic.lag = c("network", "outcome", "network x outcome"),
                           optim.method,
                           optim.control){     
  
  nt <- n*t #the total number of panel data points NT
  ###--------------------------------------------------------------------------###
  #    creating the L matrix for temporal filtering to create the lags!
  ###--------------------------------------------------------------------------###
  #creating the lag filtering block-diagonal matrix
  L <- matrix(0,nt,nt) #an empty matrix of 0s
  k<-length(W)
  start <- n+1 #the start of the lower diagonal
  for(i in start:nt){ #for all elements across sub lower diagonal
    L[i, i - n] <- 1 #add a value of 1
  }
  if(dyanmic.lag == "outcome"){ X <- as.matrix(data.frame(X,lag.y=L%*%Y))}
  if(dyanmic.lag == "network"){
    for(i in 1:length(W)){
      X <- as.matrix(data.frame(X,L%*%W[[i]]%*%Y))
      colnames(X)[ncol(X)] <- paste0("lag.net.",names(W)[i])
    }
  }
  if(dyanmic.lag == "network x outcome"){
    X <- as.matrix(data.frame(X,lag.y = L%*%Y))
    for(i in 1:length(W)){
      X <- as.matrix(data.frame(X,L%*%W[[i]]%*%Y))
      colnames(X)[ncol(X)] <- paste0("lag.net.",names(W)[i])
    }
  }
  #dropping the first time point for conditional log likelihood estimation
  X <- X[-(1:n),]
  Y <- Y[-(1:n)]
  #also, correcting the W matrix!
  W <- lapply(W,function(net){net[-(1:n),-(1:n)]})
  t <- t - 1 #correcting for the dropping of the first wave
  #luckily, given the structure, the results can be found using the fixed model formulation
  #an r side wrapper here to source to c++
  pnam_searchR <- function(par,Y,X,W,n,t){
    rho <- par #the vector of possible rho values (in case of an mstar model)
    return(pnam_ll_fixed_find(Y,X,W,n,t,rho))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      n = n,
                      t = t,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = rep(length(W),-0.99999),
                      upper = rep(length(W),0.99999),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par
  net.film <- netfilter(W,rho = rho)
  beta.est <- gls_compute_beta(net.film,omegainv=diag(10),X=X,Y=Y,random=FALSE)
  vt.est <- gls_compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2(vt.est,omegainv=diag(10),random=FALSE)
  ll.est <- pnam_ll_fixed(netA = net.film,sigma2 = sigma2.est,v=vt.est,n=n,t=t)
  rho.vcov<- solve(pnamhess)
  colnames(rho.vcov) <- rownames(rho.vcov) <-  names(W)
  names(rho) <-  names(W)
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  k <- length(rho) + length(beta.est) + 1 #plus 2 for sigma2
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  
  results <- list(rho = rho,
                  vcov.rho = rho.vcov,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = vt.est,
                  fitted.values = fitted,
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}



#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param id.index the vector of indices to map the unique units. (supplied by pnam())
#' @param time.index the vector of indices to map the time units. (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters
#' @param optim.method The user-provided optim search method
#' @param fe the type of demeaning for the fixed effects (supplied by pnam())
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.fixed <- function(X,  #the model formula
                           Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                           W,   #this should be a list of T network N x N matrices
                           n,       #a vector that indexes the unit (respondents)
                           t, #a vector that indexes the time (respondents)
                           id.index,
                           time.index,
                           fe, #fixed = c("individual", "time", "two-way")
                           optim.method,
                           optim.control){    
 
  nt <- n*t #the total number of panel data points NT
  k <- length(W)  #the number of searching parameters (length of the rho vector,theta) 
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  if(fe == "actor")  Q.matrix <- diag(nt) - kronecker(matrix(1,t,t)/t, diag(n))
  if(fe == "time")  Q.matrix <- diag(nt) -kronecker(diag(t),matrix(1, n, n) / n)
  if(fe == "two-way")  Q.matrix <- (diag(nt) - kronecker(matrix(1,t,t)/t, diag(n))) %*% (diag(nt) -kronecker(diag(t),matrix(1, n, n)/ n))
  X <- Q.matrix%*%X #demeaning the matrix of covariates
  Y <- Q.matrix%*%Y #demeaning the matrix of outcome variables
  #an r side wrapper here to source to c++
  pnam_searchR <- function(par,Y,X,W,n,t){
    rho <- par #the vector of possible rho values (in case of an mstar model)
    return(pnam_ll_fixed_find(Y,X,W,n,t,rho))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      n = n,
                      t = t,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = rep(length(W),-0.99999),
                      upper = rep(length(W),0.99999),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par
  net.film <- netfilter(W,rho = rho)
  beta.est <- gls_compute_beta(net.film,omegainv=diag(10),X=X,Y=Y,random=FALSE)
  vt.est <- gls_compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2(vt.est,omegainv=diag(10),random=FALSE)
  ll.est <- pnam_ll_fixed(netA = net.film,sigma2=sigma2.est,v=vt.est,n=n,t=t)
  rho.vcov <- solve(pnamhess)
  colnames(rho.vcov) <- rownames(rho.vcov) <- names(W)
  names(rho) <- names(W)
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  k <- length(rho) + length(beta.est) + 1 #plus 2 for sigma2
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  
  results <- list(rho = rho,
                  vcov.rho = rho.vcov,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = vt.est,
                  fitted.values = fitted,
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}


#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param id.index the vector of indices to map the unique units. (supplied by pnam())
#' @param time.index the vector of indices to map the time units. (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters
#' @param optim.method The user-provided optim search method
#' @param fe the type of demeaning for the fixed effects (supplied by pnam())
#' @param W2 the NT x NT block-diagonal matrix for the error network autocorrelation. 
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.fixed.error <- function(X,  #the model formula
                           Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                           W,   #this should be a list of T network N x N matrices
                           n,       #a vector that indexes the unit (respondents)
                           t, #a vector that indexes the time (respondents)
                           id.index,
                           time.index,
                           fe, #fixed = c("individual", "time", "two-way")
                           optim.method,
                           optim.control,
                           W2){    
  
  nt <- n*t #the total number of panel data points NT
  k <- length(W) +1  #the number of searching parameters (length of the rho vector,theta) 
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  if(fe == "actor")  Q.matrix <- diag(nt) - kronecker(matrix(1,t,t)/t, diag(n))
  if(fe == "time")  Q.matrix <- diag(nt) -kronecker(diag(t),matrix(1, n, n) / n)
  if(fe == "two-way")  Q.matrix <- (diag(nt) - kronecker(matrix(1,t,t)/t, diag(n))) %*% (diag(nt) -kronecker(diag(t),matrix(1, n, n)/ n))
  X <- Q.matrix%*%X #demeaning the matrix of covariates
  Y <- Q.matrix%*%Y #demeaning the matrix of outcome variables
  #an r side wrapper here to source to c++
  pnam_searchR <- function(par,Y,X,W,n,t,W2){
    rho <- par[1:(length(par)-1)] #the vector of possible rho values (in case of an mstar model)
    lambda <- par[length(par)]
    return(pnam_ll_fixed_error_find(Y,X,W,n,t,rho,W2,lambda))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      W2=W2,
                      n = n,
                      t = t,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = rep(length(W)+1,-0.99999),
                      upper = rep(length(W)+1,0.99999),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1:(length(mle.values$par)-1)]
  lambda <- mle.values$par[length(mle.values$par)]
  net.film <- netfilter(W,rho = rho)
  B <- create_B(W2,lambda)
  LAMBDA <- create_LAMBDA(B)
  beta.est <- gls_compute_beta_error(net.film,LAMBDA,X=X,Y=Y)
  e.est <- gls_compute_e(net.film,beta=beta.est,B,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2_error(e.est,LAMBDA)
  ll.est <- pnam_ll_fixed_error(netA = net.film,B=B,sigma2=sigma2.est,e=e.est,n=n,t=t)
  noneta.vcov<- solve(pnamhess)
  vcov.rho <- as.matrix(noneta.vcov[1:length(rho),1:length(rho)])
  colnames(vcov.rho) <- rownames(vcov.rho) <-  names(W)
  names(rho) <- names(W)
  vcov.lambda <- noneta.vcov[ncol(noneta.vcov),ncol(noneta.vcov)]
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%LAMBDA%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  k <- length(rho) + length(beta.est) + 2 #plus 2 for sigma2 and the lambda
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  
  results <- list(rho = rho,
                  vcov.rho = vcov.rho,
                  lambda=lambda,
                  vcov.lambda=vcov.lambda,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = e.est,
                  fitted.values = fitted,
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  error.network = W2,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}





#' @name  panelnam.dynam.error
#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters
#' @param optim.method The user-provided optim search method
#' @param dyanmic.lag the type of lags to add to the model.  (supplied by pnam())
#' @param W2 the NT x NT block-diagonal matrix for the error network autocorrelation. (supplied by pnam())
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.dynam.error <- function(X,  #the model formula
                           Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                           W,   #this should be a list of T network N x N matrices
                           n,       #a vector that indexes the unit (respondents)
                           t, #a vector that indexes the time (respondents)
                           dyanmic.lag = c("network", "outcome", "network x outcome"),
                           optim.method,
                           optim.control,
                           W2){     
  
  nt <- n*t #the total number of panel data points NT
  ###--------------------------------------------------------------------------###
  #    creating the L matrix for temporal filtering to create the lags!
  ###--------------------------------------------------------------------------###
  #creating the lag filtering block-diagonal matrix
  L <- matrix(0,nt,nt) #an empty matrix of 0s
  k <- length(W) +1  #the number of searching parameters (length of the rho vector,theta) 
  start <- n+1 #the start of the lower diagonal
  for(i in start:nt){ #for all elements across sub lower diagonal
    L[i, i - n] <- 1 #add a value of 1
  }
  if(dyanmic.lag == "outcome"){ X <- as.matrix(data.frame(X,lag.y=L%*%Y))}
  if(dyanmic.lag == "network"){
    for(i in 1:length(W)){
      X <- as.matrix(data.frame(X,L%*%W[[i]]%*%Y))
      colnames(X)[ncol(X)] <- paste0("lag.net.",names(W)[i])
    }
  }
  if(dyanmic.lag == "network x outcome"){
    X <- as.matrix(data.frame(X,lag.y = L%*%Y))
    for(i in 1:length(W)){
      X <- as.matrix(data.frame(X,L%*%W[[i]]%*%Y))
      colnames(X)[ncol(X)] <- paste0("lag.net.",names(W)[i])
    }
  }
  #dropping the first time point for conditional log likelihood estimation
  X <- X[-(1:n),]
  Y <- Y[-(1:n)]
  #also, correcting the W matrix!
  W <- lapply(W,function(net){net[-(1:n),-(1:n)]})
  t <- t - 1 #correcting for the dropping of the first wave
  #luckily, given the structure, the results can be found using the fixed model formulation
  #an r side wrapper here to source to c++
  pnam_searchR <- function(par,Y,X,W,n,t,W2){
    rho <- par[1:(length(par)-1)] #the vector of possible rho values (in case of an mstar model)
    lambda <- par[length(par)]
    return(pnam_ll_fixed_error_find(Y,X,W,n,t,rho,W2,lambda))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      W2=W2,
                      n = n,
                      t = t,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = rep(length(W)+1,-0.99999),
                      upper = rep(length(W)+1,0.99999),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1:(length(mle.values$par)-1)]
  lambda <- mle.values$par[length(mle.values$par)]
  net.film <- netfilter(W,rho = rho)
  B <- create_B(W2,lambda)
  LAMBDA <- create_LAMBDA(B)
  beta.est <- gls_compute_beta_error(net.film,LAMBDA,X=X,Y=Y)
  e.est <- gls_compute_e(net.film,beta=beta.est,B,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2_error(e.est,LAMBDA)
  ll.est <- pnam_ll_fixed_error(netA = net.film,B=B,sigma2=sigma2.est,e=e.est,n=n,t=t)
  noneta.vcov<- solve(pnamhess)
  vcov.rho <- as.matrix(noneta.vcov[1:length(rho),1:length(rho)])
  colnames(vcov.rho) <- rownames(vcov.rho) <-  names(W)
  names(rho) <- names(W)
  vcov.lambda <- noneta.vcov[ncol(noneta.vcov),ncol(noneta.vcov)]
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%LAMBDA%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  k <- length(rho) + length(beta.est) + 2 #plus 2 for sigma2 and the lambda
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  
  results <- list(rho = rho,
                  vcov.rho = vcov.rho,
                  lambda=lambda,
                  vcov.lambda=vcov.lambda,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = e.est,
                  fitted.values = fitted,
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  error.network = W2,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}







#' @param X the NT x K covariate matrix. (supplied by pnam())
#' @param Y The NT x 1 vector of outcomes. (supplied by pnam())
#' @param W A list of J block-diagonal NT x NT matrix (supplied by pnam())
#' @param n the number of units (supplied by pnam())
#' @param t the number of time points (supplied by pnam())
#' @param optim.control The list of user-provided optim control parameters (supplied by pnam())
#' @param optim.method The user-provided optim search method (supplied by pnam())
#' @param W2 the NT x NT block-diagonal matrix for the error network autocorrelation. 
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
panelnam.random.error <- function(X,  #the model formula
                            Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                            W,        #this should be a list of k NT x NT block-diagonal network matrices
                            n,       #a vector that indexes the unit (respondents)
                            t,  #a vector that indexes the time (respondents)
                            optim.control, #the list of optim control options
                            optim.method,
                            W2 
){    
  
  nt <- n*t #the total number of panel data points NT
  k <- length(W) + 2 #the number of searching parameters (length of the rho vector,theta) 
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  
  #an r side wrapper here to source to c++
  pnam_ll_rule_searchR <- function(par,Y,X,W,n,t,W2){
    rho <- par[1:(length(par)-2)] #the vector of possible rho values (in case of an mstar model)
    lambda <- par[(length(par)-1)] #the vector of possible rho values (in case of an mstar model)
    theta <- par[length(par)] #the last element will be the theta parameter value
    return(pnam_ll_random_error_find(Y,X,W,n,t,rho,theta,W2,lambda))
  }
  
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  starting <- abs(rnorm(k, mean=0,sd=0.001)) #ensuring positive values (for the starting theta value!)
  mle.values <- optim(par = starting,
                      fn = pnam_ll_rule_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      n = n,
                      t = t,
                      W2=W2,
                      method = optim.method,
                      #row-normalized matrices so the minimum for rho is -0.9999 and the minimum for theta is 0
                      lower = c(rep(length(W)+1,-0.99999),0),
                      upper = c(rep(length(W)+1,0.99999),Inf),
                      control = optim.control,
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1:(length(mle.values$par)-2)]
  lambda <- mle.values$par[(length(mle.values$par)-1)]
  theta <- mle.values$par[length(mle.values$par)]
  net.film <- netfilter(W,rho = rho)
  omega.est <- omega_error(theta=theta,n,t,W2,lambda)
  inv.omega.est <- solve(omega.est)
  beta.est <- gls_compute_beta(net.film,omegainv=inv.omega.est,X=X,Y=Y,random=TRUE)
  vt.est <- gls_compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- gls_compute_sigma2(vt.est,omegainv=inv.omega.est,random=TRUE)
  ll.est <- pnam_ll_random_error(netA = net.film,sigma2 = sigma2.est,
                           v=vt.est,omega=omega.est,n=n,t=t)
  noneta.vcov<- solve(pnamhess)
  vcov.rho <- as.matrix(noneta.vcov[1:length(rho),1:length(rho)])
  colnames(vcov.rho) <- rownames(vcov.rho) <-  names(W)
  names(rho) <- names(W)
  vcov.lambda <- noneta.vcov[length(rho)+1,length(rho)+1]
  vcov.theta <- noneta.vcov[length(rho)+2,length(rho)+2]
  beta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%inv.omega.est%*%X)
  colnames(beta.vcov) <-rownames(beta.vcov)  <- colnames(X)
  fitted <- solve(net.film,X%*%beta.est) #the fitted values for Y: Y = A^(-1)(XB)
  beta.est <- as.vector(beta.est)
  names(beta.est) <- colnames(X)
  k <- length(rho) + length(beta.est) + 3 #plus 2 for sigma2,lambda, and theta
  df<-nt-k
  AIC <- 2*k - 2*ll.est
  BIC <- -2*ll.est + k*log(nt)
  
  results <- list(rho = rho,
                  vcov.rho = vcov.rho,
                  lambda=lambda,
                  theta=theta,
                  vcov.theta=vcov.theta,
                  vcov.lambda=vcov.lambda,
                  coefficients = beta.est,
                  vcov = beta.vcov,
                  sigma2 = sigma2.est,
                  residuals = vt.est,
                  fitted.values = fitted,
                  logLik = ll.est,
                  convergence = mle.values$convergence,
                  error.dis = list(theta = theta,
                                   idiosyncratic = sigma2.est,
                                   unit = (theta*sigma2.est)),
                  netfilter = net.film,
                  k = k, 
                  df=df,
                  AIC=AIC,BIC=BIC,
                  error.network = W2,
                  response = Y,
                  covariates = X,
                  optim.information=mle.values)
  return(results)
}




