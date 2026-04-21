#
#
#     a wrapper function for the panel network autocorrelation models that
#     completes the scoring iterations following Millo, where the GLS 
#     steps occur within the search for the network lag and error parameters
#
#
#     last updated on 04-16-2026 (0.0.1)
#
#
#

#' @title Estimate a Panel Network Autocorrelation Model with Unit-Level Random Effects
#' @name  panelnamrule
#' @param X the covariate NT x K matrix.
#' @param Y The vector of outcomes.
#' @param W The block-diagonal NT x NT matrix
#' @param n the number of units
#' @param t the number of time points
#' @import Rcpp
#' @import RcppArmadillo
#' @return The maximum likelihood estimates.
#' @export
panelnamrule <- function(X,  #the model formula
                         Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                         W,        #this should be a list of T network N x N matrices
                         n,       #a vector that indexes the unit (respondents)
                         t){     #a vector that indexes the time (respondents)
        
  nt <- n*t #the total number of panel data points NT
  
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
   
  #an r side wrapper here to source to c++
  pnam_ll_rule_searchR <- function(par,Y,X,W,n,t){
    rho <- par[1]
    theta <- par[2]
    return(pnam_ll_rule_search(Y,X,W,n,t,rho,theta))
  }
  pnam_ll_rule_hessR <- function(par,Y,X,W,n,t){
    rho <- par[1]
    theta <- par[2]
    return(-pnam_ll_rule_search(Y,X,W,n,t,rho,theta))
  }
  
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  k <- 2 #the number of searching parameters (rho,theta)
  starting <- abs(rnorm(k, mean=0,sd=0.1)) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_ll_dynam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      n = n,
                      t = t,
                      method = "L-BFGS-B",
                      lower = c(-0.999,0),
                      upper = c(0.999, Inf),
                      control = list(maxit = 500),
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1]
  theta <- mle.values$par[2]
  net.film <- netfilter(W,rho = rho)
  

  omega.est <- omega(theta=theta,n,t)
  inv.omega.est <- solve(omega.est)
  
  beta.est <- compute_beta(net.film,omegainv=inv.omega.est,X=X,Y=Y)
  
  vt.est <- compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- compute_sigma2(vt.est,omegainv=inv.omega.est)
  ll.est <- pnam_ll_rule(netA = net.film,sigma2 = sigma2.est,
                          v=vt.est,omega=omega.est,n=n,t=t)
  
  
  noneta.vcov<- solve(pnamhess)
  colnames(noneta.vcov) <-  c("rho","theta")
  eta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%inv.omega.est%*%X)
  colnames(eta.vcov) <- paste0("beta_",1:length(beta.est))
  model <- list(eta.vcov = eta.vcov,
                noneta = hessparam,
                noneta.vcov = noneta.vcov,
                rho=rho,
                theta=theta,
                sigma2 = sigma2.est,
                beta = beta.est,
                loglike = ll.est,
                convergence = mle.values$convergence,
                type = "pnam_random_effects")
  return(model)
}













#' @title Estimate a Panel Network Autocorrelation Model with Unit-Level Random Effects and Autocorrelated Errors
#' @name  panelnamrulewae
#' @param X the covariate NT x K matrix.
#' @param Y The vector of outcomes.
#' @param W The block-diagonal NT x NT matrix
#' @param W1 The block-diagonal NT x NT matrix for the spatial error term
#' @param n the number of units
#' @param t the number of time points
#' @import Rcpp
#' @import RcppArmadillo
#' @return The maximum likelihood estimates.
#' @export

panelnamrulewae <- function(X,  #the model formula
                         Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                         W,   #this should be a list of T network N x N matrices
                         W1,
                         n,       #a vector that indexes the unit (respondents)
                         t){     #a vector that indexes the time (respondents)
  
  nt <- n*t #the total number of panel data points NT
  
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  
  #an r side wrapper here to source to c++
  pnam_ll_rule_searchR <- function(par,Y,X,W,W1,n,t){
    rho <- par[1]
    theta <- par[2]
    lambda <- par[3]
    return(pnam_ll_rulewae_search(Y,X,W,W1,n,t,rho,theta,lambda))
  }
  pnam_ll_rule_hessR <- function(par,Y,X,W,W1,n,t){
    rho <- par[1]
    theta <- par[2]
    lambda <- par[3]
    return(-pnam_ll_rulewae_search(Y,X,W,W1,n,t,rho,theta,lambda))
  }
  
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  k <- 3 #the number of searching parameters (rho,theta)
  starting <- abs(rnorm(k, mean=0,sd=0.1)) #ensuring positive values
   
  mle.values <- optim(par = starting,
                      fn = pnam_ll_rule_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      W1 = W1,
                      n = n,
                      t = t,
                      method = "L-BFGS-B",
                      lower = c(-0.999, 0, -0.999),
                      upper = c(0.999, Inf, 0.999),
                      control = list(maxit = 500),
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1]
  theta <- mle.values$par[2]
  lambda <- mle.values$par[3]
  net.film <- netfilter(W,rho = rho)

  omega.est <- omega_rulewae(theta=theta,n,t,W1=W1,lambda=lambda)
  inv.omega.est <- solve(omega.est)
  
  beta.est <- compute_beta(net.film,omegainv=inv.omega.est,X=X,Y=Y)
  
  vt.est <- compute_v(net.film,beta=beta.est,X=X,Y=Y)
  sigma2.est <- compute_sigma2(vt.est,omegainv=inv.omega.est)
  ll.est <- pnam_ll_rulewae(netA = net.film,
                            sigma2 = sigma2.est,
                            v=vt.est,omega=omega.est,n=n,t=t)

  noneta.vcov<- solve(pnamhess)
  colnames(noneta.vcov) <-  c("rho","theta","lambda")
  eta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%inv.omega.est%*%X)
  colnames(eta.vcov) <- paste0("beta_",1:length(beta.est))
  model <- list(eta.vcov = eta.vcov,
                optimization.info = mle.values,
                noneta = hessparam,
                noneta.vcov = noneta.vcov,
                rho=rho,
                lambda = lambda,
                theta=theta,
                sigma2 = sigma2.est,
                beta = beta.est,
                loglike = ll.est,
                convergence = mle.values$convergence,
                type = "pnam_random_effects_with_correlatederrors")
  return(model)
}



#' @title Estimate a Dynamic Panel Network Autocorrelation Model 
#' @name  panelnam.dynam
#' @param X the covariate NT x K matrix.
#' @param Y The vector of outcomes.
#' @param W The block-diagonal NT x NT matrix
#' @param n the number of units
#' @param t the number of time points
#' @import Rcpp
#' @import RcppArmadillo
#' @return The maximum likelihood estimates.
#' @export

panelnam.dynam <- function(X,  #the model formula
                           Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                           W,   #this should be a list of T network N x N matrices
                           n,       #a vector that indexes the unit (respondents)
                           t){     #a vector that indexes the time (respondents)
  
  nt <- n*t #the total number of panel data points NT
  #creating the lag filtering block-diagonal matrix
  L <- matrix(0,nt,nt) #an empty matrix of 0s
  start <- n+1 #the start of the lower diagonal
  for(i in start:nt){ #for all elements across sub lower diagonal
    L[i, i - n] <- 1 #add a value of 1
  }
  
  ###--------------------------------------------------------------------------###
  #    creating the L matrix
  ###--------------------------------------------------------------------------###
  
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  #an r side wrapper here to source to c++
  pnam_ll_dynam_searchR <- function(par,Y,X,W,L,n,t){
    rho <- par[1]
    delta <- par[2]
    phi <- par[3]
    return(pnam_ll_dynam_search(Y,X,W,L,n,t,rho,delta,phi))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  k <- 3 #the number of searching parameters (rho,theta)
  starting <- rnorm(k, mean=0,sd=0.001) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_ll_dynam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      L = L,
                      n = n,
                      t = t,
                      method = "L-BFGS-B",
                      lower = c(-0.999, -0.999, -0.999),
                      upper = c(0.999, 0.999, 0.999),
                      control = list(maxit = 500),
                      hessian = TRUE)
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1]
  delta <- mle.values$par[2]
  phi <- mle.values$par[3]
  hessparam <- c(rho, delta,phi)
  net.film <- netfilter_dpnam(W,L,rho,delta,phi)
  max.eigen <- max(abs(eigen((rho*W + delta*L + phi*L*W))$values))
  beta.est <- compute_beta_norandom(X,Y,net.film)
  vt.est <- net.film%*%Y - X%*%beta.est
  sigma2.est <- compute_sigma2_norandom(vt.est)
  
  ll.est <- pnam_ll_dynam(netA = net.film,
                          sigma2 = sigma2.est,
                          v=vt.est,
                          n=n,
                          t=t)
  
  noneta.vcov<- solve(pnamhess)
  
  colnames(noneta.vcov) <-  c("rho","theta","lambda")
  
  eta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%X)
  colnames(eta.vcov) <- paste0("beta_",1:length(beta.est))
  
  model <- list(eta.vcov = eta.vcov,
                optimization.info = mle.values,
                noneta = hessparam,
                noneta.vcov = noneta.vcov,
                rho=rho,
                delta = delta,
                phi=phi,
                max.eigen= max.eigen,
                sigma2 = sigma2.est,
                beta = beta.est,
                loglike = ll.est,
                convergence = mle.values$convergence,
                type = "dynamic panel network autocorrelation model")
  return(model)
}




#' @title Estimate a Dynamic Panel Network Autocorrelation Model with Autocorrelated Errors
#' @name  panelnam.dynam.wae
#' @param X the covariate NT x K matrix.
#' @param Y The vector of outcomes.
#' @param W The block-diagonal NT x NT matrix
#' @param B The block-diagonal NT x NT matrix for the spatial error term
#' @param n the number of units
#' @param t the number of time points
#' @import Rcpp
#' @import RcppArmadillo
#' @return The maximum likelihood estimates.
#' @export
panelnam.dynam.wae <- function(X,  #the model formula
                               Y,   #the dataset should be stacked as repeated cross-sections (sorted by time)
                               W,   #this should be a list of T network N x N matrices
                               B,
                               n,       #a vector that indexes the unit (respondents)
                               t){     #a vector that indexes the time (respondents)
  
  nt <- n*t #the total number of panel data points NT
  #creating the lag filtering block-diagonal matrix
  L <- matrix(0,nt,nt) #an empty matrix of 0s
  start <- n+1 #the start of the lower diagonal
  for(i in start:nt){ #for all elements across sub lower diagonal
    L[i, i - n] <- 1 #add a value of 1
  }
  
  ###--------------------------------------------------------------------------###
  #    creating the L matrix
  ###--------------------------------------------------------------------------###
  
  #this package will rely upon the nlimb function for numerical optimization of
  #the log-likelihood to find the non-(beta and sigma) parameters
  #an r side wrapper here to source to c++
  pnam_ll_dynam_searchR <- function(par,Y,X,W,L,B,n,t){
    rho <- par[1]
    delta <- par[2]
    phi <- par[3]
    lambda <- par[4]
    return(pnam_ll_dynam_wae_search(Y,X,W,L,B,n,t,rho,delta,phi,lambda))
  }
  ###--------------------------------------------------------------------------###
  #    finding the maximum likelihood estimates
  ###--------------------------------------------------------------------------###
  k <- 4 #the number of searching parameters (rho,theta)
  starting <- rnorm(k, mean=0,sd=0.001) #ensuring positive values
  mle.values <- optim(par = starting,
                      fn = pnam_ll_dynam_searchR,
                      Y = Y,
                      X = X,
                      W = W,
                      L = L,
                      B = B,
                      n = n,
                      t = t,
                      method = "L-BFGS-B",
                      lower = c(-0.999, -0.999, -0.999, -0.999),
                      upper = c(0.999, 0.999, 0.999, 0.999),
                      control = list(maxit = 500),
                      hessian = TRUE)
  
  pnamhess <- mle.values$hessian
  rho <- mle.values$par[1]
  delta <- mle.values$par[2]
  phi <- mle.values$par[3]
  lambda <- mle.values$par[4]
  hessparam <- c(rho, delta,phi,lambda)
  net.film <- netfilter_dpnam(W,L,rho,delta,phi)
  max.eigen <- max(abs(eigen((rho*W + delta*L + phi*L*W))$values))
  Bdiag <- B_dpnam(B, lambda)
  
  
  beta.est <- compute_beta_norandom_wae(X,Y,net.film,Bdiag)
  vt.est <- Bdiag%*%(net.film%*%Y - X%*%beta.est)
  sigma2.est <- compute_sigma2_norandom_wae(vt.est,Bdiag)
  
  ll.est <- pnam_ll_dynam(netA = net.film,
                          sigma2 = sigma2.est,
                          v=vt.est,
                          n=n,
                          t=t)
  
  noneta.vcov<- solve(pnamhess)
  
  colnames(noneta.vcov) <-  c("rho","theta","lambda")
  
  eta.vcov <- as.numeric(sigma2.est)*solve(t(X)%*%X)
  colnames(eta.vcov) <- paste0("beta_",1:length(beta.est))
  
  model <- list(eta.vcov = eta.vcov,
                optimization.info = mle.values,
                noneta = hessparam,
                noneta.vcov = noneta.vcov,
                rho=rho,
                delta = delta,
                phi=phi,
                max.eigen= max.eigen,
                sigma2 = sigma2.est,
                beta = beta.est,
                loglike = ll.est,
                convergence = mle.values$convergence,
                type = "dynamic panel network autocorrelation model")
  return(model)
}


