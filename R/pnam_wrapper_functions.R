#
#
#     a wrapper function for the panel network autocorrelation models 
#     in the netpanel R package
#
#     last updated on 07-01-2026 (0.0.1)
#
#


#' @title Fit a Maximum Likelihood Panel Network Autocorrelation Model (PNAM) 
#' @name mlpnam
#' @description
#' A function to estimate a Maximum Likelihood Panel Network Autocorrelation Model (MLPNAM). The function
#' allows for users to input one or multiple networks to represent the network structure. If multiple networks
#' are included the PNAM becomes an M-STAR model (see XXXX). In addition, the function allows for (1) fixed and
#' random effects and (2) the inclusion of lagged outcome and network terms. Finally, the function
#' allows for the time-varying errors to exhibit network autocorrelation, that is, a network lag and error model. 
#' 
#' @details
#' need to write this.......
#' 
#' 
#' 
#' 
#' @param reg.formula A `formula` object that is the symbolic description of the network autocorrelation model to be 
#' fitted. The dependent and exogenous regressors are extracted from the `data` argument. 
#' This is the same argument found in \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#' 
#' @param net.formula A `formula` object that is the symbolic description of the network structure for the network 
#' autocorrelation model, which allows multiple networks to be included. The dependent variable should not be defined. Each of the included covariates
#' should be a list of length t (i.e., the number of time points/panels) and each element 
#' should be an N x N matrix, where N is the number of actors/units in the network. The included networks
#' should be row-normalized. The networks are extracted from the environment.
#' 
#' @param data An object of class `data.frame` that contains the variables from the `reg.formula`, 
#' `time`, and `actor` arguments. The `data.frame` should be stacked by time, and ordered within 
#' time with respect to their position in the network adjacency matrices.  
#' 
#' @param time A `formula` object where the variable that represents the time index of the observation is on the 
#' right hand side of ~. The variable is extracted from the `data` argument.
#' 
#' @param actor A `formula` object where the variable that represents the actor/unit index of the observation is on the 
#' right hand side of ~. The variable is extracted from the `data` argument.
#' 
#' @param model The type of panel network autocorrelation model to estimate (see the Details section). "fixed" indicates
#' that a within-estimator fixed effects model should be estimated. The demeaning of the data structure is based upon the
#' `fixed.effect` argument. "random" indicates that a random-intercept model should be estimated with random effects 
#' with respect to the actors/units. "dynamic" indicates that a lagged model will be estimated, and the specific lag
#' structure is based upon the `dynamic.lag` argument.
#' 
#' @param fixed.effect **optional**. When `model` is set to "fixed", how should the data structure be demeaned? "actor" indicates
#' the data will be demeaned by the unit-specific average. "time" indicates that the data will be time demeaned (i.e.,
#' the i value at time t will be subtracted by the mean value at time t across units). "two-way" indicates that 
#' two-way fixed effects will be included (i.e., time and unit fixed effects).
#' 
#' @param dynamic.lag **optional**. When `model` is set to "dynamic", what lag structure should be included? The option 
#' "network" will add the time lag for each included network structure in the `net.formula` argument. The option 
#' "outcome" will add the time lag for the outcome variable (Y) based on the `reg.formula` argument. The "network x outcome"
#' option will add the time lags for the outcome and network variables. 
#' 
#' @param errors What is the time-varying residual structure? "idiosyncratic" indicates that the time-varying residuals are assumed to 
#' be drawn from the following distribution: N(0,sigma2), where sigma2 is estimated variance. "autocorrelated" assumes that the 
#' time-varying residuals exhibit network autocorrelation and the model become a network lag x error model. See the Details section.
#' 
#' @param autocorrelated.network **optional**. When `errors` is set to "autocorrelated", a `list` object with length t, where 
#' each element is an N x N matrix and represents the network to be used in the network error model. See the Details section.
#' 
#' @param optim.method **optional**. The numeric optimization procedure to be employed by \code{\link[stats]{optim}}. Defaults to 
#' "L-BFGS-B" and, in general, should not be changed.
#' 
#' @param optim.control **optional**. A list control options for the \code{\link[stats]{optim}} function. Please see the 
#' \code{\link[stats]{optim}} function for the set of options. 
#' 
#' @param ... Additional arguments.
#' 
#' @import stats
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom Matrix bdiag
#' @return An object of class "dream_rem" as a list containing the following components:
#' \itemize{
#'   \item \code{N} - The number of unique units in the panel network dataset.
#'   \item \code{t} - The number of panels (time periods) in the panel network dataset.
#'   \item \code{type} - The type of model estimate (e.g., random effects). 
#'   \item \code{error.structure} - The assumed structure for the time-varying errors. 
#'   \item \code{fe.type} - If a fixed effects model is estimated, the type of fixed effects demeaning.
#'   \item \code{lag.type} - If a dynamic model is estimated, the type of lag structure included in the model.
#'   \item \code{max.eigen} - The maximum absolute eigenvalue of the rho*W matrix.
#'   \item \code{lambda} - If the time-varying errors are assumed to be autocorrelated, the estimated value for lambda. 
#'   \item \code{rho} - The vector of estimated maximum likelihood rho parameters. 
#'   \item \code{vcov.rho} - The asymptotic variance-covariance of the estimated rho parameters.
#'   \item \code{coefficients} - The vector of estimated exogenous ML parameters. 
#'   \item \code{vcov} - The asymptotic variance-covariance of the exogenous ML parameters. 
#'   \item \code{sigma2} - The estimated variance of the time-varying residual. 
#'   \item \code{residuals} - The residuals of the fitted model.
#'   \item \code{fitted.values} - The fitted values from the estimated model.
#'   \item \code{logLik} - The estimated model log-likelihood.
#'   \item \code{convergence} - The converge value returned from the optim function. 0 indicates successful convergence.
#'   \item \code{netfilter} - The filtering matrix of the estimated model \eqn{I_{NT} - \sum_{j} rho_j*A_j }, where \eqn{A_j} is the jth network from the
#'   `net.formula` argument.
#'   \item \code{k} - The total number of estimated parameters.
#'   \item \code{df} - The model degrees of freedom. 
#'   \item \code{AIC} - The AIC of the estimated PNAM.
#'   \item \code{BIC} - The BIC of the estimated PNAM.
#'   \item \code{response} - The outcome vector (\eqn{Y}).
#'   \item \code{covariates} - The covariate matrix (\eqn{X}).
#'   \item \code{optim.information} - The returned list from the optim function.
#'   \item \code{error.dis} - If a random effects model is estimated, a list that contains the 
#'   residual variances for the time-varying residuals (idiosyncratic), the unit-/actor-specific
#'   residual variances (unit), and the estimated ratio of the variances (theta).
#'   \item \code{theta} - If a random effects model is estimated, the estimated value for
#'   theta (the ratio of the unit and idiosyncratic variances.)
#'   \item \code{vcov.theta} - If a random effects model is estimated, the asymptotic variance of 
#'   the estimated theta value.
#'}
#' @export
#' @seealso
#' \code{\link{coef.pnam}}, \code{\link{vcov.pnam}}, \code{\link{logLik.pnam}},
#' \code{\link{AIC}}, \code{\link{residuals.pnam}}, \code{\link{fitted.pnam}},
#' \code{\link{netimpacts}}.
#' 
#' @references 
#' add to this here....
#' 
#' 
#' @examples
#' data("network.1", package = "netpanel") #the first panel network
#' data("network.2", package = "netpanel") #the second panel network
#' data("error.network", package = "netpanel") #the third panel network
#' data("simulated.data", package = "netpanel")
#' 
#' # a random effects panel network autocorrelation model
#' re.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1 + net2,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "random")
#' summary(re.pnam)              
#' 
#' # a fixed effects panel network autocorrelation model with autocorrelated errors
#' fe.pnam.error <- mlpnam(Y~x1+x2+x3, net.formula = ~  net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor", errors = "autocorrelated",
#'                autocorrelated.network = error.net)
#' summary(fe.pnam.error)   
#' 
#' # a dynamic panel network autocorrelation model with autocorrelated errors
#' dyn.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~  net1 + net2,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "dynamic",
#'                dynamic.lag = "outcome")
#' summary(dyn.pnam)   
#' 
#' 
#' 
mlpnam <- function(reg.formula,
                 net.formula, 
                 data,
                 time,
                 actor, 
                 model = c("fixed", "random", "dynamic"),
                 fixed.effect = c("actor", "time", "two-way"),
                 dynamic.lag = c("network", "outcome", "network x outcome"),
                 errors = c("idiosyncratic", "autocorrelated"),
                 autocorrelated.network = NULL,
                 optim.method = "L-BFGS-B", 
                 optim.control = list(),
                 ...){
  
  model <- match.arg(model,c("fixed", "random", "dynamic"))  #the model estimation type
  fe.type <- match.arg(fixed.effect,  c("actor", "time", "two-way"))
  dyn.lag <- match.arg(dynamic.lag,   c("network", "outcome", "network x outcome"))
  error.structure <- match.arg(errors, c("idiosyncratic", "autocorrelated"))
  if(error.structure=="autocorrelated" & is.null(autocorrelated.network)) base::stop("Please specify the `autocorrelated.network` matrix if the errors are set to type autocorrelated.")
  if(model == "fixed") formula <- update(reg.formula, .~. - 1) #removing the intercept for the fixed effects model
  variables <- model.frame(reg.formula, data = data) #the model data frame
  Y <- model.extract(variables, "response")  #the outcome (Y) vector
  X <- model.matrix(reg.formula, data = variables)  #the X matrix
  time <- unlist(model.frame(time,data=data))#the time index vector
  t <- length(unique(time)) #the number of panels in the data frame
  unit <- unlist(model.frame(actor,data=data))#the unit index vector
  N <- length(unique(unit)) #the number of unique units in the data frame
  net.regressors <- all.vars(net.formula)
  check.nets <- all(sapply(net.regressors,FUN=exists))
  if(!isTRUE(check.nets)) base::stop("One of the networks listed in the `net.formula` argument are not in the environment. Please check these inputs.")
  network <- lapply(net.regressors,get)
  check.nets <- all(sapply(network,FUN=function(i){inherits(i,"list")}))
  if(!check.nets) base::stop("The entries in the `network` argument should be a list of N x N matrices with length t.")
  #if we only have one network in the data.frame
  if(length(net.regressors) == 1){ #there is only one network matrix
     W <- list(as.matrix(bdiag(network[[1]]))) #making the network a block diagonal matrix
     #doing a dimensional check for the network
     if(!all(dim(W[[1]]) == N*t)){
       base::stop("The block-diagonal matrix based upon the `network` argument does not have dimension equal to the 
                  the length of the response vector (Y). Please check the inputted networks and reestimate the model.") }
    names(W) <- net.regressors
  }else{ #if there are multiple networks, then each element in the list should be a block-diagonal matrix (and should be named)
    W <- lapply(network,function(net){as.matrix(bdiag(net))}) #the block-diagonal matrix
    if(is.null(names(network))) names(network) <- paste0("net.",1:length(network))
    names(W) <- net.regressors 
    dim.check <- unlist(lapply(W,function(net){all(dim(net) == N*t)}))
    if(!all(dim.check)){
      base::stop("One of the block-diagonal matrix based upon the `network` argument does not have dimension equal to the 
                  the length of the response vector (Y). Please check the inputted networks and reestimate the model.")
    }
  }
  if(!is.null(autocorrelated.network) & !inherits(autocorrelated.network,"list")) base::stop("The `autocorrelated.network` argument should be a list.")
  if(!is.null(autocorrelated.network)) errorNet <- as.matrix(bdiag(autocorrelated.network)) #making the network a block diagonal matrix
  #estimation strategy
  strategy <- paste0(model,"_",error.structure) #the estimation strategy (the modeling framework + assumed time-varying error structure)
  estimates <- switch(strategy,
                    "fixed_idiosyncratic"=panelnam.fixed(X=X,Y=Y,W=W,n=N,t=t,
                                                         id.index=unit,time.index=time,
                                                         fe=fe.type,optim.method=optim.method,
                                                         optim.control=optim.control),
                    
                    "fixed_autocorrelated"=panelnam.fixed.error(X=X,Y=Y,W=W,n=N,t=t,
                                                          id.index=unit,time.index=time,
                                                          fe=fe.type,optim.method=optim.method,
                                                          optim.control=optim.control,W2=errorNet),
                    
                    "dynamic_idiosyncratic"=panelnam.dynam(X=X,Y=Y,W=W,n=N,t=t,
                                                           dyanmic.lag=dynamic.lag,
                                                           optim.method=optim.method,
                                                           optim.control=optim.control),
                    
                    "dynamic_autocorrelated"=panelnam.dynam.error(X=X,Y=Y,W=W,n=N,t=t,
                                                            dyanmic.lag=dyanmic.lag,
                                                            optim.method=optim.method,
                                                            optim.control=optim.control,W2=errorNet),
                    
                    "random_idiosyncratic"=panelnam.random(X=X,Y=Y,W=W,n=N,t=t,
                                                          optim.method=optim.method,
                                                          optim.control=optim.control),
                    
                    "random_autocorrelated"=panelnam.random.error(X=X,Y=Y,W=W,n=N,t=t,
                                                            optim.method=optim.method,
                                                            optim.control=optim.control,W2=errorNet))
  
  
  #check the eigenvalues of the rho W matrix to check that the largest eigenvalue is less than 1 for convergence
  check.mat <- matrix(0,N*t,N*t)
  for(i in 1:length(estimates$rho)){
    check.mat <- check.mat + estimates$rho[i]*W[[i]] #if there are multiple networks used in the spatial lag model
  }
  max.value <- max(abs(eigen(check.mat, only.values = TRUE)$values))
  estimates$max.eigen <- max.value
  estimates$type <- model 
  estimates$error.structure <- error.structure 
  estimates$N <- N
  estimates$t <- t
  estimates$call <- match.call()
  estimates$fe.type <- fixed.effect
  estimates$lag.type <- dynamic.lag
  class(estimates) <- "pnam"
  return(estimates)
}

#' Print Method for `pnam` Panel Network Autocorrelation Model
#'
#' @param x An object of class "pnam".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'pnam' object.
#' @export
print.pnam  <- function(x,digits=5,...) {
  if(x$type=="random") cat("Maximum Likelihood Panel Network Autocorrelation Model with Unit-Level Random Effects\n")
  if(x$type=="fixed"){ 
    cat("Maximum Likelihood Fixed Effects Panel Network Autocorrelation Model\n")
    cat(paste0("Fixed Effect Type: ",x$fe.type,"\n"),sep="")
  }
  if(x$type=="dynamic"){ 
    cat("Maximum Likelihood Dynamic Panel Network Autocorrelation Model\n")
    cat(paste0("Lag Structure: ",x$lag.type,"\n"),sep="")
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nNetwork Autocorrelation Lag Parameters:\n")
  coefs <- round(x$rho,digits)
  print(coefs)
  cat("\nCoefficient Estimates:\n")
  coefs <- round(x$coefficients,digits)
  print(coefs)
  cat("\nError Components:\n")
  cat(paste0("Time-Varying Error Structure: ",x$error.structure,"\n"))
  
  if(x$type == "random"){
    cat("Unit-Specific SD: ", round(x$error.dis$unit,digits), " (variance: ",round(sqrt(x$error.dis$unit),digits), ")\n",sep="")
    cat("Idiosyncratic SD: ", round(x$error.dis$idiosyncratic,digits), " (variance: ",round(sqrt(x$error.dis$idiosyncratic),digits), ")\n",sep="")
    cat("Theta: ", round(x$theta,digits), " (SE: ",round(sqrt((x$vcov.theta)),digits), ")\n",sep="")
  }
  if(x$type != "random"){
    cat("Idiosyncratic SD: ", round(x$sigma2,digits), "(variance: ",round(sqrt(x$sigma2),digits), ")\n",sep="")
  }
  if(x$error.structure == "autocorrelated"){
    cat("Network Error Autocorrelation Parameter:\n")
    coefs <- round(x$lambda,digits)
    print(coefs)
  }
  invisible(x)
}


#' Summary Method for `pnam` Objects
#'
#' Summarizes the results of a panel network autocorrelation model.
#'
#' @param object An object of class "pnam".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return A list of summary statistics for the panel network autocorrelation model, such 
#' as the estimated log-likelihood and the estimated parameters.
#' @export
summary.pnam <- function(object,digits=5,...) {
  se <- sqrt(diag(object$vcov))
  z<-object$coefficients/se
  p<-2*pnorm(abs(z), lower.tail = FALSE)
  exogenous.coef_table <- cbind(
    `Estimate` = round(object$coefficients,digits),
    `Std. Error` = round(se,digits),
    `z value` = round(z,digits),
    `Pr(>|z|)` = round(p,digits)
  )
  se <- sqrt(diag(object$vcov.rho))
  z<-object$rho/se
  p<-2*pnorm(abs(z), lower.tail = FALSE)
  netlag.coef_table <- cbind(
    `Estimate` = round(object$rho,digits),
    `Std. Error` = round(se,digits),
    `z value` = round(z,digits),
    `Pr(>|z|)` = round(p,digits)
  )
  
  neterror.coef_table<-NULL
  if(object$error.structure == "autocorrelated"){
    se <- sqrt(object$vcov.lambda)
    z<-object$lambda/se
    p<-2*pnorm(abs(z), lower.tail = FALSE)
    neterror.coef_table <- cbind(
      `Estimate` = round(object$lambda,digits),
      `Std. Error` = round(se,digits),
      `z value` = round(z,digits),
      `Pr(>|z|)` = round(p,digits)
    )
    rownames(neterror.coef_table) <- "lambda"
    netlag.coef_table <- rbind(netlag.coef_table,neterror.coef_table)
  }
  
  
  if(object$type=="random"){
    
    error.comp <-cbind(sd = round(c(object$error.dis$idiosyncratic,
                                    object$error.dis$unit),digits),
                       var = round(sqrt(c(object$error.dis$idiosyncratic,
                                     object$error.dis$unit)),digits))
    rownames(error.comp)<- c("idiosyncratic:", "unit-specific:")
    se <- sqrt(object$vcov.theta)
    z <- object$theta/se 
    p<-2*pnorm(abs(z), lower.tail = FALSE)
    theta.se <- round(se,digits)
    theta.est <- round(object$theta,digits)
    theta.p <- round(p,digits)
    var<-round(sqrt(c(object$error.dis$idiosyncratic,object$error.dis$unit)),digits)
    icc <- (var[2]/(var[1]+var[2]))
    
  }else{
    error.comp <-data.frame(sd = round(c(object$sigma2),digits),
                       var = round(sqrt(c(object$sigma2)),digits))
    rownames(error.comp)<-c("idiosyncratic:")
    theta.se <- NULL
    theta.est <- NULL
    theta.p <- NULL
    icc<-NULL
  }
  
  if(object$type=="random") type <- "Maximum Likelihood Panel Network Autocorrelation Model with Unit-Level Random Effects\n"
  if(object$type=="fixed") type <- paste0("Maximum Likelihood Fixed Effects Panel Network Autocorrelation Model\n Fixed Effect Type: ",
                                          object$fe.type,"\n")
  if(object$type=="dynamic") type <- paste0("Maximum Likelihood Dynamic Panel Network Autocorrelation Model\n FLag Structure: ",
                                             object$lag.type,"\n")

  
  res <- list(
    type=type,
    call = object$call,
    exo.coef = exogenous.coef_table,
    lag.coef = netlag.coef_table,
    logLik = round(object$logLik,digits),
    k=object$k,
    df=object$df,
    nt = object$N*object$t,
    N = object$N,
    t = object$t,
    AIC = round(object$AIC,digits),
    BIC = round(object$BIC,digits),
    type = object$rem.type,
    error.comp=error.comp,
    theta.se=round(theta.se,digits),
    theta.p=round(theta.p,digits),
    theta.est=round(theta.est,digits),
    error.structure=object$error.structure,
    resid = object$residuals,
    sigma2=object$sigma2,
    convergence=object$convergence,
    iterations=object$optim.information$counts[1],
    icc=round(icc,digits),
    max.eigen=object$max.eigen
  )
  class(res) <- "summary.pnam"
  return(res)
}

#' Print Method for `summary.pnam` Panel Network Autocorrelation Model Summary
#'
#' @param x An object of class "summary.pnam".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the summary of a 'pnam' object.
#' @export
print.summary.pnam <- function(x, digits = 5, ...){
  cat(x$type)
  cat("\nCall:\n")
  print(x$call)
  cat("\nResiduals:\n")
  res_sum <- summary(as.numeric(x$resid))
  res_sum<-res_sum[names(res_sum)!="Mean"]
  names(res_sum) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(round(res_sum, digits = digits))
  
  cat("\nPanel Error Variance Components:\n")
  print(x$error.comp)
  if(!is.null(x$theta.est)) cat("theta: ",x$theta.est," (SE: ",x$theta.se,"; p= ",x$theta.p,")\n",sep = "")
  if(!is.null(x$icc)) cat("icc: ",x$icc,"\n",sep = "")
  
  if(!is.null(x$theta.coef_table)){
    cat("\nRandom Effects Theta Estimate:\n")
    printCoefmat(x$theta.coef_table, P.values = TRUE, has.Pvalue = TRUE)
  }
 
  cat("\nNetwork Autocorrelation Parameters:\n")
  printCoefmat(x$lag.coef, P.values = TRUE, has.Pvalue = TRUE)
  
  cat("\nML Regression Coefficients:\n")
  printCoefmat(x$exo.coef, P.values = TRUE, has.Pvalue = TRUE)
  
  
  cat("\nModel Fit Information:\n")
  cat("Number of Observations: ", x$nt,"; N: ",x$N, "; T: ",x$t ,"\n",sep = "")
  cat("Parameters Estimated:", x$k,"\n")
  cat("Log-Likehoood: ", x$logLik," (df=", x$df,"); AIC: ",x$AIC,"; BIC: ", x$BIC,"\n", sep = "")
  cat("Residual SD: ", round(x$sigma2,digits), "\n",sep = "")
  cat("Model Convergence: ", ifelse(x$convergence==0,"Yes","No"),"\n", sep = "")
  cat("Search Iterations: ", x$iterations, "\n",sep = "")
  cat("Largest eigenvalue of rho*W: ", x$max.eigen, "\n",sep = "")
  
}



#' Extract the ML parameter estimates from Panel Network Autocorrelation Model Fits
#'
#' This function extracts the Maximum Likelihood (ML) parameter estimates from estimated
#' panel network autocorrelation Model fits.
#'
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @export
coef.pnam <- function(object,...){
  object$coefficients
}


#' Extract variance-covariance matrix from `pnam` Fits
#'
#' This function extracts the variance-covariance matrix from estimated
#' `pnam` model fits for the exogenous covariates. For the endogenous 
#' variance matrices (such as that for rho), please extract those 
#' from the model fit. For example, the variance matrix for rho
#' can be extracted as `model$vcov.rho`. 
#'
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @export
vcov.pnam <- function(object,...){
  object$vcov
}

#' Extract the model log-likelihood from  `pnam` Fits
#'
#' This function extracts the model loglikelhood from estimated
#' panel netowrk autocorrelation model fits.
#'
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @param REML From the generic `logLik` function. Set to FALSE and does not
#' need to changed by the user.
#' @export
logLik.pnam <- function(object,...,REML = FALSE){
  val <- object$logLik
  attr(val, "df") <- object$k
  attr(val, "nobs") <- length(object$residuals)
  class(val) <- "logLik"
  val
}

#' Residuals for `pnam` model fits
#'
#' This function returns the model residuals for a `pnam` model fit.
#' 
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @export
residuals.pnam <- function(object,...){
  object$residuals
}



#' Fitted values from `pnam` model fits
#'
#' This function returns the model fitted values for a `pnam` model fit.
#' 
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @export
fitted.pnam <- function(object,...){
  object$fitted.values
}




#' The network impacts method (i.e., the matrix of impacts from the estimated model)
#'
#' This function returns the estimated network impacts from the
#' "pnam" fitted object.
#' 
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @export
netimpacts <- function(object,digits=6){
  coefs <- (object$coefficients)
  vars <- names(coefs)
  vars <- vars[vars!="(Intercept)"] #dropping the intercept
  direct<-rep(0,length(vars))
  indirect<-rep(0,length(vars))
  Ainv<-solve(object$netfilter)
  for(i in 1:length(vars)){
    impact.matrix <- Ainv*(coefs[i])
    direct[i] <- 1/nrow(Ainv)*sum(diag(impact.matrix))
    diag(impact.matrix)<-0
    indirect[i] <- 1/nrow(Ainv)*sum(impact.matrix)
  }
  total<-direct+indirect
  res <- data.frame(Effect = vars,
                    Direct = round(direct,digits),
                    Indirect = round(indirect,digits),
                    Total = round(total,digits))
  res
}







