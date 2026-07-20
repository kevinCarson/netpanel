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
#' A function that contains Maximum Likelihood estimation routines to fit panel network autocorrelation models (MLPNAM). 
#' The function allows for (1) fixed and random actor-level effects, (2) the inclusion of time lagged outcome and network terms,
#' (3) multiple networks to construct the network (spatial) filtering matrix, and (4) for the time-varying
#' errors to exhibit network autocorrelation, that is, for the model to become a network mixed effects model (i.e., error + lag model). Please 
#' see the Details section for more information.
#' 
#' @details
#' This function fits a panel network autocorrelation model with time-varying social adjacency matrices
#' via Maximum Likelihood Estimation. The estimated model takes the form:
#' \deqn{Y_t = \rho_1 W^{t}_{1}Y_t + \cdots + \rho_k W^{t}_{k}Y_t + X_t\beta + v_t, \ \ t = 1,2,\cdots,T }
#' \deqn{v_t = \mu_N + e_t,\ \ t = 1,2,\cdots,T   }
#' where \eqn{Y_t} is the outcome vector for \eqn{N} actors at time \eqn{t}, \eqn{\rho_k W^{t}_{k}Y_t} is the 
#' network autocorrelation term for the \eqn{k^{th}} social network at time \eqn{t} (taken from the `net.formula` 
#' argument), \eqn{X_t} are the set of exogenous time-varying (and invariant) regressors (taken from the `reg.formula` 
#' argument). Additionally, \eqn{v_t} is the time-varying error term that is the sum
#' of the time-invariant unit-specific unobserved factors, \eqn{\mu_N}, thought to explain \eqn{Y_t} and \eqn{e_t} are the 
#' time- and unit-varying unobserved factors. Dependent upon the `model` argument, the function will (1) estimate a fixed-effects 
#' network regression, where the within estimator demeans the data by time, actor, or both, (2) estimate a random effects
#' network regression, or (3) estimate a dynamic panel network autocorrelation where the time lags of the outcome
#' and network variables are added to the set of exogenous covariates. Please see the set of references below for more 
#' information of network (spatial) cross-sectional and panel models and their respective estimation strategies (e.g., fixed effects).
#' 
#' Moreover, the model can also estimate the network mixed effects model, where the estimated formula is:
#' \deqn{Y_t = \rho_1 W^{t}_{1}Y_t + \cdots + \rho_k W^{t}_{k}Y_t + X_t\beta + v_t, \ \ t = 1,2,\cdots,T }
#' \deqn{v_t = \mu_N + e_t,\ \ t = 1,2,\cdots,T   }
#' \deqn{e_t = \lambda W^e_t e_t + \epsilon_t, \ \ t = 1,2,\cdots,T }
#' where \eqn{\lambda} is the error autocorrelation parameter, \eqn{W^e_t} is the error social adjacency matrix (if applicable,
#' the network included in the `autocorrelated.network` argument), and \eqn{\epsilon_t} are a set of i.i.d normal errors at time \eqn{t}.
#' 
#' Internally, the function follows the two-step computation strategy discussed in Millo and Piras (2012). In the first step, the set of 
#' parameter estimates for \eqn{\rho}, \eqn{\lambda}, and \eqn{\theta} are found via numerical optimization routines. Importantly, within each search,
#' the values for \eqn{\beta} and \eqn{\sigma^2} are updated via their GLS estimators. In the second step, once the ML parameter estimates for \eqn{\rho}, 
#' \eqn{\lambda}, and \eqn{\theta} are found, the final GLS estimates for \eqn{\beta} and \eqn{\sigma^2} are computed. The asymptotic 
#' variance-covariance matrix for \eqn{\rho}, \eqn{\lambda}, and \eqn{\theta} are the inverse of the negative Hessian matrix and the 
#' variance-covariance matrix for \eqn{\beta} is found via the standard GLS variance estimator. For sake of brevity, please see the below listed
#' references for more information on these estimators.
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
#' time with respect to their position in the network adjacency matrices. If you are unsure how the 
#' stacking of observations should be constructed, please see the structure of the
#'  \code{\link{simulated.data}} `netpanel` data object.
#' 
#' @param time A `formula` object where the variable that represents the time index of the observation is on the 
#' right hand side of ~. The variable is extracted from the `data` argument.
#' 
#' @param actor A `formula` object where the variable that represents the actor/unit index of the observation is on the 
#' right hand side of ~. The variable is extracted from the `data` argument.
#' 
#' @param model The type of panel network autocorrelation model to be estimated (see the Details section). "fixed" indicates
#' that a within-estimator fixed effects model should be estimated. The demeaning of the data structure is based upon the
#' `fixed.effect` argument. "random" indicates that a random-intercept model should be estimated with random effects 
#' with respect to the actors/units. "dynamic" indicates that a lagged model will be estimated, and the specific lag
#' structure is based upon the `dynamic.lag` argument.
#' 
#' @param fixed.effect **optional**. When `model` is set to "fixed", how should the data structure be demeaned? "actor" indicates
#' the data will be demeaned by the unit-specific average. "time" indicates that the data will be time demeaned (i.e.,
#' the \eqn{x_{i1t}} value at time \eqn{t} will be subtracted by the mean of \eqn{x_{1t}}, which is the vector of 
#' \eqn{x_1} values for all units at time \eqn{t}). "two-way" indicates that 
#' two-way fixed effects will be included (i.e., time and unit fixed effects).
#' 
#' @param dynamic.lag **optional**. When `model` is set to "dynamic", what lag structure should be included? The option 
#' "network" will add the time lag for each included network structure in the `net.formula` argument. The option 
#' "outcome" will add the time lag for the outcome variable (Y) based on the `reg.formula` argument. The "network x outcome"
#' option will add the time lags for the outcome and network variables. 
#' 
#' @param errors What is the assumed structure of time-varying errors? "idiosyncratic" indicates that the time-varying errors are assumed to 
#' be drawn from the following distribution: \eqn{N(0,\sigma^2)}, where \eqn{\sigma^2} is the variance. "autocorrelated" assumes that the 
#' time-varying errors exhibit network autocorrelation and the model will become a network mixed effects model. See the Details section.
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
#' @return An object of class `pnam` as a `list` containing the following elements:
#' \itemize{
#'   \item \code{N} - The number of unique units in the panel network dataset.
#'   \item \code{t} - The number of panels (time periods) in the panel network dataset.
#'   \item \code{type} - The type of model estimated per the `model` argument.
#'   \item \code{error.structure} - The assumed structure for the time-varying errors per the `errors` argument.
#'   \item \code{fe.type} - If a fixed effects model is estimated, the type of fixed effects demeaning per the `fixed.effect` argument.
#'   \item \code{lag.type} - If a dynamic model is estimated, the type of lag structure included in the model per the `dynamic.lag` argument
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
#'   \item \code{netfilter} - The filtering matrix of the estimated model (\eqn{I_{NT} - \{\sum_{j} rho_j*W^j \}}), where \eqn{W^j} is the jth network from the
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
#' Anselin, Luc. 1988. *Spatial Econometrics: Methods and Models*. Dordrecht: Springer Netherlands.
#' 
#' Cook, Scott J., Jude C. Hays and Robert J. Franzese. 2023. “STADL Up! The Spatiotemporal
#' Autoregressive Distributed Lag Model for TSCS Data Analysis.” *American Political Science Review* 117(1):59–79.
#' 
#' Doreian, Patrick. 1982. “Maximum Likelihood Methods for Linear Models: Spatial Effect and
#' Spatial Disturbance Terms.” *Sociological Methods & Research* 10(3):243–269.
#' 
#' Doriean, Patrick, Klaus Teuter and Chi-Hsein Wang. 1984. “Network Autocorrelation Models:
#' Some Monte Carlo Results.” *Sociological Methods & Research* 13(2):155–200.
#' 
#' Duxbury, Scott. 2023. *Longitudinal Network Models*. Vol. 192 of *Quantitative Applications in the Social Sciences* SAGE Publications.
#' 
#' Hays, Jude C., Aya Kachi, and Robert J. Franzese. 2010. "A Spatial Model Incorporating 
#' Dynamic, Endogenous Network Interdependence: A Political Science Application." 
#' *Statistical Methodology* 7(3): 406-428.
#' 
#' Lee, Lung-fei and Jihai Yu. 2010. “Estimation of Spatial Autoregressive Panel Data Models with
#' Fixed Effects.” *Journal of Econometrics* 154(2):165–185.
#' 
#' Lee, Lung-fei and Jihai Yu. 2012a. “QML Estimation of Spatial Dynamic Panel Data Models with
#' Time Varying Spatial Weights Matrices.” *Spatial Economic Analysis* 7(1):31–74.
#' 
#' Lee, Lung-fei and Jihai Yu. 2012b. “Spatial Panels: Random Components Versus Fixed Effects.”
#' *International Economic Review* 53(4):1369–1412.
#' 
#' Millo, Giovanni. 2014. "Maximum Likelihood Estimation of Spatially and Serially Correlated 
#' Panels with Random Effects." *Computational Statistics & Data Analysis* 71: 914-933.
#' 
#' Millo, Giovanni and Gianfranco Piras. 2012. "splm: Spatial Panel Data Models in R." 
#' *Journal of Statistical Software* 47: 1-38.
#' 
#' Neuman, Eric J. and Mark S. Mizruchi. 2010. “Structure and Bias in the Network Autocorrelation
#' Model.” *Social Networks* 32(4):290–300.
#' 
#' Wang, Wei and Jihai Yu. 2015. "Estimation of Spatial Panel Data Models with Time Varying 
#' Spatial Weights Matrices.” *Economic Letters* 128:95-99.
#' 
#' 
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("net2", package = "netpanel") #the second panel network
#' data("net3", package = "netpanel") #the third panel network
#' data("simulated.data", package = "netpanel") #the simulated variables 
#' 
#' # a random effects panel network autocorrelation model
#' re.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "random")
#' summary(re.pnam)              
#' 
#' # a fixed effects panel network autocorrelation model 
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~  net1 + net2,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' summary(fe.pnam)   
#' 
#' # a dynamic panel network autocorrelation model with autocorrelated errors
#' dyn.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~  net1 + net2,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "dynamic",
#'                dynamic.lag = "outcome")
#' summary(dyn.pnam)   
#' 
#' 
#' \donttest{
#' # a fixed effects panel network autocorrelation model with autocorrelated errors
#' fe.pnam.error <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                errors = "autocorrelated",
#'                autocorrelated.network = net3)
#' summary(fe.pnam.error)   
#' }
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
  if(model=="fixed") X <- X[,!colnames(X) == "(Intercept)"]
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
                                                            dyanmic.lag=dynamic.lag,
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

#' Print Method for `pnam` model objects.
#'
#' @param x An object of class `pnam`.
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a `pnam` object.
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
    cat("Unit-Specific SD: ", round(sqrt(x$error.dis$unit),digits), " (variance: ",round((x$error.dis$unit),digits), ")\n",sep="")
    cat("Idiosyncratic SD: ", round(sqrt(x$error.dis$idiosyncratic),digits), " (variance: ",round((x$error.dis$idiosyncratic),digits), ")\n",sep="")
    cat("Theta: ", round(x$theta,digits), " (SE: ",round(sqrt((x$vcov.theta)),digits), ")\n",sep="")
  }
  if(x$type != "random"){
    cat("Idiosyncratic SD: ", round(sqrt(x$sigma2),digits), " (variance: ",round((x$sigma2),digits), ")\n",sep="")
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
#' Summarizes the results of a fitted `pnam` model object.
#'
#' @param object An object of class `pnam`.
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
  }
  
  
  if(object$type=="random"){
    
    error.comp <-cbind(sd = round(sqrt(c(object$error.dis$idiosyncratic,
                                    object$error.dis$unit)),digits),
                       var = round((c(object$error.dis$idiosyncratic,
                                     object$error.dis$unit)),digits))
    rownames(error.comp)<- c("idiosyncratic:", "unit-specific:")
    se <- sqrt(object$vcov.theta)
    z <- object$theta/se 
    p<-2*pnorm(abs(z), lower.tail = FALSE)
    theta.se <- round(se,digits)
    theta.est <- round(object$theta,digits)
    theta.p <- round(p,digits)
    var<-round((c(object$error.dis$idiosyncratic,object$error.dis$unit)),digits)
    icc <- round((var[2]/(var[1]+var[2])),digits)
    
  }else{
    error.comp <-data.frame(sd =  round(sqrt(c(object$sigma2)),digits),
                       var = round((c(object$sigma2)),digits))
    rownames(error.comp)<-c("idiosyncratic:")
    theta.se <- NULL
    theta.est <- NULL
    theta.p <- NULL
    icc<-NULL
  }
  
  if(object$type=="random") type <- "Maximum Likelihood Panel Network Autocorrelation Model with Unit-Level Random Effects\n"
  if(object$type=="fixed") type <- paste0("Maximum Likelihood Fixed Effects Panel Network Autocorrelation Model\nFixed effect type: ",
                                          object$fe.type,"\n")
  if(object$type=="dynamic") type <- paste0("Maximum Likelihood Dynamic Panel Network Autocorrelation Model\nLag structure: ",
                                             object$lag.type,"\n")

  
  res <- list(
    type=type,
    call = object$call,
    exo.coef = exogenous.coef_table,
    lag.coef = netlag.coef_table,
    error.coef =neterror.coef_table,
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
    theta.se=(theta.se),
    theta.p=(theta.p),
    theta.est=(theta.est),
    error.structure=object$error.structure,
    resid = object$residuals,
    sigma2=object$sigma2,
    convergence=object$convergence,
    iterations=object$optim.information$counts[1],
    icc=(icc),
    max.eigen=object$max.eigen
  )
  class(res) <- "summary.pnam"
  return(res)
}

#' Print Method for `summary.pnam` Panel Network Autocorrelation Model Summary
#'
#' @param x An object of class `summary.pnam`.
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the summary of a `pnam` object.
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
  
  if(!is.null(x$theta.est)){
  cat("\nPanel error variance components:\n")
  print(x$error.comp)
  if(!is.null(x$theta.est)) cat("theta: ",x$theta.est," (SE: ",x$theta.se,"; p= ",x$theta.p,")\n",sep = "")
  if(!is.null(x$icc)) cat("icc: ",x$icc,"\n",sep = "")
  }
  
  if(!is.null(x$theta.coef_table)){
    cat("\nRandom effects theta estimate:\n")
    printCoefmat(x$theta.coef_table, P.values = TRUE, has.Pvalue = TRUE,signif.legend = FALSE)
  }
 
  cat("\nNetwork autocorrelation parameters:\n")
  printCoefmat(x$lag.coef, P.values = TRUE, has.Pvalue = TRUE,signif.legend = FALSE)
  
  if(!is.null(x$error.coef)){
  cat("\nNetwork error parameters:\n")
  printCoefmat(x$error.coef, P.values = TRUE, has.Pvalue = TRUE,signif.legend = FALSE)
  }
  
  cat("\nML regression coefficients:\n")
  printCoefmat(x$exo.coef, P.values = TRUE, has.Pvalue = TRUE)
  
  
  cat("\nModel fit information:\n")
  cat(" -> number of observations: ", x$nt,"; N: ",x$N, "; T: ",x$t ,"\n",sep = "")
  cat(" -> parameters estimated:", x$k,"\n")
  cat(" -> log-likehoood: ", x$logLik," (df=", x$df,"); AIC: ",x$AIC,"; BIC: ", x$BIC,"\n", sep = "")
  cat(" -> residual SD: ", round(sqrt(x$sigma2),digits), "\n",sep = "")
  cat(" -> optim convergence: ", ifelse(x$convergence==0,"yes","no"),"\n", sep = "")
  cat(" -> search iterations: ", x$iterations, "\n",sep = "")
  cat(" -> largest eigenvalue of rho*W: ", x$max.eigen, "\n",sep = "")
  
}



#' Extract the ML parameter estimates from `pnam` model fits
#'
#' This function extracts the Maximum Likelihood (ML) parameter estimates for the 
#' exogenous covariates from a `pnam` model fit.
#'
#' @param object An object of class "pnam".
#' @param ... Additional arguments for other methods.
#' @return The vector of estimated parameter values for the exogenous covariates.
#' @export
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' #' a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' coef(fe.pnam) #the mle parameter values
coef.pnam <- function(object,...){
  object$coefficients
}


#' Extract the asymptotic variance-covariance matrix from `pnam` model fits
#'
#' This function extracts the asymptotic variance-covariance matrix from
#' `pnam` model fits for the exogenous covariates. For the endogenous 
#' variance matrices (such as the variance-covariance matrix for \eqn{\rho}), please extract those 
#' from the model fit. For example, the variance matrix for \eqn{\rho}
#' can be extracted as `model$vcov.rho`, where `model` is the `pnam` object name.
#'
#' @param object An object of class `pnam`.
#' @param ... Additional arguments for other methods.
#' @return The asymptotic variance-covariance matrix for the exogenous covariates.
#' @export
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' 
#' # a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' vcov(fe.pnam) #the asymptotic variance-covariance matrix
vcov.pnam <- function(object,...){
  object$vcov
}

#' Extract the estimated model log-likelihood from `pnam` model fits
#'
#' This function extracts the estimated model loglikelhood from a `pnam` model fit.
#'
#' @param object An object of class `pnam`.
#' @param ... Additional arguments for other methods.
#' @param REML From the generic `logLik` function. Set to FALSE and does not
#' need to changed by the user.
#' @return The estimated model log-likelihood.
#' @export
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' 
#' # a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' logLik(fe.pnam) #the estimated model log-likelihood.

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
#' @param object An object of class `pnam`.
#' @param ... Additional arguments for other methods.
#' @return The vector of `pnam` model fit residuals.
#' @export
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' 
#' # a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' residuals(fe.pnam) #the vector of model residuals
#' 
#' 
residuals.pnam <- function(object,...){
  object$residuals
}



#' Fitted values from `pnam` model fits
#'
#' @description
#' This function returns the model fitted values, \eqn{\hat{y}}, for a fitted `pnam` model. The function returns 
#' two types of fitted values to the user: linear regression fitted values and network regression fitted values. The
#' network regression fitted values are formulated as: \eqn{\hat{y} = (I_{NT} - \sum_i \hat{\rho}_i W^i)^{-1}(X\hat{\beta})}, where \eqn{W^i} 
#' is the \eqn{i^{th}} included network. In comparison, the linear regression fitted values are formulated as:
#' \eqn{\hat{y} = X\hat{\beta}}.
#' 
#' @param object An object of class `pnam`.
#' @param ... Additional arguments for other methods.
#' @return A `data.frame` object that stores the fitted values from a `pnam` model fit, where the column "ols.fitted" stores the
#' linear regression fitted values and "network.fitted" stores the network regression fitted values (see above).
#' @export
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' 
#' # a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' fitted(fe.pnam) #the data.frame of fitted values
#' 
fitted.pnam <- function(object,...){
  data.frame(ols.fitted = object$covariates%*%object$coefficients,
             network.fitted =  object$fitted.values)
}



#' @title Compute the direct, indirect, and total impacts from a fitted `pnam` object
#' @name netimpacts
#' @param object An object of class `pnam`.
#' @param return.impact.matrix TRUE/FALSE. If TRUE, the function will also return the estimated
#' impact matrix for each estimated covariate. Set to FALSE by default. 
#' @param digits The number of digits to round the estimates after the decimal point.
#' @export
#' @return A `data.frame` object that stores the direct, indirect, and total impact for each estimated effect in the fitted `pnam` model.
#' @description
#' This function returns the estimated average direct, indirect, and total 
#' impacts from a fitted panel network autocorrelation model. 
#' @details
#' In classical linear regression, the marginal effect of a regressor, \eqn{x_1}, on the 
#' outcome, \eqn{Y}, is \eqn{\hat{\beta_1}}. In panel network (and spatial) effects models, the 
#' marginal effect is:
#' \deqn{\frac{\partial Y}{\partial x_1} = (I_{NT} - \sum_k \hat{\rho_k}W^{k})^{-1}\hat{\beta_1}}
#' The direct impact is the average of the diagonal value of the above matrix, the
#' indirect impact is the average of the row sum of the off-diagonal values, and 
#' the total impact is sum of both. The function will also return the impact matrix for each 
#' estimated effect if `return.impact.matrix` is set to TRUE. 
#' 
#' @examples
#' data("net1", package = "netpanel") #the first panel network
#' data("simulated.data", package = "netpanel") #the simulated variables
#' 
#' # a fixed effects panel network autocorrelation model
#' fe.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1,
#'                data =  simulated.data, time = ~panel,
#'                actor = ~unit, model = "fixed",
#'                fixed.effect = "actor")
#' netimpacts(fe.pnam) #compute the impact values
#' 
netimpacts <- function(object,
                       return.impact.matrix=FALSE,
                       digits=5){
  coefs <- (object$coefficients)
  coefs <- coefs[names(coefs)!="(Intercept)"] #dropping the intercept
  direct<-rep(0,length(coefs))
  indirect<-rep(0,length(coefs))
  Ainv<-solve(object$netfilter)
  for(i in 1:length(coefs)){
    impact.matrix <- Ainv*(coefs[i])
    direct[i] <- 1/nrow(Ainv)*sum(diag(impact.matrix))
    diag(impact.matrix)<-0
    indirect[i] <- 1/nrow(Ainv)*sum(impact.matrix)
  }
  total<-direct+indirect
  res <- data.frame(Regressor = names(coefs),
                    Direct = round(direct,digits),
                    Indirect = round(indirect,digits),
                    Total = round(total,digits))
  rownames(res) <- NULL
  if(return.impact.matrix){
    impact.mat <- list()
    for(i in 1:length(coefs)){
      impact.mat[[i]] <- Ainv*(coefs[i])
    }
    names(impact.mat) <- names(coefs)
    res <- list(impacts = res,
                impact.mat=impact.mat)
  }
  res
}







