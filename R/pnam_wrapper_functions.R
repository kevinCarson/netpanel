#
#
#     a wrapper function for the panel network autocorrelation models that
#     completes the scoring iterations following Millo, where the GLS 
#     steps occur within the search for the network lag and error parameters
#
#
#     last updated on 07-01-2026 (0.0.1)
#
#
#
pnam <- function(formula,
                 data,
                 network, 
                 time,
                 actor, 
                 model = c("fixed", "random", "dynamic"),
                 fixed.effect = c("actor", "time", "two-way"),
                 dyanmic.lag = c("network", "outcome", "network x outcome"),
                 errors = c("idiosyncratic", "autocorrelated"),
                 autocorrelated.network = NULL,
                 optim.method = "L-BFGS-B", 
                 optim.control = list()){
  
  model <- match.arg(model,c("fixed", "random", "dynamic"))  #the model estimation type
  fe.type <- match.arg(fixed.effect,  c("actor", "time", "two-way"))
  dyn.lag <- match.arg(dyanmic.lag,   c("network", "outcome", "network x outcome"))
  error.structure <- match.arg(errors, c("idiosyncratic", "autocorrelated"))
  if(error.structure=="autocorrelated" & is.null(autocorrelated.network)) base::stop("Please specify the `autocorrelated.network` matrix if the errors are set to type autocorrelated.")
  if(model == "fixed") formula <- update(formula, .~. - 1) #removing the intercept for the fixed effects model
  variables <- model.frame(formula, data = data) #the model data frame
  Y <- model.extract(variables, "response")  #the outcome (Y) vector
  X <- model.matrix(formula, data = variables)  #the X matrix
  time <- unlist(model.frame(time,data=data))#the time index vector
  t <- length(unique(time)) #the number of panels in the data frame
  unit <- unlist(model.frame(actor,data=data))#the unit index vector
  N <- length(unique(unit)) #the number of unique units in the data frame
  if(!inherits(network,"list")) base::stop("The `network` argument should be a list.")
  #if we only have one network in the data.frame
  if(length(network) == t){
     W <- list(as.matrix(Matrix::bdiag(network))) #making the network a block diagonal matrix
     #doing a dimensional check for the network
     if(!all(dim(W[[1]]) == N*t)){
       base::stop("The block-diagonal matrix based upon the `network` argument does not have dimension equal to the 
                  the length of the response vector (Y). Please check the inputted networks and reestimate the model.") }
    names(W) <- "net.1"
  }else{ #if there are multiple networks, then each element in the list should be a block-diagonal matrix (and should be named)
    W <- lapply(network,function(net){as.matrix(Matrix::bdiag(net))}) #the block-diagonal matrix
    if(is.null(names(network))) names(network) <- paste0("net.",1:length(network))
    names(W) <- names(network) 
    dim.check <- unlist(lapply(W,function(net){all(dim(net) == N*t)}))
    if(!all(dim.check)){
      base::stop("One of the block-diagonal matrix based upon the `network` argument does not have dimension equal to the 
                  the length of the response vector (Y). Please check the inputted networks and reestimate the model.")
    }
  }
  if(!is.null(autocorrelated.network) & !inherits(autocorrelated.network,"list")) base::stop("The `autocorrelated.network` argument should be a list.")
  if(!is.null(autocorrelated.network)) errorNet <- as.matrix(Matrix::bdiag(autocorrelated.network)) #making the network a block diagonal matrix
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
                                                           dyanmic.lag=dyanmic.lag,
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
  estimates$lag.type <- dyanmic.lag
  class(estimates) <- "pnam"
  return(estimates)
}


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
logLik.dream_rem <- function(object,...,REML = FALSE){
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




# the spatial impacts method (i.e., the matrix of impacts from the estimated model)
netimpacts <- function(object,digits=6){
  vars <- names(object$coefficients)
  vars <- vars[vars!="(Intercept)"] #dropping the intercept
  direct<-rep(0,length(vars))
  indirect<-rep(0,length(vars))
  Ainv<-solve(object$netfilter)
  for(i in 1:length(vars)){
    impact.matrix <- Ainv*(object$coefficients[i])
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







