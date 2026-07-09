#' Maximum Likelihood Estimation Routines for Panel Network Autocorrelation Models
#'
#' @description
#' The `netpanel` package contains maximum likelihood estimation routines to fit panel 
#' network autocorrelation models (PNAMs) with time-varying social adjacency matrices. The 
#' package contains functionality for (1) fixed and random effects PNAMs, (2) dynamic PNAMs that 
#' include time lags of the outcome and networks, and (3) the estimation of network mixed 
#' effects models where, in addition to the lag of the outcome, the time-varying residuals are 
#' assumed to exhibit network autocorrelation. Importantly, the package allows for multiple 
#' network lags to be fitted simultaneously. Please see the list of below references for more information
#' on panel network (spatial) autocorrelation models.  
#'
#' @details
#' The `netpanel` package primarily works through the wrapper function \code{\link{mlpnam}}. The \code{\link{mlpnam}} 
#' function returns an S3 object of class `pnam` that stores the relevant results of the fitted PNAM. The full set of functions (including user-helper S3 methods) are:
#'  - \code{\link{mlpnam}}: fits a panel network autocorrelation model to an empirical panel network dataset via Maximum Likelihood Estimation.
#'  - \code{\link{netimpacts}}: computes the direct, indirect, and total network impacts for each included regressor from a fitted `pnam` S3 object.
#'  - \code{\link{coef.pnam}}: extracts the vector of ML parameter estimates from a fitted `pnam` S3 object.
#'  - \code{\link{vcov.pnam}}: extracts the asymptotic variance-covariance matrix from a fitted `pnam` S3 object.
#'  - \code{\link{logLik.pnam}}: extracts the estimated log-likelihood from a fitted `pnam` S3 object.
#'  - \code{\link{residuals.pnam}}: extracts the model residuals from a fitted `pnam` S3 object.
#'  - \code{\link{fitted.pnam}}: extracts the linear and network regression fitted values from a fitted `pnam` S3 object.
#'  
#' Additionally, the `netpanel` package includes a set of simulated pseudo data objects for users to start fitting panel network 
#' autocorrelation models. Please see the \code{\link{simulated.data}} `netpanel` help page for more information on the 
#' simulated dataset and respective networks. 
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
#' @author Kevin A. Carson <kacarson@arizona.edu>
#' @name netpanel
#' @useDynLib netpanel, .registration=TRUE
"_PACKAGE"
NULL
#> NULL