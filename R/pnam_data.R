# the master file for pnam R package that contains the simulated data files


#' Simulated panel network datasets
#'
#' A simulated panel network dataset that contains the outcome (Y) and three exogenous
#' regressors (x1, x3, and x3), where x1 and x3 are drawn from the standard normal 
#' distribution, and x2 is a random bernoulli (binary) regressor. The `panel`
#' variable represents the time of data collection and `unit` indexes the 
#' actor. 
#' 
#' @format ## `simulated.data`
#' A simulated panel network dataset that contains the outcome (Y) and three exogenous
#' regressors (x1, x3, and x3), where x1 and x3 are drawn from the standard normal 
#' distribution, and x2 is a random bernoulli (binary) regressor.
#'
#' @usage data(simulated.data)

"simulated.data"



#' Simulated panel networks (1)
#'
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point. The networks are generated with the
#' \code{\link[sna]{rgraph}} function in the sna package (Butts 2024).
#' 
#' @format ## `network.1`
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point.
#'
#' @usage data(network.1)
#' @source  Butts CT (2024). _sna: Tools for Social Network Analysis_. doi:10.32614/CRAN.package.sna
#' <https://doi.org/10.32614/CRAN.package.sna>, R package version 2.8, <https://CRAN.R-project.org/package=sna>.
#' 
#' 

"net1"


#' Simulated panel networks (2)
#'
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point. The networks are generated with the
#' \code{\link[sna]{rgraph}} function in the sna package (Butts 2024).
#' 
#' @format ## `network.2`
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point.
#'
#' @usage data(network.2)
#' @source  Butts CT (2024). _sna: Tools for Social Network Analysis_. doi:10.32614/CRAN.package.sna
#' <https://doi.org/10.32614/CRAN.package.sna>, R package version 2.8, <https://CRAN.R-project.org/package=sna>.
#' 
#' 
"net2"


#' Simulated panel networks (3)
#'
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point. The networks are generated with the
#' \code{\link[sna]{rgraph}} function in the sna package (Butts 2024).
#' 
#' @format ## `error.network`
#' A `list` of length 10 that contains simulated row-normalized 50 x 50 social 
#' adjacency matrices for each time point.
#'
#' @usage data(error.network)
#' @source  Butts CT (2024). _sna: Tools for Social Network Analysis_. doi:10.32614/CRAN.package.sna
#' <https://doi.org/10.32614/CRAN.package.sna>, R package version 2.8, <https://CRAN.R-project.org/package=sna>.
#' 
#' 
"error.net"


