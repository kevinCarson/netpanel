# the master file for pnam R package that contains the simulated data files


#' Simulated panel network dataset
#'
#' A simulated panel network dataset that contains a response variable (\eqn{Y}) and three
#' exogenous regressors (\eqn{x_1}, \eqn{x_2}, and \eqn{x_3}). The exogenous regresors
#' \eqn{x_1} and \eqn{x_3} vectors are drawn from standard normal distribution, and 
#' \eqn{x_2} is a vector of random bernoulli variables. The `panel`
#' variable represents the time of data collection and the `unit` variable indexes the 
#' unique actors across the panels. The outcome is generated from the following deterministic
#' equation:
#' \deqn{Y_t = \rho_1 W^1 Y_t + \rho_2 W^w Y_t + \beta_0 + \beta_1 x_t1 + \beta_2 x_t2 + \beta_3 x_t3 + v_t, \ \ t = 1,2,\cdots,10}
#' \deqn{v_t = \mu_t + e_t, \ \ t = 1,2,\cdots,10}
#' \deqn{e_t = \lambda W^3 e_t + \epsilon_t, \ \ t = 1,2,\cdots,10}
#' where \eqn{W^1} is the \code{\link{net1}} `netpanel` data object, \eqn{W^2} is the  \code{\link{net2}} `netpanel` data object, and
#' \eqn{W^3} is the  \code{\link{net3}} `netpanel` data object.
#' 
#' @format ## `simulated.data`
#' A simulated panel network dataset with size \eqn{500 \times 5} the following variables: \eqn{Y}, \eqn{x_1},
#' \eqn{x_2}, \eqn{x_3}, `panel`, and `unit.`
#' 
#' @usage data(simulated.data)

"simulated.data"



#' A `list` of simulated panel social networks (1)
#'
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#' 
#' @format ## `net1`
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#'
#' @usage data(net1)
#' 

"net1"


#' A `list` of simulated panel social networks (2)
#'
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#' 
#' @format ## `net2`
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#'
#' @usage data(net2)
#' 

"net2"


#' A `list` of simulated panel social networks (3)
#'
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#' 
#' @format ## `net3`
#' A `list` of length 10 where each element is an simulated row-normalized 50 x 50 social 
#' adjacency matrix, where a non-zero \eqn{W_{ijt}} entry denotes a directed relationship 
#' from \eqn{i} to \eqn{j} at time \eqn{t}.
#'
#' @usage data(net3)
#' 

"net3"


