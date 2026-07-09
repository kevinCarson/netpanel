
<!-- README.md is generated from README.Rmd. Please edit that file -->

# netpanel: An R Package that Contains Maximum Likelihood Estimation Routines for Panel Network Autocorrelation Models

<!-- badges: start -->

[![R-CMD-check](https://github.com/kevinCarson/netpanel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kevinCarson/netpanel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `netpanel` package contains maximum likelihood estimation routines
to fit panel network autocorrelation models (PNAMs) with time-varying
social adjacency matrices. The package contains functionality for (1)
fixed and random effects PNAMs, (2) dynamic PNAMs that include time lags
of the outcome and networks, and (3) the estimation of network mixed
effects models where, in addition to the lag of the outcome, the
time-varying residuals are assumed to exhibit network autocorrelation.
Importantly, the package allows for multiple network lags to be fitted
simultaneously. Please see the list of below references for more
information on panel network (spatial) autocorrelation models. For
generality, a panel network autocorrelation model is a linear regression
model that takes the following form:

$$Y_t = \lambda_1 W^{t}_{1}Y_t + \cdots + \lambda_k W^{t}_{k}Y_t + X_t\beta + v_t, \ \ t = 1,2,\cdots,T $$

where $Y_t$ is the outcome vector for $N$ actors at time $t$,
$\lambda_k W^{t}_{k}Y_t$ is the network autocorrelation term for the
$k^{th}$ social network at time $t$, $X_t$ are the set of exogenous
time-varying (and invariant) regressors. Moreover, the broad definition
of the time-varying error term $v_t$ is:

$$v_t = \mu_N + e_t,\ \ t = 1,2,\cdots,T $$

where $\mu_N$ are time-invariant unit-specific unobserved factors
thought to explain $Y_t$ and $e_t$ are the time- and unit-varying
unobserved factors.

## Authors

**Kevin A. Carson**  
- Author & Maintainer  
- PhD Candidate at the [University of Arizona School of
Sociology](https://sociology.arizona.edu/)  
- Email: <kacarson@arizona.edu>  
- Website: <https://kevincarson.github.io/>

## Installation

You can install the developmental version of `netpanel` from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("kevinCarson/netpanel")
```

## The `netpanel` Package API

The `netpanel` package API primarily works through the wrapper function
`mlpnam()`. The `mlpnam()` function returns an S3 object of class `pnam`
that stores the relevant results of the fitted PNAM. The full set of
functions (including user-helper S3 methods) are:

- `mlpnam()`: fits a panel network autocorrelation model to an empirical
  panel network dataset via Maximum Likelihood Estimation.
- `netimpacts()`: computes the direct, indirect, and total network
  impacts for each included regressor from a fitted `pnam` S3 object.
- `coef.pnam()`: extracts the vector of ML parameter estimates from a
  fitted `pnam` S3 object.
- `vcov.pnam()`: extracts the asymptotic variance-covariance matrix from
  a fitted `pnam` S3 object.
- `logLik.pnam()`: extracts the estimated log-likelihood from a fitted
  `pnam` S3 object.
- `residuals.pnam()`: extracts the model residuals from a fitted `pnam`
  S3 object.
- `fitted.pnam()`: extracts the linear and network regression fitted
  values from a fitted `pnam` S3 object.

## Fitting a Random Effects Panel Network Autocorrelation Model with `netpanel`

Importantly, the `netpanel` package includes a set of simulated pseudo
data objects for users to start fitting panel network autocorrelation
models. Below, we rely upon them to provide a brief example of how the
`netpanel` package fits these types of models.

``` r
library(netpanel)
data("net1", package = "netpanel") #the first panel network
data("net2", package = "netpanel") #the second panel network
data("net3", package = "netpanel") #the third panel network
data("simulated.data", package = "netpanel") #the simulated variables 

# a random effects panel network autocorrelation model
re.pnam <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1 + net2,
               data =  simulated.data, time = ~panel,
               actor = ~unit, model = "random")
```

``` r
summary(re.pnam) 
#> Maximum Likelihood Panel Network Autocorrelation Model with Unit-Level Random Effects
#> 
#> Call:
#> mlpnam(reg.formula = Y ~ x1 + x2 + x3, net.formula = ~net1 + 
#>     net2, data = simulated.data, time = ~panel, actor = ~unit, 
#>     model = "random")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -4.48611 -0.93835 -0.08678  0.89373  4.23752 
#> 
#> Panel error variance components:
#>                     sd     var
#> idiosyncratic: 0.90180 0.81324
#> unit-specific: 1.06094 1.12560
#> theta: 1.3841 (SE: 0.31323; p= 1e-05)
#> icc: 0.58055
#> 
#> Network autocorrelation parameters:
#>      Estimate Std. Error z value  Pr(>|z|)    
#> net1  0.26437    0.05388  4.9067 < 2.2e-16 ***
#> net2  0.29327    0.06024  4.8684 < 2.2e-16 ***
#> 
#> ML regression coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)  9.93222    0.16010  62.037 < 2.2e-16 ***
#> x1           0.85819    0.04277  20.068 < 2.2e-16 ***
#> x2          -1.13678    0.08582 -13.246 < 2.2e-16 ***
#> x3           1.99977    0.04052  49.356 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Model fit information:
#>  -> number of observations: 500; N: 50; T: 10
#>  -> parameters estimated: 8 
#>  -> log-likehoood: -901.1132 (df=492); AIC: 1818.226; BIC: 1851.943
#>  -> residual SD: 0.9018
#>  -> optim convergence: yes
#>  -> search iterations: 15
#>  -> largest eigenvalue of rho*W: 0.5576419
```

``` r
impacts <- netimpacts(re.pnam) #extracting the network impacts
```

``` r
print(impacts) 
#>   Regressor   Direct Indirect    Total
#> 1        x1  0.87035  1.06970  1.94004
#> 2        x2 -1.15288 -1.41694 -2.56982
#> 3        x3  2.02809  2.49262  4.52071
```

We can extend the above model by fitting a mixed network random effects
panel autocorrelation with `netpanel` that allows for the time-varying
residuals to be autocorrelated across the network ties.

``` r
# a random effects network mixed effects panel autocorrelation model
re.pnam.mixed <- mlpnam(Y~x1+x2+x3, net.formula = ~ net1 + net2,
               data =  simulated.data, time = ~panel,
               actor = ~unit, model = "random",
               errors = "autocorrelated",
               autocorrelated.network = net3)
```

``` r
summary(re.pnam.mixed) 
#> Maximum Likelihood Panel Network Autocorrelation Model with Unit-Level Random Effects
#> 
#> Call:
#> mlpnam(reg.formula = Y ~ x1 + x2 + x3, net.formula = ~net1 + 
#>     net2, data = simulated.data, time = ~panel, actor = ~unit, 
#>     model = "random", errors = "autocorrelated", autocorrelated.network = net3)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -4.48412 -0.91269 -0.07685  0.93722  4.18802 
#> 
#> Panel error variance components:
#>                     sd     var
#> idiosyncratic: 0.89532 0.80159
#> unit-specific: 1.05886 1.12118
#> theta: 1.3987 (SE: 0.31596; p= 1e-05)
#> icc: 0.58311
#> 
#> Network autocorrelation parameters:
#>      Estimate Std. Error z value Pr(>|z|)    
#> net1  0.25753    0.05519  4.6661  < 2e-16 ***
#> net2  0.25198    0.06829  3.6899  0.00022 ***
#> 
#> Network error parameters:
#>        Estimate Std. Error z value Pr(>|z|)   
#> lambda  0.35450    0.13641  2.5989  0.00935 **
#> 
#> ML regression coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 10.96107    0.16641  65.869 < 2.2e-16 ***
#> x1           0.86139    0.04237  20.331 < 2.2e-16 ***
#> x2          -1.14165    0.08430 -13.543 < 2.2e-16 ***
#> x3           1.99746    0.04015  49.756 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Model fit information:
#>  -> number of observations: 500; N: 50; T: 10
#>  -> parameters estimated: 9 
#>  -> log-likehoood: -898.1102 (df=491); AIC: 1814.22; BIC: 1852.152
#>  -> residual SD: 0.89532
#>  -> optim convergence: yes
#>  -> search iterations: 17
#>  -> largest eigenvalue of rho*W: 0.5095095
```

## Questions, Comments, or Suggestions!

If you have any questions, comments, or suggestions please feel free to
open an issue!
