# netpanel 0.0.1 (2026-07-09)

* This is the first release of the `netpanel` R package. As a suggestion, please see the 
`netpanel` help package (i.e., `help("netpanel")`) to get started. Additionally, any 
coding issues, suggestions, or errors can be made by opening an issue here: https://github.com/kevinCarson/netpanel/issues. 
**Happy netpaneling!!**
* The `mlpnam()` function includes maximum likelihood estimation routines to fit 
panel network autocorrelation models with one or more time-varying social adjacency matrices. The following
types of models can be estimated: fixed effects, random effects, and dynamic models that includes outcome and network
time lags. The function returns a S3 object of class `pnam`. 
* The `netimpacts()` function computes the direct, indirect, and total impacts for each estimated effect 
included in the `pnam` S3 object. 
* The following methods are included for the S3 `pnam` object class: (1) `logLik()`, (2) `coef()`, (3) `vcov()`, (4) `residuals()`, and (5) `fitted()`.
* The following data objects are included: (1) `net1`, (2) `net2`, (3) `net3`, and (4) `simulated.data`. Please see 
their respective help pages for information on the data objects (e.g., `?netpanel::net1`).