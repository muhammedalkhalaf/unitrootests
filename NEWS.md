# unitrootests 1.1.2

* Bug fix in `qadf()` (the same corrections as in the qadf package 1.0.2): the quantile autoregression is now estimated in levels, so `rho_tau` is rho rather than rho - 1; the statistic follows equation (9) of Koenker and Xiao (2004); the critical values of Hansen (1995) are interpolated in the estimated delta^2 instead of being indexed by tau; `delta2` is the squared correlation between the differenced series and psi_tau of the quantile residuals; and lag selection uses a common sample. For `model = "c"` the result agrees with the Stata command qadf (SSC).

# unitrootests 1.1.1

* Corrected the DOI of Hansen (1995) to 10.1017/S0266466600009993 in R and Rd files. No changes to code.

