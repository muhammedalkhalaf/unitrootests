---
title: 'unitrootests: Comprehensive Unit Root and Stationarity Tests'
tags:
  - R
  - econometrics
  - unit root
  - stationarity
  - time series
  - structural breaks
  - GARCH
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 25 March 2026
bibliography: paper.bib
---

# Summary

`unitrootests` provides a unified framework for unit root and stationarity testing in R, integrating a wide range of classical and modern tests under a single, consistent interface. The package includes: quantile ADF tests following @Koenker2004 that allow the order of integration to differ across the conditional distribution; GARCH-based unit root tests with endogenous structural breaks following @Narayan2015 that account for volatility clustering common in financial and energy price series; and a comprehensive battery of standard tests (Dickey-Fuller, Phillips-Perron, KPSS, ERS/DF-GLS, Zivot-Andrews, and Kobayashi-McAleer) guided by the Elder-Kennedy decision strategy @Elder2001. The package wraps and extends `tseries`, `urca`, and `strucchange` with higher-level wrappers that automate model specification and summarize results.

# Statement of Need

Determining the order of integration is a fundamental step in time series econometrics. While R packages such as `tseries` and `urca` provide individual tests, practitioners must manually coordinate test selection, handle structural break detection, and apply decision rules across multiple test results. No single R package offers GARCH-based unit root tests with endogenous breaks, quantile-distributional unit root tests, and classical tests within a unified workflow. Researchers in energy economics, finance, and macroeconomics frequently require all three approaches: classical tests for baseline inference, GARCH-based tests for volatile series, and quantile tests to examine tail behavior. `unitrootests` provides this integrated toolkit and implements the systematic Elder-Kennedy decision strategy to guide model selection across test specifications.

# Usage

## Classical Unit Root Battery with Elder-Kennedy Strategy

```r
library(unitrootests)

# Automated battery: ADF, PP, KPSS, DF-GLS, Zivot-Andrews
result_ur <- urstat(y, strategy = "elder_kennedy")
summary(result_ur)
print(result_ur)
```

## Quantile ADF Test

```r
# Quantile ADF: test stationarity across quantiles (Koenker & Xiao, 2004)
result_qadf <- qadf(y, taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
                    lags = 2)
print(result_qadf)
plot(result_qadf)
```

## GARCH Unit Root Test with Structural Breaks

```r
# GARCH-based unit root test with endogenous break (Narayan & Liu, 2015)
result_garch <- garchur(y, max_breaks = 1, model = "GARCH")
summary(result_garch)
```

# Implementation

`unitrootests` is built on top of the `quantreg`, `tseries`, `urca`, and `strucchange` packages. The quantile ADF test in `qadf()` follows @Koenker2004 by estimating quantile autoregressions and comparing the quantile-specific autoregressive coefficient against the null of a unit root, using critical values from the asymptotic distribution. The `garchur()` function implements @Narayan2015 by fitting a GARCH(1,1) model with structural break dummies, computing GLS-detrended data, and applying ADF-type tests on the residuals. The `urstat()` function applies the full battery of classical tests and implements the @Elder2001 decision tree to navigate model specification (trend, drift, or none) based on sequential test outcomes. All test results are returned as structured S3 objects with consistent `print()` and `summary()` methods.

# References
