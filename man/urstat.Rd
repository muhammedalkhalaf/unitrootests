\name{urstat}
\alias{urstat}
\title{Comprehensive Unified Unit Root and Stationarity Testing}
\usage{
urstat(
  x,
  tests = "ALL",
  max_lag = 12L,
  crit = "BIC",
  pp_lag = 4L,
  kpss_lags = 8L,
  za_trim = 0.15,
  level = 0.05,
  stars = TRUE,
  strategy = TRUE
)
}
\arguments{
\item{x}{A numeric vector, \code{ts} object, or a named list/data frame of
numeric vectors.}

\item{tests}{Character vector of tests to run.  Possible values:
\code{"ADF"}, \code{"PP"}, \code{"KPSS"}, \code{"ERS"}, \code{"ZA"},
\code{"KM"}.  Default \code{"ALL"} runs all supported tests.}

\item{max_lag}{Integer. Maximum lag order for ADF lag selection (default 12).}

\item{crit}{Character. Information criterion for ADF lag selection:
\code{"BIC"} (default) or \code{"AIC"}.}

\item{pp_lag}{Integer. Bandwidth (Newey-West lags) for PP test (default 4).}

\item{kpss_lags}{Integer. Lag truncation for KPSS (default 8).}

\item{za_trim}{Numeric. Trimming proportion for Zivot-Andrews test
(default 0.15).}

\item{level}{Numeric. Significance level used for decisions (default 0.05).}

\item{stars}{Logical. Print significance stars (default \code{TRUE}).}

\item{strategy}{Logical. Print Elder-Kennedy decision table
(default \code{TRUE}).}
}
\value{
A list (invisibly) containing one element per series.  Each element
is itself a named list with components:
\describe{
  \item{\code{adf}}{Results from ADF tests (list).}
  \item{\code{pp}}{Results from PP tests (list).}
  \item{\code{kpss}}{Results from KPSS tests (list).}
  \item{\code{ers}}{Results from ERS/DF-GLS tests (list).}
  \item{\code{za}}{Results from Zivot-Andrews test (list).}
  \item{\code{decision}}{Character. Inferred integration order.}
  \item{\code{process}}{Character. Suggested data transformation.}
}
}
\description{
Runs a comprehensive battery of unit root and stationarity tests on one or
more time series, producing formatted summary tables and an integration-order
decision following the Elder and Kennedy (2001) strategy.
}
\examples{
set.seed(42)
x <- cumsum(rnorm(60))
res <- urstat(x, tests = c("ADF", "PP", "KPSS"), strategy = FALSE)
}
\references{
Elder, J. and Kennedy, P. E. (2001). Testing for unit roots: What should
students be taught? \emph{Journal of Economic Education}, 32(2), 137-146.
\doi{10.1080/00220480109595179}

Dickey, D. A. and Fuller, W. A. (1979). Distribution of the estimators for
autoregressive time series with a unit root.
\emph{Journal of the American Statistical Association}, 74(366), 427-431.
\doi{10.2307/2286348}

Phillips, P. C. B. and Perron, P. (1988). Testing for a unit root in time
series regression. \emph{Biometrika}, 75(2), 335-346.
\doi{10.1093/biomet/75.2.335}

Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., and Shin, Y. (1992).
Testing the null hypothesis of stationarity against the alternative of a
unit root. \emph{Journal of Econometrics}, 54(1-3), 159-178.
\doi{10.1016/0304-4076(92)90104-Y}

Elliott, G., Rothenberg, T. J., and Stock, J. H. (1996). Efficient tests
for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
\doi{10.2307/2171846}

Zivot, E. and Andrews, D. W. K. (1992). Further evidence on the great
crash, the oil-price shock, and the unit-root hypothesis.
\emph{Journal of Business & Economic Statistics}, 10(3), 251-270.
\doi{10.1080/07350015.1992.10509904}
}
