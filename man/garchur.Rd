\name{garchur}
\alias{garchur}
\title{GARCH Unit Root Test with Endogenous Structural Breaks}
\usage{
garchur(y, breaks = 2L, model = c("ct", "c"), trim = 0.15)
}
\arguments{
\item{y}{Numeric vector. The time-series variable to be tested.}

\item{breaks}{Integer. Number of structural breaks (1, 2, or 3).
  Default is \code{2}.}

\item{model}{Character. Deterministic specification: \code{"ct"}
  (constant plus linear trend, default) or \code{"c"} (constant only).}

\item{trim}{Numeric. Trimming proportion for the break-date grid search,
  between 0.05 and 0.30. Default is \code{0.15}.}
}
\value{
A list of class \code{"garchur"} containing: \code{stat} (t-statistic),
\code{rho}, \code{kappa}, \code{alpha}, \code{beta}, \code{ab},
\code{halflife}, \code{loglik}, \code{break_dates}, \code{nobs},
\code{cv1}, \code{cv5}, \code{cv10}, \code{model}, \code{breaks},
and \code{decision}.
}
\description{
Implements the trend-GARCH(1,1) unit root test with up to three endogenous
structural breaks proposed by Narayan and Liu (2015).
}
\details{
The mean equation is:
\deqn{y_t = a_0 + a_1 t + \rho y_{t-1} + \sum_j \gamma_j DU_{jt} + \varepsilon_t}
The null hypothesis is \eqn{H_0: \rho = 1} (unit root). Break dates are
selected by sequential maximum \eqn{|t|} search. GARCH parameters are
estimated by approximate maximum likelihood via Nelder-Mead optimisation.
Critical values are interpolated from Table III of Narayan and Liu (2015).
}
\examples{
set.seed(42)
y <- cumsum(rnorm(80))
res <- garchur(y, breaks = 2, model = "ct")
print(res)
}
\references{
Narayan, P. K., & Liu, R. (2015). A unit root model for trending
time-series energy variables. \emph{Energy Economics}, 46, 1--9.
\doi{10.1016/j.eneco.2014.11.021}
}
