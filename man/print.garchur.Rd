\name{print.garchur}
\alias{print.garchur}
\title{Print Method for garchur Objects}
\usage{
\method{print}{garchur}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"garchur"}.}

\item{...}{Further arguments passed to or from other methods (unused).}
}
\value{
Invisibly returns \code{x}.
}
\description{
Prints a formatted summary of GARCH unit root test results.
}
\examples{
set.seed(1)
y <- cumsum(rnorm(60))
res <- garchur(y, breaks = 2)
print(res)
}
