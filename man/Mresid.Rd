% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvForecastExtraFunctions.R
\name{Mresid}
\alias{Mresid}
\title{Time series forecast residual analysis}
\usage{
Mresid(fit)
}
\arguments{
\item{forecast}{object of class forecast}
}
\description{
Compute six residual test on a forecast object an return its p.values. The null hypothesis that the residual don't have the tested characteristic. If p.value is > 5\% we reject null.
}
\examples{
fcnm <- forecastMethod(lynx)
fun <- get(fcnm[[1]])
fit <- fun(AirPassengers, h=10, onlyfc = FALSE)
Mresid(fit)
}
\seealso{
See \code{\link{Box.test}}, \code{\link{t.test}}, \code{\link{LB.test}}, \code{\link{jarque.bera.test}}, \code{\link{bptest}}   and \code{\link{dwtest}} for more information about the tests.
}

