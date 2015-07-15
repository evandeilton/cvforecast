\name{datasample}
\alias{datasample}
\docType{data}
\title{
Data set for testing the functions in the package
}
\description{
A data frame with 1000 columns where the first is a date vector in the format "dd/mm/aaaa hh:mm:ss" ande the last are variables simulating from 2 to 499 CPU consumption in percentage and from 500 to 1000 simulating Physical memory consumption also in percent.
}
\usage{data("datasample")}

\examples{
data(datasample)
str(datasample)
plot(datasample[,2])
}

