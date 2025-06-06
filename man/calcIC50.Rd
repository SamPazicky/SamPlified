% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcIC50.R
\name{calcIC50}
\alias{calcIC50}
\title{calculateIC50}
\usage{
calcIC50(data = NULL)
}
\arguments{
\item{data}{A data frame with columns "x" and "y".}
}
\value{
A list containing:
\item{IC50}{Calculated IC50.}
\item{fit}{Fit object.}
\item{coefs}{Fitted coefficients.}
\item{plot}{Plot with the fitted data.}
}
\description{
Calculates IC50 from dose-response data
}
\examples{
calcIC50(data)
}
