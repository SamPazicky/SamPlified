% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateIDC.R
\name{estimateIDC}
\alias{estimateIDC}
\title{estimateIDC}
\usage{
estimateIDC(
  data = NULL,
  reference = samplified::Pf_ref_IDC,
  data.hpi = 42,
  reference.hpi = 48,
  scale = TRUE,
  log2transform = FALSE,
  cor.method = "spearman"
)
}
\arguments{
\item{data}{A data frame with columns "Accession","Value" and "Sample".}

\item{reference}{The reference timepoint data frame.}

\item{data.hpi}{How long is the cycle of the Plasmodium strain in samples?}

\item{reference.hpi}{How long is the cycle of the Plasmodium strain in the reference?}

\item{scale}{Logical: Should the data be scaled against reference.hpi?}

\item{log2transform}{Logical: Should the values be log2 transformed?}

\item{cor.method}{Corelation method as in argument 'method' of the function cor().}
}
\value{
A list containing:
\item{data}{Calculated data.}
\item{plot}{Correlation plot.}
}
\description{
Estimate IDC of a sample based on a reference
}
\examples{
estimateIDC(data)
}
