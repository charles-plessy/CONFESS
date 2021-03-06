% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_fluo_NBE.R
\name{CVsampler}
\alias{CVsampler}
\title{CVsampler}
\usage{
CVsampler(data, f)
}
\arguments{
\item{data}{List. The output of Fluo_inspection() or Fluo_modeling(). It requires existence of the @GAPgroups slot.}

\item{f}{Float. The percentage of samples that are used in the CV analysis (the rest is re-estimated).}
}
\value{
An index with the data that will remain in the analysis.
}
\description{
It samples a data subset for the cross-validation analysis.
}
\keyword{internal}

