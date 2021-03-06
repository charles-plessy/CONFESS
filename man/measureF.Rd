% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_fluo_est.R
\name{measureF}
\alias{measureF}
\title{measureF}
\usage{
measureF(img, area, BFarea)
}
\arguments{
\item{img}{Data matrix. The matrix of image data from a channel.}

\item{area}{Data matrix. The bright (spot) coordinates around the spot center.}

\item{BFarea}{Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.}
}
\value{
The foreground (spot) signal estimates
}
\description{
It estimates the foreground signal for an identified spot or for a predefined area within
  the capture site.
}
\keyword{internal}

