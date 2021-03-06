% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cases.R
\name{caseof0s}
\alias{caseof0s}
\title{caseof0s}
\usage{
caseof0s(centerR, centerG, origImg, chaImgs, minDiff, despeckle, ImgLimits,
  BFarea, chip.type, separator, image.type)
}
\arguments{
\item{centerR}{Data matrix. The location statistics of one channel.}

\item{centerG}{Data matrix. The location statistics of the other channel.}

\item{origImg}{Data matrix. The original BF image to be read and processed.}

\item{chaImgs}{List. The channel image data (data matrices) of a sample.}

\item{minDiff}{Float. the mu_hat of the H0: image-to-noise ratio =
log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
implies that the identified spot is brighter than background. Default is 0.5.}

\item{despeckle}{Logical. If TRUE the BF image is despeckled.}

\item{ImgLimits}{Integer. It instructs the algorithm to find spots in the specified central image area.
For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
[cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.}

\item{BFarea}{Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
this parameter is also used as a cut-off to find matched spots across channel of the same sample image.}

\item{chip.type}{Character string. It specifies the type of Fluidigm chip to be analyzed.}

\item{separator}{Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
image file names.}

\item{image.type}{Character string. A triplet of IDs to characterize the type of images under study.}
}
\value{
A list of location estimates
}
\description{
It processes the case of 0 spots in both channels. It performs BF image modelling.
}
\keyword{internal}

