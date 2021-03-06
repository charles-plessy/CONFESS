% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_fluo_est.R
\name{spotCoords}
\alias{spotCoords}
\title{spotCoords}
\usage{
spotCoords(centerR, centerG, origImg, chaImgs, minDiff, despeckle, ImgLimits,
  BFarea, chip.type, separator, image.type, show.possible.contamination)
}
\arguments{
\item{centerR}{Data matrix. The location statistics of one channel.}

\item{centerG}{Data matrix. The location statistics of the other channel.}

\item{origImg}{Data matrix. The original BF image to be read and processed.}

\item{chaImgs}{List. A list of the red and green channel images (data matrices) of a sample}

\item{minDiff}{Float. the mu_hat of the H0: image-to-noise ratio =
log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
implies that the identified spot is brighter than background. Default is 0.5.}

\item{despeckle}{Logical. If TRUE, the BF image is descpeckled.}

\item{ImgLimits}{Integer. It instructs the algorithm to find spots in a certain central image area.
For example, for a 512 x 512 image with ImgLimits = 50, it will search for spots in
the central area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.}

\item{BFarea}{Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.}

\item{separator}{Character string. Removes the Bright Field ("BF") and channel indicators from the
image file names.}

\item{image.type}{Character string. A triplet of IDs to characterize the type of images under study.
They refer to the ImageType part of the original image or txt file names.}

\item{show.possible.contamination}{Logical. If TRUE it reports all identified unmatched spots in both channels.}

\item{type}{Character string. It specifies the type of Fluidigm chip to be analyzed.}
}
\value{
Location statistics under fluorescence-based estimation
}
\description{
It estimates the spot location statistics by fluorescence signal in each channel. Then, it integrates
  the channel-specific data into a single estimate
}
\keyword{internal}

