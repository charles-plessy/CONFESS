% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fluo_est.R
\name{spotEstimator}
\alias{spotEstimator}
\title{spotEstimator}
\usage{
spotEstimator(files, correctionAlgorithm, subset = c(),
  foregroundCut = seq(0.5, 0.7, 0.02), denoise = FALSE, despeckle = FALSE,
  chip.type = "medium/large", cutSides = 0, BFarea = 7,
  log.transform = TRUE, minDiff = 0.5, show.possible.contamination = TRUE,
  cutoff = 50, QCdata = 0, median.correction = TRUE, savePlot = getwd())
}
\arguments{
\item{files}{Character string. The file names to be read and analyzed. This is the output of readFiles()}

\item{correctionAlgorithm}{Logical. Its value specifies the estimation stage. If FALSE,
the function processes all data using the standard operations of spotCoords(), i.e. case detection and fluorescence signal
estimation. This is the first estimation stage. If TRUE, the function processes the BF image modeling estimates of outlier images
obtained by defineLocClusters(). The BF image modeling is internally applied during the first stage. Note that
correctionAlgorithm = TRUE is strictly used in the second (outliers adjustment / correction) stage of the process.}

\item{subset}{Numeric vector. It can be a series sample index numbers (a subset) that specifies the samples to be analyzed.
The index numbers are obtained from readFiles() (the position of the sample in each listed vector). By default subset = c().
The parameter is mainly used in the second estimation stage where spotEstimator() processes the outlier images (the index numbers}

\item{foregroundCut}{Numeric vector. The binary segmentation image analysis cutoffs for normalized image data. Pixels with normalized signals
higher than the cutoff belong to foreground. Default is seq(0.5,0.7,0.02).}

\item{denoise}{Logical. If TRUE it denoises the channel images with la8, universal, hard. Default is FALSE.}

\item{despeckle}{Logical. If TRUE the bf image is descpeckled in the ImageJ fashion. Default is FALSE.}

\item{chip.type}{Character string. It specifies the type of Fluidigm chip to be analyzed. Default is "medium/large". The alternative
option is "small".}

\item{cutSides}{Integer. It instructs the algorithm to find spots in a certain central image area. For example, for a 512 x 512
image with cutSides = 50, spotEstimator() will search for spots in the central area [cutSides:(512-cutSides),cutSides:(512-cutSides)]
of the image matrix. Default is 0.}

\item{BFarea}{Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is mainly used in BF image
modeling where a fluorescence spot could not be originally detected. The value of this parameter is also used as a cut-off
to find matched spots across channel of the same sample image. Default is 7.}

\item{log.transform}{Logical. If TRUE the image data are plotted in the log scale. Default is TRUE}

\item{minDiff}{Float. The mu_hat of the H0: image-to-noise ratio =
log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
implies that the identified spot is brighter than background. Default is 0.5.}

\item{show.possible.contamination}{Logical. If TRUE it reports all identified unmatched spots in both channels. Default is TRUE.}

\item{cutoff}{Integer. A cutoff of the distance between the estimated spot location of an outlier sample (X, Y) and the median
location of all non-outliers of the same run and well set (medX,medY), i.e. (X-medX, Y-medY). An outlier sample can either
have a fluorescence-based location (X, Y) or a BF-based location (X*, Y*) or both. It is re-adjusted as follows: (1) if
min{(X-medX, Y-medY)} > cutoff and min{(X*-medX, Y*-medY)} > cutoff, the sample's location is set to (medX, medY); (2) if
min{(X*-medX, Y*-medY)} <= cutoff, the sample's location is set to (X*, Y*); (3) if min{(X-medX, Y-medY)} <= cutoff and
min{(X*-medX, Y*-medY)} > cutoff, the algorithm can either produce the solution of (1) or the solution of (2) depending
on the value of median.correction parameter below. By default cutoff = 50.}

\item{QCdata}{List. The output of defineLocClusters().}

\item{median.correction}{Logical. If TRUE, the algorithm re-adjusts the location of the outlier sample as the median of all
non-outliers of the same run and well ID (if necessary).}

\item{savePlot}{Character string. Directory to store the plots. Its value can be an existing directory
or "screen" that prints the plot only on the screen. Default is the current working directory, getwd().}
}
\value{
A list of the following components:
  SpotResults: the matrix of the location and fluorescence signal estimates. It contains the index number of each sample, the X,Y
    coordinates of the spot center, the spot size, the type of estimation that have been performed (fluorescence based indicating the channels
    in which the spot has been found or BF image modelling based), the fluorescence foreground and background signals of each channel,
    the signal-to-noise ratio (logForeground - logBackground) for each channel, the associated P-value of significance of the signal-to-noise
    ratio and a column indicating the coordinates of other spots that are not matched in both images. Existence of such spots (values that are
    different from 0) indicate contaminated image or highly noisy images or images with other artefacts. If correctionAlgorithm=TRUE (second
    spotEstimator() step), there is an extra column generated indicating outlier samples (see the QCgroup column in defineLocClusters()).
  Outlier.Estimates: The estimates obtained from BF modeling (if necessary to be obtained). These are alternative location estimates that will
    be used if the original estimates of the SpotResults table are flagged as outliers.
  Processed.Files: the samples that have been processed by spotEstimator().
  BFarea: the pseudospot size.
  image.type: the image type IDs.
  dateIndex: a date index to be used in saving the output files.
}
\description{
The main function to produce the raw fluorescence signal estimation results by analysis of the Fluidigm images.
}
\details{
Triplets of images of the same sample are sequentially considered to estimate the channel-specific
  fluorescence signals (if detectable) or perform BF image modeling. The main result of this function is a table
  of location and fluorescence estimates for each sample.
}
\examples{
### set your directories
basedir<-"~/"
#data_path<-system.file("extdata",package="CONFESSdata")
#files<-readFiles(iDirectory=NULL,
#                 BFdirectory=paste(data_path,"/BF",sep=""),
#                 CHdirectory=paste(data_path,"/CH",sep=""),
#                 separator = "_",image.type = c("BF","Green","Red"),
#                 bits=2^16)

### an example where the second image produces a clear outlier!
#estimates <- spotEstimator(files=files,subset=1:3,foregroundCut=seq(0.6,0.76,0.02),
#                           correctionAlgorithm=FALSE,savePlot="screen")
}

