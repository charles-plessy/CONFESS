% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fluo_NBE.R
\name{getFluo}
\alias{getFluo}
\title{getFluo}
\usage{
getFluo(data, areacut = 0)
}
\arguments{
\item{data}{List. The output of the Fluo_adjustment().}

\item{areacut}{Integer. The "artificial" area size (BFarea^2) of the cells estimated
by BF image modelling. Default is 0, implying that the area sizes to be corrected will
by estimated automatically from the data (not recommended if prior knowledge exists).}
}
\value{
A list of estimates to be used in subsequent analysis (the slots are the same to those of getFluo_byRun()):
  index: The sample indices.
  samples: the sample IDs.
  batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
      original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
  Size: the estimated cell size.
  corrected.exprs: the background corrected channel signals (case of a single run).
  corrected.transformed.exprs: the background transformed corrected channel signals (case of a single run). The
    transformation is defined in the transformation parameter.
  correctedAreas: the log-transformed areas after correction and imputation.
  areacut: the above areacut if different from 0 or the automatically calulated one otherwise.
  transformation: the transformation applied on the fluorescence signals.
  image.type: the image type IDs as defined in readFiles(). The parameter is kept in order to enable the user to
    use this function independently of the image analysis step.
  dateIndex: the date index used.
  single.batch.analysis: the reference run of the run effect correction by flexmix.
  BGmethod: the background correction methods used.
  maxMix: the maxMix parameter used.
  prior.pi: the prior.pi parameter used.
  flex.reps: the flex.reps parameter used.
  flexmethod: the flexmethod parameter used.
  RNG: the seed that is used to generate the results.
}
\description{
It retrieves the run effect and background corrected signals.
}
\examples{
step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
package = "CONFESS"),separator="_")
step2.1 <- getFluo(data=step2)
}

