\name{estimates.2}
\alias{estimates.2}
\docType{data}
\title{
estimates.2
}
\description{
Example output from the 2nd run of the spotEstimator function
}
\usage{data("estimates.2")}
\format{
  The format is:
List of 6
 $ SpotResults      :'data.frame':	14 obs. of  15 variables:
  ..$ SampleID       : chr [1:14] "1772-062-248_A01" "1772-062-248_A02" "1772-062-248_A03" "1772-062-248_A04" ...
  ..$ X              : num [1:14] 259 491 262 261 261 258 259 189 498 194 ...
  ..$ Y              : num [1:14] 367 219 368 369 335 367 336 278 20 284 ...
  ..$ Size           : num [1:14] 31 49 19 152 141 43 59 15 49 32 ...
  ..$ Estimation.Type: chr [1:14] "Fluorescence-based" "Fluorescence-based" "Fluorescence-based" "Fluorescence-based" ...
  ..$ fore_Green     : num [1:14] 48.4 36 26.2 45.7 32.6 ...
  ..$ back_Green     : num [1:14] 17.1 17.3 16.7 16.8 17.1 ...
  ..$ fore_Red       : num [1:14] 219.1 27.6 86.5 18.4 48 ...
  ..$ back_Red       : num [1:14] 17.4 18.8 17.5 18 18.1 ...
  ..$ Green.StN      : num [1:14] 1.451 1.017 0.621 1.391 0.89 ...
  ..$ Green.Pvalue   : num [1:14] 6.03e-07 8.99e-04 9.49e-02 5.16e-27 4.57e-23 ...
  ..$ Red.StN        : num [1:14] 3.5791 0.5319 2.2409 0.0304 1.3607 ...
  ..$ Red.Pvalue     : num [1:14] 6.16e-07 3.63e-01 7.13e-05 1.00 3.33e-25 ...
  ..$ Other.Spots    : chr [1:14] "0" "0" "X = 30, Y = 204 (Green) | X = 262, Y = 368 (Red)" "0" ...
  ..$ QCgroup        : chr [1:14] "confidence" "confidence" "confidence" "confidence" ...
 $ Outlier.Estimates:List of 14
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A01"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A02"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A03"
  .. ..$ centerR: num [1:2] 263 370
  .. ..$ centerG: num [1:2] 263 370
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : chr "BF"
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A04"
  .. ..$ centerR: num [1:2] 265 370
  .. ..$ centerG: num [1:2] 265 370
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : chr "BF"
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A05"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A06"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-062-248_A07"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A01"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A02"
  .. ..$ centerR: num [1:2] 195 250
  .. ..$ centerG: num [1:2] 195 250
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : chr "BF"
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A03"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A04"
  .. ..$ centerR: num [1:2] 191 281
  .. ..$ centerG: num [1:2] 191 281
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : chr "BF"
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A05"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A06"
  .. ..$ centerR: num [1:2] 187 274
  .. ..$ centerG: num [1:2] 187 274
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : chr "BF"
  ..$ :List of 6
  .. ..$ sample : chr "1772-067-039_A07"
  .. ..$ centerR: num [1:2] 0 0
  .. ..$ centerG: num [1:2] 0 0
  .. ..$ arR    : NULL
  .. ..$ arG    : NULL
  .. ..$ warn   : NULL
 $ Processed.Files  :List of 6
  ..$ BF        : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/BF/1772-062-248_A01_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/BF/1772-062-248_A02_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/BF/1772-062-248_A03_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/BF/1772-062-248_A04_BF.txt" ...
  ..$ CH1       : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A01_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A02_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A03_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A04_Green.txt" ...
  ..$ CH2       : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A01_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A02_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A03_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/fucci/extdata/CH/1772-062-248_A04_Red.txt" ...
  ..$ separator : chr "_"
  ..$ image.type: chr [1:3] "BF" "Green" "Red"
  ..$ dateIndex : chr "WedMar2313:16:562016"
 $ BFarea           : num 7
 $ image.type       : chr [1:3] "BF" "Green" "Red"
 $ dateIndex        : chr "WedMar2313:16:562016"
}
\value{
example intermediates
}
\keyword{datasets}
