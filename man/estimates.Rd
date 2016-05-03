\name{estimates}
\alias{estimates}
\docType{data}
\title{
estimates
}
\description{
Example output of the SpotEstimator function
}
\usage{data("estimates")}
\format{
  The format is:
List of 6
 $ SpotResults      :'data.frame':	14 obs. of  14 variables:
  ..$ SampleID       : chr [1:14] "1772-062-248_A01" "1772-062-248_A02" "1772-062-248_A03" "1772-062-248_A04" ...
  ..$ X              : num [1:14] 259 491 262 261 261 258 259 189 498 194 ...
  ..$ Y              : num [1:14] 367 219 368 369 335 367 336 278 20 284 ...
  ..$ Size           : num [1:14] 31 49 19 152 141 43 59 15 49 32 ...
  ..$ Estimation.Type: chr [1:14] "Both.Channels" "Both.Channels" "One.Channel" "One.Channel" ...
  ..$ fore_Green     : num [1:14] 48.4 36 26.2 45.7 32.6 ...
  ..$ back_Green     : num [1:14] 17.2 17.3 16.6 16.9 17.1 ...
  ..$ fore_Red       : num [1:14] 219.1 27.6 86.5 18.4 48 ...
  ..$ back_Red       : num [1:14] 17.5 18.6 17.5 18.1 18 ...
  ..$ Green.StN      : num [1:14] 1.442 1.01 0.626 1.389 0.889 ...
  ..$ Green.Pvalue   : num [1:14] 6.03e-07 1.08e-03 5.55e-02 5.16e-27 4.57e-23 ...
  ..$ Red.StN        : num [1:14] 3.5689 0.5455 2.2422 0.0256 1.3664 ...
  ..$ Red.Pvalue     : num [1:14] 6.16e-07 2.68e-01 7.13e-05 1.00 3.33e-25 ...
  ..$ Other.Spots    : chr [1:14] "0" "0" "X = 30, Y = 204 (Green) | X = 262, Y = 368 (Red)" "0" ...
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
  ..$ BF        : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/BF/1772-062-248_A01_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/BF/1772-062-248_A02_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/BF/1772-062-248_A03_BF.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/BF/1772-062-248_A04_BF.txt" ...
  ..$ CH1       : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A01_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A02_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A03_Green.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A04_Green.txt" ...
  ..$ CH2       : chr [1:14] "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A01_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A02_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A03_Red.txt" "/home/diana/R/x86_64-pc-linux-gnu-library/3.2/CONFESSdata/extdata/CH/1772-062-248_A04_Red.txt" ...
  ..$ separator : chr "_"
  ..$ image.type: chr [1:3] "BF" "Green" "Red"
  ..$ dateIndex : chr "WedApr611:21:282016"
 $ BFarea           : num 7
 $ image.type       : chr [1:3] "BF" "Green" "Red"
 $ dateIndex        : chr "WedApr611:21:282016"
}
\value{
example intermediates
}
\keyword{datasets}
