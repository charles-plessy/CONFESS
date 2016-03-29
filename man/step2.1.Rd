\name{step2.1}
\alias{step2.1}
\docType{data}
\title{
step2.1
}
\description{
Example output of the getFluo function
}
\usage{data("step2.1")}
\format{
  The format is:
List of 18
 $ index                      : int [1:246] 1 2 3 4 5 6 7 8 9 10 ...
 $ samples                    : chr [1:246] "1772-062-248_A01" "1772-062-248_A03" "1772-062-248_A04" "1772-062-248_A06" ...
 $ batch                      : chr [1:246, 1:2] "1772-062-248" "1772-062-248" "1772-062-248" "1772-062-248" ...
 $ Size                       : num [1:246] 31 19 152 43 59 21 72 81 31 56 ...
 $ corrected.exprs            : num [1:246, 1:2] 34.3 11.96 31.95 19.7 9.95 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Green" "Red"
 $ corrected.transformed.exprs: num [1:246, 1:2] 3.54 2.48 3.46 2.98 2.3 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Green" "Red"
 $ correctedAreas             : num [1:246] 3.43 2.94 5.02 3.76 4.08 ...
 $ areacut                    : num 49
 $ transformation             : chr "log"
 $ image.type                 : chr [1:3] "BF" "Green" "Red"
 $ dateIndex                  : chr "WedMar2313:29:522016"
 $ single.batch.analysis      : num 5
 $ BGmethod                   : chr "normexp"
 $ maxMix                     : num 3
 $ prior.pi                   : num 0.1
 $ flex.reps                  : num 5
 $ flexmethod                 : chr "BIC"
 $ RNG                        : NULL
}
\value{
example intermediates
}
\keyword{datasets}
