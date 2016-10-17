\name{step4}
\alias{step4}
\docType{data}
\title{
step4
}
\description{
Example output of the Fluo_modeling function
}
\usage{data("step4")}
\format{
  The format is:
List of 39
 $ index                        : int [1:246] 1 2 3 4 5 6 7 8 9 10 ...
 $ samples                      : chr [1:246] "1772-062-248_A01" "1772-062-248_A03" "1772-062-248_A04" "1772-062-248_A06" ...
 $ batch                        : chr [1:246, 1:2] "1772-062-248" "1772-062-248" "1772-062-248" "1772-062-248" ...
 $ Size                         : num [1:246] 31 19 152 43 59 21 72 81 31 56 ...
 $ corrected.exprs              : num [1:246, 1:2] 34.3 11.96 31.95 19.7 9.95 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Green" "Red"
 $ corrected.transformed.exprs  : num [1:246, 1:2] 3.54 2.48 3.46 2.98 2.3 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Green" "Red"
 $ correctedAreas               : num [1:246] 3.43 2.94 5.02 3.76 4.08 ...
 $ areacut                      : num 49
 $ transformation               : chr "log"
 $ image.type                   : chr [1:3] "BF" "Green" "Red"
 $ dateIndex                    : chr "WedMar2313:29:522016"
 $ single.batch.analysis        : num 5
 $ BGmethod                     : chr "normexp"
 $ maxMix                       : num 3
 $ prior.pi                     : num 0.1
 $ flex.reps                    : num 5
 $ flexmethod                   : chr "BIC"
 $ RNG                          : NULL
 $ GAPgroups                    : num [1:246, 1:2] 3 6 2 6 6 6 4 6 5 4 ...
 $ clusterFUN                   : chr "kmeans"
 $ normal.sigma                 : num 200
 $ centroids                    : num [1:7, 1:3] 2 1 5 4 3 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:3] "Cluster" "Green" "Red"
 $ fixClusters                  : num 0
 $ Kmax                         : num 15
 $ B.kmeans                     : num 5
 $ Path                         : num [1:6] 6 3 5 4 2 1
 $ Path.type                    : chr [1:2] "circular" "clockwise"
 $ UpdatedPath                  : num [1:7] 1 2 3 4 5 6 7
 $ DataSorts                    : chr [1:246, 1:2] "0.453392432450554" "0.0087217837649943" "0.201631879709111" "-0.395060882232867" ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Distance" "Pseudotime"
 $ DDHFupdate                   : logi FALSE
 $ corrected.VStransformed.exprs: num [1:246, 1:2] 93.7 87.7 90.3 89.3 87.7 ...
 $ VSmethod                     : chr "DDHFmv"
 $ Progression                  : num [1:246, 1:2] 77 18 201 22 31 11 169 34 119 171 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Pseudotime" "transf.Difference"
 $ Updated.groups               : num [1:246] 2 1 5 1 1 1 4 1 3 4 ...
 $ CPs                          : num [1:6] 56 93 140 178 206 229
 $ CPmethod                     : chr "ECP"
 $ CPsig                        : num 0.01
 $ CPgroups                     : num 5
 $ CPmingroup                   : num 10
}
\value{
example intermediates
}
\keyword{datasets}
