\name{steps2_4}
\alias{steps2_4}
\docType{data}
\title{
steps2_4
}
\description{
example results of the Fluo_CV_modeling function
}
\usage{data("steps2_4")}
\format{
  The format is:
List of 2
 $ Batch5   :List of 39
  ..$ index                        : int [1:246] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ samples                      : chr [1:246] "1772-062-248_A01" "1772-062-248_A03" "1772-062-248_A04" "1772-062-248_A06" ...
  ..$ batch                        : chr [1:246, 1:2] "1772-062-248" "1772-062-248" "1772-062-248" "1772-062-248" ...
  ..$ Size                         : num [1:246] 31 19 152 43 59 21 72 81 31 56 ...
  ..$ corrected.exprs              : num [1:246, 1:2] 26.64 9.5 24.84 15.43 8.04 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "Green" "Red"
  ..$ corrected.transformed.exprs  : num [1:246, 1:2] 3.28 2.25 3.21 2.74 2.08 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "Green" "Red"
  ..$ correctedAreas               : num [1:246] 3.43 2.94 5.02 3.76 4.08 ...
  ..$ areacut                      : num 49
  ..$ transformation               : chr "log"
  ..$ image.type                   : chr [1:3] "BF" "Green" "Red"
  ..$ dateIndex                    : chr "ThuMar2416:02:312016"
  ..$ single.batch.analysis        : num 5
  ..$ BGmethod                     : chr "normexp"
  ..$ maxMix                       : num 3
  ..$ prior.pi                     : num 0.1
  ..$ flex.reps                    : num 5
  ..$ flexmethod                   : chr "BIC"
  ..$ RNG                          : num 999
  ..$ GAPgroups                    : num [1:246, 1:2] 6 4 5 4 4 4 1 4 2 4 ...
  ..$ clusterFUN                   : chr "kmeans"
  ..$ normal.sigma                 : num 200
  ..$ centroids                    : num [1:6, 1:3] 3 6 5 4 2 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:3] "Cluster" "Green" "Red"
  ..$ fixClusters                  : num 0
  ..$ Kmax                         : num 15
  ..$ B.kmeans                     : num 5
  ..$ Path                         : num [1:6] 3 4 6 2 1 5
  ..$ Path.type                    : chr [1:2] "circular" "clockwise"
  ..$ UpdatedPath                  : num [1:6] 1 2 3 4 5 6
  ..$ DataSorts                    : chr [1:246, 1:2] "-0.0524253515714255" "0.0573114934260645" "-0.0658947839764907" "-0.10952104926524" ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "Distance" "Pseudotime"
  ..$ DDHFupdate                   : logi FALSE
  ..$ corrected.VStransformed.exprs: num [1:246, 1:2] 90.2 86.4 88.6 88.1 86.4 ...
  ..$ VSmethod                     : chr "DDHFmv"
  ..$ Progression                  : num [1:246, 1:2] 110 58 229 64 68 51 215 76 166 85 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "Pseudotime" "transf.Difference"
  ..$ Updated.groups               : num [1:246] 3 3 6 3 3 3 5 3 4 3 ...
  ..$ CPs                          : num [1:5] 22 38 118 169 222
  ..$ CPmethod                     : chr "ECP"
  ..$ CPsig                        : num 0.01
  ..$ CPgroups                     : num 5
  ..$ CPmingroup                   : num 15
 $ init.path: chr [1:2] "bottom/left" "bottom/left"
}
\value{
cross validation modeling
}
\keyword{datasets}
