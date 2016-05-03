# Cases for the estimation of fluorescence signal from the C1 Chip (internal functions)


#' caseof0s
#'
#' It processes the case of 0 spots in both channels. It performs BF image modelling.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof0s <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type) {
    res <-
      getCsFAIL(
        centerR = centerR,
        centerG = centerG,
        origImg = origImg,
        chaImgs = chaImgs,
        minDiff =
          minDiff,
        despeckle = despeckle,
        ImgLimits = ImgLimits,
        BFarea = BFarea,
        chip.type = chip.type,
        separator = separator,
        image.type = image.type
      )
    centerR[[1]] <- res$centers
    centerG[[1]] <- res$centers
    arR <- res$areaR
    arG <- res$areaG
    warn <- res$warning
    fai <- res$fail$failParameters
    out <- res$outlier.estimates
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = 0
      )
    )
  }


#' caseof1R0G
#'
#' It processes the case of 1 spot in one channel and 0 spots in the other channel.
#'   BF image modelling is not necessarily performed.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof1R0G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type) {
    out <- list(
      sample = "",
      centerR = c(0, 0),
      centerG = c(0, 0),
      arR = c(),
      arG = c(),
      warn = c()
    )
    
    if (centerR[[1]][1, 1] == 0) {
      centerR <- centerG
      arR <- centerG$brightcoordinates
      arG <- centerG$brightcoordinates
      warn <- giveWarning(2)
      fai <- as.list(rep(0, 7))
      names(fai) <- c(
        "reducedCols",
        "reducedRows",
        "strLines",
        "strLinesFar",
        "strLinesCutparam",
        "binImgDespeckled",
        "spotBetweenRows"
      )
    }
    if (centerG[[1]][1, 1] == 0) {
      centerG <- centerR
      arR <- centerR$brightcoordinates
      arG <- centerR$brightcoordinates
      warn <- giveWarning(2)
      fai <- as.list(rep(0, 7))
      names(fai) <- c(
        "reducedCols",
        "reducedRows",
        "strLines",
        "strLinesFar",
        "strLinesCutparam",
        "binImgDespeckled",
        "spotBetweenRows"
      )
    }
    
    res <-
      getCsFAIL(
        centerR = centerR,
        centerG = centerG,
        origImg = origImg,
        chaImgs = chaImgs,
        minDiff =
          minDiff,
        despeckle = despeckle,
        ImgLimits = ImgLimits,
        BFarea = BFarea,
        chip.type = chip.type,
        separator = separator,
        image.type = image.type
      )
    out <- res$outlier.estimates
    
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = 0
      )
    )
  }

#' caseof2Rs0G
#'
#' It processes the case of >1 spots in one channel and 0 spots in the other channel.
#'   It performs BF image modelling and reports possible contamination.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof2Rs0G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type,
           show.possible.contamination) {
    # report possible contamination
    w.out <- 0
    if (centerG[[1]][1, 1] == 0) {
      w.out <-
        paste(
          "X = ",
          centerR[[1]][, 1],
          ", Y = ",
          centerR[[1]][, 2],
          " (",
          image.type[2],
          ")",
          sep = "",
          collapse = " & "
        )
    }
    if (centerR[[1]][1, 1] == 0) {
      w.out <-
        paste(
          "X = ",
          centerG[[1]][, 1],
          ", Y = ",
          centerG[[1]][, 2],
          " (",
          image.type[3],
          ")",
          sep = "",
          collapse = " & "
        )
    }
    
    res <-
      getCsFAIL(
        centerR = centerR,
        centerG = centerG,
        origImg = origImg,
        chaImgs = chaImgs,
        minDiff = minDiff,
        despeckle = despeckle,
        ImgLimits = ImgLimits,
        BFarea = BFarea,
        chip.type = chip.type,
        separator = separator,
        image.type = image.type
      )
    centerR[[1]] <- res$centers
    centerG[[1]] <- res$centers
    arR <- res$areaR
    arG <- res$areaG
    warn <- res$warning
    fai <- res$fail$failParameters
    out <- res$outlier.estimates
    
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = w.out
      )
    )
  }


#' caseof1R1G
#'
#' It processes the case of 1 spot in both channels. BF image modelling is not necessarily performed.
#' It reports possible contamination.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof1R1G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type,
           show.possible.contamination) {
    ro <- abs(centerR[[1]][1, 1] - centerG[[1]][1, 1])
    co <- abs(centerR[[1]][1, 2] - centerG[[1]][1, 2])
    
    if (ro <= BFarea & co <= BFarea) {
      centerR[[1]] <-
        c(floor(mean(c(
          centerR[[1]][1, 1], centerG[[1]][1, 1]
        ))), floor(mean(c(
          centerR[[1]][1, 2], centerG[[1]][1, 2]
        ))))
      centerG[[1]] <- centerR[[1]]
      arR <- centerR$brightcoordinates
      arG <- centerG$brightcoordinates
      warn <- giveWarning(1)
      fai <- as.list(rep(0, 7))
      names(fai) <- c(
        "reducedCols",
        "reducedRows",
        "strLines",
        "strLinesFar",
        "strLinesCutparam",
        "binImgDespeckled",
        "spotBetweenRows"
      )
      out <- list(
        sample = "",
        centerR = c(0, 0),
        centerG = c(0, 0),
        arR = c(),
        arG = c(),
        warn = c()
      )
      
    }
    
    # report possible contamination
    w.out <- 0
    
    if (ro > BFarea | co > BFarea) {
      w.out <-
        paste(
          "X = ",
          centerR[[1]][1, 1],
          ", Y = ",
          centerR[[1]][1, 2],
          " (",
          image.type[2],
          ") | ",
          "X = ",
          centerG[[1]][1, 1],
          ", Y = ",
          centerG[[1]][1, 2],
          " (",
          image.type[3],
          ")",
          sep = ""
        )
      res <-
        getCsFAIL(
          centerR = centerR,
          centerG = centerG,
          origImg = origImg,
          chaImgs = chaImgs,
          minDiff =
            minDiff,
          despeckle = despeckle,
          ImgLimits = ImgLimits,
          BFarea = BFarea,
          chip.type = chip.type,
          separator = separator,
          image.type = image.type
        )
      fai <- res$fail$failParameters
      centerR[[1]] <- res$centers
      centerG[[1]] <- res$centers
      arR <- res$areaR
      arG <- res$areaG
      warn <- res$warning
      out <- res$outlier.estimates
    }
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = w.out
      )
    )
  }


#' caseof1R2G
#'
#' It processes the case of 1 spot in the red channel and >1 spots in the green channel.
#'   BF image modelling is not necessarily performed. It reports possible contamination.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof1R2G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type,
           show.possible.contamination) {
    res <- matrix(0, nrow(centerG[[1]]), 4)
    for (j in 1:nrow(centerG[[1]])) {
      r1 <-
        max(c(abs(centerR[[1]][1, 1] - centerG[[1]][j, 1]), abs(centerR[[1]][1, 2] -
                                                                  centerG[[1]][j, 2])))
      r2 <- nrow(centerG$brightcoordinates[[j]])
      res[j,] <- c(1, j, r1, r2)
    }
    w1 <- which(res[, 3] == min(res[, 3]) & res[, 3] <= BFarea)
    
    # report possible contamination
    w.out <- 0
    if (show.possible.contamination == TRUE) {
      w.out <- which(res[, 3] > BFarea & res[, 4] > BFarea)
      if (length(w.out) > 0) {
        w.out <-
          paste(
            "X = ",
            centerG[[1]][w.out, 1],
            ", Y = ",
            centerG[[1]][w.out, 2],
            " (",
            image.type[3],
            ")",
            sep = "",
            collapse = " & "
          )
      } else {
        w.out <- 0
      }
    }
    
    w <- w1
    if (length(w1) > 0) {
      if (length(w1) > 1) {
        w <- w1[which(res[w1, 4] == max(res[w1, 4]))]
      }
      centerR[[1]] <- centerR[[1]][1,]
      centerG[[1]] <- centerG[[1]][res[w[1], 2],]
      centerR[[1]] <-
        c(floor(mean(c(
          centerR[[1]][1], centerG[[1]][1]
        ))), floor(mean(c(
          centerR[[1]][2], centerG[[1]][2]
        ))))
      centerG[[1]] <- centerR[[1]]
      arR <- centerR$brightcoordinates
      arG <- centerG$brightcoordinates[[res[w[1], 2]]]
      warn <- giveWarning(1)
      fai <- as.list(rep(0, 7))
      names(fai) <-
        c(
          "reducedCols",
          "reducedRows",
          "strLines",
          "strLinesFar",
          "strLinesCutparam",
          "binImgDespeckled",
          "spotBetweenRows"
        )
      out <- list(
        centerR = c(0, 0),
        centerG = c(0, 0),
        arR = c(),
        arG = c(),
        warn = c()
      )
      
    } else if (length(w1) == 0) {
      res <-
        getCsFAIL(
          centerR = centerR,
          centerG = centerG,
          origImg = origImg,
          chaImgs = chaImgs,
          minDiff =
            minDiff,
          despeckle = despeckle,
          ImgLimits = ImgLimits,
          BFarea = BFarea,
          chip.type = chip.type,
          separator = separator,
          image.type = image.type
        )
      fai <- res$fail$failParameters
      centerR[[1]] <- res$centers
      centerG[[1]] <- res$centers
      arR <- res$areaR
      arG <- res$areaG
      warn <- res$warning
      out <- res$outlier.estimates
    }
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = w.out
      )
    )
  }


#' caseof2R1G
#'
#' It processes the case of >1 spots in the red channel and 1 spot in the green channel.
#'   BF image modelling is not necessarily performed. It reports possible contamination.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof2R1G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type,
           show.possible.contamination) {
    res <- matrix(0, nrow(centerR[[1]]), 4)
    for (i in 1:nrow(centerR[[1]])) {
      r1 <-
        max(c(abs(centerR[[1]][i, 1] - centerG[[1]][1, 1]), abs(centerR[[1]][i, 2] -
                                                                  centerG[[1]][1, 2])))
      r2 <- nrow(centerR$brightcoordinates[[i]])
      res[i,] <- c(i, 1, r1, r2)
    }
    w1 <- which(res[, 3] == min(res[, 3]) & res[, 3] <= BFarea)
    
    # report possible contamination
    w.out <- 0
    if (show.possible.contamination == TRUE) {
      w.out <- which(res[, 3] > BFarea & res[, 4] > BFarea)
      if (length(w.out) > 0) {
        w.out <-
          paste(
            "X = ",
            centerR[[1]][w.out, 1],
            ", Y = ",
            centerR[[1]][w.out, 2],
            " (",
            image.type[2],
            ")",
            sep = "",
            collapse = " & "
          )
      } else {
        w.out <- 0
      }
    }
    
    w <- w1
    if (length(w1) > 0) {
      if (length(w1) > 1) {
        w <- w1[which(res[w1, 4] == max(res[w1, 4]))]
      }
      centerG[[1]] <- centerG[[1]][1,]
      centerR[[1]] <- centerR[[1]][res[w[1], 1],]
      centerG[[1]] <-
        c(floor(mean(c(
          centerR[[1]][1], centerG[[1]][1]
        ))), floor(mean(c(
          centerR[[1]][2], centerG[[1]][2]
        ))))
      centerR[[1]] <- centerG[[1]]
      arR <- centerR$brightcoordinates[[res[w[1], 1]]]
      arG <- centerG$brightcoordinates
      warn <- giveWarning(1)
      fai <- as.list(rep(0, 7))
      names(fai) <-
        c(
          "reducedCols",
          "reducedRows",
          "strLines",
          "strLinesFar",
          "strLinesCutparam",
          "binImgDespeckled",
          "spotBetweenRows"
        )
      out <- list(
        sample = "",
        centerR = c(0, 0),
        centerG = c(0, 0),
        arR = c(),
        arG = c(),
        warn = c()
      )
      
    } else if (length(w1) == 0) {
      res <-
        getCsFAIL(
          centerR = centerR,
          centerG = centerG,
          origImg = origImg,
          chaImgs = chaImgs,
          minDiff =
            minDiff,
          despeckle = despeckle,
          ImgLimits = ImgLimits,
          BFarea = BFarea,
          chip.type = chip.type,
          separator = separator,
          image.type = image.type
        )
      fai <- res$fail$failParameters
      centerR[[1]] <- res$centers
      centerG[[1]] <- res$centers
      arR <- res$areaR
      arG <- res$areaG
      warn <- res$warning
      out <- res$outlier.estimates
    }
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = w.out
      )
    )
  }


#' caseof2R2G
#'
#' It processes the case of >1 spots in both channels. BF image modelling is not necessarily performed
#' It implies image contamination.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. The channel image data (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in the specified central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central area
#'   [cutSides:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators (IDs) from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'
#' @return A list of location estimates
#'
#' @keywords internal
caseof2R2G <-
  function(centerR,
           centerG,
           origImg,
           chaImgs,
           minDiff,
           despeckle,
           ImgLimits,
           BFarea,
           chip.type,
           separator,
           image.type) {
    res <- matrix(0, 1, 4)
    for (i in 1:nrow(centerR[[1]])) {
      for (j in 1:nrow(centerG[[1]])) {
        r1 <-
          max(c(abs(centerR[[1]][i, 1] - centerG[[1]][j, 1]), abs(centerR[[1]][i, 2] -
                                                                    centerG[[1]][j, 2])))
        r2 <-
          nrow(centerR$brightcoordinates[[i]]) + nrow(centerG$brightcoordinates[[j]])
        res <- matrix(rbind(res, c(i, j, r1, r2)), ncol = 4)
      }
    }
    res <- matrix(res[-1,], ncol = ncol(res))
    w1 <- which(res[, 3] <= BFarea)
    if (length(w1) > 0) {
      res <- matrix(res[w1,], ncol = ncol(res))
      su <-
        apply(matrix(cbind(
          as.numeric(duplicated(res[, 1])), as.numeric(duplicated(res[, 2]))
        ), nrow = nrow(res)), 1, sum)
      w2 <- which(su == 0)
      res <- matrix(res[w2,], ncol = ncol(res))
      if (length(w2) > 1) {
        sl <- sort.list(res[, 4], decreasing = TRUE)
        w2 <- w2[sl[1:2]]
        res <- matrix(res[sl[1:2],], ncol = ncol(res))
      }
      #    c1<-matrix(centerR[[1]][res[w2,1],],ncol=ncol(centerR[[1]]))
      #    c2<-matrix(centerG[[1]][res[w2,2],],ncol=ncol(centerG[[1]]))
      c1 <-
        matrix(centerR[[1]][res[, 1],], ncol = ncol(centerR[[1]]))
      c2 <-
        matrix(centerG[[1]][res[, 2],], ncol = ncol(centerG[[1]]))
      centerR[[1]] <- c1[1,]
      centerG[[1]] <- c2[1,]
      centerR[[1]] <-
        c(floor(mean(c(
          centerR[[1]][1], centerG[[1]][1]
        ))), floor(mean(c(
          centerR[[1]][2], centerG[[1]][2]
        ))))
      centerG[[1]] <- centerR[[1]]
      a1 <- centerR$brightcoordinates[res[, 1]]
      a2 <- centerG$brightcoordinates[res[, 2]]
      arR <- a1[[1]]
      arG <- a2[[1]]
      fai <- as.list(rep(0, 7))
      names(fai) <-
        c(
          "reducedCols",
          "reducedRows",
          "strLines",
          "strLinesFar",
          "strLinesCutparam",
          "binImgDespeckled",
          "spotBetweenRows"
        )
      
      if (length(w2) > 1) {
        c.alt <-
          c(floor(mean(c(c1[2, 1], c2[2, 1]))), floor(mean(c(c1[2, 2], c2[2, 2]))))
        out <- list(
          sample = "",
          centerR = c.alt,
          centerG = c.alt,
          arR = a1[[2]],
          arG = a2[[2]],
          warn = giveWarning(5)
        )
        warn <- giveWarning(5)
      } else {
        out <- list(
          sample = "",
          centerR = c(0, 0),
          centerG = c(0, 0),
          arR = c(),
          arG = c(),
          warn = c()
        )
        warn <- giveWarning(1)
      }
      
    } else if (length(w1) == 0) {
      res <-
        getCsFAIL(
          centerR = centerR,
          centerG = centerG,
          origImg = origImg,
          chaImgs = chaImgs,
          minDiff =
            minDiff,
          despeckle = despeckle,
          ImgLimits = ImgLimits,
          BFarea = BFarea,
          chip.type = chip.type,
          separator = separator,
          image.type = image.type
        )
      fai <- res$fail$failParameters
      centerR[[1]] <- res$centers
      centerG[[1]] <- res$centers
      arR <- res$areaR
      arG <- res$areaG
      warn <- res$warning
      out <- res$outlier.estimates
    }
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out,
        possible.contamination = 0
      )
    )
  }
