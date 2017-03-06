# Estimation of fluorescence signal from the C1 Chip (internal functions)

#' unnormalizeC01
#'
#' It calculates the unnormalized signal of a normalized image. The normalization has been originally performed
#'   as x[i,j]/bits where x[i,j] is the i,j element of the signal matrix x and bits is the value of the bits.
#'
#' @param data Data matrix. The matrix of the normalized image signals.
#' @param bits Numeric. The image bits.
#'
#' @return A matrix of unnormalized image signals
#'
#' @keywords internal
unnormalizeC01 <- function(data, bits) {
  data <- t(bits * data)
  return(data)
}



#' spotStats
#'
#' It produces a table of estimated spot locations and fluorescence signals accompanied by informative plots.
#'   It can process the results of either spotCoords() for fluorescence-based estimation or forceBF() for
#'   BF image modelling estimation.
#'
#' @param img Data matrix. The original data of the BF image to be read and processed.
#' @param chaImgs List. A list of the channel images (data matrices) of a sample.
#' @param binChaImgs List. A list of binary segmented channel image data.
#' @param center Data matrix. The 2-dimensional location of the spot's representative center.
#' @param BFcoords List. A list of statistics describing the C1 chip line patterns (for BF image modelling).
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is
#'   mainly used in BF image modeling where a fluorescence spot could not be originally detected. The value of
#'   this parameter is also used as a cut-off to find matched spots across channel of the same sample image.
#' @param log.transform Logical. If TRUE the image data are plotted in the log scale.
#' @param warning Character string. An indicator of the estimation type that has been internally performed, i.e.
#'   fluorescence-based or BF image modelling.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators from the
#'   image file names (see <<separator2>> in readFiles()).
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study. They refer to
#'   the ImageType part of the original image or txt file names.
#'
#' @return A table of location and fluorescence estimated with accompanied plots
#'
#' @keywords internal
SpotStats <-
  function(img,
           chaImgs,
           binChaImgs,
           center,
           other.spots,
           BFcoords,
           BFarea,
           log.transform,
           warning,
           minDiff,
           separator,
           image.type) {
    area <- list(c(), 0)
    foreG <- foreR <- rep(0, 2)
    backG <- backR <- rep(0, 2)
    BFimg <- readOriImg(
      imgName = img,
      despeckle = FALSE,
      pix = 0,
      thresh = 0,
      separator = separator,
      image.type = image.type
    )
    ll <- as.numeric(unlist(lapply(BFcoords[1:2], length)))
    N <- c(1, 1)
    
    area <-
      joinAreas(
        areaR = binChaImgs$CH1,
        areaG = binChaImgs$CH2,
        center = center,
        chaImg =
          c(),
        areaBased = BFarea,
        warning = warning
      )
    
    if (warning == "Both.Channels" | warning == "One.Channel") {
      foreR <- measureF(img = chaImgs$CH1,
                        area = area[[1]],
                        BFarea = c())
      foreG <-
        measureF(img = chaImgs$CH2,
                 area = area[[1]],
                 BFarea = c())
      backR <-
        measureB(
          img = chaImgs$CH1,
          area = area[[1]],
          iter = 100,
          BFarea = c()
        )
      backG <-
        measureB(
          img = chaImgs$CH2,
          area = area[[1]],
          iter = 100,
          BFarea = c()
        )
      est <-
        plotImages(
          number = 1,
          origImg = BFimg[[1]],
          chaImgs = chaImgs,
          binChaImgs = binChaImgs,
          stats =
            list(backCH1 = backR[1], backCH2 = backG[1]),
          pix = area[[1]],
          log.transform =
            log.transform,
          minDiff = minDiff,
          sample = BFimg[[3]],
          image.type = image.type
        )
    }
    if (warning == "BF" |
        warning == "BFmedian" | warning == "BFmedian2") {
      foreR <-
        measureF(img = chaImgs$CH1,
                 area = area[[1]],
                 BFarea = BFarea)
      foreG <-
        measureF(img = chaImgs$CH2,
                 area = area[[1]],
                 BFarea = BFarea)
      backR <-
        measureB(
          img = chaImgs$CH1,
          area = area[[1]],
          iter = 100,
          BFarea = BFarea
        )
      backG <-
        measureB(
          img = chaImgs$CH2,
          area = area[[1]],
          iter = 100,
          BFarea = BFarea
        )
      est <-
        plotImages(
          number = 2,
          origImg = BFimg[[1]],
          chaImgs = chaImgs,
          binChaImgs = binChaImgs,
          stats = list(backCH1 = backR[1], backCH2 = backG[1]),
          pix = list(foreR[[2]], foreG[[2]]),
          log.transform = log.transform,
          minDiff =
            minDiff,
          sample = BFimg[[3]],
          image.type = image.type
        )
    }
    
    l1 <-
      log(as.numeric(foreR[[1]][1]) + 1, 2) - log(as.numeric(backR[1]) + 1, 2)
    l2 <-
      log(as.numeric(foreG[[1]][1]) + 1, 2) - log(as.numeric(backG[1]) + 1, 2)
    res <-
      list(
        stats = c(
          BFimg[[3]],
          center,
          area[[2]],
          warning,
          foreR[[1]][1],
          backR[1],
          foreG[[1]][1],
          backG[1],
          l1,
          est[1],
          l2,
          est[2],
          other.spots
        ),
        pixels = area[[1]]
      )
    return(res)
  }


#' joinAreas
#'
#' It estimates the spot area by joining the red and green bright spot location estimates.
#'
#' @param areaR Data matrix. The bright (spot) coordinates in one channel.
#' @param areaG Data matrix. The bright (spot) coordinates in the other channel.
#' @param center Data matrix. The 2-dimensional location of the spot's center.
#' @param chaImgs Data matrix. The channel binary segmented image data that are used for reference to obtain the
#'    area of interest.
#' @param areaBased Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#' @param warning Character string. An indicator of the estimation type that has been internally performed, i.e.
#'   fluorescence-based or BF image modelling.
#'
#' @return The coordinates of the spot area and its length in pixels.
#'
#' @keywords internal
joinAreas <- function(areaR,
                      areaG,
                      center,
                      chaImg,
                      areaBased,
                      warning) {
  if (warning == "Both.Channels" | warning == "One.Channel") {
    if (is.list(areaR) == TRUE) {
      areaR <- areaR[[1]]
    }
    if (is.list(areaG) == TRUE) {
      areaG <- areaG[[1]]
    }
    area <- unique(matrix(rbind(areaR, areaG), ncol = 2))
    larea <- nrow(area)
    if (larea <= areaBased) {
      a <- max(round(areaBased / 4, 0), 2)
      r1 <-
        max(1, (center[1] - a)):min(nrow(chaImg), (center[1] + a))
      r2 <-
        max(1, (center[2] - a)):min(ncol(chaImg), (center[2] + a))
      area <- as.matrix(expand.grid(r1, r2))
      larea <- areaBased ^ 2
    }
  } else {
    if (warning == "BF") {
      a <- floor(areaBased)
    }
    if (warning == "BFmedian" | warning == "BFmedian2") {
      a <- 2 * floor(areaBased)
    }
    r1 <-
      max(1, (center[1] - a)):min(nrow(chaImg), (center[1] + a))
    r2 <-
      max(1, (center[2] - a)):min(ncol(chaImg), (center[2] + a))
    area <- as.matrix(expand.grid(r1, r2))
    larea <- areaBased ^ 2
  }
  
  return(list(area, larea))
}


#' plotImages
#'
#' It generates the plotted results.
#'
#' @param number Integer. An index number that regulates the fluorescence estimation procedure.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. A list of the channel images (data matrices) of a sample.
#' @param binChaImgs List. A list of binary segmented channel image data.
#' @param stats List. A series of foreground and background statistics.
#' @param pix List. A list of bright spot coordinates in the channels.
#' @param log.transform Logical. If TRUE the image data are plotted in the log scale.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param sample Character string. The sample ID.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study. They refer to
#'   the ImageType part of the original image or txt file names.
#'
#' @return the plotted results and statistics of the signal-to-noise ratio
#'
#' @import ggplot2 reshape2 stats
#' @keywords internal
plotImages <-
  function(number,
           origImg,
           chaImgs,
           binChaImgs,
           stats,
           pix,
           log.transform,
           minDiff,
           sample,
           image.type) {
    if (number == 1) {
      dataR <- dataG <- rep(0, nrow(pix))
      for (i in 1:nrow(pix)) {
        dataR[i] <- chaImgs$CH1[pix[i, 1], pix[i, 2]]
        dataG[i] <- chaImgs$CH2[pix[i, 1], pix[i, 2]]
      }
      i1 <-
        chaImgs$CH1[max(1, (min(pix[, 1]) - 10)):min(nrow(chaImgs$CH1), (max(pix[, 1]) + 10)),
                    max(1, (min(pix[, 2]) - 10)):min(nrow(chaImgs$CH1), (max(pix[, 2]) + 10))]
      i2 <-
        chaImgs$CH2[max(1, (min(pix[, 1]) - 10)):min(nrow(chaImgs$CH2), (max(pix[, 1]) + 10)),
                    max(1, (min(pix[, 2]) - 10)):min(nrow(chaImgs$CH2), (max(pix[, 2]) + 10))]
      
      if (log.transform == TRUE) {
        i1 <- log(i1 + 1)
        i2 <- log(i2 + 1)
      }
      if (length(dataR) == 1) {
        dataR <- rep(dataR, 2)
      }
      if (length(dataG) == 1) {
        dataG <- rep(dataG, 2)
      }
      
      orimage <-
        origImg[max(1, (min(pix[, 1]) - 40)):min(nrow(origImg), (max(pix[, 1]) + 40)),
                max(1, (min(pix[, 2]) - 40)):min(ncol(origImg), (max(pix[, 2]) + 40))]
      colnames(orimage) <- NULL
      rownames(orimage) <- NULL
      
      options(warn = -1)
      res1 <-
        wilcox.test((log(dataR + 1, 2) - log(stats$backCH1 + 1, 2)),
                    mu = minDiff, alternative = "greater")$p.value
      res2 <-
        wilcox.test((log(dataG + 1, 2) - log(stats$backCH2 + 1, 2)),
                    mu = minDiff, alternative = "greater")$p.value
      
      res <- c(res1, res2)
    }
    
    if (number == 2) {
      dataR <- matrix(0, length(pix[[1]]), nrow(pix[[1]][[1]][[1]]))
      dataG <- matrix(0, length(pix[[2]]), nrow(pix[[2]][[1]][[1]]))
      rangR <- matrix(0, length(pix[[1]]), 4)
      rangG <- matrix(0, length(pix[[2]]), 4)
      for (i in 1:length(pix[[1]])) {
        for (j in 1:nrow(pix[[1]][[i]][[1]])) {
          dataR[i, j] <-
            chaImgs$CH1[pix[[1]][[i]][[1]][j, 1], pix[[1]][[i]][[1]][j, 2]]
          rangR[i,] <-
            c(min(pix[[1]][[i]][[1]][, 1]),
              max(pix[[1]][[i]][[1]][, 1]),
              min(pix[[1]][[i]][[1]][, 2]),
              max(pix[[1]][[i]][[1]][, 2]))
        }
      }
      dataR <- apply(dataR, 2, mean)
      for (i in 1:length(pix[[2]])) {
        for (j in 1:nrow(pix[[2]][[i]][[1]])) {
          dataG[i, j] <-
            chaImgs$CH2[pix[[2]][[i]][[1]][j, 1], pix[[2]][[i]][[1]][j, 2]]
          rangG[i,] <-
            c(min(pix[[2]][[i]][[1]][, 1]),
              max(pix[[2]][[i]][[1]][, 1]),
              min(pix[[2]][[i]][[1]][, 2]),
              max(pix[[2]][[i]][[1]][, 2]))
        }
      }
      dataG <- apply(dataG, 2, mean)
      
      i1 <-
        chaImgs$CH1[max(1, (min(rangR[, 1]) - 10)):min(nrow(chaImgs$CH1), (max(rangR[, 2]) + 10)),
                    max(1, (min(rangR[, 3]) - 10)):min(ncol(chaImgs$CH1), (max(rangR[, 4]) + 10))]
      i2 <-
        chaImgs$CH2[max(1, (min(rangG[, 1]) - 10)):min(nrow(chaImgs$CH2), (max(rangG[, 2]) + 10)),
                    max(1, (min(rangG[, 3]) - 10)):min(ncol(chaImgs$CH2), (max(rangG[, 4]) + 10))]
      if (log.transform == TRUE) {
        i1 <- log(i1 + 1)
        i2 <- log(i2 + 1)
      }
      if (length(dataR) == 1) {
        dataR <- rep(dataR, 2)
      }
      if (length(dataG) == 1) {
        dataG <- rep(dataG, 2)
      }
      
      imgrowsR <-
        max(1, (min(rangR[, 1]) - 40)):min(nrow(chaImgs$CH1), (max(rangR[, 2]) + 40))
      imgrowsG <-
        max(1, (min(rangG[, 1]) - 40)):min(nrow(chaImgs$CH2), (max(rangG[, 2]) + 40))
      imgcolsR <-
        max(1, (min(rangR[, 3]) - 40)):min(ncol(chaImgs$CH1), (max(rangR[, 4]) + 40))
      imgcolsG <-
        max(1, (min(rangG[, 3]) - 40)):min(ncol(chaImgs$CH2), (max(rangG[, 4]) + 40))
      
      
      orimage <-
        origImg[min(imgrowsR, imgrowsG):max(imgrowsR, imgrowsG),
                min(imgcolsR, imgcolsG):max(imgcolsR, imgcolsG)]
      colnames(orimage) <- NULL
      rownames(orimage) <- NULL
      
      #res1<-t.test(d1$x,alternative="greater")$p.value
      #res2<-t.test(d2$x,alternative="greater")$p.value
      options(warn = -1)
      res1 <-
        wilcox.test((log(dataR + 1, 2) - log(stats$backCH1 + 1, 2)), mu = minDiff, alternative =
                      "greater")$p.value
      res2 <-
        wilcox.test((log(dataG + 1, 2) - log(stats$backCH2 + 1, 2)), mu = minDiff, alternative =
                      "greater")$p.value
      
      res <- c(res1, res2)
    }
    
    if (log.transform == TRUE) {
      img.m <- melt(log(orimage + 1, 2))
    } else {
      img.m <- melt(orimage)
    }
    
    names(img.m) <- c("x", "y", "z")
    plot_original <-
      ggplot(img.m, aes_string(x = 'x', y = 'y', fill = 'z')) + geom_raster() + theme_bw() + scale_fill_distiller(palette = "Greys") +
      ggtitle(paste(image.type[1], " spot (", sample, ")", sep = "")) + theme(legend.position =
                                                                                "none")
    
    #  filled contour
    colnames(i1) <- NULL
    i1.m <- melt(i1)
    names(i1.m) <- c("x", "y", "z")
    colnames(i2) <- NULL
    i2.m <- melt(i2)
    names(i2.m) <- c("x", "y", "z")
    
    plot_contour1 <-
      ggplot(i1.m, aes_string(x = 'x', y = 'y', fill = 'z')) + geom_raster(aes_string(fill = 'z')) + coord_equal() +
      # stat_contour(geom = "polygon",
      #              aes_string(fill = '..level..'),
      #              inherit.aes = TRUE) +
      scale_fill_distiller(palette = "Spectral") + theme_bw() + theme(legend.position =
                                                                        "none") +
      ggtitle(paste(image.type[2], " spot (", sample, ")", sep = ""))
    plot_contour2 <-
      ggplot(i2.m, aes_string(x = 'x', y = 'y', fill = 'z')) + geom_raster(aes_string(fill = 'z')) + coord_equal() +
      # stat_contour(geom = "polygon",
      #              aes_string(fill = '..level..'),
      #              inherit.aes = TRUE) +
      scale_fill_distiller(palette = "Spectral") + theme_bw() + theme(legend.position =
                                                                        "none") +
      ggtitle(paste(image.type[3], " spot (", sample, ")", sep = ""))
    
    #density plot
    d12 <- data.frame(
      d1 = log(dataR + 1, 2) - log(stats$backCH1 + 1, 2),
      d2 = log(dataG + 1, 2) - log(stats$backCH2 + 1, 2)
    )
    d12.m <- melt(d12, measure.vars = c("d1", "d2"))
    plot_dens <-
      ggplot(d12.m) + geom_density(aes_string(x = 'value', linetype = 'variable')) + labs(x = NULL) +
      geom_vline(xintercept = 0,
                 colour = "grey",
                 linetype = "longdash") +
      geom_hline(yintercept = 0, colour = "grey") +
      ggtitle("BG-subtracted log-signals") +
      theme_bw() + theme(legend.position = "none")
    
    suppressWarnings(multiplot(
      plot_original,
      plot_dens,
      plot_contour1,
      plot_contour2,
      cols = 2
    ))
    
    return(res)
  }

#' getCoordinates_stats
#'
#' It finds the spot coordinates using the spot centers or BF image modeling
#'
#' @param centerR Data matrix. The location statistics in one channel.
#' @param centerG Data matrix. The location statistics in thew other channel.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param chaImgs List. A list of the channel images (data matrices) of a sample.
#' @param ll Data matrix. An internal parameter specifying the spot center.
#' @param ws List. An internal parameter specifying the spot center.
#' @param estCenter Data matrix. The estimated spot center by BF image modelling.
#'
#' @return A series of location estimates includying the channel-specific spot center and spot areas
#'
#' @import stats
#' @keywords internal
getCoordinates_stats <-
  function(centerR,
           centerG,
           minDiff,
           chaImgs,
           ll,
           ws,
           estCenter) {
    ss <-
      significantSignal(
        centerR = centerR,
        centerG = centerG,
        minDiff = minDiff,
        chaImgs = chaImgs
      )
    signifR <- ss$pvR
    signifG <- ss$pvG
    rres <- c(ll[ws$wR[1], 1], ll[ws$wG[1], 2])
    
    if (signifR[rres[1]] == 1 & signifG[rres[2]] == 0) {
      center <- centerR[[1]]
      arR <- centerR$brightcoordinates[[rres[1]]]
      arG <- arR
      sizeR <- centerR$warning[[rres[1]]]
      sizeG <- 0
      rr <- unique(ll[ws$wR1, 1])
      nspotsR <- centerR$warning[rr]
      nn <- names(nspotsR)
      nspotsR <- nspotsR[na.omit(nn)]
      nspotsG <- list(c())
    }
    
    if (signifR[rres[1]] == 0 & signifG[rres[2]] == 1) {
      center <- centerG[[1]]
      arG <- centerG$brightcoordinates[[rres[2]]]
      arR <- arG
      sizeG <- centerG$warning[[rres[2]]]
      sizeR <- 0
      rr <- unique(ll[ws$wG1, 2])
      nspotsG <- centerG$warning[rr]
      nn <- names(nspotsG)
      nspotsG <- nspotsG[na.omit(nn)]
      nspotsR <- list(c())
    }
    
    if (signifR[rres[1]] == 1 & signifG[rres[2]] == 1) {
      wei <-
        matrix(rbind(abs(centerR[[1]][rres[1],] - estCenter), abs(centerG[[1]][rres[2],] -
                                                                    estCenter)), nrow = 2)
      wei <- wei + 0.01
      weiR <-
        c((1 - wei[1, 1] / sum(wei[, 1])), (1 - wei[1, 2] / sum(wei[, 2])))
      weiG <- 1 - weiR
      wm1 <-
        weighted.mean(c(centerR[[1]][rres[1], 1], centerG[[1]][rres[2], 1]), c(weiR[1], weiG[1]))
      wm2 <-
        weighted.mean(c(centerR[[1]][rres[1], 2], centerG[[1]][rres[2], 2]), c(weiR[2], weiG[2]))
      center <- c(floor(wm1), floor(wm2))
      arR <- centerR$brightcoordinates[[rres[1]]]
      arG <- centerG$brightcoordinates[[rres[2]]]
      sizeR <- centerR$warning[[rres[1]]]
      sizeG <- centerG$warning[[rres[2]]]
      rr <- unique(ll[ws$wR1, 1])
      nspotsR <- centerR$warning[rr]
      nn <- names(nspotsR)
      nspotsR <- nspotsR[na.omit(nn)]
      rr <- unique(ll[ws$wG1, 2])
      nspotsG <- centerG$warning[rr]
      nn <- names(nspotsG)
      nspotsG <- nspotsG[na.omit(nn)]
    }
    
    if (signifR[rres[1]] == 0 & signifG[rres[2]] == 0) {
      w <- which(ss$ts == min(ss$ts))
      if (w == 1 & length(w) == 1) {
        center <- centerR[[1]][rres[1],]
        arR <- centerR$brightcoordinates[[rres[1]]]
        arG <- arR
        sizeR <- centerR$warning[[rres[1]]]
        sizeG <- 0
        rr <- unique(ll[ws$wR1, 1])
        nspotsR <- centerR$warning[rr]
        nn <- names(nspotsR)
        nspotsR <- nspotsR[na.omit(nn)]
        nspotsG <- list(c())
      }
      if (w == 2 & length(w) == 1) {
        center <- centerG[[1]][rres[2],]
        arG <- centerG$brightcoordinates[[rres[2]]]
        arR <- arG
        sizeG <- centerG$warning[[rres[2]]]
        sizeR <- 0
        rr <- unique(ll[ws$wG1, 2])
        nspotsG <- centerG$warning[rr]
        nn <- names(nspotsG)
        nspotsG <- nspotsG[na.omit(nn)]
        nspotsR <- list(c())
      }
      if (length(w) > 1) {
        center <- c(0, 0)
        arG <- 0
        arR <- 0
        sizeG <- 0
        sizeR <- 0
        nspotsG <- list(c())
        nspotsR <- list(c())
      }
    }
    return(
      list(
        center = center,
        areaR = arR,
        areaG = arG,
        sizeR = sizeR,
        sizeG = sizeG,
        nspotsR =
          nspotsR,
        nspotsG = nspotsG
      )
    )
  }


#' zoomInBF
#'
#' It estimates the chip characterisrics for BF image modelling
#'
#' @param img Data matrix. The BF image data.
#' @param pattern.search Integer. A cutoff to find horizontal and vertical lines on the chip.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central
#'   area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#'
#' @return The locations of the straight lines on the chip
#'
#' @keywords internal
zoomInBF <- function(img, pattern.search, ImgLimits, chip.type) {
  scl <-
    straightColLines(
      img = img,
      pattern.search = pattern.search,
      ImgLimits = ImgLimits,
      chip.type = chip.type
    )
  srl <-
    straightRowLines(colData = scl,
                     pattern.search = pattern.search,
                     ImgLimits =
                       ImgLimits)
  res <-
    list(
      reducedCols = scl$Cols,
      reducedRows = scl$Rows,
      strLines = scl$mainStrLines,
      strLinesFar =
        scl$farLines,
      straightRows = srl[[2]],
      straightRowsIndex = srl[[1]],
      binImgDespeckled =
        scl$spotMat
    )
  return(res)
}


#' SpotbyStrLines
#'
#' It estimates the spot location using the BF image modeling parameters
#'
#' @param binImg List. The binary segmented image data in each channel.
#' @param pattern.search Integer. A cutoff to find horizontal and vertical lines on the chip.
#' @param stats List. The estimated parameters of the BF image modelling.
#'
#' @return A series of spot location estimates
#'
#' @keywords internal
SpotbyStrLines <- function(binImg, pattern.search, stats) {
  ind <- rep(0, 2)
  centerRows1 <- c()
  centerRows2 <- c()
  
  res1 <- c(0, 0)
  if (stats$strLinesFar[1] > 0) {
    centerRows1 <- ceiling(mean(stats$reducedRows))
    ind[1] <- 1
    if (mean(stats$strLinesFar) < (0.5 * ncol(binImg))) {
      res1 <- c(centerRows1, min(stats$strLines))
    } else {
      res1 <- c(centerRows1, max(stats$strLines))
    }
  }
  
  res2 <- c(0, 0)
  if (length(stats$straightRowsIndex) > 0) {
    centerRows2 <- ceiling(mean(stats$straightRows))
    ind[2] <- 1
    if (stats$straightRowsIndex[1] == 1) {
      res2 <- c(centerRows2, min(stats$strLines))
    } else {
      res2 <- c(centerRows2, max(stats$strLines))
    }
  }
  
  if (sum(ind) == 2) {
    res <-
      c(floor(mean(c(res1[1], res2[1]))), round(mean(c(res1[2], res2[2])), 0))
  }
  if (sum(ind) == 1) {
    if (ind[1] == 1) {
      res <- res1
    }
    if (ind[2] == 1) {
      res <- res2
    }
  }
  if (sum(ind) == 0) {
    if (length(stats$strLines[stats$strLines > (0.5 * ncol(binImg))]) >= 0.5) {
      res <- c(ceiling(mean(stats$reducedRows)), max(stats$strLines))
    } else {
      res <- c(ceiling(mean(stats$reducedRows)), min(stats$strLines))
    }
  }
  return(res)
}

#' readOriImg
#'
#' It reads and processes the original BF image
#'
#' @param imgName Character string. The file name of the image.
#' @param despeckle Logical. If TRUE the BF image is despeckled.
#' @param pix Integer. A cutoff specifing the size of the area to search for speckles.
#' @param thresh Integer. A cutoff to perform the despeckle function. If pixel signal > median object signal + thresh,
#'   the object is a speckle and the median object signal is returned.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study. They refer to
#'   the ImageType part of the original image or txt file names.
#'
#' @return A list of BF image estimates
#'
#' @import utils
#'
#' @keywords internal
readOriImg <-
  function(imgName,
           despeckle,
           pix,
           thresh,
           separator,
           image.type) {
    lab <- unlist(strsplit(imgName, "/"))
    lab <-
      unlist(strsplit(lab[length(lab)], paste(separator, image.type[1], ".txt", sep =
                                                "")))
    temp <- as.matrix(read.table(imgName, sep = "\t"))
    x <- temp
    if (despeckle == TRUE) {
      x <- despecklefun(img = temp,
                        pix = pix,
                        thresh = thresh)
    }
    return(list(temp, x, lab, despeckle))
  }


#' readChaImg
#'
#' It reads and processes the channel image data
#'
#' @param imgName Character string. The file name of the images.
#' @param denoise Logical. If TRUE the channel image is denoised with 2-dimesnional la8, universal
#'   and hard thresholding.
#'
#' @return A list of channel image estimates
#'
#' @import utils
#'
#' @keywords internal
readChaImg <- function(imgNames, denoise) {
  tempR <- as.matrix(read.table(imgNames$CH1, sep = "\t"))
  tempG <- as.matrix(read.table(imgNames$CH2, sep = "\t"))
  tempRr <- processImg(img = tempR, denoise = denoise)
  tempGr <- processImg(img = tempG, denoise = denoise)
  return(list(
    CH1 = tempR,
    CH2 = tempG,
    CH1.proc = tempRr,
    CH2.proc = tempGr
  ))
}


#' forceBF
#'
#' It re-estimates the location of the outlier samples
#'
#' @param data List. The location statistics of both channels.
#' @param cutoff Integer. A cutoff do detect outlier locations.
#' @param median.correction Logical. If TRUE it performs median adjustment for outlier locations.
#' @param medians Data matrix. The estimated medians of non-outlier samples by run and well ID.
#' @param Wells Data matrix. The directionality of the well IDs.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'   They refer to the ImageType part of the original image or txt file names.
#'
#' @return A list of re-estimated locations
#'
#' @keywords internal
forceBF <-
  function(data,
           cutoff,
           median.correction,
           medians,
           Wells,
           image.type) {
    res <- as.list(rep(0, 4))
    names(res) <-
      c("center", "areaR", "areaG", "warning")
    
    spot.estimates <- data[[1]]
    QCdata <- data[[2]]
    w1 <- which(Wells[1] == medians[, 1] & Wells[3] == medians[, 2])
    
    if (spot.estimates[1] != QCdata$sample) {
      stop(paste("The location and the ", image.type[1], " data do not match!", sep =
                   ""))
    }
    
    if (spot.estimates[ncol(spot.estimates)] == "outlier") {
      BFcrit <-
        max(abs(QCdata$centerR[1] - as.numeric(medians[w1, 3])),
            abs(QCdata$centerR[2] - as.numeric(medians[w1, 4])))
      Fluocrit <-
        max(abs(as.numeric(spot.estimates[2]) - as.numeric(medians[w1, 3])),
            abs(as.numeric(spot.estimates[3]) - as.numeric(medians[w1, 4])))
      
      if (BFcrit > cutoff & Fluocrit > cutoff) {
        res$center <- round(as.numeric(medians[w1, 3:4]), 0)
        if (QCdata$centerR[1] > 0) {
          res$warning <- giveWarning(4)
        } else {
          res$warning <- giveWarning(5)
        }
        res$areaR = res$areaG = 0
      }
      
      if (BFcrit <= cutoff) {
        res$center <- QCdata$centerR
        res$areaR = res$areaG = 0
        res$warning = QCdata$warn
      }
      
      if (BFcrit > cutoff & Fluocrit <= cutoff) {
        if (median.correction == TRUE) {
          res$center <- round(as.numeric(medians[w1, 3:4]), 0)
          res$areaR = res$areaG = 0
          res$warning = giveWarning(4)
        } else {
          res$center <- QCdata$centerR
          res$areaR = res$areaG = 0
          res$warning = QCdata$warn
        }
      }
    } else {
      print(
        paste(
          "Sample ",
          spot.estimates[1],
          " is not an outlier. The correction is not performed",
          sep = ""
        )
      )
    }
    
    if (is.null(res$warning)) {
      res$warning <- "Both.Channels"
    }
    
    return(
      list(
        center = res$center,
        areaCH1 = res$areaR,
        areaCH2 = res$areaG,
        other.spots = spot.estimates$Other.Spots,
        warning =
          res$warning,
        failure = c(),
        Outlier.Estimates = c()
      )
    )
  }



#' distfromcenter
#'
#' It calculates the Eucledian distance between two 2-dimensional locations.
#'
#' @param data Data matrix. A 2-dimensional location.
#' @param center Data matrix. A second 2-dimensional location.
#'
#' @return The Eucledian distance between the two locations
#'
#' @keywords internal
distfromcenter <- function(data, center) {
  sqrt(((data[1] - center[1]) ^ 2) + ((data[2] - center[2]) ^ 2))
}


#' grubbs
#'
#' It performs the grubbs test for outliers.
#'
#' @param data Numeric vector. An 1-dimensional vector of spot distances to check for outliers.
#'
#' @return All potential outliers (indices)
#'
#' @keywords internal
grubbs <- function(data) {
  p <- 0
  data1 <- data
  result <- matrix(0, 1, 2)
  while (p <= 0.05 & length(data1) > 0) {
    res <- grubbs.test(data1)
    p <- res[[3]]
    val <-
      as.numeric(unlist(strsplit(unlist(
        strsplit(res[[2]], "value ")
      )[2], " is"))[1])
    result <- matrix(rbind(result, c(val, p)), ncol = 2)
    ww <- which(abs(data1 - val) == min(abs(data1 - val)))
    data1 <- data1[-ww]
  }
  p <- 0
  data1 <- data
  while (p <= 0.05 & length(data1) > 0) {
    res <- grubbs.test(data1, opposite = TRUE)
    p <- res[[3]]
    val <-
      as.numeric(unlist(strsplit(unlist(
        strsplit(res[[2]], "value ")
      )[2], " is"))[1])
    result <- matrix(rbind(result, c(val, p)), ncol = 2)
    ww <- which(abs(data1 - val) == min(abs(data1 - val)))
    data1 <- data1[-ww]
  }
  result <- result[-1,]
  result <- result[which(result[, 2] <= 0.05), 1]
  ww <- c()
  if (length(result) > 0) {
    for (i in 1:length(result)) {
      ww <-
        c(ww, which(abs(data - result[i]) == min(abs(data - result[i]))))
    }
  }
  return(ww)
}

#' getSpot
#'
#' Identifies one or multiple spot(s) in the image data matrix.
#'
#' @param img Data matrix. The binary segmented channel image data.
#' @param rad Integer. A cut-off to separate the spots (densely clustered 1s) from small random signals
#'   (loosely clustered 1s) on the binary image.
#'
#' @return A matrix of bright (spot) coordinates
#'
#' @keywords internal
getSpot <- function(img, rad) {
  mm <- matrix(0, nrow(img), ncol(img))
  for (i in (rad + 1):(nrow(img) - rad)) {
    ro <- (i - rad):(i + rad)
    for (j in (rad + 1):(ncol(img) - rad)) {
      co <- (j - rad):(j + rad)
      mm[i, j] <- mean(c(img[ro, co]))
    }
  }
  ww <- matrix(which(mm == max(mm), arr.ind = TRUE), ncol = 2)[1,]
  rad1 <- seq(2, rad, 1)
  stat <- c()
  for (i in 1:length(rad1)) {
    mat <-
      c(img[max(1, (ww[1] - rad1[i])):min(nrow(img), (ww[1] + rad1[i])), max(1, (ww[2] -
                                                                                   rad1[i])):min(ncol(img), (ww[2] + rad1[i]))])
    stat <- c(stat, length(mat[mat > 0]) / length(mat))
  }
  #return(c(ww,stat))
  return(ww)
}


#' spotCenter
#'
#' It estimates a series of spot statistics on each channel.
#'
#' @param img Data matrix. The channel image data.
#' @param foregroundCut Numeric vector. The binary segmentation image analysis cutoffs
#'   for normalized image data. Pixels with normalized signals higher than the cutoff
#'   belong to foreground.
#' @param howbig Integer. An user defined value of the minimum expected spot size.
#' @param ImgLimits a cutoff that determines where the spot is supposed to be found.
#'
#' @import stats
#' @importFrom EBImage Image bwlabel fillHull paintObjects
#' @importFrom raster quantile
#'
#' @return A list of spot coordinate estimates
#'
#' @keywords internal
spotCenter <- function(img, foregroundCut, howbig, ImgLimits) {
  howbig <- 1 / howbig
  IMG <- Image(img, c(nrow(img), ncol(img)))
  cut <- foregroundCut
  i <- 1
  numbers <- 1
  spots <- matrix(1, nrow(img), ncol(img))
  listedCritSeries <- matrix(0, 2, 1)
  crit <- 1
  while (numbers > 1 &
         numbers != 0 &
         i <= length(cut) |
         max(crit) > howbig  & i <= length(cut)) {
    img1 <- IMG > cut[i]
    imglabel = bwlabel(img1)
    imglabel = fillHull(imglabel)
    spots.suggested <- paintObjects(img1, imglabel)@.Data
    spots.suggested[spots.suggested < 1] <- 0
    s1 <- sort(unique(c(spots.suggested)))
    s1 <- s1[s1 > 0.5]
    if (length(s1) > 0) {
      spots <- spots.suggested
      crit <- rep(1, length(s1))
      for (j in 1:length(s1)) {
        spots[spots == s1[j]] <- j
        crit[j] <- length(spots[spots == j]) / length(spots)
      }
      critSeries <-
        matrix(rbind(matrix(1:j, nrow = 1), matrix(crit, nrow = 1)), nrow = 2)
      listedCritSeries <-
        matrix(cbind(listedCritSeries, critSeries), nrow = 2)
    }
    if (length(s1) == 0) {
      critSeries <- matrix(0, 2, 0)
    }
    numbers <- max(spots)
    i <- i + 1
  }
  if (numbers == 0) {
    numbers <- ncol(critSeries)
  }
  listedCritSeries <- matrix(listedCritSeries[,-1], nrow = 2)
  critSeries <-
    matrix(rbind(matrix(0, 1, ncol(critSeries)), critSeries), nrow = 3)
  listedCritSeries <-
    matrix(rbind(matrix(0, 1, ncol(listedCritSeries)), listedCritSeries), nrow =
             3)
  
  
  if (ncol(critSeries) > 0) {
    w <- which(critSeries[3,] > howbig)
    if (length(w) > 0 & length(w) < ncol(critSeries)) {
      critSeries <- matrix(critSeries[,-w], nrow = 3)
      for (i in 1:length(w)) {
        spots[which(spots == w[i])] <- 0
      }
      for (i in 1:ncol(critSeries)) {
        critSeries[1, i] <- i
        spots[which(spots == critSeries[2, i])] <- critSeries[1, i]
        w1 <- which(listedCritSeries[2,] == critSeries[2, i])
        listedCritSeries[1, w1] <- critSeries[1, i]
      }
      w2 <- which(listedCritSeries[1,] == 0)
      if (length(w2) > 0) {
        listedCritSeries <- matrix(listedCritSeries[,-w2], nrow = 3)
      }
      numbers <- ncol(critSeries)
    }
    if (length(w) == 0 | length(w) == ncol(critSeries)) {
      critSeries[1,] <- critSeries[2,]
      for (i in 1:ncol(critSeries)) {
        w1 <- which(listedCritSeries[2,] == critSeries[2, i])
        listedCritSeries[1, w1] <- critSeries[1, i]
      }
      w2 <- which(listedCritSeries[1,] == 0)
      if (length(w2) > 0) {
        listedCritSeries <- matrix(listedCritSeries[,-w2], nrow = 3)
      }
    }
    numbers <- ncol(critSeries)
    critSeries <- matrix(critSeries[-2,], nrow = 2)
    listedCritSeries <- matrix(listedCritSeries[-2,], nrow = 2)
    listedCritSeries <-
      as.list(aggregate(listedCritSeries[2,], list(listedCritSeries[1,]), list))
    names(listedCritSeries$x) <- listedCritSeries$Group.1
    listedCritSeries <- listedCritSeries$x
  }
  
  if (ncol(critSeries) == 0) {
    spot <- matrix(0, nrow(img), ncol(img))
    numbers <- 0
  }
  
  if (numbers == 0) {
    res <- rep(0, 2)
    spots <- c()
    ind <- c()
    critSeries <- c()
    listedCritSeries <- c()
  }
  
  if (numbers > 0) {
    N <- numbers
    if (N == 1) {
      ind <- which(spots == 1, arr.ind = TRUE)
      #here: it was +/-10
      indR <- (min(ind[, 1]) - 4):(max(ind[, 1]) + 4)
      indC <- (min(ind[, 2]) - 4):(max(ind[, 2]) + 4)
      if (max(indR) <= (nrow(spots) - ImgLimits) &
          min(indR) >= ImgLimits &
          max(indC) <= (ncol(spots) - ImgLimits) &
          min(indC) >= ImgLimits) {
        ss <- spots[indR, indC]
        i1 <- img[indR, indC]
        res <- getSpot(img = (i1 * ss), rad = 3)
        res[1] <- indR[res[1]]
        res[2] <- indC[res[2]]
      } else {
        res <- c(0, 0)
        ind <- c()
      }
    }
    if (N > 1) {
      ind <- as.list(rep(0, N))
      ss <- as.list(rep(0, N))
      res <- matrix(0, N, 2)
      for (i in 1:N) {
        ind[[i]] <- which(spots == i, arr.ind = TRUE)
        for (b in 1:2) {
          ind[[i]] <-
            unique(matrix(rbind(
              ind[[i]], (ind[[i]] + c(0, 1)), (ind[[i]] - c(0, 1)), (ind[[i]] + c(1, 0)), (ind[[i]] -
                                                                                             c(1, 0))
            ), ncol = 2))
        }
        #here: it was +/- 10
        indR <- (min(ind[[i]][, 1]) - 4):(max(ind[[i]][, 1]) + 4)
        indC <- (min(ind[[i]][, 2]) - 4):(max(ind[[i]][, 2]) + 4)
        if (max(indR) <= nrow(spots) &
            min(indR) >= 1 &
            max(indC) <= ncol(spots) & min(indC) >= 1) {
          ss[[i]] <- spots[indR, indC]
          i1 <- img[indR, indC]
          res[i,] <-
            getSpot(img = (i1 * ss[[i]]), rad = 3)
          res[i, 1] <- indR[res[i, 1]]
          res[i, 2] <- indC[res[i, 2]]
        }
      }
      ww <- which(res[, 1] > 0)
      ww1 <- which(res[, 1] == 0)
      if (length(ww) > 0) {
        res <- res[ww,]
        ind <- ind[ww]
      } else {
        res <- c(0, 0)
        ind <- c()
      }
      if (length(ww1) > 0) {
        for (i in 1:length(ww1)) {
          spots[which(spots == ww[i])] <- 0
        }
        ss <- sort(unique(c(spots)))
        ss <- ss[ss > 0]
        if (length(ss) > 0) {
          for (i in 1:length(ss)) {
            spots[which(spots == ss[i])] <- i
          }
        }
      }
    }
  }
  if (is.list(ind) == FALSE) {
    ind <- list(ind)
  }
  res <- matrix(res, ncol = 2)
  if (nrow(res) > 1) {
    w <- as.numeric(duplicated(res))
    res <- matrix(res[which(w == 0),], ncol = ncol(res))
    ind <- ind[which(w == 0)]
    listedCritSeries <- listedCritSeries[which(w == 0)]
  }
  
  return(
    list(
      center = res,
      binSpotMat = spots,
      brightcoordinates = ind,
      warning = listedCritSeries
    )
  )
}

#' giveWarning
#'
#' It generates a coded message of the estimation type that is being performed.
#'
#' @param number Integer. An internally defined number that produces a message .
#'
#' @return A coded message
#'
#' @keywords internal
giveWarning <- function(number) {
  ll <-
    c(# estimated from both Ch2/Ch3 with or without location adjustment
      # an area is estimated
      "Both.Channels",
      
      # estimated from Ch2 or Ch3 but not both with or without location adjustment
      # an area is estimated
      "One.Channel",
      
      # Estimated from BF images only
      # the area is not estimated
      "BF",
      
      # Estimated from BF images only via median correction
      # the area is not estimated
      "BFmedian",
      
      # Estimated from BF images only via forced median correction
      # (the BF center was (0,0) indicating initial estimation from both channels)
      # the area is not estimated
      "BFmedian2")
  
  return(ll[number])
}


#' subsetAnalysis
#'
#' It reads a subset of the original file names.
#'
#' @param files List. A list of all BF and channel file names.
#' @param sub Character string. A vector of the file names to be read.
#'
#' @return A list of files to be analyzed
#'
#' @keywords internal
subsetAnalysis <- function(files, sub) {
  if (length(sub) == 1) {
    ss <- sample(1:length(ls$BF), sub)
  }
  if (length(sub) > 1) {
    ss <- sub
  }
  res <- as.list(rep(0, 3))
  names(res) <- names(files)
  for (i in 1:3) {
    res[[i]] <- files[[i]][ss]
  }
  return(res)
}

#' processImg
#'
#' Denoises the image data to be used for spot location estimation.
#'
#' @param img Data matrix. The matrix of image data from a channel.
#' @param  denoise Logical. If TRUE the channel image is denoised with 2-dimensional la8, universal
#'   and hard thresholding.
#'
#' @importFrom waveslim denoise.dwt.2d
#'
#' @return A denoised image
#'
#' @keywords internal
processImg <- function(img, denoise) {
  if (denoise == TRUE) {
    tempr <-
      denoise.dwt.2d(
        img,
        wf = "la8",
        J = 4,
        method = "universal",
        H = 0.5,
        noise.dir = 3,
        rule = "hard"
      )
    tempr <- tempr / max(tempr)
  } else {
    tempr <- img / max(img)
  }
  if (max(tempr) > 1) {
    stop("The maximum normalized signal is greater than 1")
  }
  return(tempr)
}

#' myt
#'
#' A helper to test whether the foreground signal is statistically higher than the background.
#'
#' @param data Numeric vector. The 1-dimensinal signal of log(foreground) - log(background).
#' @param minDiff Float. The mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background.
#'
#' @return A test P-value
#'
#' @import stats
#'
#' @keywords internal
myt <- function(data, minDiff) {
  u <- unique(data)
  res <- 1
  if (length(u) == 1 & u[1] > minDiff) {
    res <- 0
  }
  if (length(u) > 1) {
    res <- t.test(data, mu = minDiff, alternative = "greater")$p.value
  }
  return(res)
}


#' significantSignal
#'
#' It tests whether the foreground signal is statistically higher than the background.
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param chaImgs List. A list of the channel images (data matrices) of a sample.
#'
#' @return P-values and test statistics for both channels
#'
#' @keywords internal
significantSignal <- function(centerR, centerG, minDiff, chaImgs) {
  signifR <- rep(0, nrow(centerR[[1]]))
  tR <- 1
  if (centerR[[1]][1] > 0) {
    for (i in nrow(centerR[[1]])) {
      ar <- centerR$brightcoordinates[i]
      fore <-
        measureF(img = chaImgs$CH1,
                 area = ar[[1]],
                 BFarea = c())$intensity
      back <-
        measureB(
          img = chaImgs$CH1,
          area = ar[[1]],
          iter = 10,
          BFarea = c()
        )
      signal <- log(fore + 1, 2) - log(back[1] + 1, 2)
      tR <- myt(data = signal, minDiff = minDiff)
      if (tR < 0.01) {
        signifR[i] <- 1
      }
    }
  }
  signifG <- rep(0, nrow(centerG[[1]]))
  tG <- 1
  if (centerG[[1]][1] > 0) {
    for (i in nrow(centerG[[1]])) {
      ar <- centerG$brightcoordinates[i]
      fore <-
        measureF(img = chaImgs$CH2,
                 area = ar[[1]],
                 BFarea = c())$intensity
      back <-
        measureB(
          img = chaImgs$CH2,
          area = ar[[1]],
          iter = 10,
          BFarea = c()
        )
      signal <- log(fore + 1, 2) - log(back[1] + 1, 2)
      tG <- myt(data = signal, minDiff = minDiff)
      if (tG < 0.01) {
        signifG[i] <- 1
      }
    }
  }
  return(list(
    pvR = signifR,
    pvG = signifG,
    ts = c(tR, tG)
  ))
}

#' getCsFAIL
#'
#' A helper to re-estimate the spot location or find the capture site location
#'   (BF image modelling).
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. A list of the channel images (data matrices) of a sample.
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background.
#' @param despeckle Logical. If TRUE, the BF image is descpeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with ImgLimits = 50, it will search for spots in
#'   the central area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'   They refer to the ImageType part of the original image or txt file names.
#'
#' @return Location statistics under BF image modelling
#'
#' @keywords internal
getCsFAIL <-
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
    fai <-
      failurecase(
        img = origImg,
        pattern.search = 5,
        despeckle = despeckle,
        ImgLimits =
          ImgLimits,
        chip.type = chip.type,
        separator = separator,
        image.type = image.type
      )
    if (fai[[1]][1] == 0) {
      centerR[[1]] <- c(0, 0)
      centerG[[1]] <- c(0, 0)
      arR <- c()
      arG <- c()
      warn <- fai$warning
      out <- list(
        sample = "---",
        centerR = c(0, 0),
        centerG = c(0, 0),
        arR = c(),
        arG = c(),
        warn = c()
      )
    }
    
    if (fai[[1]][1] > 0) {
      centerR[[1]] <- matrix(centerR[[1]], ncol = 2)
      centerG[[1]] <- matrix(centerG[[1]], ncol = 2)
      
      res <- matrix(0, 1, 6)
      for (i in 1:nrow(centerR[[1]])) {
        for (j in 1:nrow(centerG[[1]])) {
          pixR <- 0
          pixG <- 0
          if (centerR[[1]][1] > 0 & nrow(centerR[[1]]) == 1) {
            pixR <- nrow(centerR$brightcoordinates[[1]])
          }
          if (centerR[[1]][1] > 0 & nrow(centerR[[1]]) > 1) {
            pixR <- nrow(centerR$brightcoordinates[[i]])
          }
          if (centerG[[1]][1] > 0 & nrow(centerG[[1]]) == 1) {
            pixG <- nrow(centerG$brightcoordinates[[1]])
          }
          if (centerG[[1]][1] > 0 & nrow(centerG[[1]]) > 1) {
            pixG <- nrow(centerG$brightcoordinates[[j]])
          }
          diR <-
            c(abs(centerR[[1]][i, 1] - fai[[1]][1]), abs(centerR[[1]][i, 2] - fai[[1]][2]))
          diG <-
            c(abs(centerG[[1]][j, 1] - fai[[1]][1]), abs(centerG[[1]][j, 2] - fai[[1]][2]))
          res <-
            matrix(rbind(res, c(
              i, j, max(diR), max(diG), pixR, pixG
            )), ncol = 6)
        }
      }
      res <- matrix(res[-1,], ncol = ncol(res))
      
      wR1 <- which(res[, 3] <= BFarea)
      wG1 <- which(res[, 4] <= BFarea)
      wR <- wR1
      wG <- wG1
      if (length(wR1) > 1) {
        wR <- wR1[which(res[wR1, 5] == max(res[wR1, 5]))]
      }
      if (length(wG1) > 1) {
        wG <- wG1[which(res[wG1, 6] == max(res[wG1, 6]))]
      }
      
      if (length(wR1) == 0 & length(wG1) == 0) {
        ss <-
          significantSignal(
            centerR = centerR,
            centerG = centerG,
            minDiff = minDiff,
            chaImgs = chaImgs
          )
        signifR <- ss$pvR
        signifG <- ss$pvG
        
        if (sum(signifR) == 1 & sum(signifG) == 0) {
          wh <- which(signifR == 1)
          centerR[[1]] <- matrix(centerR[[1]][wh,], ncol = 2)
          centerG[[1]] <- centerR[[1]]
          arR <- centerR$brightcoordinates[[wh]]
          arG <- centerR$brightcoordinates[[wh]]
          warn <- giveWarning(2)
        }
        if (sum(signifR) == 0 & sum(signifG) == 1) {
          wh <- which(signifG == 1)
          centerG[[1]] <- matrix(centerG[[1]][wh,], ncol = 2)
          centerR[[1]] <- centerG[[1]]
          arR <- centerG$brightcoordinates[[wh]]
          arG <- centerG$brightcoordinates[[wh]]
          warn <- giveWarning(2)
        }
        
        if (sum(signifR) == 0 &
            sum(signifG) == 0 |
            sum(signifR) > 1 |
            sum(signifG) > 1 |
            sum(signifR) >= 1 & sum(signifG) >= 1) {
          centerR[[1]] <- fai[[1]]
          centerG[[1]] <- fai[[1]]
          arR <- c()
          arG <- c()
          warn <- fai$warning
        }
        
      }
      
      if (length(wR1) > 0 & length(wG1) > 0) {
        estim <-
          getCoordinates_stats(
            centerR = centerR,
            centerG = centerG,
            minDiff = minDiff,
            chaImgs = chaImgs,
            ll =
              res,
            ws = list(
              wR = wR,
              wG = wG,
              wR1 = wR1,
              wG1 = wG1
            ),
            estCenter = fai[[1]]
          )
        if (estim$center[1] > 0) {
          #        fai[[1]] <- estim$center
          centerR[[1]] <- matrix(estim$center, ncol = 2)
          centerG[[1]] <- matrix(estim$center, ncol = 2)
          arR <- estim$areaR
          arG <- estim$areaR
          warn <- giveWarning(1)
        }
        if (estim$center[1] == 0) {
          centerR[[1]] <- fai[[1]]
          centerG[[1]] <- fai[[1]]
          arR <- c()
          arG <- c()
          warn <- fai$warning
        }
      }
      
      if (length(wR1) > 0 & length(wG1) == 0) {
        #      fai[[1]] <- centerR[[1]][res[wR[1],1],]
        centerR[[1]] <-
          matrix(centerR[[1]][res[wR[1], 1],], ncol = 2)
        centerG[[1]] <- centerR[[1]]
        arR <- centerR$brightcoordinates[[res[wR[1], 1]]]
        arG <- centerR$brightcoordinates[[res[wR[1], 1]]]
        w <- unlist(strsplit(fai$warning, ":"))[1]
        warn <- giveWarning(2)
      }
      if (length(wR1) == 0 & length(wG1) > 0) {
        #      fai[[1]] <- centerG[[1]][res[wG[1],1],]
        centerG[[1]] <-
          matrix(centerG[[1]][res[wG[1], 1],], ncol = 2)
        centerR[[1]] <- centerG[[1]]
        arR <- centerG$brightcoordinates[[res[wG[1], 1]]]
        arG <- centerG$brightcoordinates[[res[wG[1], 1]]]
        w <- unlist(strsplit(fai$warning, ":"))[1]
        warn <- giveWarning(2)
      }
      out <- list(
        sample = "---",
        centerR = fai[[1]],
        centerG = fai[[1]],
        arR = c(),
        arG = c(),
        warn = fai$warning
      )
    }
    return(
      list(
        centers = centerR[[1]],
        areaR = arR,
        areaG = arG,
        warning = warn,
        fail = fai,
        outlier.estimates = out
      )
    )
  }

#' spotCoords
#'
#' It estimates the spot location statistics by fluorescence signal in each channel. Then, it integrates
#'   the channel-specific data into a single estimate
#'
#' @param centerR Data matrix. The location statistics of one channel.
#' @param centerG Data matrix. The location statistics of the other channel.
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param chaImgs List. A list of the red and green channel images (data matrices) of a sample
#' @param minDiff Float. the mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param despeckle Logical. If TRUE, the BF image is descpeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with ImgLimits = 50, it will search for spots in
#'   the central area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#' @param type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'   They refer to the ImageType part of the original image or txt file names.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels.
#'
#' @return Location statistics under fluorescence-based estimation
#'
#' @keywords internal
spotCoords <-
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
    centerR[[1]] <- matrix(centerR[[1]], ncol = 2)
    centerG[[1]] <- matrix(centerG[[1]], ncol = 2)
    ros <- c(nrow(centerR[[1]]), nrow(centerG[[1]]))
    
    #case of 0s
    if (centerR[[1]][1, 1] == 0 & centerG[[1]][1, 1] == 0) {
      res <-
        caseof0s(
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
    }
    
    #case of 1R and 0G
    if (max(ros) == 1 &
        centerR[[1]][1, 1] == 0 &
        centerG[[1]][1, 1] != 0 |
        max(ros) == 1 &
        centerG[[1]][1, 1] == 0 & centerR[[1]][1, 1] != 0) {
      res <-
        caseof1R0G(
          centerR,
          centerG,
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
    }
    
    #case of 2R and 0G or opposite
    if (max(ros) > 1 &
        centerR[[1]][1, 1] == 0 |
        max(ros) > 1 & centerG[[1]][1, 1] == 0) {
      res <-
        caseof2Rs0G(
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
          image.type = image.type,
          show.possible.contamination = show.possible.contamination
        )
    }
    
    # case of 1R and 1G
    if (max(ros) == 1 & centerR[[1]][1, 1] != 0 &
        centerG[[1]][1, 1] != 0) {
      res <-
        caseof1R1G(
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
          image.type = image.type,
          show.possible.contamination = show.possible.contamination
        )
    }
    
    #case of 1R and 2G
    if (ros[1] == 1 & ros[2] > 1 & centerR[[1]][1, 1] != 0) {
      res <-
        caseof1R2G(
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
          image.type = image.type,
          show.possible.contamination = show.possible.contamination
        )
    }
    
    #case of 2R and 1G
    if (ros[1] > 1 & ros[2] == 1 & centerG[[1]][1, 1] != 0) {
      res <-
        caseof2R1G(
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
          image.type = image.type,
          show.possible.contamination = show.possible.contamination
        )
    }
    
    #case of many R/G
    if (ros[1] > 1 & ros[2] > 1) {
      res <-
        caseof2R2G(
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
    }
    
    return(
      list(
        center = res$centers,
        areaCH1 = res$areaR,
        areaCH2 = res$areaG,
        other.spots = res$possible.contamination,
        warning =
          res$warning,
        failure = res$fail,
        outlier.estimates = res$outlier.estimates
      )
    )
  }

#' extractBFarea
#'
#' It estimates the spot or capture site location by BF image modelling.
#'
#' @param cc Data matrix. An estimated spot center.
#' @param img Data matrix. The matrix of image data from a channel.
#' @param area Data matrix. The bright (spot) coordinates around the spot center.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#'
#' @return The area identified by BF image modelling
#'
#' @keywords internal
extractBFArea <- function(cc, img, area, BFarea) {
  b1 <- (cc[1] - floor(BFarea / 2)):(cc[1] + floor(BFarea / 2))
  b2 <- (cc[2] - floor(BFarea / 2)):(cc[2] + floor(BFarea / 2))
  newarea <- expand.grid(b1, b2)
  newarea[, 1] <- area[[1]][newarea[, 1]]
  newarea[, 2] <- area[[2]][newarea[, 2]]
  return(list(newarea))
}


#' measureF
#'
#' It estimates the foreground signal for an identified spot or for a predefined area within
#'   the capture site.
#'
#' @param img Data matrix. The matrix of image data from a channel.
#' @param area Data matrix. The bright (spot) coordinates around the spot center.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#'
#' @importFrom raster raster focal
#'
#' @return The foreground (spot) signal estimates
#'
#' @keywords internal
measureF <- function(img, area, BFarea) {
  if (length(BFarea) == 0) {
    int <- rep(0, nrow(area))
    for (i in 1:nrow(area)) {
      int[i] <- img[area[i, 1], area[i, 2]]
    }
    meanint <- 2 ^ mean(log(int + 1, 2))
    totint <- sum(int)
    NewArea <- area
  }
  if (length(BFarea) > 0) {
    if (BFarea %% 2 == 0) {
      BFarea <- BFarea + 1
    }
    a1 <- min(area[, 1]):max(area[, 1])
    a2 <- min(area[, 2]):max(area[, 2])
    x <- img[a1, a2]
    xx <- log(x + 1, 2)
    xr <- raster(xx)
    int1 <-
      2 ^ raster::as.matrix(focal(
        xr,
        matrix(1, BFarea, BFarea),
        mean,
        pad = TRUE,
        padvalue = 0
      ))
    int1[is.na(int1)] <- 0
    meanint <- as.numeric(quantile(int1[int1 > 0], 1))
    diff1 <- int1 - meanint
    w1 <-
      matrix(which(abs(diff1) == min(abs(diff1)), arr.ind = TRUE), ncol = 2)
    NewArea <-
      as.list(apply(
        w1,
        1,
        extractBFArea,
        img = x,
        area = list(a1, a2),
        BFarea = BFarea
      ))
    totint <- meanint * (BFarea ^ 2)
    int <- int1[int1 > 0]
  }
  return(list(c(meanint, totint), NewArea, intensity = int))
}

#' measureB
#'
#' It estimates the background signal for an image.
#'
#' @param img Data matrix. The matrix of image data from a channel.
#' @param area Data matrix. The bright (spot) coordinates around the spot center.
#' @param iter Integer. A number of iterations (pseudo-spots) to be summarized.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated.
#'
#' @importFrom raster raster focal
#'
#' @return The background image estimates
#'
#' @keywords internal
measureB <- function(img, area, iter, BFarea) {
  if (length(BFarea) > 0) {
    area <-
      expand.grid((floor(mean(area[, 1])) - floor(BFarea / 2)):(floor(mean(area[, 1])) +
                                                                  floor(BFarea / 2)),
                  (floor(mean(area[, 2])) - floor(BFarea / 2)):(floor(mean(area[, 2])) +
                                                                  floor(BFarea / 2)))
  }
  rowdif <-
    c(50, min((min(area[, 1]) - 1), (nrow(img) - max(area[, 1]) - 50)))
  coldif <-
    c(50, min((min(area[, 2]) - 1), (ncol(img) - max(area[, 2]) - 50)))
  s <-
    matrix(cbind(
      sample(rowdif[1]:rowdif[2], iter, replace = TRUE),
      sample(coldif[1]:coldif[2], iter, replace =
               TRUE)
    ), ncol = 2)
  s[, 1] <- sample(c(-1, 1), nrow(s), replace = TRUE) * s[, 1]
  s[, 2] <- sample(c(-1, 1), nrow(s), replace = TRUE) * s[, 2]
  s <- unique(s)
  bint <- matrix(0, nrow(s), 2)
  for (b in 1:nrow(s)) {
    m1 <- matrix(area[, 1] + s[b, 1], ncol = 1)
    m2 <- matrix(area[, 2] + s[b, 2], ncol = 1)
    if (length(m1[m1 < 1]) > 0 | length(m1[m1 > nrow(img)]) > 0) {
      m1 <- matrix(area[, 1] - s[b, 1], ncol = 1)
    }
    if (length(m2[m2 < 1]) > 0 | length(m2[m2 > ncol(img)]) > 0) {
      m2 <- matrix(area[, 2] - s[b, 2], ncol = 1)
    }
    barea <- matrix(cbind(m1, m2), ncol = 2)
    int <- rep(0, nrow(barea))
    for (i in 1:nrow(barea)) {
      int[i] <- img[barea[i, 1], barea[i, 2]]
    }
    bint[b, 1] <- 2 ^ mean(log(int + 1, 2))
    bint[b, 2] <- sum(int)
  }
  res <- apply(bint, 2, mean)
  return(res)
}


#' findPattern
#'
#' A helper to find the characteristic straight lines of the BF image .
#'
#' @param imgline Integer. A row or a column number of the BF image data matrix.
#' @param bpattern Integer. A user-defined cutoff that specifies whether a given row
#'   or column contains a characteristic straight line.
#'
#' @return An estimate for the existence of a characteristic line
#'
#' @keywords internal
findPattern <- function(imgline, bpattern) {
  m <- rollmean(imgline, k = bpattern)
  length(m[m == 0])
}


#' despecklefun
#'
#' It despeckles the BF image data matrix (similar to the despeckle function of ImageJ)
#'
#' @param img Data matrix. The BF image data matrix.
#' @param pix Integer. A cutoff specifing the area to be examined for speckles.
#' @param thresh Integer. A cutoff to perform the despeckle function. If pixel signal > median object signal + thresh,
#'   the object is a speckle and the median object signal is returned.
#'
#' @return A despeckled image
#'
#' @import stats
#'
#' @keywords internal
despecklefun <- function(img, pix, thresh) {
  for (i in 1:nrow(img)) {
    if (i == 1) {
      ro <- 1:(i + pix)
    }
    if (i == nrow(img)) {
      ro <- (i - pix):nrow(img)
    }
    if (i > 1 & i < nrow(img)) {
      ro <- (i - pix):(i + pix)
    }
    for (j in 1:ncol(img)) {
      if (j == 1) {
        co <- 1:(j + pix)
      }
      if (j == ncol(img)) {
        co <- (j - pix):ncol(img)
      }
      if (j > 1 & j < ncol(img)) {
        co <- (j - pix):(j + pix)
      }
      m <- median(img[ro, co])
      if (img[i, j] > (m + thresh)) {
        img[i, j] <- m
      }
    }
  }
  return(img)
}


#' grouplines
#'
#' It reconstructs the BF image characteristic lines.
#'
#' @param data Integer. A row or a column of the BF image data matrix.
#'
#' @return An estimate for the existence of a characteristic line
#'
#' @keywords internal
grouplines <- function(data) {
  data[data > 0] <- 1
  ww <- which(data == 1)
  tt <- split(ww, cumsum(c(1, diff(ww) != 1)))
  tt <- as.numeric(unlist(lapply(tt, length)))
  return(tt)
}


#' straightRowLines
#'
#' It identifies the horizontal BF image characteristic lines.
#'
#' @param colData List. The data to be analyzed.
#' @param pattern.search Integer. A cutoff to find horizontal and vertical lines on the chip
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central
#'   area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#'
#' @import utils stats
#' @importFrom zoo rollmean
#'
#' @return The estimated horizontal BF characteristic lines
#'
#' @keywords internal
straightRowLines <- function(colData, pattern.search, ImgLimits) {
  fac <- ncol(colData$spotMat) / 512
  line.lims <- fac * c(34, 51)
  
  if (colData$Rows[1] > 0) {
    xx <- colData$spotMat[colData$Rows,]
  } else {
    rr <- ImgLimits:(nrow(colData$spotMat) - ImgLimits)
    xx <- colData$spotMat[rr,]
  }
  f1 <- as.numeric(apply(xx, 1, findPattern, pattern.search))
  ff <- c(0, round(rollmean(f1, 3) * 3))
  ff <- sort.list(ff, decreasing = TRUE)[1:10]
  if (colData$Rows[1] > 0) {
    ff <- colData$Rows[ff]
  } else {
    ff <- rr[ff]
  }
  ind <- c()
  ro <- c()
  if (length(ff) > 0) {
    cmbn <- matrix(combn(ff, 2), nrow = 2)
    di <- abs(apply(cmbn, 2, diff))
    ww <- which(di > line.lims[1] & di < line.lims[2])
    if (length(ww) > 0) {
      cmbn <-
        matrix(cmbn[, which(di > line.lims[1] &
                              di < line.lims[2])], nrow = 2)
      cmbn <- apply(cmbn, 2, sort)
      ro <- round(apply(cmbn, 1, median), 0)
      ro <- ro[1]:ro[2]
      cc1 <-
        matrix(colData$spotMat[unique(cmbn[1,]),], nrow = length(unique(cmbn[1,])))
      cc2 <-
        matrix(colData$spotMat[unique(cmbn[2,]),], nrow = length(unique(cmbn[2,])))
      cc <- apply(matrix(rbind(cc1, cc2), ncol = ncol(cc1)), 2, sum)
      res <- cumsum(cc) / sum(cc)
      ind <- 1
      if (which(res > 0.5)[1] > (0.5 * ncol(colData$spotMat))) {
        ind <- 0
      }
    }
  }
  return(list(ind, ro))
}

#' straightColLines
#'
#' It identifies the vertical BF image characteristic lines.
#'
#' @param img Data matrix. The BF image data matrix.
#' @param pattern.search Integer. A cutoff to find horizontal and vertical lines on the chip.
#' @param cut Integer. A normalized signal cutoff above which a pixel is considered as foreground.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with cutSides = 50, it will search for spots in the central
#'   area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#'
#' @import stats
#' @importFrom zoo rollmean rollmedian
#' @importFrom raster raster focal
#'
#' @return Estimated vertical BF characteristic lines.
#'
#' @keywords internal
straightColLines <-
  function(img,
           pattern.search,
           cut = seq(0.3, 1.5, 0.02),
           ImgLimits,
           chip.type) {
    ffCols <- as.list(rep(0, 2))
    ffRows <- 0
    anothercut <- seq(0.85, 0.65,-0.05)
    j <- 1
    while (length(ffCols[[1]]) < 2 &
           j <= length(cut) | ffRows[1] == 0 & j <= length(cut)) {
      xx <- round(img / (cut[j] * max(img)), 0)
      #    xx <- round(img / cut[j],0)
      f1 <- as.numeric(apply(xx, 2, findPattern, pattern.search))
      f11 <- length(f1[f1 == 0])
      k <- 1
      while (length(ffCols[[1]]) < 2 & k <= 5 | ffRows[1] == 0 &
             k <= 5) {
        if (f11 / ncol(img) < anothercut[k]) {
          ff <- as.numeric(apply(xx, 2, findPattern, pattern.search))
          ff <- c(0, round(rollmean(ff, 3) * 3))
          ff <- sort.list(ff, decreasing = TRUE)[1:15]
          keep <- rep(0, length(ff))
          ro <- matrix(0, length(ff), 2)
          for (i in 1:length(ff)) {
            mcols <- max(1, (ff[i] - 2)):min(ncol(img), (ff[i] + 2))
            if (length(mcols) < 5 & mcols[1] == 1) {
              mcols <- c(mcols, 5)
            }
            if (length(mcols) < 5 &
                mcols[length(mcols)] == ncol(img)) {
              mcols <- c((ncol(img) - 4), mcols)
            }
            fr <- raster(xx[, mcols])
            xx1 <-
              raster::as.matrix(focal(
                fr,
                matrix(1, 5, 5),
                mean,
                pad = TRUE,
                padvalue = 0
              ))
            xx1[is.na(xx1)] <- 0
            xx1[xx1 == 1] <- 0
            rr1 <- rollmedian(xx1[1:nrow(xx1), 3], 51)
            rr2 <- rollmedian(xx1[nrow(xx1):1, 3], 51)
            t1 <- grouplines(data = rr1)
            if (max(t1) >= 70) {
              keep[i] <- ff[i]
              rr1 <- range(which(rr1 > 0))
              rr2 <- nrow(xx) - range(which(rr2 > 0))
              ro[i,] <- c(rr1[1], rr2[1])
            }
          }
          ff <- ff[which(ro[, 1] > 0)]
          ro <- ro[which(ro[, 1] > 0),]
          if (length(ff) > 1) {
            ffCols <-
              highlight.cols(ff, chip.type = chip.type, fac = ncol(img) / 512)
            ro1 <- matrix(0, 1, 2)
            if (length(ffCols[[1]]) > 0 & ffCols[[1]][1] > 0) {
              mm <- match(ffCols[[1]], ff)
              ro1 <- matrix(ro[mm,], ncol = 2)
            }
            ro2 <- matrix(0, 1, 2)
            if (ffCols[[2]][1] > 0) {
              mm <- match(ffCols[[2]], ff)
              ro2 <- matrix(ro[mm,], ncol = 2)
            }
            ffRows <- matrix(rbind(ro1, ro2), ncol = 2)
            if (sum(ro[, 1]) == 0) {
              ffRows <- 0
            } else {
              ffRows <- round(apply(ffRows, 2, median), 0)
              ffRows <- ffRows[1]:ffRows[2]
            }
          }
        }
        k <- k + 1
      }
      j <- j + 1
    }
    if (ffCols[[1]][1] > ImgLimits &
        ffRows[1] > ImgLimits |
        ffCols[[1]][length(ffCols[[1]])] < (ncol(img) - ImgLimits) &
        ffRows[length(ffRows)] < (nrow(img) - ImgLimits)) {
      res <-
        max(1, (min(ffCols[[1]]) - 35)):min(ncol(img), (max(ffCols[[1]]) + 35))
    } else {
      res <- 0
      ffRows <- 0
      ffCols[[1]] <- 0
      ffCols[[2]] <- 0
      xx <- matrix(0, nrow(xx), ncol(xx))
    }
    return(list(
      Rows = ffRows,
      Cols = res,
      mainStrLines = ffCols[[1]],
      farLines = ffCols[[2]],
      spotMat =
        xx
    ))
  }

#' highlight.cols
#'
#' A helper to identify the vertical BF image characteristic lines.
#'
#' @param data Numeric. A column of the BF image data matrix.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param fac Float. An internally estimated factor that adjusts for different chip types.
#'
#' @return Estimated vertical BF characteristic lines.
#'
#' @import utils
#'
#' @keywords internal
highlight.cols <- function(data, chip.type, fac) {
  if (chip.type == "medium/large") {
    line.lims <- fac * c(15, 35, 45, 219, 276)
  } else {
    line.lims <- fac * c(25, 35, 190, 231)
  }
  res <- matrix(0, length(data), length(data))
  indices <- combn(1:length(data), 2)
  for (i in 1:ncol(indices)) {
    d <- abs(data[indices[1, i]] - data[indices[2, i]])
    
    if (chip.type == "medium/large") {
      if (d > line.lims[1] &
          d < line.lims[2] |
          d > line.lims[2] &
          d < line.lims[3] | d > line.lims[4] & d < line.lims[5]) {
        res[indices[1, i], indices[2, i]] <- d
      } else {
        res[indices[1, i], indices[2, i]] <- 0
      }
    }
    
    if (chip.type == "small") {
      if (d > line.lims[1] &
          d < line.lims[2] | d > line.lims[3] & d < line.lims[4]) {
        res[indices[1, i], indices[2, i]] <- d
      } else {
        res[indices[1, i], indices[2, i]] <- 0
      }
    }
  }
  ll <- rep(0, nrow(res))
  for (i in 1:nrow(res)) {
    ll[i] <- length(res[i,][res[i,] > 0])
  }
  w1 <- which(ll == max(ll))[1]
  w2 <- which(res[w1,] > 0, arr.ind = TRUE)
  data <- data[unique(c(w1, w2))]
  if (length(data) > 1) {
    res <- matrix(0, length(data), length(data))
    indices <- combn(1:length(data), 2)
    for (i in 1:ncol(indices)) {
      d <- abs(data[indices[1, i]] - data[indices[2, i]])
      
      if (chip.type == "medium/large") {
        if (d > line.lims[1] &
            d < line.lims[2] |
            d > line.lims[2] & d < line.lims[3]) {
          res[indices[1, i], indices[2, i]] <- d
        } else {
          res[indices[1, i], indices[2, i]] <- 0
        }
      }
      
      if (chip.type == "small") {
        if (d > line.lims[1] & d < line.lims[2]) {
          res[indices[1, i], indices[2, i]] <- d
        } else {
          res[indices[1, i], indices[2, i]] <- 0
        }
      }
      
    }
    result <- data[unique(c(which(res > 0, arr.ind = TRUE)))]
  } else {
    result <- data
  }
  ss <- setdiff(data, result)
  if (length(ss) == 0) {
    ss <- 0
  }
  return(list(result, ss))
}


#' failurecase
#'
#' Another helper to re-estimate the spot location or find the capture site location
#'   (BF image modelling).
#'
#' @param origImg Data matrix. The original BF image to be read and processed.
#' @param pattern.search Integer. A cutoff to find horizontal and vertical lines on the chip.
#' @param despeckle Logical. If TRUE, the BF image is descpeckled.
#' @param ImgLimits Integer. It instructs the algorithm to find spots in a certain central image area.
#'   For example, for a 512 x 512 image with ImgLimits = 50, it will search for spots in
#'   the central area [ImgLimits:(512-ImgLimits),ImgLimits:(512-ImgLimits)] of the image matrix.
#' @param type Character string. It specifies the type of Fluidigm chip to be analyzed.
#' @param separator Character string. Removes the Bright Field ("BF") and channel indicators from the
#'   image file names.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#'
#' @importFrom EBImage Image bwlabel fillHull paintObjects
#' @importFrom zoo rollmean
#'
#' @return Location statistics and characteristic lines of BF image modelling
#'
#' @keywords internal
failurecase <-
  function(img,
           pattern.search,
           despeckle,
           ImgLimits,
           chip.type,
           separator,
           image.type) {
    images <- readOriImg(
      imgName = img,
      despeckle = despeckle,
      pix = 1,
      thresh = 50,
      separator = separator,
      image.type = image.type
    )
    parZoom <-
      zoomInBF(
        img = images[[2]],
        pattern.search = pattern.search,
        ImgLimits = ImgLimits,
        chip.type = chip.type
      )
    
    if (parZoom$reducedRows[1] > 0 & parZoom$reducedCols[1] > 0) {
      tt1 <- images[[2]][parZoom$reducedRows, parZoom$reducedCols]
    } else {
      rr <- ImgLimits:(nrow(images[[1]]) - ImgLimits)
      cc <- ImgLimits:(ncol(images[[1]]) - ImgLimits)
      tt1 <- images[[2]][rr, cc]
    }
    
    qq <- as.numeric(quantile(tt1, 0.95))
    tt2 <- tt1 - qq
    img1 <- Image(tt2, c(nrow(tt2), ncol(tt2)))
    img2 <- img1 > 0
    imglabel = bwlabel(img2)
    imglabel = fillHull(imglabel)
    spots <- paintObjects(img2, imglabel)@.Data
    spots[spots > 0] <- 1
    
    res.lines <- c(0, 0)
    if (parZoom$reducedRows[1] > 0 & parZoom$strLines[1] > 0) {
      res.lines <-
        SpotbyStrLines(
          binImg = parZoom$binImgDespeckled,
          pattern.search = pattern.search,
          stats =
            parZoom
        )
    }
    
    res.img <- getSpot(img = (tt1 * spots), rad = 3)
    if (parZoom$reducedRows[1] > 0 &
        parZoom$strLines[1] > 0 & res.img[1] > 0) {
      res.img[1] <- parZoom$reducedRows[res.img[1]]
      res.img[2] <- parZoom$reducedCols[res.img[2]]
      bf <- min(parZoom$strLines):max(parZoom$strLines)
      fi <-
        min(findInterval(res.img[2], bf), (length(bf) - findInterval(res.img[2], bf)))
      spotsNew <- matrix(0, nrow(images[[1]]), ncol(images[[1]]))
      spotsNew[parZoom$reducedRows, parZoom$reducedCols] <- spots
      if (sum(fi) > 7 |
          res.img[1] < min(parZoom$reducedRows) &
          res.img[1] > max(parZoom$reducedRows)) {
        res.img <- rep(0, 2)
        spotsNew <- matrix(0, nrow(images[[1]]), ncol(images[[1]]))
      }
    } else {
      res.img[1] <- rr[res.img[1]]
      res.img[2] <- cc[res.img[2]]
      spotsNew <- matrix(0, nrow(images[[1]]), ncol(images[[1]]))
      spotsNew[rr, cc] <- spots
    }
    
    di <-
      c(abs(res.lines[1] - res.img[1]), abs(res.lines[2] - res.img[2]))
    warn <- giveWarning(3)
    
    res <- c(0, 0)
    if (res.lines[1] == 0 & res.img[1] > 0) {
      res <- res.img
    }
    if (res.lines[1] > 0) {
      res <- res.lines
    }
    return(list(
      res,
      spotsNew,
      warning = warn,
      failParameters = parZoom
    ))
  }

#' filled.contour3
#'
#' It generates a contour plot of a channel specific spot. It is a modification by Ian Taylor of the filled.contour()
#'   to remove the key and facilitate overplotting with contour(). It has been further modified by Carey McGilliard and
#'   Bridget Ferris to allow multiple plots on one page http://wiki.cbr.washington.edu/qerm/sites/qerm/images/1/16/Filled.contour3.R
#'   We have added some extra parameters to adapt the function to our application.
#'
#' @param x,y,z Numeric vectors. Some plot coordinates.
#' @param ... Other parameters of the function.
#'
#' @return A plotted spot
#'
#' @import grDevices graphics
#' @keywords internal
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)),
            z,
            xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE),
            zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels),
            nlevels = 20,
            color.palette = cm.colors,
            col = color.palette(length(levels) - 1),
            joinedPlots,
            plot.title,
            plot.axes,
            key.title,
            key.axes,
            asp = NA,
            xaxs = "i",
            yaxs = "i",
            las = 1,
            axes = TRUE,
            frame.plot = axes,
            mar,
            ...)
  {
    if (length(joinedPlots) > 0) {
      levels = pretty(joinedPlots, nlevels)
      col = color.palette(length(levels) - 1)
    }
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else
        stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim,
                ylim,
                "",
                xaxs = xaxs,
                yaxs = yaxs,
                asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "",
              xlab = "",
              ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else
      plot.axes
    if (frame.plot)
      graphics::box()
    if (missing(plot.title))
      title(...)
    else
      plot.title
    invisible()
  }

#' multiplot
#'
#' Multiple plot function.
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#'  cols:   Number of columns in layout.
#'  layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#'  then plot 1 will go in the upper left, 2 will go in the upper right, and
#'  3 will go all the way across the bottom.
#'
#' @return ggplot2 multiplot
#'
#' @importFrom grid grid.newpage grid.layout pushViewport viewport
#'
#' @keywords internal
multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
    }
    
    if (numPlots == 1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <-
          as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]],
              vp = viewport(
                layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col
              ))
      }
    }
  }
