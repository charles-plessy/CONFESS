# Estimation of fluorescence signal from the C1 Chip

#' readFiles
#'
#' Reads the image data that are going to be analyzed. It converts the images into txt files. The images should
#' be in .C01 (high resolution) or .BMP, or .JPG or .PNG format. The file names should be of the form:
#'
#'                   "RunID(separator1)WellID(separator2)ImageType.ImageFormat
#'
#' For example in "1772-062-248_A01@BF.C01", RunID = 1772-062-248", separator1 = _, WellID = A01, separator2 = @
#'   ImageType = BF, ImageFormat = C01. The function expects to see both Bright Field and channel images. It will
#'   store them in different directories. It will return a list of the respective .txt file names. Note that separator1
#'   and separator2 CAN BE the same character (e.g. "_").
#'
#' If the images have been already converted, then the txt files should be stored in the above form with ImageFormat = txt.
#'
#' readFiles() will take the minimum overlapping sets. Converted images not present in any
#' of the channels or the Bright Field list will be reported and discarded.
#'
#' @param iDirectory Character string. The directory where all images are stored. The images should be in
#'   the same format. Available choices are: C01, BMP, JPEG and PNG. The function recognizes the format automatically.
#'   If ommitted, the function assumes that the txt data already exist at the predefined folders.
#' @param BFdirectory Character string. The directory to store the .txt converted Bright Field images.
#' @param CHdirectory Character string. The directory to store the .txt converted channel (e.g. Red/Green) images
#' @param separator Character string. This is the <<separator2>> parameter that removes the Bright Field ("BF") and
#' channel indicators (IDs) from the image file names. Default is "_".
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study. They refer to
#'   the ImageType part of the original image or txt file names. Default is c("BF","Red","Green").
#' @param bits Numeric. The image bits. It is used to unnormalize the C01 signals from readCellomics(). It does
#'   not affect the signals of other image types. Default is 2^16.
#'
#' @import readbitmap
#'
#' @return A list with the followign components:
#'   BF: the files names of the Bright Field converted data matrices.
#'   CH1: the files names of the converted data matrices of one channel.
#'   CH2: the files names of the converted data matrices of the other channel.
#'   separator: the separator being used.
#'   image.type: the image type IDs.
#'   dateIndex: a date index to be used in saving the output files.
#'
#' @export
#'
#' @examples
#' library(fucci)
#'
#' ### set your directories
#' basedir<-"~/"
#' fucci_path<-system.file("extdata",package="fucci")
#' 
#' ## to read txt files
#' files<-readFiles(iDirectory=NULL,
#'                  BFdirectory=paste(fucci_path,"/BF",sep=""),
#'                  CHdirectory=paste(fucci_path,"/CH",sep=""),
#'                  separator = "_",image.type = c("BF","Green","Red"),
#'                  bits=2^16)
#'                  
#' ## to convert from BMP/JPEG images
#' #write_dir<-"~/converted_images/"
#' #files<-readFiles(iDirectory=fucci_path,
#' #                 BFdirectory=paste(write_dir,"/BF",sep=""),
#' #                 CHdirectory=paste(write_dir,"/CH",sep=""),
#' #                 separator = "_",image.type = c("BF","Green","Red"),
#' #                 bits=2^16)                  
readFiles <-function(iDirectory,BFdirectory,CHdirectory,separator = "_",
                     image.type = c("BF", "Red", "Green"),
                     bits = 2 ^ 16) {
    if (length(iDirectory) > 0) {
      if (as.numeric(file.access(iDirectory)) < 0) {
        stop("The iDirectory cannot be found in the system")
      }
    }
    if (as.numeric(file.access(BFdirectory)) < 0) {
      dir.create(BFdirectory)
      message("The BFdirectory could not be found in the system and has been automatically created")
    }
    if (as.numeric(file.access(CHdirectory)) < 0) {
      dir.create(CHdirectory)
      message("The CHdirectory could not be found in the system and has been automatically created")
    }
    
    datin <- paste(unlist(strsplit(date(), " ")), collapse = "")
    
    if (length(iDirectory) > 0) {
      # read the image data (if they exist)
      all.ls <- as.list(list.files(iDirectory, full.names = TRUE))
      test <- unlist(strsplit(all.ls[[1]], ".", fixed = TRUE))
      image.format <- test[length(test)]
      
      if (image.format != "BMP" & image.format != tolower("BMP") &
          image.format != "PNG" &
          image.format != tolower("PNG") &
          image.format != "JPEG" &
          image.format != tolower("JPEG") &
          image.format != "JPG" &
          image.format != tolower("JPG")) {
        stop("This image format is not supported")
      }
      
      txtdata <- lapply(all.ls, read.bitmap)
      # if (image.format == "C01" | image.format == tolower("C01")) {
      #   txtdata <- lapply(all.ls, readCellomics)
      #   txtdata <- lapply(txtdata, unnormalizeC01, bits = bits)
      # } else {
      #   txtdata <- lapply(all.ls, read.bitmap)
      # }
      names(txtdata) <- all.ls
      # store the image data
      for (i in 1:length(txtdata)) {
        u <- unlist(strsplit(names(txtdata)[i], "/"))
        #print(txtdata[i])
        u <- unlist(strsplit(u[length(u)], ".", fixed = TRUE))
        #print(u[1])
        if (length(grep(paste(separator, image.type[1], sep = ""), u[1])) >
            0) {
          write(
            t(txtdata[[i]]),
            paste(BFdirectory, "/", u[1], ".txt", sep = ""),
            ncolumns = ncol(txtdata[[1]]),
            sep = "\t"
          )
        } else {
          write(
            t(txtdata[[i]]),
            paste(CHdirectory, "/", u[1], ".txt", sep = ""),
            ncolumns = ncol(txtdata[[1]]),
            sep = "\t"
          )
        }
        message(paste(
          "Storing the txt data of image ",
          i,
          " of ",
          length(txtdata),
          sep = ""
        ))
      }
    }
    
    # keep the triplets
    lsBF <- list.files(BFdirectory, pattern = image.type[1])
    lsR <- list.files(CHdirectory, pattern = image.type[2])
    lsG <- list.files(CHdirectory, pattern = image.type[3])
    
    if (length(lsBF) == 0 | length(lsR) == 0 | length(lsG) == 0) {
      stop(
        "There are no txt data available. Convert the images of iDirectory or specify the correct txt folders and insert the appropriate separator and image.type!"
      )
    }
    
    l1 <-
      unlist(strsplit(lsBF, paste(separator, image.type[1], sep = "")))[seq(1, 2 *
                                                                              length(lsBF), by = 2)]
    l2 <-
      unlist(strsplit(lsR, paste(separator, image.type[2], sep = "")))[seq(1, 2 *
                                                                             length(lsR), by = 2)]
    l3 <-
      unlist(strsplit(lsG, paste(separator, image.type[3], sep = "")))[seq(1, 2 *
                                                                             length(lsG), by = 2)]
    
    complete_sets <- Reduce(intersect, list(l1, l2, l3))
    full_sets <- Reduce(union, list(l1, l2, l3))
    missed1 <- full_sets[!(full_sets %in% l1)]
    missed2 <- full_sets[!(full_sets %in% l2)]
    missed3 <- full_sets[!(full_sets %in% l3)]
    
    ## Report missing files
    if (length(missed1) > 0) {
      message(
        paste(
          "Warning: Some ",
          image.type[1],
          " data are missing! Verify the separator and image.type parameters.",
          sep = ""
        )
      )
      for (i in missed1)
        print(i)
    }
    if (length(missed2) > 0) {
      message(
        paste(
          "Warning: Some ",
          image.type[2],
          " data are missing! Verify the separator and image.type parameters.",
          sep = ""
        )
      )
      for (i in missed2)
        print(i)
    }
    if (length(missed3) > 0) {
      message(
        paste(
          "Warning: Some ",
          image.type[3],
          " data are missing! Verify the separator and image.type parameters.",
          sep = ""
        )
      )
      for (i in missed3)
        print(i)
    }
    
    return(
      list(
        BF = file.path(BFdirectory, lsBF[l1 %in% complete_sets]),
        CH1 = file.path(CHdirectory, lsR[l2 %in% complete_sets]),
        CH2 = file.path(CHdirectory, lsG[l3 %in% complete_sets]),
        separator = separator,
        image.type = image.type,
        dateIndex = datin
      )
    )
  }

#' defineLocClusters
#'
#' It performs quality check on the estimated location of spotEstimator() in order to flag possible outliers.
#' The flagging is done both visually and statistically using the Grubbs test.
#'
#' The outlier locations will be re-estimated by BF image modelling or adjusted as the 2-dimensional median of all
#'   non-outlying locations.
#'
#' @param LocData The table of the location estimates obtained by spotEstimator().
#' @param dims Numeric vector. The dimensions of the image data. Default is rep(512,2).
#' @param out.method Character string. The method by which to flag outliers: "interactive.clustering" or "interactive.manual"
#'   or "manual". Default is "interactive.clustering".
#'   The interactive options work through interactive plots: "interactive.clustering" enables
#'   the user to highight the outliers via co-centric circles in the plot while "interactive.manual"
#'   asks the user to click on the plot to highlight the outliers (to confirm & finalize the picks
#'   in each plot the user has to select the "stop" command (Windows) or press the right click in
#'   Linux/Mac. Note that 'interactive.clustering' works when one has more then or equal to 15 samples IN EACH CATEGORY
#'   (Run/Well combination).
#'   The "manual" option simply gives back the original table of location estimates with the last column
#'   being a series of "confidence". The outliers should be manually annotated by inserting "outlier" in the
#'   appropriate rows of the last column.
#' @param subset List. It allows the user to run the algorithm for a subset of data (run ids and wells).
#'   Default c() using all data. Otherwise put the run IDs and the wells (left and/or right) in a list,
#'   e.g. list(c("1772-115-xxx","1772-115-yyy"),"left").
#' @param separator Character string. It refers to <<separator1>> parameter described in readFiles() that separates the
#'   run ID from the Well ID in the original image (converted) file names. Default is "_".
#' @param savePlot Character string. Directory to store the plots if out.method = manual. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen. Default is the current working directory, getwd().
#'
#' @return A list of components summarizing the location estimates and their quality control statistics:
#'   Results: The table of the location estimates from spotEstimator() with an extra "QCgroup" labelled column that flags the samples
#'     either by "confidence" or by "outlier" (the locations that have been selected as outliers from the interactive plots). If
#'     out.method = "manual" the column includes a series of "confidence" entries. The outliers should be manually labelled.
#'   BFdata: the outlier estimates of spotEstimator(). They are kept here for processing in the second spotEstimator() step. See spotEstimator()
#'     for more details.
#'   Processed.Files: the samples that have been processed by spotEstimator(). Also kept from the first spotEstimator() step. They
#'     will be processed in the second spotEstimator() step.
#'   Outlier.indices: a vector of outlier sample indices. They are generated from the flagging of the outliers via interactive plots.
#'     They have to be manually specified if out.method = "manual".
#'   Medians: the 2-dimenional medians by run ID and wellID sets.
#'   Wellsets: a matrix showing the directionality of the well IDs.
#'   BFarea: the size of the pseudospot.
#'   image.type: the image type IDs.
#'   dateIndex: a date index to be used in saving the output files.
#'
#' @import stats graphics grDevices
#' @importFrom plotrix draw.circle
#'
#' @export
#'
#' @examples
#' library(fucci)
#' ### set your directories
#' basedir<-"~/"
#' fucci_path<-system.file("extdata",package="fucci")
#' files<-readFiles(iDirectory=NULL,
#'                  BFdirectory=paste(fucci_path,"/BF",sep=""),
#'                  CHdirectory=paste(fucci_path,"/CH",sep=""),
#'                  separator = "_",image.type = c("BF","Green","Red"),
#'                  bits=2^16)
#'
#' #this example is run using out.method="manual" (not interactive) 
#' clu <- defineLocClusters(LocData=estimates,out.method="manual",savePlot="screen")
defineLocClusters <- function(LocData,
                              dims = rep(512, 2),
                              out.method = "interactive.clustering",
                              subset = c(),
                              separator = "_",
                              savePlot = "screen") {
  if (savePlot != "screen") {
    if (as.numeric(file.access(savePlot)) < 0) {
      stop("The savePlot directory cannot be found in the system or invalid value for savePlot")
    }
  }
  
  bf <- LocData$BFarea
  it <- LocData$image.type
  datin <- LocData$dateIndex
  ff <- LocData$Processed.Files
  extra.info <- LocData$Outlier.Estimates
  if (!is.data.frame(LocData)) {
    LocData <- LocData$SpotResults
  }
  
  if (colnames(LocData)[ncol(LocData)] == "QCgroup") {
    legend <- c(colnames(LocData), "QCgroup2")
  } else {
    legend <- c(colnames(LocData), "QCgroup")
  }
  
  b <-
    t(matrix(unlist(strsplit(
      LocData$SampleID, separator
    )), nrow = 2))
  ffRun <- b[, 1]
  ffWell <- substr(b[, 2], 2, 3)
  WellID <- rep("left", length(ffWell))
  WellID[ffWell == "02" |
           ffWell == "05" |
           ffWell == "07" |
           ffWell == "09" |
           ffWell == "10" | ffWell == "12"] <- "right"
  Wellsets <- cbind(b, WellID)
  uR <- unique(ffRun)
  uW <- unique(WellID)
  res <- c()
  
  if (length(subset) > 0) {
    ssRun <- setdiff(subset[[1]], uR)
    if (length(ssRun) > 0) {
      stop("Some subset runs do not exist")
    }
    ssWell <- setdiff(subset[[2]], uW)
    if (length(ssWell) > 0) {
      stop("Some subset well directions do not exist")
    }
    uR <- subset[[1]]
    uW <- subset[[2]]
  }
  
  if (!(out.method %in% c("interactive.clustering", "interactive.manual", "manual"))) {
    stop("This QC method is not available")
  }
  
  test.interactive <- table(Wellsets[, 1], Wellsets[, 3])
  if (out.method == "interactive.clustering") {
    if (min(as.numeric(test.interactive)) < 15) {
      stop(
        "The interactive clustering method cannot be used in small datasets (below 15 samples). Consider one between interactive.manual or manual methods"
      )
    }
  }
  median.ests <- matrix(0, 1, 4)
  if (out.method != "manual") {
    for (i in 1:length(uR)) {
      for (j in 1:length(uW)) {
        w <- which(ffRun == uR[i] & WellID == uW[j])
        ulocs <- LocData[w,]
        plot(
          ulocs[, 2],
          ulocs[, 3],
          cex = 0.5,
          main = paste("Run=", uR[i], " / Well=", uW[j], sep = ""),
          xlim = c(0, dims[1]),
          ylim = c(0, dims[2]),
          xlab = "X coordinate (in pixels)",
          ylab = "Y coordinate (in pixels)"
        )
        mX <- median(as.numeric(ulocs[, 2]))
        mY <- median(as.numeric(ulocs[, 3]))
        median.ests <-
          matrix(rbind(median.ests, c(uR[i], uW[j], mX, mY)), ncol = ncol(median.ests))
        points(mX, mY, col = 3, pch = 3)
        
        if (out.method == "interactive.clustering") {
          umatrix <- data.matrix(ulocs[, 2:3])
          dists <-
            apply(umatrix, 1, distfromcenter, center = c(mX, mY))
          sl <- sort.list(dists)
          umatrix <- umatrix[sl,]
          ulocs <- ulocs[sl,]
          dists <- dists[sl]
          dd <- diff(sort(dists))
          cp <- grubbs(dd)
          if (length(cp) > 0.5 * length(dd)) {
            stop(
              "This out.method estimated an unreliable number of outliers. Try out.method = interactive.manual"
            )
          }
          if (length(cp) > 0) {
            estDiff <- rep(0, length(cp))
            for (k in 1:length(cp)) {
              estDiff[k] <-
                max(abs(mX - umatrix[cp[k], 1]), abs(mY - umatrix[cp[k], 2]))
              text(
                c(dims[1] * 0.8),
                c(dims[2] * 0.8) - (20 * k),
                paste(
                  "The estimated radius is ",
                  estDiff[k] - 0.5,
                  " pixels",
                  sep = ""
                ),
                col = (k + 1),
                cex = 0.7
              )
              draw.circle(mX,
                          mY,
                          radius = estDiff[k],
                          border = (k + 1))
            }
          } else {
            text(
              c(dims[1] * 0.8),
              c(dims[2] * 0.8) - 20,
              paste(
                "All locations are reliable. Set radius to ",
                dims[1],
                sep = ""
              ),
              col = 1,
              cex = 0.7
            )
          }
          
          qu <-
            readline(
              "Enter the radius of the area containing the reliable coordinates only (suggested values are on the figure):"
            )
          
          if (as.numeric(qu) >= 0 | qu != "") {
            g <- rep(1, nrow(umatrix))
            g[dists > as.numeric(qu)] <- 2
            plot(
              as.numeric(ulocs[, 2]),
              as.numeric(ulocs[, 3]),
              cex = 0.5,
              main = paste("Run=", uR[i], " / Well=", uW[j], sep = ""),
              xlim = c(0, dims[1]),
              ylim = c(0, dims[2]),
              xlab = "X coordinate (in pixels)",
              ylab = "Y coordinate (in pixels)"
            )
            points(umatrix, col = g, cex = 0.5)
            ab <-
              readline("Hit Enter to move to the next image or A + Enter to Abort: ")
            if (ab == "") {
              tres <- cbind(ulocs, g)
              names(tres)[ncol(tres)] <- legend[length(legend)]
              if (is.null(res)) {
                res <- tres
              } else {
                res <- rbind(res, tres)
              }
            } else {
              stop("The analysis has been stopped")
            }
          }
          
          if (qu == "") {
            dev.off()
            stop("Error in entering the number of clusters: the analysis has been stopped")
          }
        }
        
        
        if (out.method == "interactive.manual") {
          resout <- matrix(as.numeric(unlist(locator())), nrow = 1)
          if (length(resout) > 2) {
            resout <- matrix(resout, ncol = 2)
          }
          g <- rep(1, nrow(ulocs))
          if (length(resout) > 0) {
            for (b in 1:nrow(resout)) {
              ss <-
                abs(resout[b, 1] - as.numeric(ulocs[, 2])) + abs(resout[b, 2] - as.numeric(ulocs[, 3]))
              ss <- which(ss == min(ss))
              g[ss] <- 2
            }
          }
          plot(
            ulocs[, 2],
            ulocs[, 3],
            cex = 0.5,
            col = g,
            main = paste("Run=", uR[i], " / Well=", uW[j], sep = ""),
            xlim = c(0, dims[1]),
            ylim = c(0, dims[2]),
            xlab = "X coordinate (in pixels)",
            ylab = "Y coordinate (in pixels)"
          )
          ab <-
            readline("Hit Enter to move to the next image or A + Enter to Abort: ")
          if (ab == "") {
            tres <- cbind(ulocs, g)
            names(tres)[ncol(tres)] <- legend[length(legend)]
            if (is.null(res)) {
              res <- tres
            } else {
              res <- rbind(res, tres)
            }
          } else {
            stop("The analysis has been stopped")
          }
        }
      }
    }
    median.ests <-
      matrix(median.ests[-1,], ncol = ncol(median.ests))
  }
  
  if (out.method == "manual") {
    res <- cbind(LocData, rep("confidence", nrow(LocData)))
    names(res)[ncol(res)] <- "QCgroup"
    if(savePlot!="screen") pdf(paste(savePlot, "/defineLocClusters_", 
                                     datin, ".pdf", sep = ""),paper = "a4r")
    for (i in 1:length(uR)) {
      for (j in 1:length(uW)) {
        w <- which(ffRun == uR[i] & WellID == uW[j])
        ulocs <- LocData[w,]
        plot(
          ulocs[, 2],
          ulocs[, 3],
          cex = 0.5,
          main = paste("Run=", uR[i], " / Well=", uW[j], sep = ""),
          xlim = c(0, dims[1]),
          ylim = c(0, dims[2]),
          xlab = "X coordinate (in pixels)",
          ylab = "Y coordinate (in pixels)"
        )
        mX <- median(as.numeric(ulocs[, 2]))
        mY <- median(as.numeric(ulocs[, 3]))
        median.ests <-
          matrix(rbind(median.ests, c(uR[i], uW[j], mX, mY)), ncol = ncol(median.ests))
        points(mX, mY, col = 3, pch = 3)
        if (uW[j] == "left") {
          text((dims[1] * 0.85),
               (dims[2] * 0.1),
               "xxx 01/03/04/06/08/11",
               cex = 0.6)
        } else {
          text((dims[1] * 0.85),
               (dims[2] * 0.1),
               "xxx 02/05/07/09/10/12",
               cex = 0.6)
        }
      }
    }
    dev.off()
    print(
      "Identify outlier locations manually: Separate cells into confidence or outlier at the last column of the output table"
    )
    median.ests <-
      matrix(median.ests[-1,], ncol = ncol(median.ests))
  }
  
  if (legend[length(legend)] == "QCgroup2") {
    res <- res[,-c(ncol(res) - 1)]
    names(res)[ncol(res)] <- "QCgroup"
  }
  w1 <- which(as.numeric(res[, ncol(res)]) == 1)
  w2 <- which(as.numeric(res[, ncol(res)]) == 2)
  res[w1, ncol(res)] <- "confidence"
  res[w2, ncol(res)] <- "outlier"
  res <- res[sort.list(as.numeric(as.character(rownames(res)))),]
  
  return(
    list(
      Results = res,
      BFdata = extra.info,
      Processed.Files = ff,
      Outlier.indices = which(res[, ncol(res)] == "outlier"),
      Medians = median.ests,
      Wellsets = Wellsets,
      BFarea = bf,
      image.type = it,
      dateIndex = datin
      
    )
  )
}


#' spotEstimator
#'
#' The main function to produce the raw fluorescence signal estimation results by analysis of the Fluidigm images.
#'
#' Triplets of images of the same sample are sequentially considered to estimate the channel-specific
#'   fluorescence signals (if detectable) or perform BF image modeling. The main result of this function is a table
#'   of location and fluorescence estimates for each sample.
#'
#' @param files Character string. The file names to be read and analyzed. This is the output of readFiles()
#' @param correctionAlgorithm Logical. Its value specifies the estimation stage. If FALSE,
#'   the function processes all data using the standard operations of spotCoords(), i.e. case detection and fluorescence signal
#'   estimation. This is the first estimation stage. If TRUE, the function processes the BF image modeling estimates of outlier images
#'   obtained by defineLocClusters(). The BF image modeling is internally applied during the first stage. Note that
#'   correctionAlgorithm = TRUE is strictly used in the second (outliers adjustment / correction) stage of the process.
#' @param subset Numeric vector. It can be a series sample index numbers (a subset) that specifies the samples to be analyzed.
#'   The index numbers are obtained from readFiles() (the position of the sample in each listed vector). By default subset = c().
#'   The parameter is mainly used in the second estimation stage where spotEstimator() processes the outlier images (the index numbers
#    are automatically specified).
#' @param foregroundCut Numeric vector. The binary segmentation image analysis cutoffs for normalized image data. Pixels with normalized signals
#'   higher than the cutoff belong to foreground. Default is seq(0.5,0.7,0.02).
#' @param denoise Logical. If TRUE it denoises the channel images with la8, universal, hard. Default is FALSE.
#' @param despeckle Logical. If TRUE the bf image is descpeckled in the ImageJ fashion. Default is FALSE.
#' @param chip.type Character string. It specifies the type of Fluidigm chip to be analyzed. Default is "medium/large". The alternative
#'   option is "small".
#' @param cutSides Integer. It instructs the algorithm to find spots in a certain central image area. For example, for a 512 x 512
#'   image with cutSides = 50, spotEstimator() will search for spots in the central area [cutSides:(512-cutSides),cutSides:(512-cutSides)]
#'   of the image matrix. Default is 0.
#' @param BFarea Integer. Defines a rectangular pseudo-spot size whose fluorescence will be estimated. This is mainly used in BF image
#'   modeling where a fluorescence spot could not be originally detected. The value of this parameter is also used as a cut-off
#'   to find matched spots across channel of the same sample image. Default is 7.
#' @param log.transform Logical. If TRUE the image data are plotted in the log scale. Default is TRUE
#' @param minDiff Float. The mu_hat of the H0: image-to-noise ratio =
#'   log(foreground_signal) - log(background_signal) = mu_hat. Rejection of H0
#'   implies that the identified spot is brighter than background. Default is 0.5.
#' @param show.possible.contamination Logical. If TRUE it reports all identified unmatched spots in both channels. Default is TRUE.
#' @param cutoff Integer. A cutoff of the distance between the estimated spot location of an outlier sample (X, Y) and the median
#'   location of all non-outliers of the same run and well set (medX,medY), i.e. (X-medX, Y-medY). An outlier sample can either
#'   have a fluorescence-based location (X, Y) or a BF-based location (X*, Y*) or both. It is re-adjusted as follows: (1) if
#'   min{(X-medX, Y-medY)} > cutoff and min{(X*-medX, Y*-medY)} > cutoff, the sample's location is set to (medX, medY); (2) if
#'   min{(X*-medX, Y*-medY)} <= cutoff, the sample's location is set to (X*, Y*); (3) if min{(X-medX, Y-medY)} <= cutoff and
#'   min{(X*-medX, Y*-medY)} > cutoff, the algorithm can either produce the solution of (1) or the solution of (2) depending
#'   on the value of median.correction parameter below. By default cutoff = 50.
#' @param QCdata List. The output of defineLocClusters().
#' @param median.correction Logical. If TRUE, the algorithm re-adjusts the location of the outlier sample as the median of all
#'   non-outliers of the same run and well ID (if necessary).
#' @param savePlot Character string. Directory to store the plots. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen. Default is the current working directory, getwd().
#'
#' @return A list of the following components:
#'   SpotResults: the matrix of the location and fluorescence signal estimates. It contains the index number of each sample, the X,Y
#'     coordinates of the spot center, the spot size, the type of estimation that have been performed (fluorescence based indicating the channels
#'     in which the spot has been found or BF image modelling based), the fluorescence foreground and background signals of each channel,
#'     the signal-to-noise ratio (logForeground - logBackground) for each channel, the associated P-value of significance of the signal-to-noise
#'     ratio and a column indicating the coordinates of other spots that are not matched in both images. Existence of such spots (values that are
#'     different from 0) indicate contaminated image or highly noisy images or images with other artefacts. If correctionAlgorithm=TRUE (second
#'     spotEstimator() step), there is an extra column generated indicating outlier samples (see the QCgroup column in defineLocClusters()).
#'   Outlier.Estimates: The estimates obtained from BF modeling (if necessary to be obtained). These are alternative location estimates that will
#'     be used if the original estimates of the SpotResults table are flagged as outliers.
#'   Processed.Files: the samples that have been processed by spotEstimator().
#'   BFarea: the pseudospot size.
#'   image.type: the image type IDs.
#'   dateIndex: a date index to be used in saving the output files.
#'
#' @import grDevices
#' @export
#'
#' @examples
#' ### set your directories
#' basedir<-"~/"
#' #fucci_path<-system.file("extdata",package="fucci")
#' #files<-readFiles(iDirectory=NULL,
#' #                 BFdirectory=paste(fucci_path,"/BF",sep=""),
#' #                 CHdirectory=paste(fucci_path,"/CH",sep=""),
#' #                 separator = "_",image.type = c("BF","Green","Red"),
#' #                 bits=2^16)
#'
#' ### an example where the second image produces a clear outlier!
#' #estimates <- spotEstimator(files=files,subset=1:3,foregroundCut=seq(0.6,0.76,0.02),
#' #                           correctionAlgorithm=FALSE,savePlot="screen")
spotEstimator <-
  function(files,
           correctionAlgorithm,
           subset = c(),
           foregroundCut = seq(0.5, 0.7, 0.02),
           denoise = FALSE,
           despeckle = FALSE,
           chip.type = "medium/large",
           cutSides = 0,
           BFarea = 7,
           log.transform = TRUE,
           minDiff = 0.5,
           show.possible.contamination = TRUE,
           cutoff = 50,
           QCdata = 0,
           median.correction = TRUE,
           savePlot = getwd()) {
    # do not stop the process if there are no outliers. Just store the correct output
    if (correctionAlgorithm == TRUE & length(subset) == 0) {
      result <- QCdata$Results
      out <- QCdata$BFdata
      missed.outs <-
        which(as.numeric(result$X) == 0 & as.numeric(result$Y) == 0)
      result$QCgroup[missed.outs] <- "outlier"
      
      # fix the contamination column
      new <- as.character(result$QCgroup)
      ww <-
        which(
          as.character(result$Estimation.Type) == "Both.Channels" &
            as.character(result$QCgroup) == "outlier" |
            as.character(result$Estimation.Type) == "BFmedian2"
        )
      new[ww] <- "contamination"
      result$QCgroup <- new
      
      # fix the Estimation.Type
      result$Estimation.Type[grep("Channel", as.character(result$Estimation.Type))] <-
        "Fluorescence-based"
      result$Estimation.Type[grep("BF", as.character(result$Estimation.Type))] <-
        "Chip.Pattern-based"
      
      message(
        "There are no outliers to be corrected.
        Use LocationMatrix() directly with the rules of filtering (vignette)!"
      )
    } else {
      if (length(subset) > 0) {
        subset <- subset[which(subset <= length(files[[1]]) & subset > 0)]
        
        if (length(subset) > 0 &
            length(which(subset > length(files[[1]]) |
                         subset <= 0)) > 0) {
          print(paste("Some subset numbers have been ignored", sep = ""))
        }
        if (length(subset) == 0) {
          stop("The subset numbers are not valid")
        }
      } else if (length(subset) == 0) {
        subset <- 1:length(files[[1]])
      }
      
      if (correctionAlgorithm == TRUE) {
        BFarea <- QCdata$BFarea
        if (length(QCdata$Processed.Files) == 0) {
          stop(
            "The QCdata parameter is not correctly specified. Use an appropriate value for QCdata parameter"
          )
        }
        if (any(files$BF != QCdata$Processed.Files$BF)) {
          files <- QCdata$Processed.Files
        }
      }
      
      if(savePlot != "screen") {
        if (as.numeric(file.access(savePlot)) < 0) {
          stop(
            "The savePlot directory cannot be found in the system or invalid value for savePlot"
          )
        } else if (correctionAlgorithm)
          pdf(
            paste(
              savePlot,
              "/Step2_spotEstimator_",
              files$dateIndex,
              ".pdf",
              sep = ""
            ),
            paper = "a4r"
          )
        else
          pdf(
            paste(
              savePlot,
              "/Step1_spotEstimator_",
              files$dateIndex,
              ".pdf",
              sep = ""
            ),
            paper = "a4r"
          )
      }
      
      for (i in 1:3) {
        files[[i]] <- files[[i]][subset]
      }
      
      
      result <-
        data.frame(matrix("", nrow = length(files$BF), ncol = 14), stringsAsFactors = FALSE)
      names(result) <-
        c(
          "SampleID",
          "X",
          "Y",
          "Size",
          "Estimation.Type",
          paste("fore_", files$image.type[2], sep = ""),
          paste("back_", files$image.type[2], sep = ""),
          paste("fore_", files$image.type[3], sep = ""),
          paste("back_", files$image.type[3], sep = ""),
          paste(files$image.type[2], ".StN", sep = ""),
          paste(files$image.type[2], ".Pvalue", sep = ""),
          paste(files$image.type[3], ".StN", sep = ""),
          paste(files$image.type[3], ".Pvalue", sep = ""),
          "Other.Spots"
        )
      out <- as.list(rep(0, length(files$BF)))
      #par(mfrow = c(2,2))
      for (k in 1:length(files$BF)) {
        if (correctionAlgorithm == FALSE) {
          message(paste("Analyzing image", k, "of", length(files$BF)))
        } else {
          message(paste("Analyzing outlier image", k, "of", length(files$BF)))
        }
        
        chaImages <-
          readChaImg(imgNames = list(CH1 = files$CH1[k], CH2 = files$CH2[k]),
                     denoise = denoise)
        
        if (correctionAlgorithm == FALSE) {
          R1 <-
            spotCenter(
              img = chaImages$CH1.proc,
              foregroundCut = foregroundCut,
              howbig = BFarea,
              ImgLimits = cutSides
            )
          G1 <-
            spotCenter(
              img = chaImages$CH2.proc,
              foregroundCut = foregroundCut,
              howbig = BFarea,
              ImgLimits = cutSides
            )
          
          spot <-
            spotCoords(
              centerR = R1,
              centerG = G1,
              origImg = files$BF[k],
              chaImgs = list(CH1 = chaImages$CH1, CH2 = chaImages$CH2),
              minDiff = minDiff,
              despeckle = despeckle,
              ImgLimits = cutSides,
              BFarea = BFarea,
              chip.type = chip.type,
              separator = files$separator,
              image.type = files$image.type,
              show.possible.contamination = show.possible.contamination
            )
        } else {
          spot <-
            forceBF(
              data = list(QCdata$Results[subset[k],], QCdata$BFdata[[subset[k]]]),
              cutoff = cutoff,
              median.correction = median.correction,
              medians = QCdata$Medians,
              Wells = QCdata$Wellsets[subset[k],],
              image.type = files$image.type
            )
        }
        
        result1 <-
          SpotStats(
            img = files$BF[k],
            chaImgs = list(CH1 = chaImages$CH1, CH2 = chaImages$CH2),
            binChaImgs =
              list(CH1 = spot$areaCH1, CH2 = spot$areaCH2),
            center = spot$center,
            BFcoords =
              spot$failure,
            BFarea = BFarea,
            warning = spot$warning,
            minDiff = minDiff,
            other.spots = spot$other.spots,
            log.transform = log.transform,
            separator = files$separator,
            image.type = files$image.type
          )
        
        result[k,] <- result1$stats
        if (correctionAlgorithm == FALSE) {
          out[[k]] <- spot$outlier.estimates
          out[[k]]$sample <- as.character(result1$stats[1])
        }
      }
      if (savePlot != "screen") {
        dev.off()
      }
      
      if (correctionAlgorithm == TRUE) {
        mm <- match(result[, 1], QCdata$Results[, 1])
        for (i in 1:length(mm)) {
          QCdata$Results[mm[i],] <- c(result[i,], "confidence")
        }
        result <- QCdata$Results
        out <- QCdata$BFdata
      }
      
      # fixing columns
      for (i in c(2, 3, 4, 6:13)) {
        result[, i] <- as.numeric(result[, i])
      }
      
      if (correctionAlgorithm == TRUE) {
        # fix the contamination column
        new <- as.character(result$QCgroup)
        ww <-
          which(
            as.character(result$Estimation.Type) == "Both.Channels" &
              as.character(result$QCgroup) == "outlier" |
              as.character(result$Estimation.Type) == "BFmedian2"
          )
        new[ww] <- "contamination"
        result$QCgroup <- new
        
        # fix the Estimation.Type
        result$Estimation.Type[grep("Channel", as.character(result$Estimation.Type))] <-
          "Fluorescence-based"
        result$Estimation.Type[grep("BF", as.character(result$Estimation.Type))] <-
          "Chip.Pattern-based"
        
        message("")
        message("")
        message("See the rules of filtering in the vignette")
      }
    }
    
    return(
      list(
        SpotResults = result,
        Outlier.Estimates = out,
        Processed.Files = files,
        BFarea = BFarea,
        image.type = files$image.type,
        dateIndex = files$dateIndex
      )
    )
    }



#' LocationMatrix
#'
#' It generates the final cell location and fluorescnece signal estimates and summarizes the quality
#'   control statistics.
#'
#' @param data Data matrix. The matrix of the location and fluorescence signal estimates after two rounds (maximum)
#'   of spotEstimator().
#' @param filter.by Data matrix. A series of filtering criteria and cut-offs that specify which samples are KEPT for
#'   further analysis (see vignette). By default it flags by FDR (alpha = 0.005) and outlier index (keeps only the 'confident'
#'   estimates).
#' @param report.by.signif Character string. It returns the pre-defined channel-specific signal-to-noise ratio and test
#'   statistics for each sample. If "min", the algorithm only reports the P-values/FDRs and signal-to-noise of the channel
#'   with the minimum signal-to-noise ratio. If "max", the algorithm only reports the P-values/FDRs and signal-to-noise of
#'   the channel with the maximum signal-to-noise ratio.  Default is "max".
#'
#' @return List. The first component is a data matrix of the final table of estimates. The main body of this table has been generated by spotEstimator().
#'   It summarizes the location, the raw fluorescence signal estimates (foreground and background) and the quality control statistics. It keeps
#'   only the signal-to-noise ratio and the associated P-value/FDR of a predefined channel (see parameter report.by.signif). The last column
#'   ("Cells") consists of 1s for the samples that pass the filtering step (filter.by) and are used for further analysis. The rest of the samples
#'   are assigned 0s. The user should always inspect them along with the images to obtain the final list of samples to be used for further analysis.
#'   The second component is the date index for storing the output files. It is transfered to the next step.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' ### the results matrix (column 'Cells') indicates three empty capture chambers
#' ### (thus not only outliers were associated with the absense of a cell!)
#' Results <- LocationMatrix(data=estimates.2,
#'            filter.by = matrix(c("FDR","Out.Index",0.005,"confidence"),ncol=2))
LocationMatrix <-
  function(data,
           filter.by = matrix(c("FDR", "Out.Index", 0.005, "confidence"), ncol = 2),
           report.by.signif = "max") {
    it <- data$image.type
    datin <- data$dateIndex
    data <- data[[1]]
    
    # fix the contamination column (in case the second spotEstimator was followed by outlier plots to get more outliers)
    new <- as.character(data$QCgroup)
    ww <-
      which(
        as.character(data$Estimation.Type) == "Both.Channels" &
          as.character(data$QCgroup) == "outlier" |
          as.character(data$Estimation.Type) == "BFmedian2"
      )
    new[ww] <- "contamination"
    data$QCgroup <- new
    
    p <- matrix(1, nrow(data), 2)
    ww <- which(data$QCgroup == "confidence")
    vec <-
      matrix(c(as.numeric(data[ww, 11]), as.numeric(data[ww, 13])), ncol = 2)
    vec <- t(matrix(p.adjust(c(t(
      vec
    )), "BH"), nrow = 2))
    for (i in 1:length(ww)) {
      p[ww[i],] <- vec[i,]
    }
    colnames(p) <-
      c(paste(it[2], ".FDR", sep = ""), paste(it[3], ".FDR", sep = ""))
    
    pv <-
      data.matrix(data[, c(paste(it[2], ".Pvalue", sep = ""), paste(it[3], ".Pvalue", sep =
                                                                      ""))])
    f <-
      data.matrix(data[, c(paste(it[2], ".StN", sep = ""), paste(it[3], ".StN", sep =
                                                                   ""))])
    ww <- which(data$QCgroup == "outlier")
    if (length(ww) > 0) {
      pv[ww,] <- c(1, 1)
      f[ww,] <- c(0, 0)
      p[ww,] <- c(1, 1)
    }
    
    
    # report the appropriate FC,PV,FDR (min/max)
    if (report.by.signif == "min") {
      ww <- apply(f, 1, which.min)
    }
    if (report.by.signif == "max") {
      ww <- apply(f, 1, which.max)
    }
    res <- matrix(0, length(ww), 4)
    for (i in 1:length(ww)) {
      res[i, 1:4] <-
        c(f[i, ww[i]], pv[i, ww[i]], p[i, ww[i]], data$QCgroup[i])
    }
    
    
    # set the filters
    if (length(filter.by) == 1) {
      stop("Parameter filter.by should contain the filter type and the filter rule")
    }
    if (length(filter.by) == 2) {
      filter.by <- matrix(filter.by, nrow = 1)
    }
    
    cells <- c()
    if (length(filter.by) == 0) {
      message("All data are kept for further analysis!")
      cells <- rep(1, nrow(data))
    }
    
    if (length(filter.by) > 0) {
      cells <- matrix(0, nrow(data), nrow(filter.by))
      
      # stop the process if filtering is misspelled
      match.criteria <- match(
        filter.by[, 1],
        c(
          "Size",
          "Estimation.Type",
          "Pvalue",
          "StN",
          "FDR",
          paste("Pvalue/", it[2], sep =
                  ""),
          paste("StN/", it[2], sep = ""),
          paste("FDR/", it[2], sep = ""),
          paste("Pvalue/", it[3], sep =
                  ""),
          paste("StN/", it[3], sep = ""),
          paste("FDR/", it[3], sep = ""),
          "Out.Index",
          "Other.Spots"
        ),
        nomatch = 0
      )
      if (length(match.criteria[match.criteria == 0]) > 0) {
        stop("Some filter.by types are misspelled or not available. Check the vignette")
      }
      
      for (i in 1:nrow(filter.by)) {
        if (filter.by[i, 1] == "Size") {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(as.numeric(data[, 4]) >= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == "Estimation.Type") {
          cut <- unlist(strsplit(filter.by[i, 2], "/"))
          cells1 <- matrix(0, nrow(data), length(cut))
          for (j in 1:length(cut)) {
            cells1[which(as.character(data[, 5]) == cut[j]), j] <- 1
          }
          cells[, i] <-
            apply(matrix(cells1, ncol = length(cut)), 1, min)
        }
        
        if (filter.by[i, 1] == "Pvalue") {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          pp <- as.numeric(apply(data[, c(11, 13)], 1, min))
          cells[which(pp <= as.numeric(filter.by[i, 2])), i] <- 1
        }
        
        if (filter.by[i, 1] == "StN") {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          pp <- as.numeric(apply(abs(data[, c(10, 12)]), 1, max))
          cells[which(pp >= as.numeric(filter.by[i, 2])), i] <- 1
        }
        
        if (filter.by[i, 1] == "FDR") {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          pp <- as.numeric(apply(p, 1, min))
          cells[which(pp <= as.numeric(filter.by[i, 2])), i] <- 1
        }
        
        if (filter.by[i, 1] == paste("Pvalue/", it[2], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(as.numeric(data[, 11]) <= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == paste("StN/", it[2], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(as.numeric(data[, 10]) >= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == paste("FDR/", it[2], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(p[, 1] <= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == paste("Pvalue/", it[3], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(as.numeric(data[, 13]) <= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == paste("StN/", it[3], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(as.numeric(data[, 12]) >= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == paste("FDR/", it[3], sep = "")) {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          cells[which(p[, 2] <= as.numeric(filter.by[i, 2])), i] <-
            1
        }
        
        if (filter.by[i, 1] == "Out.Index") {
          cut <- unlist(strsplit(filter.by[i, 2], "/"))
          cells1 <- matrix(0, nrow(data), length(cut))
          for (j in 1:length(cut)) {
            cells1[which(as.character(data$QCgroup) == cut[j]), j] <- 1
          }
          cells[, i] <-
            apply(matrix(cells1, ncol = length(cut)), 1, min)
        }
        
        if (filter.by[i, 1] == "Other.Spots") {
          if (is.na(as.numeric(filter.by[i, 2]))) {
            stop("The second element of filter.by should be numeric")
          }
          bb <- rep(0, nrow(data))
          for (j in 1:nrow(data)) {
            bb1 <-
              as.character(unlist(strsplit(
                as.character(data$Other.Spots)[j], " / "
              )))
            bb[j] <- length(bb1[bb1 != "0"])
          }
          cells[which(bb <= as.numeric(filter.by[i, 2])), i] <- 1
        }
      }
      cells <- apply(matrix(cells, nrow = nrow(data)), 1, min)
      if (sum(cells) == 0) {
        message(
          "There are no data matching the criteria of filter.by. All samples are excluded from further analysis!"
        )
      }
    }
    result <- cbind(data[, c(1:9)], res[, 1:4], data[, 14], cells)
    colnames(result) <-
      c(
        colnames(data)[c(1:9)],
        "Signal-to-Noise",
        "Pvalue",
        "FDR",
        "Out.Index",
        "Other.Spots",
        "Cells"
      )
    return(list(Output = result, dateIndex = datin))
  }
