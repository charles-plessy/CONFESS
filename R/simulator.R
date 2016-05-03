# The data simulation functions including the helpers

#' spot.simulator
#'
#' A helper to simulate spots in a image.
#'
#' @param location Numeric vector. A vector of 2d coordinates for the spot center.
#' @param size Integer. The number of pixels the spot consists of.
#' @param average.signal Float. The parameter of the poisson distribution to generate the
#    pixel signals for a spot. This is the average of all pixels.
#' @param dimension Numeric vector. The image dimensions.
#'
#' @return The image with the generated spot(s)
#'
#' @import stats
#'
#' @keywords internal
spot.simulator <-
  function(location, size, average.signal, dimension) {
    testimg <- matrix(0, dimension[1], dimension[2])
    ff <- ceiling(sqrt(size) / 2)
    coords <-
      expand.grid((location[1] - ff):(location[1] + ff),
                  (location[2] - ff):(location[2] + ff))
    coords <- as.matrix(coords[sample(1:nrow(coords)),], ncol = 2)
    dcoords <-
      matrix(cbind(matrix(abs(
        coords[, 1] - location[1]
      ), ncol = 1), matrix(abs(
        coords[, 2] - location[2]
      ), ncol = 1)), ncol = 2)
    coords <- coords[sort.list(as.numeric(apply(dcoords, 1, sum))),]
    signal <- sort(rpois(size, average.signal), decreasing = TRUE)
    for (i in 1:length(signal)) {
      testimg[coords[i, 1], coords[i, 2]] <- signal[i]
    }
    return(list(mat = testimg, coo = coords[1:length(signal),]))
  }

#' signal.from.matrix
#'
#' A helper to simulate spots in a image.
#'
#' @param coords Numeric vector. A vector of 2d coordinates for the spot center.
#' @param mat Data Matrix. The data matrix which generates the image of interest.
#'
#' @return Some values of interest
#'
#' @keywords internal
signal.from.matrix <- function(coords, mat) {
  return(mat[coords[1], coords[2]])
}

#' mean_signal
#'
#' A helper to simulate the spot signal.
#'
#' @param data Numeric vector. An 1-dimensional vector of spot signals.
#' @param noise.level Float. The noise level of the image.
#'
#' @return The average spot's noisy signal
#'
#' @keywords internal
mean_signal <- function(data, noise.level) {
  ans <- 2 ^ mean(log(data + noise.level, 2))
  return(ans)
}


#' simcells
#'
#' The main function to simulate spots of various numbers, sizes, signals in one or multiple images
#'   of a given dimension.
#'
#' @param channels Integer. The number of channels for each sample. Default is 2.
#' @param spots.per.image Numeric vector. The number of spots in each image (channel). The length
#'   of the vector equals to the number of channels. Default is one spot per channel.
#' @param one.location Numeric vector. The central location of the matched spots across the channels
#'   (in pixel) coordinates. Default is (X,Y) = (50,50).
#' @param image.dimension Numeric vector. The image dimension (in pixels). Default is 100 x 100.
#' @param signal.level List. The lambda parameter of the Poisson distribution that generates the
#'   true spot (pixel) signals. The list has as many components (length) as the number of channels.
#'   The number of elements of each component equals to the number of spots in each particular channel.
#'   Default is list(700,700).
#' @param noise.level Numeric vector. The sigma parameter of the Normal distribution that generates the
#'   image noise level. The length of the vector equals to the number of channels. Default is c(200,200).
#' @param spot.size List. The size of each spot on each channel (in pixels). The list has as many components
#'   (length) as the number of channels. The number of elements of each component equals to the number of
#'   spots in each particular channel. Default is list(30,30).
#' @param agreement.number Integer. It defines how many spot pairs are matched, i.e. they are located in the
#'   same coordinates across channels. These reflect true cells. Default is 1 corresponding to a single-cell
#'   case study.
#'
#' @return The image(s) with the generated spot(s). It consists of the data matrices and the location of the spot centers.
#'
#' @import stats graphics
#' @export
#'
#' @examples
#' r<-simcells(channels = 2, spots.per.image = c(2, 3), one.location = c(50, 50),
#' image.dimension = rep(200, 2), signal.level = list(c(1000, 1000), c(1000, 700, 300)),
#' noise.level = c(100, 100),spot.size = list(c(81, 100), c(26, 29, 50)), agreement.number = 1)
#'
#' r<-simcells(channels = 2, spots.per.image = c(0, 0), image.dimension = rep(200, 2),
#' signal.level = list(c(),c()),noise.level = c(0, 0), spot.size = list(c(), c()))
simcells <-
  function(channels = 2,
           spots.per.image = c(1, 1),
           one.location = c(50, 50),
           image.dimension = rep(100, 2),
           signal.level = list(700, 700),
           noise.level = c(200, 200),
           spot.size = list(30, 30),
           agreement.number = 1) {
    res <- as.list(rep(0, channels))
    ff <- floor(0.2 * image.dimension[1])
    
    signal.vec <- as.list(rep(0, channels))
    for (i in 1:channels) {
      signal.vec[[i]] <- as.list(rep(0, max(spots.per.image[i], 1)))
    }
    
    loc1 <- c()
    loc2 <- c()
    par(mfrow = c(1, 2))
    
    if (channels > 2) {
      stop("This function supports maximum 2 channels")
    }
    if (channels != length(spots.per.image)) {
      stop(
        "Parameter spots.per.image should be a vector whose size equals to the number of channels"
      )
    }
    if (length(one.location) != 2) {
      stop(
        "Parameter one.location is a vector of size two containing the 2d coordinates of a single spot"
      )
    }
    if (length(image.dimension) == 1) {
      stop("Fix the image dimensions as (number of rows) (number of columns)")
    }
    if (channels != length(signal.level)) {
      stop("Parameter signal.level should be a list whose size equals the number of channels")
    }
    if (channels != length(signal.level)) {
      stop("Parameter signal.level should be a list whose size equals the number of channels")
    }
    if (any(unlist(lapply(signal.level, length)) != spots.per.image)) {
      stop(
        "Parameter signal.level should be a list whose size equals the number of channels. The size of each component equals to the number of spots for each image"
      )
    }
    if (channels != length(noise.level)) {
      stop("Parameter noise.level is a vector with size equal the number of channels")
    }
    if (channels != length(spot.size)) {
      stop("Parameter spot.size should be a list whose size equals the number of channels")
    }
    if (any(unlist(lapply(spot.size, length)) != spots.per.image)) {
      stop(
        "Parameter spot.size should be a list whose size equals the number of channels. The size of each component equals to the number of spots for each image"
      )
    }
    if (agreement.number > min(spots.per.image)) {
      agreement.number <- min(spots.per.image)
    }
    
    for (i in 1:channels) {
      if (i == 1) {
        res[[i]] <-
          matrix(
            rnorm(image.dimension[1] ^ 2, 0, noise.level[i]),
            image.dimension[1],
            image.dimension[2]
          )
        set.row <- ff:(image.dimension[1] - ff)
        set.col <- ff:(image.dimension[2] - ff)
        if (spots.per.image[i] > 0) {
          loc1 <- matrix(0, spots.per.image[i], 2)
          for (j in 1:spots.per.image[i]) {
            if (length(one.location) > 0 &
                j == 1 &
                min(spots.per.image) > 0 &
                agreement.number > 0 & channels == 2) {
              if (one.location[1] <= min(set.row)) {
                one.location[1] <- min(set.row)
                print(
                  paste(
                    "The row number of the one.location parameter has been truncated to ",
                    min(set.row),
                    sep = ""
                  )
                )
              }
              if (one.location[1] >= max(set.row)) {
                one.location[1] <- max(set.row)
                print(
                  paste(
                    "The row number of the one.location parameter has been truncated to ",
                    max(set.row),
                    sep = ""
                  )
                )
              }
              if (one.location[2] <= min(set.col)) {
                one.location[2] <- min(set.col)
                print(
                  paste(
                    "The column number of the one.location parameter has been truncated to ",
                    min(set.col),
                    sep = ""
                  )
                )
              }
              if (one.location[2] >= max(set.col)) {
                one.location[2] <- max(set.col)
                print(
                  paste(
                    "The column number of the one.location parameter has been truncated to ",
                    max(set.col),
                    sep = ""
                  )
                )
              }
              
              loc1[j,] <- one.location
            } else {
              loc1[j,] <- c(sample(set.row, 1), sample(set.col, 1))
            }
            wrow <-
              which(abs(set.row - loc1[j, 1]) > (1.5 * sqrt(spot.size[[i]][j])))
            wcol <-
              which(abs(set.col - loc1[j, 2]) > (1.5 * sqrt(spot.size[[i]][j])))
            set.row <- set.row[wrow]
            set.col <- set.col[wcol]
            if (length(set.row) < 3 & length(set.col) < 3) {
              stop(
                "The image is to dense with spots. Consider increasing the image dimension or decreasing the number of spots"
              )
            }
            res1 <-
              spot.simulator(
                location = loc1[j,],
                size = spot.size[[i]][j],
                average.signal = signal.level[[i]][j],
                dimension = image.dimension
              )
            res[[i]] <- res[[i]] + res1$mat
            signal.vec[[i]][[j]] <-
              as.numeric(apply(res1$coo, 1, signal.from.matrix, mat = res[[i]]))
          }
        }
        noi <- abs(min(res[[i]]))
        res[[i]] <- round(res[[i]] + noi, 0)
        img <- res[[i]][nrow(res[[i]]):1,] / max(res[[i]])
        img.m <- melt(img)
        names(img.m) <- c("x", "y", "z")
        sim1 <-
          ggplot(img.m, aes_string(x = 'x', y = 'y', fill = 'z')) + geom_raster() + theme_bw() + scale_fill_distiller(palette = "Greys") +
          ggtitle(paste(
            "The ",
            spots.per.image[i],
            " simulated spot(s) on channel ",
            i,
            sep = ""
          )) + theme(legend.position = "none")
        #image(t(res[[i]][nrow(res[[i]]):1,]/max(res[[i]])),xaxt="n",yaxt="n",main=paste("The ",spots.per.image[i]," simulated spot(s) on channel ",i,sep=""))
      }
      
      if (i > 1) {
        res[[i]] <-
          matrix(
            rnorm(image.dimension[1] ^ 2, 0, noise.level[i]),
            image.dimension[1],
            image.dimension[2]
          )
        
        if (spots.per.image[i] > 0) {
          loc2 <- matrix(0, spots.per.image[i], 2)
          
          if (length(loc1) == 0) {
            agreement.number <- 0
          }
          if (agreement.number == 0) {
            for (j in 1:spots.per.image[i]) {
              loc2[j,] <- c(sample(set.row, 1), sample(set.col, 1))
              wrow <-
                which(abs(set.row - loc2[j, 1]) > (1.5 * sqrt(spot.size[[i]][j])))
              wcol <-
                which(abs(set.col - loc2[j, 2]) > (1.5 * sqrt(spot.size[[i]][j])))
              set.row <- set.row[wrow]
              set.col <- set.col[wcol]
              if (length(set.row) < 3 & length(set.col) < 3) {
                stop(
                  "The image is to dense with spots. Consider increasing the image dimension or decreasing the number of spots"
                )
              }
              res1 <-
                spot.simulator(
                  location = loc2[j,],
                  size = spot.size[[i]][j],
                  average.signal = signal.level[[i]][j],
                  dimension = image.dimension
                )
              res[[i]] <- res[[i]] + res1[[1]]
              signal.vec[[i]][[j]] <-
                as.numeric(apply(res1$coo, 1, signal.from.matrix, mat = res[[i]]))
            }
          }
          
          if (agreement.number > 0) {
            loc2[1:agreement.number,] <- loc1[1:agreement.number,]
            for (j in 1:agreement.number) {
              res1 <-
                spot.simulator(
                  location = loc2[j,],
                  size = spot.size[[i]][j],
                  average.signal = signal.level[[i]][j],
                  dimension = image.dimension
                )
              res[[i]] <- res[[i]] + res1$mat
              signal.vec[[i]][[j]] <-
                as.numeric(apply(res1$coo, 1, signal.from.matrix, mat = res[[i]]))
            }
            
            if (agreement.number < spots.per.image[i]) {
              for (j in (1 + agreement.number):spots.per.image[i]) {
                if (length(set.row) < 3 & length(set.col) < 3) {
                  stop(
                    "The image is to dense with spots. Consider increasing the image dimension or decreasing the number of spots"
                  )
                }
                loc2[j,] <-
                  c(sample(set.row, 1), sample(set.col, 1))
                wrow <-
                  which(abs(set.row - loc2[j, 1]) > (1.5 * sqrt(spot.size[[i]][j])))
                wcol <-
                  which(abs(set.col - loc2[j, 2]) > (1.5 * sqrt(spot.size[[i]][j])))
                set.row <- set.row[wrow]
                set.col <- set.col[wcol]
                if (length(set.row) < 3 & length(set.col) < 3) {
                  stop(
                    "The image is to dense with spots. Consider increasing the image dimension or decreasing the number of spots"
                  )
                }
                res1 <-
                  spot.simulator(
                    location = loc2[j,],
                    size = spot.size[[i]][j],
                    average.signal = signal.level[[i]][j],
                    dimension = image.dimension
                  )
                res[[i]] <- res[[i]] + res1[[1]]
                signal.vec[[i]][[j]] <-
                  as.numeric(apply(res1$coo, 1, signal.from.matrix, mat = res[[i]]))
              }
            }
          }
        }
        noi <- abs(min(res[[i]]))
        res[[i]] <- round(res[[i]] + noi, 0)
        img <- res[[i]][nrow(res[[i]]):1,] / max(res[[i]])
        img.m <- melt(img)
        names(img.m) <- c("x", "y", "z")
        sim2 <-
          ggplot(img.m, aes_string(x = 'x', y = 'y', fill = 'z')) + geom_raster() + theme_bw() + scale_fill_distiller(palette = "Greys") +
          ggtitle(paste(
            "The ",
            spots.per.image[i],
            " simulated spot(s) on channel ",
            i,
            sep = ""
          )) + theme(legend.position = "none")
        #image(t(res[[i]][nrow(res[[i]]):1,]/max(res[[i]])),xaxt="n",yaxt="n",main=paste("The ",spots.per.image[i]," simulated spot(s) on channel ",i,sep=""))
      }
    }
    
    multiplot(sim1, sim2, cols = 2)
    spot.center = list(Ch1 = loc1, Ch2 = loc2)
    for (i in 1:channels) {
      if (spots.per.image[i] > 0) {
        ll <- lapply(signal.vec[[i]], mean_signal, noise.level = noi)
        spot.center[[i]] <-
          matrix(cbind(spot.center[[i]], ll), nrow = nrow(spot.center[[i]]))
      }
    }
    return(list(Img = res, Spots = spot.center))
  }
