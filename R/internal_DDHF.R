# Data Driven Haar Fisz for multivariate data (internal functions)

#' DDHFfit
#'
#' An internal function to produce the DDHFmv clustering and other model estimates.
#'
#' @param data List. A list of fluorescence signals with their change-points and clusters.
#' @param den.method Character string. A method to perform the denoising. One of "wavelets", "splines"
#'   or "lregr" (for linear regression).
#' @param savePlot Character string. The directory to store the plots.
#'
#' @import grDevices stats utils
#' @return The DDHFmv clusters and model estimates
#'
#' @keywords internal
DDHFfit <- function(data, den.method, savePlot) {
  ccc <-
    colors()[c(24,
               26,
               32,
               81,
               630,
               68,
               152,
               254,
               382,
               498,
               536,
               552,
               450,
               547,
               652,
               635,
               617,
               541)]
  newColors <- ccc[data$Updated.groups]
  dd <-
    matrix(cbind(data$corrected.VStransformed.exprs, newColors),
           ncol = 3)
  
  infl.stats <- matrix(c("index", "Cook"), nrow = 1)
  xout1 <-
    matrix(cbind(data$index, data$Progression, data$Updated.groups),
           ncol = 4)
  for (i in 1:length(sort(unique(xout1[, 4])))) {
    xx <- xout1[which(xout1[, 4] == i), c(1, 3, 2)]
    ols <- lm(xx[, 2] ~ xx[, 3])
    ss <-
      capture.output(summary(influence.measures(ols)), file = NULL)
    if (length(rownames(ss)) > 0) {
      a <- xx[as.numeric(rownames(ss)), 1]
      b <- as.numeric(ss[, 5])
      infl.stats <- matrix(rbind(infl.stats, cbind(a, b)), ncol = 2)
    }
  }
  
  xout <- denoiser(data = xout1, method = "lregr")
  fdata <- matrix(0, 1, ncol(xout))
  fouts <- c()
  for (i in 1:length(sort(unique(xout[, 4])))) {
    ww <- which(xout[, 4] == i)
    xx <- xout[ww,]
    xxout <- GrubbsOutliers(data = xx, alpha = 0.01)
    fdata <- matrix(rbind(fdata, xxout[[1]]), ncol = ncol(fdata))
    fouts <- c(fouts, xxout[[2]])
  }
  
  
  outflags <- rep(1, nrow(dd))
  infstat <- as.numeric(infl.stats[-1, 1])
  if (length(infstat) > 0) {
    outflags[unique(c(fouts, infstat))] <- 4
  } else {
    outflags[unique(fouts)] <- 4
  }
  
  estimated_residuals <-
    denoiser(data = xout1, method = den.method)[, 5]
  d1 <-
    data.frame(matrix(as.numeric(dd[, 1:2]), ncol = 2))
  colnames(d1) <- c("x", "y")
  centroid <- data$centroids
  centroid <-
    data.frame(matrix(matrix(centroid[order(centroid[, 1]),], ncol = 3)[, 2:3], ncol =
                        2))
  colnames(centroid) <- c("x", "y")
  
  corrdata1 <-
    matrix(cbind(data$Progression, data$Updated.groups),
           nrow = nrow(data$Progression))
  corrdata1 <-
    matrix(corrdata1[sort.list(corrdata1[, 1]),], ncol = ncol(corrdata1))
  corrdata1 <- diff(corrdata1[, 3])
  data$CPs <- which(corrdata1 != 0) + 1
  
  if (savePlot != "OFF" & savePlot != "screen") {
    pdf(paste(
      savePlot,
      "/Fluo_ordering_progression_",
      data$dateIndex,
      ".pdf",
      sep = ""
    ))
  }
  
  d1.plot <-
    ggplot(d1) + geom_point(aes_string(
      x = 'x',
      y = 'y',
      colour = as.factor(data$Updated.groups)
    )) +
    geom_point(
      data = centroid,
      aes_string(x = 'x', y = 'y'),
      shape = 3,
      size = 5
    ) +
    geom_text(
      data = centroid,
      aes_string(x = 'x', y = 'y'),
      label = c(1:nrow(centroid)),
      size = 7,
      hjust = -0.1,
      vjust = -0.1
    ) +
    labs(
      x = paste(data$VSmethod, " corrected ", data$image.type[2], " signal", sep =
                  ""),
      y = paste(data$VSmethod, " corrected ", data$image.type[3], " signal", sep =
                  "")
    ) + theme(legend.position = "none") +
    ggtitle(
      paste(
        "Estimated cell progression groups by ",
        data$CPmethod,
        " at a = ",
        data$CPsig,
        sep = ""
      )
    ) +
    theme(plot.title = element_text(size = 20))
  if (sum(outflags) > length(outflags)) {
    outpts <-
      data.frame(x = d1[which(outflags > 1), 1], y = d1[which(outflags > 1), 2])
    d1.plot <-
      d1.plot + geom_text(data = outpts,
                          aes_string(x = 'x', y = 'y'),
                          label = which(outflags > 1))
  }
  dlogratio.plot <-
    ggplot(data.frame(x = data$Progression[, 1], y = data$Progression[, 2])) +
    geom_point(aes_string(
      x = 'x',
      y = 'y',
      colour = as.factor(data$Updated.groups)
    )) +
    geom_vline(xintercept = data$CPs, colour = "black") +
    geom_hline(yintercept = 0,
               colour = "black",
               linetype = "dotted") +
    theme(legend.position = "none") +
    ggtitle(paste(
      "Estimated cell progression by ",
      data$CPmethod,
      " at a = ",
      data$CPsig,
      sep = ""
    )) +
    labs(
      x = "reconstructed time index",
      y = paste(
        data$VSmethod,
        "(",
        data$image.type[2],
        ") - ",
        data$VSmethod,
        "(",
        data$image.type[3],
        ") (corrected)",
        sep = ""
      )
    ) + theme(legend.position = "none") +
    theme(plot.title = element_text(size = 20))
  
  if (sum(outflags) > length(outflags)) {
    outpts2 <-
      data.frame(x = data$Progression[which(outflags > 1), 1], y = data$Progression[which(outflags >
                                                                                            1), 2])
    dlogratio.plot <-
      dlogratio.plot + geom_text(data = outpts2,
                                 aes_string(x = 'x', y = 'y'),
                                 label = which(outflags > 1))
  }
  multiplot(d1.plot, dlogratio.plot)
  if (savePlot != "OFF" & savePlot != "screen") {
    dev.off()
  }
  
  outflags[outflags == 4] <- "outlier"
  outflags[outflags != "outlier"] <- "normal"
  
  return(c(
    data,
    list(Outliers = outflags, Residuals = estimated_residuals)
  ))
}


#' diagnoseResiduals
#'
#' An internal function to perform residual diagnostic tests.
#'
#' @param data List. A list of fluorescence signals with their model estimates.
#' @param savePlot Character string. The directory to store the plots.
#'
#' @import ggplot2 grDevices stats
#' @importFrom raster density
#'
#' @return The residual diagnostics and plots
#'
#' @keywords internal
diagnoseResiduals <- function(data, savePlot = "OFF") {
  if (savePlot != "OFF" & savePlot != "screen") {
    pdf(paste(
      savePlot,
      "/Fluo_ordering_diagnostics_",
      data$dateIndex,
      ".pdf",
      sep = ""
    ))
  }
  
  hist.plot <-
    ggplot(data.frame(x = data$Residuals), aes_string(x = 'x')) +
    geom_histogram(bins = 50,
                   colour = "black",
                   fill = "grey") +
    labs(x = "standardized residuals") + ggtitle("Residual diagnostics") +
    theme(plot.title = element_text(size = 20))
  
  res.plot <-
    ggplot(data.frame(
      x = 1:length(data$Residuals),
      y = data$Residuals[sort.list(data$Progression[, 1])]
    ),
    aes_string(x = 'x', y = 'y')) +
    geom_point(shape = 1) + geom_hline(yintercept = 0,
                                       colour = "black",
                                       linetype = "dotted") +
    labs(x = "Reconstructed time index", y = "standardized residuals") + ggtitle("Residual diagnostics") +
    theme(plot.title = element_text(size = 20))
  
  multiplot(hist.plot, res.plot)
  if (savePlot != "OFF" & savePlot != "screen") {
    dev.off()
  }
  
  ago <- round(agostino.test(data$Residuals)$p.value, 3)
  bon <- round(bonett.test(data$Residuals)$p.value, 3)
  jar <- round(jarque.test(data$Residuals)$p.value, 3)
  ks <-
    round(ks.test(data$Residuals, "pnorm", 0, sqrt(var(data$Residuals)))$p.value, 3)
  sk <- round(skewness(data$Residuals), 3)
  ku <- round(kurtosis(data$Residuals), 3)
  legR <-
    c(
      "Skewness (ideal=0)",
      "Kurtosis (ideal=3)",
      "Agostino test P-value for Skewness",
      "Bonett test P-value for Kurtosis",
      "Jarque test P-value for Normality",
      "KS-test P-value for Normality"
    )
  report <-
    matrix(cbind(legR, c(sk, ku, ago, bon, jar, ks)), ncol = 2)
  return(c(data, list(Residuals_diagnostics = report)))
}

#' ddhft.np.2
#'
#' The original DDHF function (Motakis et al, 2006).
#'
#' @param data Numeric vector. A vector of data exhibiting monotonically increasing
#'   mean-variance relationship. The data will be transformed.
#'
#' @return The DDHF transformed data
#'
#' @keywords internal
ddhft.np.2 <- function (data) {
  n <- length(data)
  nhalf <- n / 2
  J <- logb(n, 2)
  hft <- data
  factors <- rep(0, n)
  sm <- rep(0, nhalf)
  det <- sm
  odd <- data[seq(from = 1,
                  by = 2,
                  length = nhalf)]
  even <- data[seq(from = 2,
                   by = 2,
                   length = nhalf)]
  det <- odd - even
  sigma2 <- 1 / 2 * det ^ 2
  mu <- (odd + even) / 2
  ord.mu <- order(mu)
  mu <- mu[ord.mu]
  sigma2 <- sigma2[ord.mu]
  sigma <- sqrt(isotone(sigma2, increasing = TRUE))
  vv <- ll <- list()
  for (i in 1:J) {
    sm[1:nhalf] <- (hft[2 * (1:nhalf) - 1] + hft[2 * (1:nhalf)]) / 2
    det[1:nhalf] <-
      (hft[2 * (1:nhalf) - 1] - hft[2 * (1:nhalf)]) / 2
    v <- function.from.vector(mu, sigma, sm[1:nhalf])
    ll[[i]] <- sm[1:nhalf]
    vv[[i]] <- v
    det[v > 0] <- det[v > 0] / v[v > 0]
    hft[1:nhalf] <- sm[1:nhalf]
    hft[(nhalf + 1):n] <- det[1:nhalf]
    factors[(nhalf + 1):n] <- v
    n <- n / 2
    nhalf <- n / 2
    sm <- 0
    det <- 0
  }
  nhalf <- 1
  n <- 2
  for (i in 1:J) {
    sm[1:nhalf] <- hft[1:nhalf]
    det[1:nhalf] <- hft[(nhalf + 1):n]
    hft[2 * (1:nhalf) - 1] <- sm[1:nhalf] + det[1:nhalf]
    hft[2 * (1:nhalf)] <- sm[1:nhalf] - det[1:nhalf]
    nhalf <- n
    n <- 2 * n
  }
  return(hft)
}


#' revDDHFinput
#'
#' It reverts the DDHF sorted fluorescence signals into the original sorting.
#'
#' @param data Numeric Data matrix. A data matrix of DDHF transformed data.
#'
#' @return The reverted data
#'
#' @keywords internal
revDDHFinput <- function(data, hft) {
  new <- matrix(0, nrow(data), 2)
  for (i in 1:nrow(hft[[2]])) {
    wh <- which(data[, 5] == hft[[2]][i, 1])
    if (as.numeric(hft[[2]][i, 2]) == 2) {
      new[wh, 1] <-
        hft[[1]][as.numeric(hft[[2]][i, 3]):as.numeric(hft[[2]][i, 4])]
    }
    if (as.numeric(hft[[2]][i, 2]) == 3) {
      new[wh, 2] <-
        hft[[1]][as.numeric(hft[[2]][i, 3]):as.numeric(hft[[2]][i, 4])]
    }
  }
  return(new)
}


#' DDHFinput
#'
#' It sorts the fluorescence signals of both channels data for DDHF.
#'
#' @param data Data matrix. A data matrix of fluorescence signals.
#' @param ms Data matrix. A matrix of estimated cluster centroids.
#'
#' @return The sorted fluorescence signals
#'
#' @keywords internal
DDHFinput <- function(data, ms) {
  sl <- sort.list(ms[, 2:3])
  datax <- c()
  slIndex <- matrix(0, length(sl), 4)
  for (i in 1:length(sl)) {
    llstart <- length(datax) + 1
    if (sl[i] <= max(as.numeric(ms[, 1]))) {
      datax <-
        c(datax, data[which(as.numeric(data[, 5]) == ms[sl[i], 1]), 3])
      llend <- length(datax)
      slIndex[i,] <- c(ms[sl[i], 1], 2, llstart, llend)
    }
    if (sl[i] > max(as.numeric(ms[, 1]))) {
      datax <-
        c(datax, data[which(as.numeric(data[, 5]) == ms[(sl[i] - max(as.numeric(ms[, 1]))), 1]), 4])
      llend <- length(datax)
      slIndex[i,] <-
        c(ms[(sl[i] - max(as.numeric(ms[, 1]))), 1], 3, llstart, llend)
    }
  }
  return(list(as.numeric(datax), slIndex))
}

#' isotone
#'
#' The original function to perform isotone regression (Motakis et al 2006).
#'
#' @param x Numeric vector. A vector of appropriately sorted data.
#' @param ... Other parameters.
#'
#' @return The isotone regression model estimates
#'
#' @keywords internal
isotone <- function (x,
                     wt = rep(1, length(x)),
                     increasing = TRUE) {
  nn <- length(x)
  if (nn == 1) {
    return(x)
  }
  if (!increasing) {
    x <- -x
  }
  ip <- (1:nn)
  dx <- diff(x)
  nx <- length(x)
  while ((nx > 1) && (min(dx) < 0)) {
    #cat("diff: ", dx, "\n")
    jmax <- (1:nx)[c(dx <= 0, FALSE) & c(TRUE, dx > 0)]
    jmin <- (1:nx)[c(dx > 0, TRUE) & c(FALSE, dx <= 0)]
    #cat("jmax: ", jmax, "\n")
    #cat("jmin: ", jmin, "\n")
    for (jb in (1:length(jmax))) {
      ind <- (jmax[jb]:jmin[jb])
      wtn <- sum(wt[ind])
      x[jmax[jb]] <- sum(wt[ind] * x[ind]) / wtn
      wt[jmax[jb]] <- wtn
      x[(jmax[jb] + 1):jmin[jb]] <- NA
    }
    ind <- !is.na(x)
    x <- x[ind]
    wt <- wt[ind]
    ip <- ip[ind]
    dx <- diff(x)
    nx <- length(x)
  }
  #print(x)
  #print(ip)
  jj <- rep(0, nn)
  jj[ip] <- 1
  z <- x[cumsum(jj)]
  
  if (!increasing) {
    z <- -z
  }
  return(z)
}

#' function.from.vector
#'
#' A helper for DDHF.
#'
#' @param x,y,argument.vector Appropriate vectors for analysis.
#'
#' @return Preliminary DDHF results
#'
#' @keywords internal
function.from.vector <- function(x, y, argument.vect) {
  indices <- sapply(argument.vect, which.min.diff, x)
  return(y[indices])
}


#' which.min.diff
#'
#' A helper for DDHF
#'
#' @param a,vector Appropriate vectors for analysis
#'
#' @return Preliminary DDHF results
#'
#' @keywords internal
which.min.diff <- function(a, vect) {
  return(which.min(abs(a - vect)))
}

#' GrubbsOutliers
#'
#' It identifies potential outliers by the Grubbs test.
#'
#' @param data Data matrix. A data matrix of fluorescence signals and model residuals.
#' @param alpha Float. A significance level for the grubbs test .
#'
#' @import outliers stats
#'
#' @return The fluorescence signals and the potential outliers
#'
#' @keywords internal
GrubbsOutliers <- function(data, alpha) {
  pp <- 0
  outliers <- matrix(0, 1, ncol(data))
  while (pp < alpha) {
    res <- grubbs.test(data[, 5], opposite = FALSE)
    pp <- res$p.value
    val <-
      as.numeric(unlist(strsplit(unlist(
        strsplit(res[[2]], "value ")
      )[2], " is"))[1])
    ww <- which(abs(val - data[, 5]) == min(abs(val - data[, 5])))
    if (pp < alpha) {
      outliers <-
        matrix(rbind(outliers, data[ww,]), ncol = ncol(outliers))
      data[ww, 3] <- median(data[, 3])
      data[ww, 5] <- 0
    }
  }
  pp <- 0
  while (pp < alpha) {
    res <- grubbs.test(data[, 5], opposite = TRUE)
    pp <- res$p.value
    val <-
      as.numeric(unlist(strsplit(unlist(
        strsplit(res[[2]], "value ")
      )[2], " is"))[1])
    ww <- which(abs(val - data[, 5]) == min(abs(val - data[, 5])))
    if (pp < alpha) {
      outliers <-
        matrix(rbind(outliers, data[ww,]), ncol = ncol(outliers))
      data[ww, 3] <- median(data[, 3])
      data[ww, 5] <- 0
    }
  }
  data <- data[sort.list(data[, 2]),]
  if (nrow(outliers) > 1) {
    outliers <- outliers[-1, 1]
  } else {
    outliers <- c()
  }
  return(list(data, outliers))
}
