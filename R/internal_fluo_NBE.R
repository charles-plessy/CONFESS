# Fluorescent signal normalization, background correction & cell progression dynamics estimation (internal functions)

#' denoisefun
#'
#' An internal function to run data denoising.
#'
#' @param data Numeric vector. An 1-dimensional vector of pseudotime sorted transformed signal differences.
#' @param method Character string. A method to perform the denoising. One of "wavelets", "splines" or "lregr"
#'   (linear rergression).
#'
#' @return The denoised data
#'
#' @import stats
#' @importFrom wavethresh wd threshold wr
#'
#' @keywords internal
denoisefun <- function(data, method) {
  if (method != "wavelets" &
      method != "splines" & method != "lregr") {
    stop("Not an appropriate denoising method. Choose one among: wavelets, splines and lregr")
  }
  if (method == "wavelets") {
    datanew <-
      c(data, data[(length(data) - 1):1])[1:2 ^ ceiling(log(length(data), 2))]
    while (length(datanew) < 9) {
      datanew <- c(datanew, rev(datanew))
    }
    wdy1 <-
      wd(
        datanew,
        filter.number = 1,
        family = "DaubExPhase",
        type = "wavelet",
        bc =
          "symmetric",
        verbose = FALSE,
        min.scale = 0,
        precond = TRUE
      )
    filty1 <- threshold(wdy1, policy = "sure")
    filtered <- wr(filty1)
    filtered <- filtered[1:length(data)]
  }
  if (method == "splines") {
    filtered <- smooth.spline(data)[[2]]
  }
  if (method == "lregr") {
    filtered <- lm(data ~ c(1:length(data)))$fitted.values
  }
  return(filtered)
}

#' denoiser
#'
#' An internal function to run transformed data denoising that estimates the model residuals.
#'
#' @param data Data matrix. A matrix contains the sample indices, the pseudotimes, the Ch2-Ch3 transformed data
#'   and the cluster numbers.
#' @param method Character string. A method to perform the denoising. One of "wavelets", "splines" or "lregr"
#'   (linear regression).
#'
#' @return The denoised data and the model residuals
#'
#' @keywords internal
denoiser <- function(data, method) {
  ug <- unique(data[, 4])
  res <- matrix(0, 1, 4)
  for (i in 1:length(ug)) {
    dd <- data[which(data[, 4] == i),]
    dd <- dd[sort.list(dd[, 2]),]
    fdd <- denoisefun(dd[, 3], method = method)
    fdd <- matrix(cbind(dd[, 1:2], fdd, dd[, 4]), ncol = ncol(res))
    res <- matrix(rbind(res, fdd), ncol = ncol(res))
  }
  res <- res[sort.list(res[, 1]),]
  res <- res[-1,]
  resids <- data[, 3] - res[, 3]
  res <- matrix(cbind(res, resids), nrow = nrow(res))
  return(res)
}

#' BGcorrectFluo
#'
#' It performs the background correction of cell fluorescence signals.
#'
#' @param data Data matrix. A data matrix with the foreground and background raw signals
#'   from each channel.
#' @param method Character string. The type of background correction to be performed. One
#'   of "normexp" or "subtract".
#' @param old.offset Float. An offset for the background correction method.
#' @param bg Logical. If TRUE foreground - background is performed.
#'
#' @import methods
#' @importFrom limma backgroundCorrect
#'
#' @return The background corrected data
#'
#' @keywords internal
BGcorrectFluo <- function(data, method, old.offset, bg) {
  if (bg == TRUE) {
    RG <-
      new("RGList", list(
        R = data[, 1],
        G = data[, 3],
        Rb = data[, 2],
        Gb = data[, 4]
      ))
    s <-
      matrix(cbind(apply(as.matrix(data[, 2:1], ncol = 2), 1, diff), apply(as.matrix(data[, 4:3], ncol =
                                                                                       2), 1, diff)), ncol = 2)
    offset <- abs(min(apply(matrix(s, ncol = ncol(
      s
    )), 2, min))) + 1
    if (method == "subtract") {
      core <- backgroundCorrect(RG, method = method, offset = offset)
    } else {
      core <- backgroundCorrect(RG, method = method)
    }
  }
  if (bg == FALSE) {
    RG <- new("RGList", list(R = data[, 1], G = data[, 2]))
    core <-
      backgroundCorrect(RG, method = method, offset = old.offset)
  }
  core <- matrix(cbind(core$R, core$G), ncol = 2)
  return(list(core, offset))
}

#' boxcoxMatrix
#'
#' A helper to run the Box Cox transformation on a data matrix of adjusted fluorescence signals.
#'
#' @param data Numeric vector. The adjusted signals of one channel.
#'
#' @return The transformed data
#'
#' @import stats
#' @importFrom MASS boxcox
#'
#' @keywords internal
boxcoxMatrix <- function(data) {
  est <- boxcox(data ~ 1, lambda = seq(-4, 4, 0.01), plotit = FALSE)
  lmb <- est$x[which(est$y == max(est$y))]
  rang <- range(est$x[est$y > max(est$y) - qchisq(0.95, 1) / 2])
  if (rang[1] * rang[2] < 0) {
    data <- log(data)
    lmb <- 0
  } else {
    lmb <- mean(lmb)
    data <- ((data ^ lmb) - 1) / lmb
  }
  return(list(data, lmb))
}



#' boxcoxMatrixEst
#'
#' A helper to run the Box Cox transformation on a data matrix of adjusted fluorescence signals.
#'
#' @param data Numeric vector. The adjusted signals of one channel.
#'
#' @return The transformed data
#'
#' @keywords internal
boxcoxMatrixEst <- function(data) {
  if (data[1] == 0) {
    res <- log(data[2:length(data)])
  } else {
    res <- ((data[2:length(data)] ^ data[1]) - 1) / data[1]
  }
  return(res)
}


#' doTransform
#'
#' It transforms the adjusted fluorescence signals of a matrix.
#'
#' @param data Data matrix. The adjusted signals of both channels.
#' @param transform Character string. The type of transformation to be performed. One
#'   of "bc" (Box-Cox), "log", "log10" or "asinh".
#' @param lpar Float. The lambda parameter of the Box-Cox.
#'
#' @return The transformed data matrix
#'
#' @keywords internal
doTransform <- function(data, transformation, lpar = NULL) {
  dat <- c()
  if (transformation == "bc") {
    if (is.null(lpar)) {
      est <- apply(matrix(data, ncol = ncol(data)), 2, boxcoxMatrix)
      dat <-
        matrix(cbind(matrix(est[[1]][[1]], ncol = 1), matrix(est[[2]][[1]], ncol =
                                                               1)), ncol = 2)
      lpar <- c(est[[1]][[2]], est[[2]][[2]])
    } else {
      dat <-
        apply(matrix(rbind(lpar, data), ncol = 2), 2, boxcoxMatrixEst)
    }
  }
  if (transformation == "log") {
    dat <- matrix(log(data), ncol = 2)
    lpar <- c()
  }
  if (transformation == "log10") {
    dat <- matrix(log(data, 10), ncol = 2)
    lpar <- c()
  }
  if (transformation == "asinh") {
    dat <- matrix(asinh(data), ncol = 2)
    lpar <- c()
  }
  if (transformation == "none") {
    dat <- data
    lpar <- c()
  }
  return(list(dat, lpar))
}


#' invTransform
#'
#' It back-transforms the transformed adjusted cell fluorescence signals of a matrix.
#'
#' @param data Data matrix. The adjusted signals of both channels.
#' @param transform Character string. The type of transformation to be performed. One
#'   of "bc" (Box-Cox), "log", "log10" or "asinh".
#' @param lambda Float. The lambda parameter of the Box-Cox.
#'
#' @return The back-transformed data matrix
#'
#' @keywords internal
invTransform <- function(data, lambda, transformation) {
  res <- c()
  if (transformation == "bc") {
    if (lambda == 0) {
      res <- exp(data)
    } else {
      res <- (lambda * data + 1) ^ (1 / lambda)
    }
  }
  if (transformation == "log10") {
    res <- 10 ^ data
  }
  if (transformation == "log") {
    res <- exp(data)
  }
  if (transformation == "asinh") {
    res <- sinh(data)
  }
  return(res)
}

#' boxFluo
#'
#' It generates the density plots of the uncorrected and corrected cell fluorescence signals.
#'
#' @param data List. A list with the fluorescence data. Typically, the output of createFluo().
#' @param transform Character string. The type of transformation to be performed. One
#'   of "bc" (Box-Cox), "log", "log10" or "asinh".
#' @param legends Character vector. Puts the "uncorrected" or "corrected" legends on the signal density
#'   plots.
#' @param reference Integer. The run number to be used as baseline for the run correction.
#' @param batchnames original run name.
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param savePlot Character. A switch to generate the density plots.
#'
#' @return The density plots of the fluorescence data
#'
#' @import ggplot2 stats
#'
#' @importFrom raster hist
#'
#'
#' @keywords internal
boxFluo <-
  function(data,
           transformation,
           reference,
           legends,
           batchnames,
           image.type,
           savePlot) {
    edata <- data$exprs
    f <- factor(data$batch)
    f <- as.numeric(levels(relevel(f, ref = reference)))
    d <- as.list(rep(0, (2 * max(f))))
    est <-
      doTransform(data = edata, transformation = transformation)
    dat <- est[[1]]
    lpar <- est[[2]]
    d[[f[1]]] <- dat[which(data$batch == f[1]), 1]
    for (i in 2:(2 * max(f))) {
      if (i <= max(f)) {
        d[[f[i]]] <- dat[which(data$batch == f[i]), 1]
      } else {
        d[[(f[(i - max(f))] + max(f))]] <-
          dat[which(data$batch == f[(i - max(f))]), 2]
      }
    }
    legs <-
      c(
        paste(image.type[2], ": ", max(f), " Runs", sep = ""),
        paste(image.type[3], ": ", max(f), " Runs", sep = "")
      )
    densX <- as.list(rep(0, length(d)))
    densY <- as.list(rep(0, length(d)))
    for (i in 1:length(d)) {
      densX[[i]] <- density(d[[i]])$x
      densY[[i]] <- density(d[[i]])$y
    }
    #  runs<-unique(unlist(strsplit(batchnames,separator))[c(TRUE,FALSE)])
    runs <- unique(batchnames)
    adX <- unlist(densX)
    adY <- unlist(densY)
    dens.m <- melt(d)
    dens.ch2 <-
      dens.m[dens.m$L1 <= max(f),]
    dens.ch2$L1 <- as.factor(dens.ch2$L1)
    dens.ch3 <-
      dens.m[dens.m$L1 > max(f),]
    dens.ch3$L1 <- as.factor(dens.ch3$L1)
    
    plot1 <- plot2 <- "Plots are not produced when savePlot = OFF!"
    
    
    if (savePlot != "OFF") {
      plot1 <- ggplot(dens.ch2, aes_string(x = "value", colour = 'L1')) +
        geom_line(stat = "density") +
        scale_colour_discrete(name = "Runs", labels = runs) +
        theme(legend.title = element_text(size = 12, face = "bold")) + theme_bw() +
        labs(title = paste(legends, legs[1]))
      
      plot2 <-
        ggplot(dens.ch3, aes_string(x = "value", colour = 'L1')) +
        geom_line(stat = "density") +
        scale_colour_discrete(name = "Runs", labels = runs) +
        theme(legend.title = element_text(size = 12, face = "bold")) + theme_bw() +
        labs(title = paste(legends, legs[2]))
    }
    
    
    return(list(dat, data$batch, lpar, p1 = plot1, p2 = plot2))
  }


#' refineMixes
#'
#' An helper internal function to generate the results of flexmix.
#'
#' @param data Numeric vector. The background corrected fluorescence signals of
#'   a single channel.
#' @param batch Integer. The run number.
#' @param model Model. The output of getModel() from flexmix.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @return The flexmix mixture components
#' @importFrom stats aggregate
#' @importFrom flexmix clusters
#'
#' @keywords internal
refineMixes <- function(data, batch, model, seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cc1 <- clusters(model)
  cc <- rep(0, length(cc1))
  for (i in 1:max(batch)) {
    w <- which(batch == i)
    cl <- cc1[w]
    d <- data[w]
    sl <- sort.list(aggregate(d, list(cl), mean)[, 2])
    cl1 <- rep(0, length(cl))
    for (j in 1:length(sl)) {
      cl1[which(cl == sl[j])] <- j
    }
    cc[w] <- cl1
  }
  return(cc)
}


#' lmFluo
#'
#' It estimates the optimal number of mixtures for the flexmix model on data from multiple runs.
#'
#' @param data Numeric vector. An 1-dimensional vector of adjusted data from a single channel.
#' @param batch Integer. The run number.
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model.
#' @param reference Numeric vector. Specifies the runs to be used as baseline (iteratively).
#' @param prior.pi Float. The prior probability to accept a component.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate
#'   the flexmix model.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL".
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @return The flexmix mixture components and other statistics
#'
#' @importFrom flexmix getModel stepFlexmix
#'
#' @keywords internal
lmFluo <-
  function(data,
           batch,
           maxMix,
           reference,
           prior.pi,
           flex.reps,
           flexmethod,
           seed) {
    Batch <- factor(batch)
    if (length(levels(Batch)) < 2) {
      stop("This function is used only for data with multiple runs")
    }
    Batch <- relevel(Batch, ref = reference)
    options(warn = -1)
    try <-
      capture.output(stepFlexmix(
        data ~ 1,
        k = 1,
        nrep = 2,
        control = list(minprior = prior.pi)
      ))
    
    if (maxMix < 1) {
      maxMix <- 1
      print(paste("maxMix < 1 is not supported. maxMix was set to 1", sep =
                    ""))
    }
    if (flexmethod != "BIC" & flexmethod != "AIC" &
        flexmethod != "ICL") {
      stop("This flexmix method is not supported")
    }
    
    if (maxMix > 1) {
      if (!is.null(seed)) {
        set.seed(min(seed * 37, 100000000))
      }
      fmod <-
        stepFlexmix(
          data ~ Batch,
          k = 1:maxMix,
          nrep = flex.reps,
          control = list(minprior = prior.pi)
        )
      if (!is.null(seed)) {
        set.seed(min(seed * 3, 100000000))
      }
      
      gg <- flexmix::getModel(fmod, which = flexmethod)
      Comp <- refineMixes(data, batch, gg, seed)
      w1 <- c()
      
      if (length(unique(batch)) > 1 & length(unique(Comp)) > 1) {
        fmod <- lm(data ~ factor(Comp) * factor(Batch))
        desMat <- model.matrix( ~ factor(Comp) * factor(Batch))
        rr <- rownames(summary(fmod)[[4]])
        w1 <-
          which(rr == paste("factor(Comp)", max(Comp), sep = ""))
      }
      if (length(unique(batch)) > 1 & length(unique(Comp)) == 1) {
        fmod <- lm(data ~ factor(Batch))
        desMat <- model.matrix( ~ factor(Batch))
        rr <- rownames(summary(fmod)[[4]])
        w1 <- which(rr == "(Intercept)")
      }
    }
    if (maxMix == 1) {
      fmod <- lm(data ~ factor(Batch))
      desMat <- model.matrix( ~ factor(Batch))
      rr <- rownames(summary(fmod)[[4]])
      w1 <- which(rr == "(Intercept)")
      Comp <- rep(1, length(data))
    }
    
    rr <- rownames(summary(fmod)[[4]])
    cc <- c(colnames(summary(fmod)[[4]]), "FDR")
    compEsts <-
      matrix(as.numeric(summary(fmod)[[4]][, c(1, 4)]), ncol = 2)
    compEsts <-
      matrix(cbind(compEsts, matrix(p.adjust(compEsts[, ncol(compEsts)], "BH"), ncol =
                                      1)), nrow = nrow(compEsts))
    resids <- fmod$residuals
    stdresids <- rstandard(fmod)
    fitted <- fmod$fitted.values
    compEsts <-
      matrix(cbind(rr, compEsts), ncol = (ncol(compEsts) + 1))
    compEsts <-
      matrix(rbind(c("", cc[c(1, 4, 5)]), compEsts), ncol = ncol(compEsts))
    ww <- as.list(rep(0, max(Comp)))
    for (i in 1:max(Comp)) {
      ww[[i]] <- which(Comp == i)
    }
    wall <- length(rr)
    subtInd <- c()
    if (length(w1) > 0) {
      subtInd <- (w1 + 2):(wall + 1)
    }
    return(list(compEsts, Comp, ww, resids, stdresids, fitted, subtInd, desMat))
  }


#' BatchFluo
#'
#' It performs the run effect correction of the cell fluorescence signals by flexmix or 2-way ANOVA.
#'
#' @param data List. A list of transformed and adjusted fluorescence signals.
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model. If maxMix=1 or if the the optimal number of the estimated components
#'   is 1, the model reduces to the classical 2-way ANOVA.
#' @param reference Numeric vector. Specifies the runs to be used as baseline (iteratively).
#' @param prior.pi Float. The prior probability to accept a component.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate
#'   the flexmix model.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL".
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @return A list of fluorescence signals, mixture components and flexmix model estimates
#'
#' @keywords internal
BatchFluo <-
  function(data,
           maxMix,
           reference,
           prior.pi,
           flex.reps,
           flexmethod,
           seed) {
    edata <- data[[1]]
    om <- apply(matrix(edata, ncol = ncol(edata)), 2, mean)
    ests <-
      apply(
        matrix(edata, ncol = ncol(edata)),
        2,
        lmFluo,
        batch = data[[2]],
        maxMix = maxMix,
        reference = reference,
        prior.pi =
          prior.pi,
        flex.reps = flex.reps,
        flexmethod = flexmethod,
        seed = seed
      )
    coefs <-
      list(as.numeric(ests[[1]][[1]][ests[[1]][[7]], 2]), as.numeric(ests[[2]][[1]][ests[[2]][[7]], 2]))
    new <- resids <- stdresids <- fitted <- edata
    for (i in 1:2) {
      resids[, i] <- as.numeric(ests[[i]][[4]])
      stdresids[, i] <- as.numeric(ests[[i]][[5]])
      fitted[, i] <- as.numeric(ests[[i]][[6]])
      if (length(ests[[i]][[7]]) > 0) {
        for (b in 1:length(ests[[i]][[7]])) {
          new[, i] <-
            new[, i] - coefs[[i]][b] * ests[[i]][[8]][, (ests[[i]][[7]][b] - 1)]
        }
      }
      new[, i] <- new[, i] - mean(new[, i]) + om[i]
    }
    return(
      list(
        exprs = new,
        batch = data[[2]],
        mixesCh2 = ests[[1]][[2]],
        mixesCh3 = ests[[2]][[2]],
        residuals =
          resids,
        standardized.residuals = stdresids,
        fitted.values = fitted,
        estCh2 =
          ests[[1]][[1]],
        estCh3 = ests[[2]][[1]],
        designCh2 = ests[[1]][[8]],
        designCh3 =
          ests[[2]][[8]]
      )
    )
  }

#' fclust
#'
#' The main flowClust function used in this application.
#'
#' @param data Data matrix. A matrix of run effect and background corrected fluorescence signals
#'   in both channels.
#' @param k Integer. The maximum number of clusters to be generated.
#' @param nstart Integer. A flowClust parameter specifing the number of random sets to
#'   be chosen for the clustering estimation.
#'
#' @return The flowClust results
#'
#' @keywords internal
fclust <- function(data, k, nstart = 1) {
  ff <- flowclust_step1(data, k, nstart)
  return(flowclust_step2(ff))
}

#' flowclust_step1
#'
#' A helper function for flowClust analysis.
#'
#' @param data Data matrix. A matrix of run effect and background corrected fluorescence signals
#'   in both channels.
#' @param k Integer. The maximum number of clusters to be generated.
#' @param nstart Integer. A flowClust parameter specifing the number of random sets to
#'   be chosen for the clustering estimation.
#'
#' @return Preliminary flowClust results
#'
#' @keywords internal
flowclust_step1 <- function(data, k, nstart) {
  dd <- data.frame(data)
  colnames(dd) <- c("Ch2", "Ch3")
  
  r1 <-
    tryCatch({
      res1 <-
        suppressWarnings(flowClust(dd, varNames = colnames(dd), K = k))
    },
    error = function(e) {
      stop("flowClust failed to estimate the clusters. Use a smaller k.max parameter")
    })
  r2 <- tryCatch({
    res2 <- suppressWarnings(split(dd, r1))
  },
  error = function(e) {
    stop("flowClust failed to estimate the clusters. Use a smaller k.max parameter")
  })
  return(r2)
}


#' flowclust_step2
#'
#' Another helper function for flowClust analysis.
#'
#' @param data Data matrix. A matrix of flowclust results from flowclust_step1().
#'
#' @return Preliminary flowClust results
#'
#' @keywords internal
flowclust_step2 <- function(data) {
  n <- lapply(data, rownames)
  gg <- rep(0, length(unlist(n)))
  for (i in 1:length(data)) {
    gg[as.numeric(n[[i]])] <- i
  }
  return(list(cluster = gg))
}


#' adjustFluo
#'
#' It performs the run effect correction and normalization of the raw fluorescence signals.
#'
#' @param data List. A list with the fluorescence signal information of both channels.
#' @param transform Character string. The type of transformation to be performed. One
#'   of "bc" (Box-Cox), "log", "log10" or "asinh".
#' @param BGmethod Character string. The type of image background correction to be performed.
#'   One of "normexp" or "subtract".
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model. If maxMix=1 or if the the optimal number of the estimated components
#'   is 1, the model reduces to the classical 2-way ANOVA.
#' @param reference Numeric vector. Specifies the runs to be used as baseline (iteratively).
#' @param prior.pi Float. The prior probability to accept a component.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate
#'   the flexmix model.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL".
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param savePlot Character string. The directory to store the plots or an option to print them on the screen.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#' @param dateIndex Character string. a date index to be used in saving the output files.
#'
#' @import moments ggplot2 grDevices graphics stats
#'
#' @return A list with the fluorescence signals, mixture components and flexmix model estimates
#'
#' @keywords internal
adjustFluo <-
  function(data,
           transformation,
           BGmethod,
           maxMix,
           reference,
           prior.pi,
           flex.reps,
           flexmethod,
           image.type,
           savePlot,
           seed,
           dateIndex) {
    if (reference > max(as.numeric(data$batch[, 2]))) {
      reference <- max(as.numeric(data$batch[, 2]))
      message(
        paste(
          "the baseline dataset for batch effect correction has changed to ",
          reference,
          sep = ""
        )
      )
    }
    
    if (savePlot != "OFF" & savePlot != "screen") {
      pdf(
        paste(
          savePlot,
          "/Fluo_adjustment_densities_batch=",
          reference,
          "_",
          dateIndex,
          ".pdf",
          sep = ""
        ),
        paper = "a4r"
      )
    }
    
    dd <-
      BGcorrectFluo(
        data$RGexprs,
        method = "subtract",
        old.offset = c(),
        bg = TRUE
      )
    dd <-
      list(exprs = dd[[1]],
           batch = as.numeric(data$batch[, 2]),
           offset = dd[[2]])
    res1 <-
      boxFluo(
        dd,
        transformation = transformation,
        reference = reference,
        legends = "uncorrected",
        batchnames = data$batch[, 1],
        image.type = image.type,
        savePlot = savePlot
      )
    res2 <-
      BatchFluo(
        res1,
        maxMix = maxMix,
        reference = reference,
        prior.pi = prior.pi,
        flexmethod =
          flexmethod,
        flex.reps = flex.reps,
        seed = seed
      )
    
    if (BGmethod == "normexp") {
      res2$exprs <-
        BGcorrectFluo(
          matrix(cbind(
            invTransform(
              res2$exprs[, 1],
              lambda = res1[[3]][1],
              transformation = transformation
            ),
            invTransform(
              res2$exprs[, 2],
              lambda =
                res1[[3]][2],
              transformation = transformation
            )
          ), ncol = 2),
          method = "normexp",
          old.offset = dd$offset,
          bg = FALSE
        )[[1]]
    } else {
      res2$exprs <- res2$exprs + dd$offset
    }
    
    res3 <-
      boxFluo(
        res2,
        transformation = transformation,
        reference = reference,
        legends = "corrected",
        batchnames = data$batch[, 1],
        image.type = image.type,
        savePlot = savePlot
      )
    
    if (savePlot != "OFF") {
      suppressWarnings(multiplot(
        res1$p1,
        res1$p2,
        res3$p1,
        res3$p2,
        layout = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
      ))
    }
    
    if (savePlot != "OFF" & savePlot != "screen") {
      dev.off()
    } else if (savePlot == "screen") {
      dev.new()
    }
    
    r <- res2$residuals
    sr <- res2$standardized.residuals
    f <- res2$fitted.values
    
    if (savePlot != "OFF" & savePlot != "screen") {
      pdf(
        paste(
          savePlot,
          "/Fluo_adjustment_diagnostics_batch=",
          reference,
          "_",
          dateIndex,
          ".pdf",
          sep = ""
        ),
        paper = "a4r"
      )
    }
    
    
    if (savePlot != "OFF") {
      par(mfrow = c(3, 3))
      options(warn = -1)
      hist(
        c(r),
        breaks = 50,
        xlab = "standardized residuals",
        main = "Residual diagnostics",
        sub = paste("KS-test for normality: ", round(
          ks.test(c(r), "pnorm", 0, sqrt(var(c(
            r
          ))))$p.value, 3
        ), sep = "")
      )
      hist(
        r[, 1],
        breaks = 50,
        xlab = "standardized residuals",
        main = paste("Residuals of ", image.type[2], " (corrected)", sep = ""),
        sub = paste("KS-test for normality: ", round(
          ks.test(r[, 1], "pnorm", 0, sqrt(var(r[, 1])))$p.value, 3
        ), sep = "")
      )
      hist(
        r[, 2],
        breaks = 50,
        xlab = "standardized residuals",
        main = paste("Residuals of ", image.type[3], " (corrected)", sep = ""),
        sub = paste("KS-test for normality: ", round(
          ks.test(r[, 2], "pnorm", 0, sqrt(var(r[, 2])))$p.value, 3
        ), sep = "")
      )
      plot(c(f),
           c(sr),
           xlab = "Model fitted values",
           ylab = "Model Standardized Residuals",
           main = "")
      plot(
        f[, 1],
        sr[, 1],
        xlab = paste(image.type[2], " fitted values", sep = ""),
        ylab = paste(image.type[2], " Standardized Residuals", sep = ""),
        main = ""
      )
      plot(
        f[, 2],
        sr[, 2],
        xlab = paste(image.type[2], " fitted values", sep = ""),
        ylab = paste(image.type[2], " Standardized Residuals", sep = ""),
        main = ""
      )
      
      cr <- c(r)
      cf <- c(f)
      acf(cr[sort.list(cf)], main = "Autocorrelation of model standardized residuals")
      acf(r[sort.list(f[, 1]), 1],
          main = paste(
            "Autocorrelation of ",
            image.type[2],
            " standardized residuals",
            sep = ""
          ))
      acf(r[sort.list(f[, 2]), 2],
          main = paste(
            "Autocorrelation of ",
            image.type[3],
            " standardized residuals",
            sep = ""
          ))
      
    }
    
    
    if (savePlot != "OFF" & savePlot != "screen") {
      dev.off()
    } else if (savePlot == "screen") {
      dev.new()
    }
    
    if (savePlot != "OFF" & savePlot != "screen") {
      pdf(
        paste(
          savePlot,
          "/Fluo_adjustment_corrected2D_batch=",
          reference,
          "_",
          dateIndex,
          ".pdf",
          sep = ""
        ),
        paper = "a4r"
      )
    }
    
    
    
    if (savePlot != "OFF") {
      plot(
        res3[[1]],
        cex = 0.7,
        xlab = paste(transformation, " corrected ", image.type[2], " signal", sep =
                       ""),
        ylab = paste(transformation, " corrected ", image.type[3], " signal", sep =
                       ""),
        main = ""
      )
      if (savePlot != "screen")
        dev.off()
    }
    
    ago <- c(
      round(agostino.test(r[, 1])$p.value, 3),
      round(agostino.test(r[, 2])$p.value, 3),
      round(agostino.test(c(r))$p.value, 3)
    )
    bon <- c(
      round(bonett.test(r[, 1])$p.value, 3),
      round(bonett.test(r[, 2])$p.value, 3),
      round(bonett.test(c(r))$p.value, 3)
    )
    jar <- c(
      round(jarque.test(r[, 1])$p.value, 3),
      round(jarque.test(r[, 2])$p.value, 3),
      round(jarque.test(c(r))$p.value, 3)
    )
    ks <- c(
      round(ks.test(r[, 1], "pnorm", 0, sqrt(var(
        r[, 1]
      )))$p.value, 3),
      round(ks.test(r[, 2], "pnorm", 0, sqrt(var(
        r[, 2]
      )))$p.value, 3),
      round(ks.test(c(r), "pnorm", 0, sqrt(var(
        c(r)
      )))$p.value, 3)
    )
    sk <-
      c(round(skewness(r[, 1]), 3), round(skewness(r[, 2]), 3), round(skewness(c(r)), 3))
    ku <-
      c(round(kurtosis(r[, 1]), 3), round(kurtosis(r[, 2]), 3), round(kurtosis(c(r)), 3))
    legR <- c(
      "",
      "Skewness (ideal=0)",
      "Kurtosis (ideal=3)",
      "Agostino test for Skewness",
      "Bonett test for Kurtosis",
      "Jarque test for Normality",
      "KS-test for Normality"
    )
    legC <-
      c(image.type[2],
        image.type[3],
        paste(image.type[2], " & ", image.type[3], sep = ""))
    report <-
      matrix(cbind(legC, sk, ku, ago, bon, jar, ks), nrow = 3)
    report <- matrix(rbind(legR, report), ncol = ncol(report))
    
    
    colnames(dd$exprs) <-
      colnames(res2$exprs) <- colnames(res3[[1]]) <-
      colnames(f) <- colnames(r) <- colnames(sr) <- image.type[2:3]
    extraData <- list(
      dd$exprs,
      res2$exprs,
      res3[[1]],
      res2$mixesCh2,
      res2$mixesCh3,
      res2$estCh2,
      res2$estCh3,
      f,
      transformation,
      r,
      sr,
      report,
      res1[[3]],
      res2$designCh2,
      res2$designCh3,
      reference
    )
    names(extraData) <-
      c(
        "exprs",
        "corrected.exprs",
        "corrected.transformed.exprs",
        paste("mixes.", image.type[2], sep = ""),
        paste("mixes.", image.type[3], sep = ""),
        paste("Batch.", image.type[2], ".est", sep = ""),
        paste("Batch", image.type[3], ".est", sep = ""),
        "fitted.values",
        "transformation",
        "model.residuals",
        "model.standardized.residuals",
        "residual.statistics",
        "lpar",
        paste("design.", image.type[2], sep = ""),
        paste("design.", image.type[3], sep = ""),
        "reference"
      )
    
    resALL <- c(data, extraData)
    
    return(resALL)
  }


#' contrastFluo
#'
#' It estimates the contrasts comparisons across runs and runs*component in each channel.
#'
#' @param data List. A list with the fluorescence signal information of both channels. Typically,
#'   the output of adjustFluo().
#' @param channel Character string. An identifier for the channel to be analyzed.
#'
#' @return A list with the fluorescence signals, mixture components, flexmix model estimates and
#'   contrast results
#'
#' @import stats utils
#' @importFrom contrast contrast
#'
#' @keywords internal
contrastFluo <- function(data, channel, legends) {
  dd <- data$exprs
  dd <- doTransform(dd, data$transform)[[1]]
  if (channel == "CH1") {
    d1 <-
      data.frame(
        Int = dd[, 1],
        mixes = factor(data[[11]]),
        batch = relevel(factor(as.numeric(data$batch[, 2])), ref =
                          data$reference)
      )
    a <- as.character(sort(unique(data[[11]])))
  }
  if (channel == "CH2") {
    d1 <-
      data.frame(
        Int = dd[, 2],
        mixes = factor(data[[12]]),
        batch = relevel(factor(as.numeric(data$batch[, 2])), ref =
                          data$reference)
      )
    a <- as.character(sort(unique(data[[12]])))
  }
  
  if (length(levels(d1$mixes)) > 1) {
    mod <- lm(Int ~ mixes * batch, data = d1)
    b <- c()
    for (i in 1:max(as.numeric(data$batch[, 2]))) {
      b <- c(b, as.character(i))
    }
    b <- t(combn(b, 2))
    res <-
      matrix(c(
        "Channel",
        "Component",
        "Run1",
        "Run2",
        "Contrast",
        "Pvalue"
      ),
      1,
      6)
    for (i in 1:length(a)) {
      for (j in 1:nrow(b)) {
        rr <-
          contrast(
            mod,
            a = list(mixes = a[i], batch = b[j, 1]),
            b = list(mixes = a[i], batch =
                       b[j, 2]),
            type = "average"
          )
        res <-
          matrix(rbind(res, c(
            channel, a[i], b[j,], as.numeric(unlist(rr)[c(1, 7)])
          )), ncol = ncol(res))
      }
    }
  }
  if (length(levels(d1$mixes)) == 1) {
    mod <- lm(Int ~ batch, data = d1)
    b <- c()
    for (i in 1:max(as.numeric(data$batch[, 2]))) {
      b <- c(b, as.character(i))
    }
    b <- t(combn(b, 2))
    res <-
      matrix(c(
        "Channel",
        "Component",
        "Run1",
        "Run2",
        "Contrast",
        "Pvalue"
      ),
      1,
      6)
    for (j in 1:nrow(b)) {
      rr <-
        contrast(
          mod,
          a = list(batch = b[j, 1]),
          b = list(batch = b[j, 2]),
          type = "average"
        )
      res <-
        matrix(rbind(res, c(
          channel, 1, b[j,], as.numeric(unlist(rr)[c(1, 7)])
        )), ncol = ncol(res))
    }
  }
  
  p <- p.adjust(as.numeric(res[2:nrow(res), ncol(res)]), "BH")
  res <- matrix(cbind(res, c("FDR", p)), nrow = nrow(res))
  ll <- list(res)
  names(ll) <- paste(legends, ".contrasts", sep = "")
  return(c(data, ll))
}

#' summarizeAdjFluo
#'
#' A wrapper of the functions used for run effect and background correction. It gives the the
#'   corrected, transformed corrected and mixture groups of each baseline run.
#'
#' @param data List. A list with the fluorescence signal information of both channels.
#' @param transform Character string. The type of transformation to be performed. One
#'   of "bc" (Box-Cox), "log", "log10" or "asinh".
#' @param BGmethod Character string. The type of image background correction to be performed.
#'   One of "normexp" or "subtract".
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model. If maxMix=1 or if the the optimal number of the estimated components
#'   is 1, the model reduces to the classical 2-way ANOVA.
#' @param reference Numeric vector. Specifies the runs to be used as baseline (iteratively).
#' @param prior.pi Float. The prior probability to accept a component.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate
#'   the flexmix model.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL".
#' @param image.type Character string. A triplet of IDs to characterize the type of images under study.
#' @param savePlot Character string. The directory to store the plots or an option to print them on the screen.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @return A list with the adjusted fluorescence signals
#'
#' @keywords internal
summarizeAdjFluo <-
  function(data,
           transformation,
           BGmethod,
           maxMix,
           reference,
           prior.pi,
           flex.reps,
           flexmethod,
           image.type,
           savePlot,
           seed) {
    res <- as.list(rep(0, length(reference)))
    datin <- data$dateIndex
    for (i in 1:length(reference)) {
      res[[i]] <-
        adjustFluo(
          data = data,
          transformation = transformation,
          BGmethod = BGmethod,
          maxMix = maxMix,
          reference = reference[i],
          prior.pi =
            prior.pi,
          flex.reps = flex.reps,
          flexmethod = flexmethod,
          image.type = image.type,
          savePlot = savePlot,
          seed = seed,
          dateIndex = datin
        )
      if (i == 1) {
        result <-
          list(
            corrected.exprs = res[[i]]$corrected.exprs,
            corrected.transformed.exprs = res[[i]]$corrected.transformed.exprs
          )
      }
      if (i > 1) {
        for (j in 1:length(result)) {
          if (j == 1) {
            result[[j]] <-
              matrix(cbind(result[[j]], res[[i]]$corrected.exprs),
                     nrow = nrow(res[[i]]$corrected.exprs))
          }
          if (j == 2) {
            result[[j]] <-
              matrix(
                cbind(result[[j]], res[[i]]$corrected.transformed.exprs),
                nrow = nrow(res[[i]]$corrected.transformed.exprs)
              )
          }
        }
      }
    }
    means2 <- matrix(0, nrow(result[[2]]), 2)
    means2[, 1] <-
      apply(matrix(result[[2]][, seq(1, (2 * length(reference)), 2)], nrow =
                     nrow(result[[2]])), 1, mean)
    means2[, 2] <-
      apply(matrix(result[[2]][, seq(2, (2 * length(reference)), 2)], nrow =
                     nrow(result[[2]])), 1, mean)
    for (j in 1:length(result)) {
      result[[j]] <-
        matrix(rbind(matrix(paste(
          "Ref=", rep(sort(reference), each = 2), sep = ""
        ), nrow = 1), result[[j]]), ncol = ncol(result[[j]]))
    }
    means1 <-
      invTransform(data = means2,
                   lambda = c(0, 0),
                   transformation = transformation)
    return(
      list(
        index = data$index,
        samples = data$samples,
        batch = data$batch,
        Size = data$Size,
        allResults =
          result,
        summ.corrected = means1,
        summ.corrected.transformed = means2,
        One.batch = res
      )
    )
  }



#' updateCentroids
#'
#' It updates the centroids of the clusters that are re-estimated by change-point analysis.
#'
#' @param data Data matrix. A matrix of appropriately transformed fluorescence signals .
#' @param centroidTable Data matrix. A previously estimated centroids table to be updated.
#'
#' @return The updated centroids table
#'
#' @import stats
#'
#' @keywords internal
updateCentroids <- function(data, centroidTable) {
  ce <- apply(matrix(data[, 1:2], ncol = 2), 2, median)
  ww <- which(centroidTable[, 1] == unique(data[, 3]))
  centroidTable[ww, 2:3] <- ce
  return(centroidTable)
}

#' FluoInspection
#'
#' It generates the initial clusters, their centroids and plots the results.
#'
#' @param data List. A list of fluorescence signal information for both channels.
#' @param dateIndex Character string. A date index to be used in saving the output files.
#' @param savePlot Character string. The directory to store the plots or an option to print
#'   them on the screen.
#'
#' @import grDevices stats
#' @importFrom flowClust flowClust getEstimates plot
#' @importFrom flowPeaks flowPeaks assign.flowPeaks
#'
#' @return A list with the adjusted fluorescence signals and the centroids
#'
#' @keywords internal
FluoInspection <- function(data, dateIndex, savePlot) {
  it <- colnames(data$corrected.transformed.exprs)
  dd <- data.frame(data$corrected.transformed.exprs)
  colnames(dd) <- c("Ch2", "Ch3")
  measure <- data$clusterFUN
  
  if (savePlot != "OFF" & savePlot != "screen") {
    pdf(
      paste(
        savePlot,
        "/Fluo_inspection_",
        measure,
        "_",
        dateIndex,
        ".pdf",
        sep = ""
      ),
      paper = "a4r"
    )
  }
  
  
  aa1 <- aggregate(dd$Ch2, list(data$GAPgroups[, 1]), median)
  aa2 <- aggregate(dd$Ch3, list(data$GAPgroups[, 1]), median)
  aa <- matrix(cbind(as.matrix(aa1), as.matrix(aa2)[, 2]), ncol = 3)
  colnames(aa) <- c("Cluster", it[1:2])
  
  cellprog <- data.frame(x = dd$Ch2, y = dd$Ch3)
  cptitle <-
    paste(data$clusterFUN, " cell progression clusters", sep = "")
  cellprog.plot <-
    ggplot(cellprog) + geom_point(aes_string(x = 'x', y = 'y'), shape = 1) +
    geom_text(
      data = data.frame(x = aa[, 2], y = aa[, 3]),
      aes_string(x = 'x', y = 'y'),
      label = aa[, 1],
      size = 7,
      hjust = -0.1,
      vjust = -0.1
    ) +
    ggtitle(bquote(atop(.(cptitle), atop(
      italic("(the numbers label the estimated groups)"), ""
    )))) +
    labs(
      x = paste(data$transformation, " corrected ", it[1], " signal", sep = ""),
      y = paste(data$transformation, " corrected ", it[2], " signal", sep = "")
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 20,
        face = "bold",
        colour = "black",
        vjust = -1
      )
    )
  
  
  if (savePlot != "OFF") {
    print(cellprog.plot)
  }
  
  
  if (savePlot != "OFF" & savePlot != "screen") {
    dev.off()
  }
  
  
  return(c(data, list(centroids = aa)))
}


#' FixPath
#'
#' It tests whether the path has been appropariately defined and produces an error if not.
#'
#' @param data List. A list of fluorescence signal information for both channels.
#' @param groups Numeric vector. A vector of cluster indices.
#'
#' @importFrom flowClust flowClust getEstimates plot
#' @importFrom flowPeaks flowPeaks assign.flowPeaks
#'
#' @return A list with the adjusted fluorescence signals and the clusters
#'
#' @keywords internal
fixPath <- function(data, groups) {
  if (length(groups) == 0) {
    stop("Insert the cell progression path using the labels of Fluo_inspection")
  }
  if (length(unique(groups)) != length(unique(data$GAPgroups[, 1]))) {
    stop(
      "Error in cell progression path: you have specified different number of groups than the one estimated"
    )
  }
  return(c(data, list(UpdatedPath = groups)))
}


#' orderFluo
#'
#' It sort the adjusted (and transfomed) fluorescence signals according to the path progression.
#'
#' @param data List. A list of adjusted fluorescence signals.
#' @param path.start Integer. A cluster number indicating the starting cluster that algorithm should use to
#'   build the path. The cluster numbers refer to the plot generated by Fluo_inspection(). Default is 1.
#'   If path.type = "circular" the number does not matter. If path.type = "A2Z" the user should inspect the
#'   Fluo_inspection() plot to detect the beginning of the path. If path.type = "other", the function will
#'   not estimate a path. The user has to manually insert the path progression (the cluster numbers) in
#'   Fluo_modeling()  .
#' @param updater Logical. An indicator of the estimation stage. If FALSE the initial groups are analyzed,
#'    otherwise the change-point based groups are analyzed.
#'
#' @return A list with the adjusted fluorescence signals, the centroids, the clusters and the pseudotimes
#'
#' @import stats
#' @importFrom raster text
#'
#' @keywords internal
orderFluo <- function(data, path.type, updater = FALSE) {
  if (path.type[1] != "circular" & path.type[1] != "A2Z") {
    stop("The path.type is not correctly specified")
  }
  if (updater == FALSE) {
    path <- c(data$UpdatedPath, data$UpdatedPath[1])
    ms <- data$centroids
    groups <- data$GAPgroups[, 1]
  }
  if (updater == TRUE) {
    path <- c(1:max(data$Updated.groups), 1)
    ddVS <- data.frame(data$corrected.VStransformed.exprs)
    colnames(ddVS) <- c("Ch2", "Ch3")
    aa1 <- aggregate(ddVS$Ch2, list(data$Updated.groups), median)
    aa2 <- aggregate(ddVS$Ch3, list(data$Updated.groups), median)
    ms <-
      matrix(cbind(as.matrix(aa1), as.matrix(aa2)[, 2]), ncol = 3)
    groups <- data$Updated.groups
  }
  
  mydata <- matrix(0, 1, 5)
  for (i in 1:(length(path) - 1)) {
    center1 <- as.numeric(ms[which(ms[, 1] == path[i]), 2:3])
    center2 <- as.numeric(ms[which(ms[, 1] == path[(i + 1)]), 2:3])
    ww <- which(groups == path[(i + 1)])
    if (updater == FALSE) {
      dots <- data$corrected.transformed.exprs[ww,]
      rgroups <- data$GAPgroups[ww]
    }
    if (updater == TRUE) {
      dots <- data$corrected.VStransformed.exprs[ww,]
      rgroups <- groups[ww]
    }
    rsamples <- data$samples[ww]
    res <-
      apply(matrix(dots, ncol = ncol(dots)),
            1,
            project,
            centers = matrix(rbind(c(center1), c(center2)), ncol = 2))
    res <-
      matrix(cbind(rsamples, dots, rgroups, res), ncol = ncol(mydata))
    res <-
      res[sort.list(as.numeric(res[, ncol(res)]), decreasing = TRUE),]
    mydata <- matrix(rbind(mydata, res), ncol = ncol(mydata))
  }
  mydata <- mydata[-1,]
  
  mydatab <- matrix(0, 1, 5)
  for (i in 1:(length(path) - 1)) {
    center1 <-
      as.numeric(ms[which(ms[, 1] == path[(length(path) - i + 1)]), 2:3])
    center2 <-
      as.numeric(ms[which(ms[, 1] == path[(length(path) - i)]), 2:3])
    ww <- which(groups == path[(length(path) - i)])
    if (updater == FALSE) {
      dots <- data$corrected.transformed.exprs[ww,]
      rgroups <- data$GAPgroups[ww]
    }
    if (updater == TRUE) {
      dots <- data$corrected.VStransformed.exprs[ww,]
      rgroups <- groups[ww]
    }
    rsamples <- data$samples[ww]
    res <-
      -apply(matrix(dots, ncol = ncol(dots)),
             1,
             project,
             centers = matrix(rbind(c(center1), c(center2)), ncol = 2))
    res <-
      matrix(cbind(rsamples, dots, rgroups, res), ncol = ncol(mydatab))
    res <-
      res[sort.list(as.numeric(res[, ncol(res)]), decreasing = TRUE),]
    mydatab <- matrix(rbind(mydatab, res), ncol = ncol(mydatab))
  }
  mydatab <- mydatab[-1,]
  
  mm <- match(mydata[, 1], mydatab[, 1])
  all <- matrix(cbind(mydata, mydatab[mm,]), nrow = nrow(mydata))
  if (path.type[1] == "A2Z") {
    path.update <- path[1:(length(path) - 1)]
    ww <- which(as.numeric(all[, 4]) == path.update[1])
    all[ww, 5] <- all[ww, 10]
    ww <-
      which(as.numeric(all[, 4]) == path.update[length(path.update)])
    all[ww, 10] <- all[ww, 5]
  }
  all <- matrix(as.numeric(all[, c(5, 10)]), ncol = 2)
  means <- apply(matrix(all, ncol = ncol(all)), 1, mean)
  mydata <-
    matrix(cbind(mydata[, 1:4], matrix(means, ncol = 1)), nrow = nrow(mydata))
  
  wh <- which(mydata[, 4] == path[length(path)])
  mydata <-
    matrix(rbind(mydata[wh,], mydata[-wh,]), ncol = ncol(mydata))
  mydata <-
    matrix(cbind(mydata, matrix(1:nrow(mydata), ncol = 1)), nrow = nrow(mydata))
  wh <- match(data$samples, mydata[, 1])
  projEsts <- mydata[wh, (ncol(mydata) - 1):ncol(mydata)]
  if (updater == FALSE) {
    result <- c(data, list(DataSorts = projEsts, DDHFupdate = updater))
  }
  if (updater == TRUE) {
    result <- data
    result$DataSorts <- projEsts
    result$DDHFupdate <- updater
  }
  return(result)
}

#' project
#'
#' A helper to estimate the pseudotimes by orthogonal projection.
#'
#' @param data Data matrix. The fluorescence signals of a particular cluster in both channels.
#' @param centers Numeric vector. A vector of the 2-dimensional cluster centroids.
#'
#' @return The pseudotimes and cell progression
#'
#' @keywords internal
project <- function(data, centers) {
  htenuse <- data - centers[2,]
  centersdiff <- centers[1,] - centers[2,]
  (htenuse[1] * centersdiff[1] + htenuse[2] * centersdiff[2]) / sqrt((centersdiff[1] ^
                                                                        2) + (centersdiff[2] ^ 2))
}

#' transformFluo
#'
#' It performs the variance stabilizing transformation of the flurescence signals.
#'
#' @param data List. A list of adjusted fluorescence signals in both channels. Typically,
#'   the output of orderFluo().
#' @param method Character string. The variance stabilizing method. One of "DDHFmv" or "log".
#'
#' @return A list with the adjusted fluorescence signals and the transformed adjusted fluorescence signals
#'
#' @keywords internal
transformFluo <- function(data, method) {
  d1 <-
    matrix(
      cbind(
        data$index,
        data$samples,
        data$corrected.exprs,
        data$GAPgroups[, 1],
        data$DataSorts
      ),
      nrow = length(data$index)
    )
  nn <- data$UpdatedPath
  dataHFT <- matrix(0, 1, ncol(d1))
  for (i in 1:length(nn)) {
    ww <- which(as.numeric(d1[, 5]) == nn[i])
    dataHFT <-
      matrix(rbind(dataHFT, d1[ww,]), ncol = ncol(dataHFT))
  }
  dataHFT <- dataHFT[-1,]
  ms <- data$centroids
  # run DDHF for each R/G & kmeans cluster group
  if (method == "DDHFmv") {
    datax <- DDHFinput(dataHFT, ms = ms)
    datax1 <-
      c(datax[[1]], datax[[1]][(length(datax[[1]]) - 1):1])[1:2 ^ ceiling(log(length(datax[[1]]), 2))]
    hft <- ddhft.np.2(datax1)[1:length(datax[[1]])]
    back <- revDDHFinput(dataHFT, list(hft, datax[[2]]))
    dataHFT[, 3:4] <- back
    dataHFT <- dataHFT[sort.list(as.numeric(dataHFT[, 1])),]
  }
  if (method != "DDHFmv") {
    dataHFT <- dataHFT[sort.list(as.numeric(dataHFT[, 1])),]
    dataHFT[, 3] <- data$corrected.transformed.exprs[, 1]
    dataHFT[, 4] <- data$corrected.transformed.exprs[, 2]
  }
  return(c(
    data,
    list(
      corrected.VStransformed.exprs = matrix(as.numeric(dataHFT[, 3:4]), ncol =
                                               2),
      VSmethod = method
    )
  ))
}


#' cpoints
#'
#' It performs the change-point analysis of the variance stabilized adjusted fluorescence signals.
#'
#' @param data List. A list of adjusted fluorescence signals in both channels. Typically,
#'   the output of transformFluo().
#' @param thresh Integer. The minimum number of values for a cluster re-estimated by the
#'   change-point analysis.
#' @param cmethod Character string. The change point method to be used. It can be one of "ECP",
#'   (non-parametric) "manualECP" (non-parametric with user-defined numner of change-points) or
#'   "PELT" (Pruned Exact Linear Time; parametric).
#' @param sig.level Float. The significance level below which we do not reject a change point.
#' @param Q Integer. The number of change-points to be kept if CPmethod = "manualECP".
#' @param path.type Character vector. A user-defined vector that characterizes the cell progression dynamics.
#'   The first element can be either "circular" or "A2Z" or "other". If "circular" the path progression is
#'   assummed to exhibit a circle-like behavior. If "A2Z" the path is assumed to have a well-defined start
#'   and a well-defined end point (e.g. a linear progression). If "other" the progression is assumed to be
#'   arbitrary without an obvious directionality.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @importFrom ecp e.divisive
#' @import ggplot2
#'
#' @return A list with the adjusted fluorescence signals and their change-points
#'
#' @keywords internal
cpoints <- function(data,
                    thresh,
                    cmethod,
                    sig.level,
                    Q,
                    path.type,
                    seed) {
  res <-
    cutpointsEstimator(
      data = data,
      thresh = thresh,
      cmethod = cmethod,
      sig.level = sig.level,
      Q = Q,
      seed = seed
    )
  res1 <-
    updateCentroidsPaths(data = data,
                         estimates = res,
                         path.type = path.type)
  if (path.type[1] != "other") {
    res <-
      cutpointsEstimator(
        data = res1[[1]],
        thresh = thresh,
        cmethod = cmethod,
        sig.level = sig.level,
        Q = Q,
        seed = seed
      )
    res1 <-
      updateCentroidsPaths(data = res1[[1]],
                           estimates = res,
                           path.type = path.type)
  }
  result <-
    c(
      res1[[1]],
      list(
        Progression = res1[[2]][[1]][, 2:3],
        Updated.groups = res1[[2]][[1]][, 4],
        CPs = res1[[2]][[2]],
        CPmethod = cmethod,
        CPsig = sig.level,
        CPgroups = Q,
        CPmingroup = thresh
      )
    )
  
  return(result)
}

#' listSorter
#'
#' A helper that sorts the data of a list variable.
#'
#' @param data List. A list variable.
#' @param sorter Numeric vector. An appropriate sorter.
#'
#' @return A list of appropriately sorted data
#'
#' @keywords internal
listSorter <- function(data, sorter) {
  d1 <-
    c("index",
      "samples",
      "Size",
      "correctedAreas",
      "Updated.groups",
      "Residuals")
  d2 <-
    c(
      "RGexprs",
      "batch",
      "corrected.transformed.exprs",
      "corrected.exprs",
      "GAPgroups",
      "DataSorts",
      "corrected.VStransformed.exprs",
      "Progression"
    )
  
  w1 <- match(d1, names(data), nomatch = 0)
  w1 <- w1[-which(w1 == 0)]
  w2 <- match(d2, names(data), nomatch = 0)
  w2 <- w2[-which(w2 == 0)]
  for (i in 1:length(w1)) {
    data[[w1[i]]] <- data[[w1[i]]][sorter]
  }
  for (i in 1:length(w2)) {
    data[[w2[i]]] <- data[[w2[i]]][sorter,]
  }
  return(data)
}

#' cpPELT
#'
#' It performs the change-point analysis of the variance stabilized adjusted fluorescence signals by PELT.
#'
#' @param data List. A list of adjusted fluorescence signals in both channels. Typically,
#'   the output of transformFluo().
#' @param sig.level Float. The significance level below which we do not reject a change point.
#' @param thresh Integer. The minimum number of values for a cluster re-estimated by the
#'   change-point analysis.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @importFrom changepoint cpt.mean cpts
#'
#' @return A list of change-points and the associated change-point clusters
#'
#' @keywords internal
cpPELT <- function(data, sig.level, thresh, seed) {
  if (!is.null(seed)) {
    set.seed(min(seed * 73, 100000000))
  }
  cc <-
    cpt.mean(data,
             penalty = "Asymptotic",
             pen.value = sig.level,
             method = "PELT")
  cc <- cpts(cc)
  cc <- sort(cc[cc != 0])
  ww <- which(cc == 1 | cc == length(data))
  if (length(ww) > 0) {
    cc <- cc[-ww]
  }
  cc1 <- c(1, cc, length(data))
  di <- diff(cc1)
  ww <- which(di < thresh & di == min(di))
  while (length(ww) > 0) {
    di <- di[-ww]
    cc1 <- cc1[-c(ww + 1)]
    di <- diff(cc1)
    ww <- which(di < thresh & di == min(di))
  }
  ww <- which(cc1 > (length(data) - thresh))
  if (length(ww) > 0) {
    cc1 <- cc1[-ww]
  }
  if (cc1[length(cc1)] != length(data)) {
    cc1 <- c(cc1, length(data))
  }
  gg <- rep(0, max(cc1))
  for (i in 2:length(cc1)) {
    gg[cc1[(i - 1)]:cc1[i]] <- (i - 1)
  }
  cc <- cc1[-which(cc1 == 1 | cc1 == length(data))]
  return(list(cluster = gg, estimates = cc))
}


#' GAPanalysis
#'
#' It perfomrs GAP analysis using different methods. It generates the cluster numbers an an indicator of outliers.
#'
#' @param data List. A list of adjusted fluorescence signals. Typically, the output of summarizeAdjFluo().
#' @param fixClusters Integer. A number that defines the number of k-mean clusters to be initially generated.
#'   If 0, the function runs GAP analysis to estimate the optimal number of clusters.
#' @param sigma Integer. A value for the sigma parameter of samSPECTRAL algorithm.
#' @param altFUN Character string. A user-defined method to generate the initial clusters. It can be one of
#'   kmeans, samSpec, fmeans,fmerge or fpeaks.
#' @param k.max Integer. This is the maximum number of clusters that can be generated by k-means (if
#'   fixClusters = 0).
#' @param B.kmeans Integer. The number of bootstrap samples for GAP analysis in altFUN = kmeans.
#' @param savePlot Character string. The directory to store the plots or an option to print them on the screen.
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @importFrom flowClust flowClust
#' @importFrom SamSPECTRAL SamSPECTRAL
#' @importFrom flowMeans flowMeans
#' @importFrom flowPeaks flowPeaks assign.flowPeaks
#' @importFrom flowCore flowFrame
#' @importFrom flowMerge fitPiecewiseLinreg flowObj BIC merge
#' @import cluster parallel stats grDevices foreach
#'
#' @return A list of adjusted fluorescence signals with cluster indices and outlier indicators (the 2s in
#'   the second column of GAPgroups).
#'
#' @keywords internal
GAPanalysis <-
  function(data,
           fixClusters,
           sigma,
           altFUN,
           k.max,
           B.kmeans,
           savePlot,
           seed) {
    if (altFUN != "kmeans" &
        altFUN != "samSpec" &
        altFUN != "fmeans" &
        altFUN != "fmerge" & altFUN != "fpeaks") {
      stop("THe clustering method is not supported. Change the value of the altFUN parameter")
    }
    
    if (altFUN == "kmeans" |
        altFUN == "samSpec" | altFUN == "fmeans") {
      dd <- data$corrected.transformed.exprs
    }
    if (altFUN == "fmerge" | altFUN == "fpeaks") {
      dd <- data$corrected.exprs
      if (altFUN == "fpeaks") {
        dd <- doTransform(dd, "asinh")[[1]]
      }
    }
    
    if (altFUN == "kmeans") {
      if (fixClusters == 0) {
        if (!is.null(seed)) {
          set.seed(min(seed * 7, 100000000))
        }
        gs <-
          clusGap(
            dd,
            FUNcluster = kmeans,
            nstart = 100,
            K.max = k.max,
            B = B.kmeans
          )
        if (!is.null(seed)) {
          set.seed(min(seed * 59, 100000000))
        }
        whmax <- order(gs[[1]][, 3], decreasing = TRUE)
        whmax <- whmax[-which(whmax < whmax[1])]
        whmax <- ceiling(mean(whmax[1:2]))
        g <- kmeans(dd, whmax, nstart = 200)[[1]]
        #      g<-RSKC(dd,whmax,alpha=0.1,nstart=100)$lab
        if (savePlot != "OFF" & savePlot != "screen") {
          pdf(
            paste(
              savePlot,
              "/Fluo_inspection_GAP_",
              data$dateIndex,
              ".pdf",
              sep = ""
            ),
            paper = "a4r"
          )
        }
        
        if (savePlot != "OFF") {
          print(plot_clusgap(
            gs,
            paste(
              altFUN,
              "-based GAP statistics: maxClusters = ",
              k.max,
              sep = ""
            )
          ))
        }
        
        
        if (savePlot != "OFF" & savePlot != "screen") {
          dev.off()
        }
        
      } else {
        if (!is.null(seed)) {
          set.seed(min(seed * 11, 100000000))
        }
        g <- kmeans(dd, fixClusters, nstart = 100)[[1]]
      }
      g <-
        matrix(cbind(matrix(g, ncol = 1), matrix(rep(1, length(
          g
        )), ncol = 1)), ncol = 2)
    }
    
    if (altFUN == "samSpec") {
      g <-
        SamSPECTRAL(
          data.points = dd,
          dimensions = c(1, 2),
          normal.sigma = sigma,
          separation.factor = 0.39,
          maximum.number.of.clusters = k.max
        )
      g1 <- rep(1, length(g))
      w <- as.numeric(is.na(g))
      g1[which(w == 1)] <- 2
      g[which(w == 1)] <- -999
      g <-
        matrix(cbind(matrix(g, ncol = 1), matrix(g1, ncol = 1)), ncol =
                 2)
    }
    
    if (altFUN == "fmeans") {
      colnames(dd) <- c("Ch2", "Ch3")
      res <- flowMeans(dd, colnames(dd))
      g <- res@Label
      g1 <- rep(1, length(g))
      w <- as.numeric(is.na(g))
      g1[which(w == 1)] <- 2
      g1[which(g == 0)] <- 2
      g[which(w == 1)] <- -999
      g[which(g == 0)] <- -999
      g <-
        matrix(cbind(matrix(g, ncol = 1), matrix(g1, ncol = 1)), ncol =
                 2)
    }
    
    if (altFUN == "fmerge") {
      options(warn = -1)
      d1 <- as.matrix(dd)
      colnames(d1) <- c("Ch2", "Ch3")
      res1 <- flowClust(d1, varNames = colnames(d1), K = 1:k.max)
      flowClust.maxBIC <- res1[[which.max(flowMerge::BIC(res1))]]
      flowClust.flowobj <- flowObj(flowClust.maxBIC, flowFrame(d1))
      flowClust.merge <-
        flowMerge::merge(flowClust.flowobj, metric = "entropy")
      i <- fitPiecewiseLinreg(flowClust.merge)
      flowClust.mergeopt <- flowClust.merge[[i]]
      g <- flowClust.mergeopt@label
      g1 <- rep(1, length(g))
      w <- as.numeric(is.na(g))
      g1[which(w == 1)] <- 2
      g[which(w == 1)] <- -999
      g <-
        matrix(cbind(matrix(g, ncol = 1), matrix(g1, ncol = 1)), ncol =
                 2)
    }
    
    if (altFUN == "fpeaks") {
      colnames(dd) <- c("Ch2", "Ch3")
      fp <- flowPeaks(dd)
      g <- fp$peaks.cluster
      g1 <- rep(1, length(g))
      fpc <- assign.flowPeaks(fp, fp$x)
      g1[which(fpc <= 0)] <- 2
      g[which(fpc <= 0)] <- -999
      g <-
        matrix(cbind(matrix(g, ncol = 1), matrix(g1, ncol = 1)), ncol =
                 2)
    }
    if (sum(g[, 2]) > nrow(g)) {
      message(
        "The model identified potential outliers:
        Use FluoSelection_byRun() to remove them from the data or
        relabel them at slot $GAPgroups and use them in the analysis"
      )
    }
    return(c(
      data,
      list(
        GAPgroups = g,
        clusterFUN = altFUN,
        normal.sigma = sigma
      )
    ))
    }


#' trigofun
#'
#' It computes the sin,cos and tan trigonometric functions for a number of points (2-dimensional fluorescence centoids)
#'   relative to the start of the axes.
#'
#' @param data Data matrix. A matrix of fluorescence centroids.
#'
#' @return The centroids and their trigonometric function values
#'
#' @keywords internal
trigofun <- function(data) {
  hypo <- sqrt(sum(data ^ 2))
  sin.est <- data[2] / hypo
  cos.est <- data[1] / hypo
  if (data[1] != 0) {
    tan.est <- data[2] / data[1]
  } else {
    tan.est <- sign(data[2]) * 1e+100
  }
  est <- tan.est
  return(c(sin.est, cos.est, est))
}

#' sortCentroids
#'
#' A helper funcion to sort the centroids .
#'
#' @param data Data matrix. A matrix of centroids with their trigonometric function values.
#' @param type Character string. A user-defined value that characterizes the cell progression dynamics.
#'   It can be either "clockwise" or "anticlockwise" depending on how the path is expected
#'   to proceed.
#'
#' @return A matrix of sorted centroids
#'
#' @keywords internal
sortCentroids <- function(data, type) {
  if (type == "clockwise") {
    data <-
      matrix(data[sort.list(data[, ncol(data)], decreasing = TRUE),], ncol =
               ncol(data))
  } else {
    data <-
      matrix(data[sort.list(data[, ncol(data)], decreasing = FALSE),], ncol =
               ncol(data))
  }
  return(data)
}


#' estimatePath
#'
#' The main function for automatic path estimation .
#'
#' @param data Data matrix. A matrix of centroids with their trigonometric function values.
#' @param type Character string. A user-defined value that characterizes the cell progression dynamics.
#'   It can be either "clockwise" or "anticlockwise" depending on how the path is expected
#'   to proceed.
#' @param start Integer. The cluster number that is assigned as the path starting point
#'
#' @return The sorted cluster indices (path)
#'
#' @keywords internal
estimatePath <- function(data, type, start) {
  res <- as.list(rep(0, 4))
  res[[1]] <-
    matrix(data[which(data[, 2] < 0 &
                        data[, 3] <= 0),], ncol = ncol(data))
  res[[1]] <- sortCentroids(res[[1]], type)
  if (type == "clockwise") {
    res[[2]] <-
      matrix(data[which(data[, 2] <= 0 &
                          data[, 3] > 0),], ncol = ncol(data))
    res[[2]] <- sortCentroids(res[[2]], type)
  } else {
    res[[4]] <-
      matrix(data[which(data[, 2] <= 0 &
                          data[, 3] > 0),], ncol = ncol(data))
    res[[4]] <- sortCentroids(res[[4]], type)
  }
  res[[3]] <-
    matrix(data[which(data[, 2] > 0 &
                        data[, 3] >= 0),], ncol = ncol(data))
  res[[3]] <- sortCentroids(res[[3]], type)
  if (type == "clockwise") {
    res[[4]] <-
      matrix(data[which(data[, 2] >= 0 &
                          data[, 3] < 0),], ncol = ncol(data))
    res[[4]] <- sortCentroids(res[[4]], type)
  } else {
    res[[2]] <-
      matrix(data[which(data[, 2] >= 0 &
                          data[, 3] < 0),], ncol = ncol(data))
    res[[2]] <- sortCentroids(res[[2]], type)
  }
  path <- do.call(rbind, res)[, 1]
  path <- rep(path, 2)
  if (length(which(path == start)) == 0) {
    stop("The defined starting point does not exist. Check the start parameter")
  }
  if (start != path[1]) {
    w <- which(path == start)
    path <- path[w[1]:(w[2] - 1)]
  } else {
    path <- path[1:(length(path) / 2)]
  }
  return(path)
}

#' pathUpdater
#'
#' A helper that updates the path sorted clusters after re-estimation by change-point analysis.
#'
#' @param data Data matrix. A matrix of centroids with their path progression indices.
#' @param path Numeric vector. The path progression indices.
#'
#' @return The sorted cluster indices (path)
#'
#' @keywords internal
pathUpdater <- function(data, path) {
  res <- matrix(0, 1, 4)
  for (i in 1:length(path)) {
    res <-
      matrix(rbind(res, data[which(data[, 4] == path[i]),]), ncol = ncol(res))
  }
  res <- res[-1,]
  #res[,2]<-1:nrow(res)
  return(res)
}


#' cutpointsEstimator
#'
#' A helper that estimates the change-points by various methods.
#'
#' @param data List. A list of adjusted fluorescence signals for both channels.
#' @param thresh Integer. The minimum number of values for a cluster re-estimated by the
#'   change-point analysis.
#' @param cmethod Character string. The change point method to be used. It can be one of "ECP",
#'   (non-parametric) "manualECP" (non-parametric with user-defined numner of change-points) or
#'   "PELT" (Pruned Exact Linear Time; parametric).
#' @param sig.level Float. The significance level below which we do not reject a change point.
#' @param Q Integer. The number of change-points to be kept if CPmethod = "manualECP".
#' @param seed Integer. An optional seed number for the Random Number Generator.
#'
#' @return The sorted transformed signal differences (path) and the associated change-points
#'
#' @keywords internal
cutpointsEstimator <-
  function(data, thresh, cmethod, sig.level, Q, seed) {
    sl <- sort.list(as.numeric(data$DataSorts[, 2]))
    newdata <- listSorter(data = data, sorter = sl)
    dd <-
      -apply(matrix(newdata$corrected.VStransformed.exprs, ncol = 2),
             1,
             diff)
    if (cmethod == "PELT") {
      cps <- cpPELT(
        data = dd,
        thresh = thresh,
        sig.level = sig.level,
        seed = seed
      )
      if (length(unique(cps$cluster)) == 1) {
        stop(
          "The progression estimation failed because PELT did not detect change points. Try ECP method."
        )
      }
      res <-
        matrix(cbind(
          matrix(newdata$index, ncol = 1),
          matrix(1:length(dd), ncol = 1),
          matrix(dd, ncol = 1),
          matrix(cps$cluster, ncol = 1)
        ), nrow = length(dd))
      cps <- cps$estimates
    }
    if (cmethod == "ECP") {
      if (!is.null(seed)) {
        set.seed(min(seed * 5 + 121, 100000000))
      }
      v <-
        e.divisive(
          X = matrix(dd, ncol = 1),
          sig.lvl = sig.level,
          R = 199,
          k = NULL,
          min.size = thresh,
          alpha = 1
        )
      if (length(unique(v$cluster)) == 1) {
        stop(
          "The progression estimation failed because ECP did not detect change points. Try PELT method."
        )
      }
      res <-
        matrix(cbind(newdata$index, 1:length(dd), dd, v$cluster), ncol = 4)
      cps <- v$estimates[2:(length(v$estimates) - 1)]
    }
    if (cmethod == "manualECP") {
      if (!is.null(seed)) {
        set.seed(min(seed * 5 + 121, 100000000))
      }
      v <-
        e.divisive(
          X = matrix(dd, ncol = 1),
          sig.lvl = sig.level,
          R = 199,
          k = NULL,
          min.size = thresh,
          alpha = 1
        )
      if (length(unique(v$cluster)) == 1) {
        stop(
          "The progression estimation failed because ECP did not detect change points. Try PELT method."
        )
      }
      cps <- v$order.found[3:length(v$order.found)]
      cps <- cps[1:min(Q, length(cps))]
      cps <- sort(cps)
      cc1 <- c(1, cps, length(dd))
      clu <- rep(1, length(dd))
      for (i in 1:(length(cc1) - 1)) {
        clu[cc1[i]:cc1[(i + 1)]] <- i
      }
      res <-
        matrix(cbind(newdata$index, 1:length(dd), dd, clu), ncol = 4)
    }
    res <- res[sort.list(res[, 1]),]
    return(list(res, cps))
  }

#' updateCentroidsPaths
#'
#' It updates the path sorted clusters after re-estimation by change-point analysis.
#'
#' @param data List. A list of adjusted fluorescence signals.
#' @param estimates List. A list of sorted Ch2-Ch3 transformed fluorescence signals with
#'   their assocuated change-points.
#' @param path.type Character vector. A user-defined vector that characterizes the cell progression dynamics.
#'   The first element can be either "circular" or "A2Z" or "other". If "circular" the path progression is
#'   assummed to exhibit a circle-like behavior. If "A2Z" the path is assumed to have a well-defined start
#'   and a well-defined end point (e.g. a linear progression). If "other" the progression is assumed to be
#'   arbitrary without an obvious directionality.
#'
#' @return A list of adjusted fluorescence signals and the updated path after the change-point analysis
#'
#' @import stats
#' @keywords internal
updateCentroidsPaths <- function(data, estimates, path.type) {
  cc <- colnames(data$centroids)
  mat <-
    matrix(cbind(data$corrected.VStransformed.exprs, estimates[[1]][sort.list(estimates[[1]][, 1]), 4]),
           ncol = 3)
  u <- unique(mat[, 3])
  data$centroids <- matrix(0, length(u), 3)
  data$centroids[, 1] <- u
  for (i in 1:length(u)) {
    data$centroids[i, 2:3] <-
      apply(matrix(mat[which(mat[, 3] == u[i]), 1:2], ncol = 2), 2, median)
  }
  if (path.type[1] != "other" & length(estimates[[2]]) > 0) {
    cdata <- data$centroids
    oc <- apply(matrix(cdata[, 2:3], ncol = 2), 2, mean)
    cdata[, 2] <- cdata[, 2] - oc[1]
    cdata[, 3] <- cdata[, 3] - oc[2]
    trigs <-
      matrix(cbind(cdata, t(apply(
        matrix(cdata[, 2:3], ncol = 2), 1, trigofun
      ))), nrow = nrow(cdata))
    data$UpdatedPath <-
      estimatePath(trigs, type = path.type[2], start =
                     1)
    estimates[[1]] <-
      pathUpdater(data = estimates[[1]], path = data$UpdatedPath)
    sl <- sort.list(estimates[[1]][, 1])
    estimates[[1]] <- estimates[[1]][sl,]
    
    #check!
    data$DataSorts[, 1] <- data$DataSorts[sl, 1]
    
    data$DataSorts[, 2] <- estimates[[1]][, 2]
  }
  colnames(data$centroids) <- cc
  return(list(data, estimates))
}

#' plot_clusgap
#'
#' adapted from plot_clusgap from phyloseq.
#'
#' @return plot gap stats
#'
#' @import ggplot2
#' @keywords internal
plot_clusgap <-
  function (clusgap, title = "Gap Statistic results") {
    SE.sim <- gap <- k <- NULL
    gstab = data.frame(clusgap$Tab, k = 1:nrow(clusgap$Tab))
    p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point()
    p = p + geom_errorbar(
      width = .8,
      size = 1,
      aes(ymax = gap + SE.sim, ymin = gap - SE.sim),
      colour = 'red'
    )
    p = p + ggtitle(title) + theme_bw() + theme(legend.position = "none")
    return(p)
  }



#' estimate.new.pseudotimes
#'
#' It estimates a unique pseudotime vector to be used for analysis with Fluo_ordering(). It is based on the cross-validation
#'   estimates for each sample and a particular method defined by pseudo.est.method parameter at Fluo_CV_modeling().
#'
#' @param data Data matrix. The estimated pseudotimes for all samples with the original data and the CV data undet two methods,
#'   i.e. "median/original" and median/null". For details see parameter pseudo.est.method at Fluo_CV_modeling().
#'
#' @return It summarizes the CV-estimated pseudotimes into a single value. There are three possible methods that may produce
#'   different results. For details see parameter pseudo.est.method at Fluo_CV_modeling().
#'
#' @keywords internal
estimate.new.pseudotimes <- function(data) {
  new.pseudos <- matrix(0, 2, ncol(data))
  for (k in 2:nrow(data)) {
    dat <- matrix(cbind(1:ncol(data), data[1,], data[k,]), ncol = 3)
    rank.pseudo <- rank(as.numeric(dat[, 3]))
    sl <- sort.list(rank.pseudo)
    rank.data <- dat[sl,]
    u <- unique(as.numeric(rank.data[, 3]))
    
    result <- matrix(0, 1, ncol(rank.data))
    for (i in 1:length(u)) {
      d1 <-
        matrix(rank.data[which(as.numeric(rank.data[, 3]) == u[i]),], ncol = ncol(rank.data))
      if (nrow(d1) == 1 | u[i] == 0) {
        result <- matrix(rbind(result, d1), ncol = ncol(result))
      } else {
        d1 <- d1[sort.list(as.numeric(d1[, 2])),]
        d1[, 3] <-
          as.numeric(d1[, 3]) + as.numeric(d1[, 3]) * (1:nrow(d1)) / 1000000
        result <- matrix(rbind(result, d1), ncol = ncol(result))
      }
    }
    
    result <- result[-1,]
    result[which(as.numeric(result[, 3]) > 0), 3] <-
      rank(as.numeric(result[which(as.numeric(result[, 3]) > 0), 3]))
    new.pseudos[(k - 1),] <- result[sort.list(result[, 1]), 3]
  }
  
  return(new.pseudos)
}


#' reestimate.pseudos.byCV
#'
#' It estimates a new pseudotime for each sample based on its cross-validation estimates.
#'
#' @param data Numeric. A vector of estimated pseudotimes for one sample. The first value corresponds to the estimate
#'   of the whole data. The second value is the difference between the estimate of the original data and the median
#'   of the CV-estimated pseudotimes and it is used as a dissimilarity measure. The rest of the values are the CV estimated
#'   pseudotimes themselves.
#' @param diff.quantile Float. The qth quantile of the distribution of the difference between the original and the CV-estimated
#'   pseudotimes. The q parameter is defined in Fluo_CV_modeling().
#' @param perc.cutoff Float. The percentage of CV-estimated pseudotimes that are similar (clustered together by k-means)
#' @param pseudotime.cutoff Integer. A user-defined value to define outlier samples, i.e. samples with Pseudotime(original) -
#'   median{Pseudotime(CV)} > pseudotime.cutoff.
#'
#' @return It summarizes the CV-estimated pseudotimes into a single value. There are three possible methods that may produce
#'   different results. For details see parameter pseudo.est.method at Fluo_CV().
#'
#' @import stats
#'
#' @keywords internal
reestimate.pseudos.byCV <-
  function(data,
           diff.quantile,
           perc.cutoff,
           pseudotime.cutoff) {
    perc.cutoff <- max(perc.cutoff, 1 - perc.cutoff)
    
    f <- data[1]
    qu <- data[2]
    data1 <- data[3:length(data)]
    
    if (length(unique(data1)) > 1) {
      k <- kmeans(data1, 2)[[1]]
      tab <- table(k) / length(data1)
      tab <- tab[sort.list(tab)]
    } else {
      k <- rep(1, length(data1))
      tab <- c("2" = 0, "1" = 1)
    }
    
    pseu.diff <- rep(0, 2)
    for (i in 1:2) {
      pseu.diff[i] <-
        abs(median(data1[which(k == as.numeric(names(tab[i])))]) - f)
      pseu.diff[is.na(pseu.diff)] <- 1e+6
    }
    
    # which has the highest frequency
    w1 <- which(k == as.numeric(names(tab))[2])
    
    #which is closest to the original
    w2 <- which(k == names(tab)[which(pseu.diff == min(pseu.diff))])
    est <- c(median(data1[w1]), median(data1[w2]))
    
    est1 <- est2 <- c()
    if (max(tab) <= perc.cutoff &
        qu > diff.quantile & qu > pseudotime.cutoff) {
      est1 <- est[2]
      est2 <- 0
    } else {
      est1 <- est2 <- est[1]
    }
    
    c(f, est1, est2)
  }


#' path.initiator
#'
#' It finds the cluster that initiates the progression path.
#'
#' @param data List. The output of Fluo_inspection().
#' @param where Character. One of "bottom/left", "bottom/right", "top/left", "top/right" that specify
#'   the starting point of the progression path.
#'
#' @return A starting point for the progression path
#'
#' @keywords internal
path.initiator <- function(data, where) {
  w <- c()
  if (where != "random") {
    w <- where
  } else {
    w <-
      sample(c("bottom/left", "bottom/right", "top/left", "top/right"),
             1)
  }
  
  if (w == "bottom/left") {
    sl1 <- sort.list(data$centroids[, 2])
    sl2 <- sl1[sort.list(data$centroids[sl1[1:2], 3])]
  } else if (w == "bottom/right") {
    sl1 <- sort.list(data$centroids[, 2], decreasing = TRUE)
    sl2 <- sl1[sort.list(data$centroids[sl1[1:2], 3])]
  } else if (w == "top/left") {
    sl1 <- sort.list(data$centroids[, 2])
    sl2 <-
      sl1[sort.list(data$centroids[sl1[1:2], 3], decreasing = TRUE)]
  } else if (w == "top/right") {
    sl1 <- sort.list(data$centroids[, 2], decreasing = TRUE)
    sl2 <-
      sl1[sort.list(data$centroids[sl1[1:2], 3], decreasing = TRUE)]
  } else {
    stop("Invalid starting point. Revise init.path parameter!")
  }
  
  return(c(sl2[1], w))
}


#' CVsampler
#'
#' It samples a data subset for the cross-validation analysis.
#'
#' @param data List. The output of Fluo_inspection() or Fluo_modeling(). It requires existence of the @GAPgroups slot.
#' @param f Float. The percentage of samples that are used in the CV analysis (the rest is re-estimated).
#'
#' @return An index with the data that will remain in the analysis.
#'
#' @keywords internal
CVsampler <- function(data, f) {
  if (f >= 1 | f <= 0) {
    stop("Parameter fold is a percentage and it should be in (0,1)")
  }
  
  u <- unique(data$GAPgroups[, 1])
  res <- c()
  
  for (i in 1:length(u)) {
    w <- which(data$GAPgroups[, 1] == u[i])
    res <- c(res, sample(w, floor(f * length(w))))
  }
  
  return(res)
}


#' predict.kmeans
#'
#' Takes a training sample and predicts the k-mean clusters of a new dataset (minimizing the Eucledian distance).
#'
#' @param data Data matrix. A 2-dimensional matrix of corrected fluorescence signals that are clustered by k-means.
#' @param centroid Data matrix. A 2-dimensional matrix of the k-means centroids.
#'
#' @return The predicted k-mean clusters.
#'
#' @import stats
#'
#' @keywords internal
predict.kmeans <- function(data, centroid) {
  di <- rep(0, nrow(centroid))
  for (i in 1:nrow(centroid)) {
    di[i] <- dist(rbind(data, centroid[i,]))
  }
  which.min(di)[1]
}


#' addKpredictions
#'
#' Adds the predicted k-mean clusters to the existing set (CV-estimated).
#'
#' @param whole List. The output of Fluo_inspection() on the original data.
#' @param cv List. The output of Fluo_inspection() on the cross-validated data.
#' @param new.clusters Numeric. The predicted k-mean clusters.
#' @param indices List. The index numbers of the old and the new data.
#'
#' @return The new output of Fluo_inspection() with the original data where @GAPgroups has been replaced with
#'   the CV estimates and the new predictions and the @centroids with the CV-estimated centroids.
#'
#' @keywords internal
addKpredictions <- function(whole, cv, new.clusters, indices) {
  for (i in 1:length(indices[[1]])) {
    whole$GAPgroups[indices[[1]][i], 1] <- cv$GAPgroups[i, 1]
  }
  for (i in 1:length(indices[[2]])) {
    whole$GAPgroups[indices[[2]][i], 1] <- new.clusters[i]
  }
  whole$centroids <- cv$centroids
  
  return(whole)
}



#' aveDiff
#'
#' It estimates the average difference between the original and the CV estimated pseudotimes. For circular path
#'   types, the difference is defined as min(diff,max(pseudotime)-diff). For example assuming 300 cells (thus the
#'   maximum pseudotime is 300) in a circular path, two pseudotimes 1 and 300 differ only by 1 and not by 299.
#'
#' @param data Numeric. A vector of pseudotimes whose first element is the originally estimated one (fromt he full data).
#' @param path.type Character. The input of path.type parameter in pathEstimator().
#' @param maxPseudo Numeric. The maximum possible pseudotime.
#'
#' @return The average difference between the original and the CV estimated pseudotimes for a sample.
#'
#' @import stats
#'
#' @keywords internal
aveDiff <- function(data, path.type, maxPseudo) {
  if (path.type != "circular") {
    median(data[1] - data[2:length(data)])
  }
  
  if (path.type == "circular") {
    res <- c()
    for (i in 2:length(data)) {
      res <-
        c(res, min(abs(data[1] - data[i]), maxPseudo - abs(data[1] - data[i])))
    }
    median(res)
  }
}
