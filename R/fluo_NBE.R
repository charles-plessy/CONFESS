
#' createFluo
#'
#' The data format creator function for the signal normalization step.
#'
#' @param data Data matrix. The output data matrix of LocationMatrix().
#' @param dateIndex a date index to be used for storing the output files. It is either transfered from LocationMatrix() or
#'     it is generated here for the first time (e.g. if image analysis was not run by CONFESS or if the analysis has
#'     been repeated many times).
#' @param from.file Logical. If TRUE the data is read from a file whose format should be the same to
#'   the output of LocationMatrix(). Default is FALSE.
#' @param separator Character string. It separates the run ID from the Well ID in the image filenames
#'   (the <<separator1>> of readFiles()). It is also used here to enable the user perform the analysis
#'   independently of the previous step (cell recognition via imaging). Default is "_".
#'
#' @return A list of reformed data to be used in subsequent analysis:
#'   index: The sample indices.
#'   RGexprs: the foreground (columns 1 and 3) and background (columns 2 and 4) signals of each channel that
#'     have been estimated by spotEstimator() and filtered in LocationMatrix().
#'   samples: the sample IDs.
#'   batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
#'     original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
#'   size: the estimated cell size.
#'   image.type: the image type IDs as defined in readFiles(). The parameter is kept in ordeer to enable the user to
#'     use this function independently of the image analysis step.
#'   dateIndex: a date index to be used for storing the output files. It is either transfered from LocationMatrix() or
#'     it is generated here for the first time (e.g. if image analysis was not run by CONFESS or if the analysis has
#'     been repeated many times).
#'
#' @import utils
#'
#' @export
#'
#' @examples
#' step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#' package = "CONFESS"),separator="_")
createFluo <-
  function(data,
           dateIndex = c(),
           from.file = FALSE,
           separator = "_") {
    if (from.file != FALSE) {
      data <-
        read.table(
          from.file,
          sep = "\t",
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      
      if (length(dateIndex) == 0) {
        datin <- paste(unlist(strsplit(date(), " ")), collapse = "")
      } else {
        datin <- dateIndex
      }
      
    }
    
    if (from.file == FALSE) {
      data <- data$Output
      datin <- data$dateIndex
    }
    
    data <- data[grep(1, data$Cells), c(1, 6:9, 4)]
    data <-
      cbind(data, Runs = unlist(strsplit(data[, 1], separator))[c(TRUE, FALSE)])
    
    options(warn = -1)
    vec <- as.numeric(is.na(as.numeric(names(data))))
    if (sum(vec) != length(vec)) {
      message("Warning: The data column names may be wrong! Check the values of the $RGexprs slot")
    }
    
    dd <- data[, 2:(ncol(data) - 2)]
    image.type = c("BF", unique(t(matrix(
      unlist(strsplit(names(dd), "_")), nrow = 2
    ))[, 2]))
    
    #samples<-data[,1]
    options(warn = -1)
    if (is.na(as.numeric(data[1, ncol(data)]))) {
      batches <- as.numeric(factor(data[, ncol(data)]))
    } else {
      batches <- as.numeric(data[, ncol(data)])
    }
    ub <- unique(batches)
    for (i in min(ub):max(ub)) {
      print(paste("Run ", unique(data[which(batches == i), ncol(data)]), " was converted into ", i, sep =
                    ""))
    }
    area <- as.numeric(data[, (ncol(data) - 1)])
    return(
      list(
        index = 1:nrow(data),
        RGexprs = dd,
        samples = data$SampleID,
        batch = matrix(cbind(as.character(data$Runs), batches), ncol =
                         2),
        Size = area,
        image.type = image.type,
        dateIndex = datin
      )
    )
  }



#' Fluo_adjustment
#'
#' A summary of the signal adjustment algorithms into a single function. It corrects the
#'   run effect (if any) and performs background adjustment for appropriately transformed data.
#'
#' @param data List. The output of createFluo().
#' @param BGmethod Character string. The type of image background correction to be performed.
#'   One of "normexp" or "subtract". Default is "normexp".
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model. If maxMix=1 or if the the optimal number of the estimated components
#'   is 1, the model reduces to the classical 2-way ANOVA. Default is 3.
#' @param single.batch.analysis Integer. The baseline run against with the run effect
#'   correction is perfomred. Default is 1. If 0, each run is used as baseline iteratively
#'   and the final corrected data are obtained as the average of all corrections.
#' @param transformation Character string. One of bc (Box-Cox), log, log10, asinh transforms
#'   applied to the data. Default is "log".
#' @param prior.pi Float. The prior probability to accept a component. Default is 0.1.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate the flexmix
#'   model. Default is 50.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL". Default is "BIC".
#' @param savePlot Character string. Directory to store the plots. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen or "OFF" that does not generate a plot (suggested
#'   only during cross-validations). Default is the current working directory, getwd().
#' @param seed Integer. An optional seed number for the Random Number Generator. Note that this seed is a 'reference'
#'   value of the actual seed used in sampling. CONFESS is using various random sampling methods. Each method's
#'   actual seed is factor*seed. The factors vary across methods. Default is NULL.
#'
#' @return A list with the data description, the normalized and corrected estimates over all runs by averaging
#'   (Summarized_estimates) AND for a particular "reference" run (Batch_estimates). Analytically, the components
#'   are:
#'   General
#'     index: The sample indices.
#'     samples: the sample IDs.
#'     batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
#'         original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
#'     Size: the estimated cell size.
#'     RGexprs: the foreground (columns 1 and 3) and background (columns 2 and 4) signals of each channel that
#'       have been estimated by spotEstimator() and filtered in LocationMatrix().
#'     exprs: the background corrected (only) signals of each channel. These data are fed into the flexmix model.
#'     image.type: the image type IDs as defined in readFiles().
#'     dateIndex: the date index used.
#'     single.batch.analysis: the reference run used for run effect correction with flexmix.
#'     BGmethod: the background correction method used.
#'     maxMix: the maxMix parameter used.
#'     prior.pi: the prior.pi parameter used.
#'     flex.reps: the flex.reps parameter used.
#'     flexmethod: the flexmethod parameter used.
#'     RNG: the seed that is used to generate the results.
#'
#'   Summarized_estimates:
#'     corrected.exprs: the background and run effect corrected channel signals (by averaging the estimates of all runs).
#'     corrected.transformed.exprs: the background and run effect transformed corrected channel signals (by averaging the
#'       estimates of all runs). The transformation is defined in the transformation parameter (see above).
#'     allResults: the background and run effect corrected and transformed corrected channel signals (two different slots) for
#'       all runs.
#'
#'   Batch_estimates: it contains the analytical results for each batch in different slots. Each slot includes:
#'     corrected.exprs: the background and run effect corrected channel signals (for a run).
#'     corrected.transformed.exprs: the background and run effect transformed corrected channel signals (for a run).
#'       The transformation is defined in the transformation parameter (see above).
#'     mixes.(image.type 1): the estimated components of the flexmix model for one channel.
#'     mixes.(image.type 2): the estimated components of the flexmix model for the other channel.
#'     Batch.(image.type 1).est: the run effects of one channel. It contains the model estimates and significance P-values/FDRs.
#'       "Comp" corresponds to the factor of flexmix components (mixes) and "Batch" to the factor of runs.
#'     Batch.(image.type 2).est: the run effects of the other channel. It contains the model estimates and significance P-values/FDRs.
#'       "Comp" corresponds to the factor of flexmix components (mixes) and "Batch" to the factor of runs.
#'     fitted.values: the fitted values of the flexmix model for each channel.
#'     transformation: the transformation applied on the fluorescence signals (it stores the value of transformation parameter).
#'     model.residuals: the flexmix residuals for each channel.
#'     model.standardized.residuals: the flexmix standardized residuals for each channel.
#'     residual.statistics: the result of various normality tests for the residuals.
#'     lpar: the lambda parameter of the Box-Cox transformation (if used).
#'     design.(image.type 1): the design matrix of one channel.
#'     design.(image.type 2): the design matrix of the other channel.
#'     reference: the run that has been used as reference.
#'     (image.type 1).contrasts: the contrasts matrix for the differences across flexmix components and runs for one channel (only for
#'       the reference batch if any).
#'     (image.type 2).contrasts: the contrasts matrix for the differences across flexmix components and runs for the other channel (only for
#'       the reference batch if any).
#'
#' @export
#'
#' @examples
#' step2 <- Fluo_adjustment(data=step1,flex.reps = 5,single.batch.analysis=5,savePlot="OFF")
Fluo_adjustment <-
  function(data,
           BGmethod = "normexp",
           maxMix = 3,
           single.batch.analysis = 1,
           transformation = "log",
           prior.pi = 0.1,
           flex.reps = 50,
           flexmethod = "BIC",
           savePlot = getwd(),
           seed = NULL) {
    if (length(unique(as.numeric(data$batch[, 2]))) == 1) {
      stop("The data come from a single run. Use directly getFluo_byRun()")
    }
    
    reference = min(as.numeric(data$batch[, 2])):max(as.numeric(data$batch[, 2]))
    if (!transformation %in% c("log", "log10", "asinh", "bc")) {
      stop("This transformation is not available")
    }
    #    reference <- reference[which(reference >= min(as.numeric(data$batch[,2])) &
    #                                  reference <= max(as.numeric(data$batch[,2])))]
    
    if (savePlot != "OFF" & savePlot != "screen") {
      if (as.numeric(file.access(savePlot)) < 0) {
        stop("The savePlot directory cannot be found in the system or invalid value for savePlot")
      }
    }
    
    step <-
      summarizeAdjFluo(
        data = data,
        transformation = transformation,
        BGmethod = BGmethod,
        maxMix = maxMix,
        reference = reference,
        prior.pi = prior.pi,
        flex.reps = flex.reps,
        flexmethod = flexmethod,
        image.type = data$image.type,
        savePlot = savePlot,
        seed = seed
      )
    
    colnames(step$summ.corrected) <-
      colnames(step$summ.corrected.transformed) <-
      data$image.type[2:3]
    stepB <- c()
    
    if (single.batch.analysis > 0) {
      if (length(step$One.batch) < single.batch.analysis) {
        single.batch.analysis <- length(step$One.batch)
        message(paste(
          "the baseline batch has changed to ",
          length(step$One.batch),
          sep = ""
        ))
      }
      stepA <- step$One.batch[[single.batch.analysis]]
      stepB <-
        contrastFluo(
          data = stepA,
          channel = "CH1",
          legends = data$image.type[2]
        )
      stepB <-
        contrastFluo(
          data = stepB,
          channel = "CH2",
          legends = data$image.type[3]
        )
    } else {
      stepA <- step$One.batch[[1]]
      stepB <-
        contrastFluo(
          data = stepA,
          channel = "CH1",
          legends = data$image.type[2]
        )
      stepB <-
        contrastFluo(
          data = stepB,
          channel = "CH2",
          legends = data$image.type[3]
        )
    }
    
    list1 <-
      list(
        index = step$index,
        samples = step$samples,
        batch = step$batch,
        Size = step$Size,
        RGexprs = stepB$RGexprs,
        exprs = stepB$exprs,
        image.type = stepB$image.type,
        dateIndex = stepB$dateIndex,
        single.batch.analysis = single.batch.analysis,
        BGmethod = BGmethod,
        maxMix = maxMix,
        prior.pi = prior.pi,
        flex.reps = flex.reps,
        flexmethod = flexmethod,
        RNG = seed
      )
    list2 <-
      list(
        corrected.exprs = step$summ.corrected,
        corrected.transformed.exprs = step$summ.corrected.transformed,
        allResults = step$allResults
      )
    
    #### change February 2016! single.batch.analysis =1 failed to give the stepB estimates
    list3.b <- list(0)
    for (i in 1:length(step$One.batch)) {
      if (i != single.batch.analysis) {
        list3.b1 <- list(step$One.batch[[i]][9:23])
        names(list3.b1[[1]]) <- names(stepB)[9:23]
      } else {
        list3.b1 <- list(c(step$One.batch[[i]][9:23], stepB[24:25]))
        names(list3.b1[[1]]) <- names(stepB)[9:25]
      }
      list3.b <- c(list3.b, list3.b1)
    }
    list3.b <- list3.b[-1]
    #### change February 2016! single.batch.analysis =1 failed to give the stepB estimates
    
    
    names(list3.b) <- paste("Batch", reference, sep = "")
    
    return(c(
      General = list(list1),
      Summarized_estimates = list(list2),
      Batch_estimates = list(list3.b)
    ))
  }


#' getFluo
#'
#' It retrieves the run effect and background corrected signals.
#'
#' @param data List. The output of the Fluo_adjustment().
#' @param areacut Integer. The "artificial" area size (BFarea^2) of the cells estimated
#'   by BF image modelling. Default is 0, implying that the area sizes to be corrected will
#'   by estimated automatically from the data (not recommended if prior knowledge exists).
#'
#' @return A list of estimates to be used in subsequent analysis (the slots are the same to those of getFluo_byRun()):
#'   index: The sample indices.
#'   samples: the sample IDs.
#'   batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
#'       original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
#'   Size: the estimated cell size.
#'   corrected.exprs: the background corrected channel signals (case of a single run).
#'   corrected.transformed.exprs: the background transformed corrected channel signals (case of a single run). The
#'     transformation is defined in the transformation parameter.
#'   correctedAreas: the log-transformed areas after correction and imputation.
#'   areacut: the above areacut if different from 0 or the automatically calulated one otherwise.
#'   transformation: the transformation applied on the fluorescence signals.
#'   image.type: the image type IDs as defined in readFiles(). The parameter is kept in order to enable the user to
#'     use this function independently of the image analysis step.
#'   dateIndex: the date index used.
#'   single.batch.analysis: the reference run of the run effect correction by flexmix.
#'   BGmethod: the background correction methods used.
#'   maxMix: the maxMix parameter used.
#'   prior.pi: the prior.pi parameter used.
#'   flex.reps: the flex.reps parameter used.
#'   flexmethod: the flexmethod parameter used.
#'   RNG: the seed that is used to generate the results.
#'
#' @import stats
#' @export
#'
#' @examples
#' step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#' package = "CONFESS"),separator="_")
#' step2.1 <- getFluo(data=step2)
getFluo <- function(data, areacut = 0) {
  testarea <- data$General$Size
  
  if (length(which(testarea == 0)) > 0) {
    message("Warning: Some sizes were 0 and have been truncated to 1")
  }
  testarea[which(testarea == 0)] <- 1
  
  
  if (areacut == 0) {
    tt <- as.numeric(table(testarea))
    w <- which(tt == max(tt))
    if (length(w) > 1 | length(tt) == 1) {
      message("Warning: Parameter areacut can not be estimated. The areas are not corrected!")
    } else {
      areacut <- as.numeric(names(which.max(table(testarea))))
    }
  }
  
  
  # correct the areas
  if (data$General$single.batch.analysis == 0) {
    signals <- data$Summarized_estimates$corrected.transformed.exprs
  } else {
    signals <-
      data$Batch_estimates[[data$General$single.batch.analysis]]$corrected.transformed.exprs
  }
  existingAreas <- which(testarea != areacut)
  nonexistingAreas <- which(testarea == areacut)
  u.areas <- log(testarea)
  bb <- as.numeric(data$General$batch[, 2])
  mod <-
    model.matrix( ~ signals[existingAreas, 1] + signals[existingAreas, 2] +
                    factor(bb[existingAreas]))
  res <- lm(u.areas[existingAreas] ~ signals[existingAreas, 1] +
              signals[existingAreas, 2] + factor(bb[existingAreas]))
  
  for (i in 1:length(existingAreas)) {
    for (j in 4:ncol(mod)) {
      u.areas[existingAreas[i]] <-
        u.areas[existingAreas[i]] - res[[1]][[j]] * mod[i, j]
    }
  }
  
  c.area <- u.areas
  
  if (length(nonexistingAreas) > 0) {
    for (i in 1:length(nonexistingAreas)) {
      c.area[nonexistingAreas[i]] <-
        res[[1]][[1]] + res[[1]][[2]] * signals[nonexistingAreas[i], 1] + res[[1]][[3]] *
        signals[nonexistingAreas[i], 2]
    }
  }
  
  # collect the Signals Of Interest
  SOI <- c(data$General[1:3], list(Size = testarea))
  if (data$General$single.batch.analysis == 0) {
    SOI <- c(SOI, data$Summarized_estimates[1:2])
  } else {
    SOI <-
      c(SOI, data$Batch_estimates[[data$General$single.batch.analysis]][1:2])
  }
  SOI <-
    c(
      SOI,
      list(correctedAreas = c.area, areacut = areacut),
      data$Batch_estimates[[1]][8],
      data$General[7:15]
    )
  
  
  return(SOI)
}


#' getFluo_byRun
#'
#' It produces the background corrected data when run correction is not needed. It can
#'   be used for data coming from a single run instead of Fluo_adjustment().
#'   Alternatively, this function can be used to visualize the fluorescence
#'   densities of a single batch before deciding the form of the normalization
#'   model.
#'
#' @param data List. The output of createFluo().
#' @param BGmethod Character string. The type of image background correction to be performed.
#'   One of "normexp" or "subtract". Default is "normexp".
#' @param transformation Character string. One of bc (Box-Cox), log, log10, asinh transforms
#'   applied to the data. Default is "log".
#' @param areacut Integer. The "artificial" area size (BFarea^2) of the cells estimated
#'   by BF image modelling. Default is 0, implying that the area sizes to be corrected will
#'   by estimated automatically from the data (not recommended if prior knowledge exists).
#' @param savePlot Character string. Directory to store the plots. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen or "OFF" that does not generate a plot (suggested
#'   only during cross-validations). Default is the current working directory, getwd().
#'
#' @return A list of corrected signal estimates. The slots are the same to those of getFluo():
#'   index: The sample indices.
#'   samples: the sample IDs.
#'   batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
#'       original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
#'   Size: the estimated cell size.
#'   corrected.exprs: the background corrected channel signals (case of a single run).
#'   corrected.transformed.exprs: the background transformed corrected channel signals (case of a single run). The
#'     transformation is defined in the transformation parameter.
#'   correctedAreas: the log-transformed areas after correction and imputation.
#'   areacut: the above areacut if different from 0 or the automatically calulated one otherwise.
#'   transformation: the transformation applied on the fluorescence signals.
#'   image.type: the image type IDs as defined in readFiles(). The parameter is kept in order to enable the user to
#'     use this function independently of the image analysis step.
#'   dateIndex: the date index used.
#'   single.batch.analysis: it returns 0 because there is no run effect correction done.
#'   BGmethod: the background correction methods used.
#'   maxMix: it returns NULL because there is no flexmix run effect correction done.
#'   prior.pi: it returns NULL because there is no flexmix run effect correction done.
#'   flex.reps: it returns NULL because there is no flexmix run effect correction done.
#'   flexmethod: it returns NULL because there is no flexmix run effect correction done.
#'   RNG: the seed that is used to generate the results.
#'
#' @import ggplot2 reshape2 stats grDevices
#' @export
#'
#' @examples
#' step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#' package = "CONFESS"),separator="_")
#'
#' ### select the samples of a single run and correct them
#' step2a <- FluoSelection_byRun(data = step1, batch = 5)
#' step2.1 <- getFluo_byRun(data=step2a,savePlot="OFF")
getFluo_byRun <-
  function(data,
           BGmethod = "normexp",
           areacut = 0,
           transformation = "log",
           savePlot = getwd()) {
    testarea <- data$Size
    
    if (length(which(testarea == 0)) > 0) {
      message("Warning: Some sizes were 0 and have been truncated to 1")
    }
    testarea[which(testarea == 0)] <- 1
    
    if (areacut == 0) {
      tt <- as.numeric(table(testarea))
      w <- which(tt == max(tt))
      if (length(w) > 1 | length(tt) == 1) {
        message("Warning: Parameter areacut can not be estimated. The areas are not corrected!")
      } else {
        areacut <- as.numeric(names(which.max(table(testarea))))
      }
    }
    
    if (transformation != "log" &
        transformation != "log10" &
        transformation != "asinh" & transformation != "bc") {
      stop("This transformation is not available")
    }
    
    if (savePlot != "OFF" & savePlot != "screen") {
      if (as.numeric(file.access(savePlot)) < 0) {
        stop("The savePlot directory cannot be found in the system or invalid value for savePlot")
      }
    }
    
    if (length(unique(as.numeric(data$batch[, 2]))) > 1) {
      stop("The data come from different runs. Use Fluo_adjustment()")
    }
    RG <-
      new(
        "RGList",
        list(
          R = data$RGexprs[, 1],
          G = data$RGexprs[, 3],
          Rb = data$RGexprs[, 2],
          Gb = data$RGexprs[, 4]
        )
      )
    offset <-
      min(apply(as.matrix(data$RGexprs[, c(2, 4)], ncol = 2), 2, min)) + 1
    core <-
      backgroundCorrect(RG, method = "subtract", offset = offset)
    
    if (BGmethod == "normexp") {
      RG <- new("RGList", list(R = core[[1]], G = core[[2]]))
      core <-
        backgroundCorrect(RG, method = "normexp", offset = offset)
    }
    
    new <-
      list(
        index = data$index,
        samples = data$samples,
        batch = data$batch,
        Size = testarea,
        corrected.exprs = matrix(cbind(core[[1]], core[[2]]), ncol = 2),
        corrected.transformed.exprs = doTransform(matrix(cbind(core[[1]], core[[2]]), ncol = 2), transformation)[[1]],
        correctedAreas = testarea,
        areacut = areacut,
        transformation = transformation,
        image.type = data$image.type,
        dateIndex = data$dateIndex,
        single.batch.analysis = 0,
        BGmethod = BGmethod,
        maxMix = NULL,
        prior.pi = NULL,
        flex.reps = NULL,
        flexmethod = NULL,
        RNG = NULL
        
      )
    
    existingAreas <- which(new$Size != areacut)
    nonexistingAreas <- which(new$Size == areacut)
    signals <- new$corrected.transformed.exprs
    u.areas <- log(new$Size)
    
    res <-
      lm(u.areas[existingAreas] ~ signals[existingAreas, 1] + signals[existingAreas, 2])
    c.area <- u.areas
    if (length(nonexistingAreas) > 0) {
      for (i in 1:length(nonexistingAreas)) {
        c.area[nonexistingAreas[i]] <-
          res[[1]][[1]] + res[[1]][[2]] * signals[nonexistingAreas[i], 1] + res[[1]][[3]] *
          signals[nonexistingAreas[i], 2]
      }
    }
    new$correctedAreas <- c.area
    
    d12 <-
      data.frame(new$corrected.transformed.exprs[, 1],
                 new$corrected.transformed.exprs[, 2])
    names(d12) <- data$image.type[2:3]
    d12.m <- melt(d12, measure.vars = data$image.type[2:3])
    plot_dens <-
      ggplot(d12.m, aes_string(x = 'value', colour = 'variable')) + geom_line(stat =
                                                                                "density", size = 1.5) +
      labs(x = paste(transformation, " corrected signal", sep = "")) +
      ggtitle(paste("Corrected signals: Run ", unique(data$batch[, 1]), sep =
                      "")) +
      theme_bw() + guides(fill = guide_legend(override.aes = list(colour = NULL))) +
      scale_fill_discrete(guide = "legend", labels = data$image.type[2:3])
    
    logsig <- data.frame(new$corrected.transformed.exprs)
    plot_signal <-
      ggplot(logsig, aes_string(x = 'X1', y = 'X2')) + geom_point(size = 5) +
      labs(
        x = paste(
          transformation,
          " corrected ",
          data$image.type[2],
          " signal",
          sep = ""
        ),
        y = paste(
          transformation,
          " corrected ",
          data$image.type[3],
          " signal",
          sep = ""
        )
      ) +
      theme_bw() + theme(legend.position = "none")
    
    if (savePlot != "OFF" & savePlot != "screen") {
      pdf(
        paste(
          savePlot,
          "/getFluo_byRun",
          unique(as.numeric(data$batch[, 2])),
          "_",
          data$dateIndex,
          ".pdf",
          sep = ""
        ),
        paper = "a4r"
      )
    }
    
    if (savePlot != "OFF") {
      suppressWarnings(print(plot_dens))
      suppressWarnings(print(plot_signal))
      if (savePlot != "screen")
        dev.off()
    }
    
    colnames(new$corrected.exprs) <-
      colnames(new$corrected.transformed.exprs) <-
      data$image.type[2:3]
    
    return(new)
  }


#' pathEstimator
#'
#' It reads the generated groups of Fluo_inspection() and estimates the path cell progression
#'   given a user-defined expected pattern. It can also join some of the groups into a single one (manual
#'   selection is required).
#'
#' @param data List. The output of Fluo_inspection().
#' @param path.start Integer. A cluster number indicating the starting cluster that algorithm should use to
#'   build the path. The cluster numbers refer to the plot generated by Fluo_inspection(). Default is 1.
#'   If path.type = "circular" the number does not matter. If path.type = "A2Z" the user should inspect the
#'   Fluo_inspection() plot to detect the beginning of the path. If path.type = "other", the function will
#'   not estimate a path. The user has to manually insert the path progression (the cluster numbers) in
#'   Fluo_modeling().
#' @param path.type Character vector. A user-defined vector that characterizes the cell progression dynamics.
#'   The first element can be either "circular" or "A2Z" or "other". If "circular" the path progression is
#'   assummed to exhibit a circle-like behavior. If "A2Z" the path is assumed to have a well-defined start
#'   and a well-defined end point (e.g. a linear progression). If "other" the progression is assumed to be
#'   arbitrary without an obvious directionality. Default is "circular".
#
#'   The second element can be either "clockwise" or "anticlockwise" depending on how the path is expected
#'   to proceed. Default is "clockwise". If the first element is "other" the second element can be ommited.
#'
#'   If path.type = "other", the function does not estimate a path. The exact path has to be manually inserted
#'   in Fluo_modeling().
#' @param joinedGroups List. A list of cluster numbers to join. E.g. list(c(2,4)) joins cluster 2 and 4 as depicted
#'   in the Fluo_inspection() plot. Alternatively, list(c(2,4),c(1,6)) joins cluster 2 and 4 and clusters 1 and 6
#'   as depicted in the Fluo_inspection() plot.Each list entry should contain 2 groups. Default is NULL.
#'
#' @return The list of adjusted signal estimates, a progression path and the defined path type. The output is
#'   essentially the output of Fluo_inspection() with the addition of the following components:
#'   Path: the estimated path (visualized in the Fluo_Inspection() helper plot).
#'   path.type: the path.type that has been used to estimate the path.
#'
#' @export
#'
#' @examples
#' step3.1 <- pathEstimator(step3,path.start=6,path.type=c("circular","clockwise"))
pathEstimator <-
  function(data,
           path.start = 1,
           path.type = c("circular", "clockwise"),
           joinedGroups = NULL) {
    if (!is.null(joinedGroups)) {
      if (!is.list(joinedGroups)) {
        joinedGroups <- list(joinedGroups)
      }
      if (min(unlist(lapply(joinedGroups, length))) <= 1) {
        stop("No groups to join: use more than one group in the analysis")
      }
      for (k in 1:length(joinedGroups)) {
        if (k > 1) {
          for (i in 1:length(joinedGroups[[k]])) {
            joinedGroups[[k]][i] <-
              o[which(o[, 1] == joinedGroups[[k]][i]), 2] - 100
          }
        }
        ww <- c()
        for (i in 1:length(joinedGroups[[k]])) {
          ww <- c(ww, which(data$GAPgroups[, 1] == joinedGroups[[k]][i]))
        }
        data$GAPgroups[ww, 1] <- max(data$GAPgroups[, 1]) + 1
        newcentr <-
          c(max(data$GAPgroups[, 1]), apply(matrix(data$centroids[joinedGroups[[k]], 2:3], ncol =
                                                     2), 2, mean))
        data$centroids <- data$centroids[-joinedGroups[[k]],]
        data$centroids <-
          matrix(rbind(data$centroids, newcentr), ncol = ncol(data$centroids))
        
        o <-
          matrix(cbind(data$centroids[, 1], (order(
            data$centroids[, 1]
          ) + 100)), ncol = 2)
        for (i in 1:nrow(o)) {
          w <- which(data$GAPgroups[, 1] == o[i, 1])
          data$GAPgroups[w, 1] <- o[i, 2]
        }
        data$GAPgroups[, 1] <- data$GAPgroups[, 1] - 100
        data$centroids[, 1] <- o[, 2] - 100
      }
      
      plot(
        data$corrected.transformed.exprs,
        cex = 0.5,
        xlim = c(min(
          c(
            data$corrected.transformed.exprs[, 1],
            data$corrected.transformed.exprs[, 2]
          )
        ),
        max(
          c(
            data$corrected.transformed.exprs[, 1],
            data$corrected.transformed.exprs[, 2]
          )
        )),
        ylim = c(min(
          c(
            data$corrected.transformed.exprs[, 1],
            data$corrected.transformed.exprs[, 2]
          )
        ),
        max(
          c(
            data$corrected.transformed.exprs[, 1],
            data$corrected.transformed.exprs[, 2]
          )
        )),
        xlab = paste(
          data$transformation,
          " corrected ",
          data$image.type[3],
          " signal",
          sep = ""
        ),
        ylab = paste(
          data$transformation,
          " corrected ",
          data$image.type[2],
          " signal",
          sep = ""
        ),
        main = paste(
          data$clusterFUN,
          " cell progression phases with joined groups (inspection)",
          sep = ""
        ),
        sub = "the numbers label the estimated groups"
      )
      text(data$centroids[, 2], data$centroids[, 3], data$centroids[, 1])
      
    }
    path = c()
    if (path.type[1] != "other") {
      if (length(path.start) > 0) {
        cdata <- matrix(data$centroids, ncol = ncol(data$centroids))
        oc <- apply(matrix(cdata[, 2:3], ncol = 2), 2, mean)
        #      ocall<-apply(matrix(cdata[,2:3],ncol=2),2,quantile,seq(0.25,0.75,0.1))
        cdata[, 2] <- cdata[, 2] - oc[1]
        cdata[, 3] <- cdata[, 3] - oc[2]
        trigs <-
          matrix(cbind(cdata, t(apply(
            matrix(cdata[, 2:3], ncol = 2), 1, trigofun
          ))), nrow = nrow(cdata))
        path <-
          estimatePath(trigs, type = path.type[2], start = path.start)
      }
    } else {
      print(
        "The path has not been estimated. Input the correct path at init.path parameter of Fluo_modeling()"
      )
    }
    return(c(data, list(Path = path, Path.type = path.type)))
  }


#' cluster2outlier
#'
#' It turns one or more selected clusters to outlier clusters, i.e. clusters consisting of outlying corrected
#'   signals.
#'
#' @param data List. The output of Fluo_inspection().
#' @param out.cluster Numeric vector. The cluster number(s) to be turned into outlier clusters.
#'
#' @return A list of corrected fluorescence signal estimates with the selected clusters turned into outlier
#'   clusters.
#'
#' @export
#' 
#' @examples
#' ### here we (erroneously) assume that cluster 1 is an outlier and we flag it so below
#' step3.withoutliers <- cluster2outlier(step3,out.cluster=1)
#'
#' ### the outlier samples can be removed by FluoSelection_byRun()
#' step3.withoutliers <- FluoSelection_byRun(step3.withoutliers,
#'                                  other=which(step3.withoutliers$GAPgroups[,1]!=-999))
cluster2outlier <- function(data, out.cluster) {
  for (i in 1:length(out.cluster)) {
    w1 <- which(data$centroids[, 1] == out.cluster[i])
    data$centroids[w1, 1] <- -999
    w2 <- which(data$GAPgroups[, 1] == out.cluster[i])
    data$GAPgroups[w2, 1] <- -999
    data$GAPgroups[w2, 2] <- 2
  }
  return(data)
}


#' Fluo_inspection
#'
#' It generates the initial cell clusters as defined by their corrected fluorescence signals. The clusters
#'   can be generated by k-means (with GAP statistic estimated number of clusters) or by flow
#'   cytometry based approaches. This function shows the number and the characteristics of the
#'   initial groups and help us inspect cells' progression type for pathEstimator().
#'
#' @param data List. The output of getFluo() or getFluo_byRun().
#' @param altFUN Character string. A user-defined method to generate the initial clusters. It can be one of
#'   kmeans, samSpec, fmeans,fmerge or fpeaks. Default is "kmeans".
#' @param fixClusters Integer. A number that defines the number of k-mean clusters to be initially generated.
#'   If 0, the function runs GAP analysis to estimate the optimal number of clusters. Default is 0.
#' @param SAM.sigma Integer. A value for the sigma parameter of SamSPECTRAL algorithm. Default is 200.
#' @param k.max Integer. This is the maximum number of clusters that can be generated by k-means (if
#'   fixClusters = 0). Default is 15.
#' @param B.kmeans Integer. The number of bootstrap samples for the calculation of the GAP statistic. Default is 50.
#' @param savePlot Character string. Directory to store the plots. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen or "OFF" that does not generate a plot (suggested
#'   only during cross-validations). Default is the current working directory, getwd().
#' @param seed Integer. An optional seed number for the Random Number Generator. Note that this seed is a 'reference'
#'   value of the actual seed used in sampling. CONFESS is using various random sampling methods. Each method's
#'   actual seed is factor*seed. The factors vary across methods. Default is NULL.
#'
#' @return A list of corrected fluorescence signal estimates and a helper plot for deciding the number of groups and
#'   the cell progression path. The output is essentially the output of getFluo() or getFluo_byRun() with the addition
#'   of the following components:
#'   GAPgroups: the groups estimated by one of the altFUN methods are depicted in the first column. The second column contains
#'     1s for non-outlier signals and 2s for outlier signals (as estimated by each of the methods).
#'   clusterFUN: the altFUN method that has been used for clustering.
#'   normal.sigma: the sigma parameter of samSpec method.
#'   centroids: the 2 dimensional medians (centroids) of the estimated clusters.
#'   fixClusters: the fixClusters parameter used.
#'   Kmax: the k.meax parameter used.
#'   B.kmeans: the B.kmeans parameter used
#'
#' @import parallel
#'
#' @export
#'
#' @examples
#' step3 <- Fluo_inspection(data=step2.1,altFUN="kmeans",B.kmeans=5,savePlot="OFF")
Fluo_inspection <-
  function(data,
           altFUN = "kmeans",
           fixClusters = 0,
           SAM.sigma = 200,
           k.max = 15,
           B.kmeans = 50,
           savePlot = getwd(),
           seed = NULL) {
    if (is.null(data$RNG)) {
      seed <- seed
    } else {
      seed <- data$RNG
    }
    if (fixClusters > k.max) {
      fixClusters <- k.max
    }
    if (fixClusters < 0) {
      fixClusters <- 0
      print("fixClusters was set to 0")
    }
    if (savePlot != "OFF" & savePlot != "screen") {
      if (as.numeric(file.access(savePlot)) < 0) {
        stop("The savePlot directory cannot be found in the system or invalid value for savePlot")
      }
    }
    
    step <-
      GAPanalysis(
        data = data,
        fixClusters = fixClusters,
        sigma = SAM.sigma,
        altFUN = altFUN,
        k.max =
          k.max,
        B.kmeans = B.kmeans,
        savePlot = savePlot,
        seed = seed
      )
    step <-
      FluoInspection(data = step,
                     savePlot = savePlot,
                     dateIndex = data$dateIndex)
    return(c(
      step,
      list(
        fixClusters = fixClusters,
        Kmax = k.max,
        B.kmeans = B.kmeans
      )
    ))
  }

#' Fluo_modeling
#'
#' It takes the initial groups and the path progression and estimates the pseudotimes of cell
#'   progression and the associated change-points (updated cell clusters).
#'
#' @param data List. The output of pathEstimator().
#' @param init.path Numeric vector. The cell path progression as it has been estimated by
#'   pathEstimator() or a user-defined path that can be deduced from Fluo_inspection(). The latter
#'   is suggested only when path.type = "other" in pathEstimator().
#' @param VSmethod Character string. The variance stabilization transformation method to be applied
#'   to the corrected fluorescence data prior to the change point analysis. IT can be one of "log"
#'   or "DDHFmv". Default is "DDHFmv".
#' @param CPmethod Character string. The change point method to be used. It can be one of "ECP",
#'   (non-parametric) "manualECP" (non-parametric with user-defined numner of change-points) or
#'   "PELT" (Pruned Exact Linear Time; parametric). Default is ECP.
#' @param CPgroups Integer. The number of change-points to be kept if CPmethod = "manualECP".
#'   Default is 5.
#' @param CPpvalue Float. The significance level below which we do not reject a change point.
#'   Default is 0.05.
#' @param CPmingroup Integer. The minimum number of values for a cluster re-estimated by the
#'   change-point analysis. Default is 10.
#' @param seed Integer. An optional seed number for the Random Number Generator. Note that this seed is a 'reference'
#'   value of the actual seed used in sampling. CONFESS is using various random sampling methods. Each method's
#'   actual seed is factor*seed. The factors vary across methods. Default is NULL.
#'
#' @return A list of corrected fluorescence signal estimates, the pseudotimes and the cell progression clusters.
#'   The output is essentially the output of pathEstimator() with the addition of the following components:
#'   UpdatedPath: the updated progression path after re-estimation by change points and clustering.
#'   DataSorts: a matrix contains the calculated distances by orthogonal projection and the pseudotimes.
#'   DDHFupdate: it takes TRUE or FALSE to signify whether the clustering/pseudotime estimation has been updated by the re-estimation
#'     procedure.
#'   corrected.VStransformed.exprs: the background and run effect transformed corrected channel signals (by one of "log" or "DDHFmv"). The
#'     transformation is defined in the VSmethod parameter.
#'   VSmethod: the transformation that has been applied to the channel signals.
#'   Progression: it describes the estimated progression by the pseudotimes (first column) and the differences between the transformed channel
#'     signals.
#'   Updated.groups: the final clusters.
#'   CPs: the final change points detected.
#'   CPmethod: the CPmethod parameter used.
#'   CPsig: the CPpvalue parameter used.
#'   CPgroups: the CPgroups parameter used.
#'   CPmingroup: the CPmingroup parameter used.
#'
#' @export
#'
#' @examples
#' step4<-Fluo_modeling(data=step3.1,init.path=step3.1$Path,VSmethod="DDHFmv",
#'                      CPmethod="ECP",CPpvalue=0.01)
Fluo_modeling <-
  function(data,
           init.path,
           VSmethod = "DDHFmv",
           CPmethod = "ECP",
           CPgroups = 5,
           CPpvalue = 0.05,
           CPmingroup = 10,
           seed = NULL) {
    type <- data$Path.type
    
    if (is.null(data$RNG)) {
      seed <- seed
    } else {
      seed <- data$RNG
    }
    
    if (length(data$Path) == 0) {
      type <- c("A2Z", "clockwise")
    }
    
    #give priority to init.path
    if (!identical(data$Path, init.path) & length(init.path) > 0) {
      data$Path <- init.path
    }
    if (length(init.path) == 0) {
      init.path <- data$Path
    }
    
    if (length(data$Path) == 0 & length(init.path) == 0) {
      stop("The path has not been specified. Input the correct path at init.path parameter! ")
    }
    
    
    
    if (VSmethod != "log" & VSmethod != "DDHFmv") {
      stop("This variance stabilization method is not supported")
    }
    if (CPmethod != "PELT" & CPmethod != "ECP" &
        CPmethod != "manualECP") {
      stop("This change-point method is not supported")
    }
    
    step <- fixPath(data = data, groups = init.path)
    step <- orderFluo(data = step, path.type = type)
    step <- transformFluo(data = step, method = VSmethod)
    step <-
      cpoints(
        data = step,
        thresh = CPmingroup,
        cmethod = CPmethod,
        sig.level = CPpvalue,
        Q = CPgroups,
        path.type = type,
        seed = seed
      )
    colnames(step$DataSorts) <- c("Distance", "Pseudotime")
    colnames(step$Progression) <-
      c("Pseudotime", "transf.Difference")
    colnames(step$corrected.VStransformed.exprs) <-
      data$Data$image.type[2:3]
    return(step)
  }


#' Fluo_ordering
#'
#' It produces the final output table of CONFESS. It includes the Sample IDs, the Run IDs, the estimated cell areas
#'   (image analysis), the corrected fluorescence signals of both channels (run and background adjustED), the
#'   pseudotimes of cell progression, the final cell clusters and other statistics of cell progression analysis.
#'
#' @param data List. The outut of Fluo_modeling().
#' @param den.method Character string. A method to denoise the transformed channel signal differences (used for change-point analysis).
#'   The denoising obtains the residuals that can be subjected to statistical testing (model assumptions). It is one of
#'   "splines", "wavelets" or "lregr" (linear regression). Default is "wavelets".
#' @param savePlot Character string. Directory to store the plots. Its value can be an existing directory
#'   or "screen" that prints the plot only on the screen or "OFF" that does not generate a plot (suggested
#'   only during cross-validations). Default is the current working directory, getwd().
#'
#' @return The list of final results in two components:
#'   Summary_results:
#'     It contains a matrix that summarizes the findings of CONFESS. It has the index number of each sample, the sample IDs, the run IDs,
#'     the estimated cell size, the estimated run corrected cell size, the estimated pseudotime, the log, and if specified, DDHFmv transformed
#'     channel signals, the log or DDHFmv transformed channel differences, the estimated clusters, the residuals and a column flagginf outlier samples.
#'
#'   Analytical results:
#'     It contains all the components of Fluo_modeling() with the addition of:
#'     Outliers: a vector having "normal" for non-outlier samples and "outlier" for outlier samples. The outliers are estimated by Grubbs
#'       statistic based on their distance from the bulk of the clustered samples.
#'     Residuals: the residuals of the fitted model for the denoising of the corrected transformed channel differences (see parameter den.method).
#'     Residuals_diagnostics: various normality tests for the estimated residuals.
#'
#'   The component of
#'
#' @export
#'
#' @examples
#' step5<-Fluo_ordering(data=step4,savePlot="OFF")
Fluo_ordering <- function(data,
                          den.method = "wavelets",
                          savePlot = "OFF") {
  if (savePlot != "OFF" & savePlot != "screen") {
    if (as.numeric(file.access(savePlot)) < 0) {
      stop("The savePlot directory cannot be found in the system or invalid value for savePlot")
    }
  }
  step <- DDHFfit(data = data,
                  den.method = den.method,
                  savePlot = savePlot)
  step <- diagnoseResiduals(data = step, savePlot = savePlot)
  if (data$VSmethod == "DDHFmv") {
    result <-
      cbind(
        step$index,
        step$samples,
        step$batch[, 1],
        step$Size,
        step$correctedArea,
        step$Progression[, 1],
        step$corrected.transformed.exprs,
        step$corrected.VStransformed.exprs,
        step$Progression[, 2],
        step$Updated.groups,
        step$Residuals,
        step$Outliers
      )
    colnames(result) <-
      c(
        "Index",
        "Samples",
        "Runs",
        "Size",
        "Corr.Size",
        "Pseudotime",
        paste("log(", data$image.type[3], ")", sep = ""),
        paste("log(", data$image.type[2], ")", sep = ""),
        paste("DDHFmv(", data$image.type[3], ")", sep = ""),
        paste("DDHFmv(", data$image.type[2], ")", sep = ""),
        paste(
          "DDHFmv(",
          data$image.type[3],
          ")-",
          "DDHFmv(",
          data$image.type[2],
          ")",
          sep = ""
        ),
        "Clusters",
        "Residuals",
        "Out.Index"
      )
  } else {
    result <-
      cbind(
        step$index,
        step$samples,
        step$batch[, 1],
        step$Size,
        step$correctedArea,
        step$Progression[, 1],
        step$corrected.VStransformed.exprs,
        step$Progression[, 2],
        step$Updated.groups,
        step$Residuals,
        step$Outliers
      )
    colnames(result) <-
      c(
        "Index",
        "Samples",
        "Runs",
        "Size",
        "Corr.Size",
        "Pseudotime",
        paste("log(", data$image.type[3], ")", sep = ""),
        paste("log(", data$image.type[2], ")", sep = ""),
        paste(
          "log(",
          data$image.type[3],
          ")-",
          "log(",
          data$image.type[2],
          ")",
          sep = ""
        ),
        "Clusters",
        "Residuals",
        "Out.Index"
      )
  }
  return(list(Summary_results = result, Analytical_results = step))
}


#' FluoSelection_byRun
#'
#' It accepts a subset of data to inspect their background corrected fluorescence signal characteristics.
#'   Typically it one can inout the data from a single run to identify an appropriate mixture model for
#'   run effect correction. Any other arbitrary subset of the data can also be used. For example, it can
#'   be used to keep certain samples and filter out outliers.
#'
#' @param data List. The output of createFluo().
#' @param batch Integer. A selected run. If it is c() then the "other" parameter should be activated. Default
#'   is 1.
#' @param other Numeric vector. It accepts the sample numbers indicating the samples to be kept for analysis,
#'   e.g. other = c(1:10, 101:110) to keep samples 1:10 and 100:110. Default is c().
#'
#' @return A list of reformed data to be used in subsequent analysis. It is essentially the same slots of createFluo()
#'  with only a subset of data included (as defined by the batch and other parameters):
#'   index: The sample indices.
#'   RGexprs: the foreground (columns 1 and 3) and background (columns 2 and 4) signals of each channel that
#'     have been estimated by spotEstimator() and filtered in LocationMatrix().
#'   samples: the sample IDs.
#'   batch: a matrix of the run IDs. The first column contains the original run IDs. The second column is the converted
#'     original IDs into numeric values (to be used in the statistical modeling step of Fluo_adjustment()).
#'   size: the estimated cell size.
#'   image.type: the image type IDs as defined in readFiles(). The parameter is kept in ordeer to enable the user to
#'     use this function independently of the image analysis step.
#'   dateIndex: a date index to be used for storing the output files. It is either transfered from LocationMatrix() or
#'     it is generated here for the first time (e.g. if image analysis was not run by CONFESS or if the analysis has
#'     been repeated many times).
#'
#' @export
#'
#' @examples
#' step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#'                    package = "CONFESS"),separator="_")
#' step2a <- FluoSelection_byRun(data = step1, batch = 4:5)
FluoSelection_byRun <- function(data,
                                batch = c(),
                                other = c()) {
  if (length(batch) > 0) {
    if (sum(batch %in% as.numeric(data$batch[, 2])) == 0) {
      stop("The batch number does not exist")
    } else if (sum(batch %in% as.numeric(data$batch[, 2])) != length(batch)) {
      message("One or more batch numbers do not exist.")
    }
  }
  
  if (length(batch) > 0 & length(other) > 0) {
    stop("Filtering by both batch and other parameters is not currently supported")
  }
  if (length(batch) == 0 & length(other) == 0) {
    stop("All filtering variables are NULL. Filtering failed!")
  }
  
  #   n <- length(data$index)
  #   l <- rep(0,n)
  if (length(batch) > 0) {
    w <- which(as.numeric(data$batch[, 2]) %in% batch)
  } else {
    w <- other
  }
  
  d1 <- c("index",
          "samples",
          "Size",
          "correctedAreas",
          "Updated.groups")
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
  for (i in 1:length(data)) {
    if (names(data)[i] %in% d1)
      data[[i]] <- data[[i]][w]
    else if (names(data)[i] %in% d2)
      data[[i]] <- data[[i]][w,]
  }
  data$index <- 1:length(data$index)
  
  cent <- as.numeric(names(data) == "centroids")
  if (sum(cent) > 0) {
    w <- which(data[[which(cent == 1)]][, 1] == -999)
    if (length(w) > 0) {
      data[[which(cent == 1)]] <- data[[which(cent == 1)]][-w,]
    }
  }
  
  #fix the batch numbers
  uu1 <- unique(as.numeric(data$batch[, 2]))
  data$batch[, 2] <- as.numeric(factor(data$batch[, 2]))
  uu2 <- unique(as.numeric(data$batch[, 2]))
  if (!identical(uu1, uu2)) {
    for (i in 1:length(uu2)) {
      print(paste("Run ", uu1[i], " was converted into ", uu2[i], sep = ""))
    }
  }
  
  
  return(data)
}


#' Fluo_CV_prep
#'
#' It generates the data that will be used in the cross-validation analysis. Essentialy, it analyzes and stores the original (full)
#'   dataset for different reference runs, seeds, starting clusters etc. It estimates the progression path automatically that is
#'   feasible only for standard paths (path.type parameter different than 'other'). For this reason this function is useful only
#'   in these cases. If otherwise, it should be ommitted from the analysis and the user is should generate it manually, i.e. run
#'   Fluo_adjustment() - Fluo_modeling() series as many times as the cases to be studied with manual init.path input in Fluo_modeling().
#'
#' The function can also be used to generate all pseudotime/clustering results up to the function of Fluo_modeling() but the starting
#'   cluster has to be defined in general terms (see init.path parameter below). For this reason, its parameters are essentially the same
#'   to the ones defined previously at the Fluo_adjustment() - Fluo_modeling() functions.
#'
#' @param data List. The output of crearteFluo(), i.e. the image analysis estimates.
#' @param init.path Character vector. It defines the starting cluster of the progression path in general terms.
#'   It can be one of "top/right", "top/left", "bottom/right" or "bottom/left" indicating the cluster of interest
#'   on the 2d scatterplot of Fluo_inspection(). Default is rep("bottom/left",2), i.e. in Fucci an EM/earlyG1 like
#'   cluster.
#' @param path.type Character vector. A user-defined vector that characterizes the cell progression dynamics.
#'   The first element can be either "circular" or "A2Z" or "other". If "circular" the path progression is
#'   assummed to exhibit a circle-like behavior. If "A2Z" the path is assumed to have a well-defined start
#'   and a well-defined end point (e.g. a linear progression). If "other" the progression is assumed to be
#'   arbitrary without an obvious directionality. Default is "circular".
#
#'   The second element can be either "clockwise" or "anticlockwise" depending on how the path is expected
#'   to proceed. Default is "clockwise". If the first element is "other" the second element can be ommited.
#'
#'   If path.type = "other", the function does not estimate a path. The cross-validation algorithm will probably
#'   fail for this kind of path.type values because it will not be able to automatically guess the progression path.
#'   It is suggested that the user runs the cross-validation manually (each time specifying the path in Fluo_modeling()),
#'   collect the data in a list similar to the one produced here and input them into Fluo_CV_modeling() to get the results.
#' @param BGmethod Character string. The type of image background correction to be performed.
#'   One of "normexp" or "subtract". Default is "normexp".
#' @param maxMix Integer. The maximum number of components to fit into the mixture of
#'   regressions model. If maxMix=1 or if the the optimal number of the estimated components
#'   is 1, the model reduces to the classical 2-way ANOVA. Default is 3.
#' @param single.batch.analysis Numeric. The baseline run(s) to perform run effect correction with flexmix. Due to iterative
#'   nature of this function it can be a series of values includying 0 (averaging of run correction estimates). Default is 1:5.
#' @param transformation Character string. One of bc (Box-Cox), log, log10, asinh transforms applied to the data. Default is "log".
#' @param prior.pi Float. The prior probability to accept a component. Default is 0.1.
#' @param flex.reps Integer. The iterations of the Expectation-Maximization algorithm to estimate the flexmix
#'   model. Default is 50.
#' @param flexmethod Character string. A method to estimate the optimal number of flexmix
#'   components. One of "BIC", "AIC", "ICL". Default is "BIC".
#' @param areacut Integer. The "artificial" area size (BFarea^2) of the cells estimated
#'   by BF image modelling. Default is 0, implying that the area sizes to be corrected will
#'   by estimated automatically from the data (not recommended if prior knowledge exists).
#' @param fixClusters Integer. A number that defines the number of k-mean clusters to be initially generated.
#'   If 0, the function runs GAP analysis to estimate the optimal number of clusters. Default is 0.
#' @param altFUN Character string. A user-defined method to generate the initial clusters. It can be one of
#'   kmeans, samSpec, fmeans,fmerge or fpeaks. Default is "kmeans".
#' @param k.max Integer. This is the maximum number of clusters that can be generated by k-means (if
#'   fixClusters = 0). Default is 15.
#' @param VSmethod Character string. The variance stabilization transformation method to be applied
#'   to the corrected fluorescence data prior to the change point analysis. IT can be one of "log"
#'   or "DDHFmv". Default is "DDHFmv".
#' @param CPmethod Character string. The change point method to be used. It can be one of "ECP",
#'   (non-parametric) "manualECP" (non-parametric with user-defined numner of change-points) or
#'   "PELT" (Pruned Exact Linear Time; parametric). Default is ECP.
#' @param CPgroups Integer. The number of change-points to be kept if CPmethod = "manualECP".
#'   Default is 5.
#' @param B.kmeans Integer. The number of bootstrap samples for the calculation of the GAP statistic. Default is 50.
#' @param CPpvalue Float. The significance level below which we do not reject a change point.
#'   Default is 0.05.
#' @param CPmingroup Integer. The minimum number of values for a cluster re-estimated by the
#'   change-point analysis. Default is 10.
#' @param savePlot Character string. Directory to store the plots of the analysis of the whole data. Its
#'   value can be an existing directory or "screen" that prints the plot only on the screen. The "OFF"
#'   option is permanently used in cross-validations). Default is the current working directory, getwd().
#' @param seed Integer. An optional seed number for the Random Number Generator. Note that this seed is a 'reference'
#'   value of the actual seed used in sampling. CONFESS is using various random sampling methods. Each method's
#'   actual seed is factor*seed. The factors vary across methods. Default is NULL.
#'
#' @return The results of Fluo_modeling() for difference reference runs (batches) are stored in different slots. An additional slot
#'   @init.path exists that stores the init.path parameter (its value to be used in the CV automatically).
#'
#'   One can directly use the run components in Fluo_ordering() to finalize the data analysis. The main purpose of this function,
#'   though, is to prepare the data for cross-validation.
#'
#' @export
#'
#' @examples
#' step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#' package = "CONFESS"),separator="_")
#' steps2_4 <- Fluo_CV_prep(data=step1,init.path = "bottom/left",path.type=c("circular","clockwise"),
#' single.batch.analysis = 5,flex.reps=5,altFUN="kmeans",VSmethod="DDHFmv",CPmethod="ECP",
#' B.kmeans=5,CPpvalue=0.01,savePlot="OFF")
Fluo_CV_prep <-
  function(data,
           init.path = "bottom/left",
           path.type = c("circular", "clockwise"),
           BGmethod = "normexp",
           maxMix = 3,
           single.batch.analysis = 1:5,
           transformation = "log",
           prior.pi = 0.1,
           flex.reps = 50,
           flexmethod = "BIC",
           areacut = 0,
           fixClusters = 0,
           altFUN = "kmeans",
           k.max = 15,
           VSmethod = "DDHFmv",
           CPmethod = "ECP",
           CPgroups = 5,
           B.kmeans =
             50,
           CPpvalue = 0.05,
           CPmingroup = 15,
           savePlot = getwd(),
           seed = NULL) {
    mm <-
      match(single.batch.analysis,
            as.numeric(data$batch[, 2]),
            nomatch = -999)
    single.batch.analysis <-
      single.batch.analysis[which(mm != -999)]
    if (length(mm[mm == -999]) > 0) {
      message(
        "Warning: Some reference runs did not exist and have been removed! Check if the run indexing has changed."
      )
    }
    
    CVresults.full <- as.list(rep(0, length(single.batch.analysis)))
    names(CVresults.full) <-
      paste("Batch", single.batch.analysis, sep = "")
    
    for (i in 1:length(single.batch.analysis)) {
      olddate <- data$dateIndex
      dt.indx <-
        paste(data$dateIndex, "_ref", single.batch.analysis[i], sep = "")
      data$dateIndex <- dt.indx
      
      mes <-
        paste(
          "| Running the full data analysis for reference run ",
          single.batch.analysis[i],
          "... |",
          sep = ""
        )
      message(rep("-", length(unlist(
        strsplit(mes, "")
      ))))
      message(mes)
      message(rep("-", length(unlist(
        strsplit(mes, "")
      ))))
      if (length(unique(data$batch[, 2])) > 1) {
        step2 <-
          Fluo_adjustment(
            data = data,
            BGmethod = BGmethod,
            maxMix = maxMix,
            single.batch.analysis = single.batch.analysis[i],
            transformation = transformation,
            prior.pi = prior.pi,
            flex.reps = flex.reps,
            flexmethod = flexmethod,
            savePlot = savePlot,
            seed = seed
          )
        step2.1 <- getFluo(data = step2, areacut = areacut)
      } else {
        step2.1 <-
          getFluo_byRun(
            data = data,
            BGmethod = BGmethod,
            areacut = areacut,
            transformation = transformation,
            savePlot = savePlot
          )
      }
      
      step3 <-
        Fluo_inspection(
          data = step2.1,
          altFUN = altFUN,
          fixClusters = fixClusters,
          SAM.sigma = 200,
          k.max = k.max,
          B.kmeans = B.kmeans,
          savePlot = savePlot,
          seed = seed
        )
      path.start.full <-
        path.initiator(data = step3, where = init.path)
      step3.1 <-
        pathEstimator(
          step3,
          path.start = as.numeric(path.start.full[1]),
          path.type = path.type,
          joinedGroups = NULL
        )
      
      # keep this!
      CVresults.full[[i]] <-
        Fluo_modeling(
          data = step3.1,
          init.path = step3.1$Path,
          VSmethod = VSmethod,
          CPmethod = CPmethod,
          CPgroups = CPgroups,
          CPpvalue = CPpvalue,
          CPmingroup = CPmingroup,
          seed = seed
        )
      CVresults.full[[i]]$dateIndex <- olddate
    }
    
    return(c(CVresults.full, list(init.path = init.path)))
  }


#' Fluo_CV_modeling
#'
#' It performs the cross-validation analysis on the estimated pseudotimes and clusters of the previous step, i.e. Fluo_CV_prep() or
#'   a manually generated list based on Fluo_modeling(). This function will evaluate the change in the estimated obtained (i) from a subset of data
#'   by f-fold cross-validation where f is the percentage of the samples from a specific group (@GAPgroups) that stay in the analysis at each
#'   CV iteration, or (ii) from a subset of runs that stay in the analysis at each CV iteration. It produces informative plots for the differences
#'   in the estimates between each iteration and the original estimates. It also summarizes the CV-estimated pseudotimes into a new set of estimates.
#'
#' @param data List. The output of Fluo_CV_prep() or any other manually retrieved list with the components of Fluo_CV_prep().
#' @param B Integer. The number of cross-validation to be performed. Default is 20.
#' @param batch Numeric. A vector of runs to remain in the cross-validation. The rest are temporarily removed. The algorithm estimates the
#'   centroids of the reduced data and then calls the out-of-bag samples and re-estimates their k-mean clusters.
#' @param perc.cutoff Float. The percentage of similar CV-estimated pseudotimes for each sample. The similarity is assessed by k-means
#'   with k = 2. It serves as a cut-off to identify outlying CV-estimated pseudotimes (along with q and pseudotime.cutoff). Default is 0.6.
#' @param q Float. The q-th quantile of the difference between the original data estimated pseudotimes and the CV-estimated pseudotimes
#'   for each sample. It serves as a cut-off to identify outlying CV-estimated pseudotimes (along with perc.cutoff and pseudotime.cutoff).
#'   Default is 0.9.
#' @param f Float. The percentage of samples from each estimated cluster (@GAPgroups) to remain in the cross-validation analysis. The rest are
#'   temporarily removed. The algorithm estimates the centroids of the reduced data and then calls the out-of-bag samples and re-estimates their
#'   k-mean clusters.
#' @param seed.it Logical. If TRUE it performs cross-validation with the seed used in the analysis of the original data, i.e. in
#'   Fluo_CV_prep(). Default is TRUE.
#' @param pseudotime.cutoff Integer. A user-defined value to define outlier samples (along with perc.cutoff and q), i.e. samples with
#'   Pseudotime(original) - median{Pseudotime(CV)} > pseudotime.cutoff. Default is 20.
#' @param savePlot Character string. Directory to store the plots of the analysis of the whole data. Its
#'   value can be an existing directory or "screen" that prints the plot only on the screen. The "OFF"
#'   option is permanently used in cross-validations). Default is the current working directory, getwd().
#'
#' @return The output of Fluo_modeling() with the original estimates and the CV-based estimated pseudotimes/clusters in different slots of component CV results.
#'   The results are categorized by run number. Each run contains the original estimates (@Original Pseudotimes), the CV-based estimates by the "median/original"
#'   method (@Reest.Pseudotimes_median/original) and the CV-based estimates by the "median/null" method (@Reest.Pseudotimes_median/null).
#'
#'   1. "median/original"
#'   It integrates the information of the CV and the originally estimated pseudotimes. It build kmean clusters of the B CV estimates for each sample
#'   and defines pseudotime(i) = median(pseudotime(set1,i)) where set1 is a subset of the B pseudotimes that exhibit some similarity. The similarity
#'   is assessed by k-means clustering. This subset should contain a large percentage of the B data (>perc.cutoff) and it's median should be lower than
#'   the q-th quantile of the average differences between the original and the CV-estimated pseudotimes across all samples. If the CV estimated pseudotimes
#'   do not satisfy the above then the algorithm returns pseudotime(i) = median(pseudotime(set2,i)) where set2 is the cluster of B pseudotimes that minimizes
#'   |median(pseudotimes(set2,i))-original.pseudotimes|.
#'
#'   2. "median/null"
#'   if set1 with similar pseudotimes that satisfies the above rules exists, it returns the pseudotime(i) = median(pseudotime(set1,i)). Otherwise it returns
#'   NULL, i.e. the sample CV-estimated pseudotimes are not similar and the algorithm cannot estimate reliably the pseudotime of interest.
#'
#'   Both solutions are then going under a final round of change-point analysis that uses the CV-estimated pseudotimes and produce the final results of
#'   Fluo_CV_modeling(). All results canbe subsequently used in Fluo_ordering().
#'   The output also includes a second component, @All.Progressions, with the original and the CV estimated pseudotimes. This information is kept for comparison
#'   reasons and it is not used further.
#'
#' @import grDevices stats graphics
#' @export
#'
#' @examples
#' print("Not run because takes a long time")
#' #step1 <- createFluo(from.file=system.file("extdata", "Results_of_image_analysis.txt",
#' #package = "CONFESS"),separator="_")
#' #steps2_4 <- Fluo_CV_prep(data=step1,init.path = "bottom/left",path.type=c("circular","clockwise"),
#' #single.batch.analysis = 5,flex.reps=5,altFUN="kmeans",VSmethod="DDHFmv",CPmethod="ECP",
#' #B.kmeans=5,CPpvalue=0.01,savePlot="OFF")
#' #steps2_4cv<-Fluo_CV_modeling(data=steps2_4,B=5,f=0.99,savePlot="OFF")
Fluo_CV_modeling <-
  function(data,
           B = 20,
           batch = 1,
           perc.cutoff = 0.6,
           q = 0.9,
           f = 0.90,
           seed.it = TRUE,
           pseudotime.cutoff = 20,
           savePlot = getwd()) {
    # error if invalid perc.cutoff
    if (perc.cutoff < 0 | perc.cutoff > 1) {
      stop("The perc.cutoff should be a value between 0.5 and 1")
    }
    
    # error if invalid q
    if (q > 1 | q < 0) {
      stop("The q should be a value between 0 and 1")
    }
    
    # fix seed
    if (seed.it) {
      seed <- data[[1]]$RNG
    } else {
      seed <- NULL
    }
    
    if (B == 2) {
      stop(
        "Setting B = 2 does not generate enough data for cross-validation with random leave-outs. Consider higher B values!"
      )
    }
    if (B == 1 & length(batch) == 0) {
      stop(
        "Setting B = 1 requires the specification of leave-out runs. Modify parameter batch accordingly!"
      )
    }
    
    
    if (B == 1) {
      message("")
      message(
        paste(
          "Setting B = 1 performs run-based cross-validation only (keeps run(s) ",
          paste(batch, collapse = ","),
          ")",
          sep = ""
        )
      )
      message("")
      
      
      ### changed here Feb 2016
      ss <-
        matrix(which(match(
          as.numeric(data[[1]]$batch[, 2]), batch, nomatch = 0
        ) > 0), nrow = 1)
      ### changed here Feb 2016
      
    } else {
      message("")
      message(
        "Setting B > 2 removes (1-f)% of the samples from each group at random to perform cross-validation"
      )
      message("")
      
      # take random sampling
      ss <-
        matrix(sort(CVsampler(data = data[[1]], f = f)), nrow = 1)
      for (b in 2:B) {
        ss <-
          matrix(rbind(ss, matrix(sort(
            CVsampler(data = data[[1]], f = f)
          ), nrow = 1)), ncol = ncol(ss))
      }
    }
    
    Results <- CVresults <-
      progressions <- as.list(rep(0, (length(data) - 1)))
    names(Results) <- names(data)[1:(length(data) - 1)]
    
    if (savePlot != "OFF" & savePlot != "screen") {
      pdf(paste(
        savePlot,
        "/Fluo_CV_modeling_",
        data[[1]]$dateIndex,
        ".pdf",
        sep = ""
      ))
    }
    
    for (i in 1:(length(data) - 1)) {
      data.step2 <- data[[i]][1:18]
      
      # define the final results variable
      Results[[i]] <- as.list(rep(0, 3))
      names(Results[[i]]) <-
        c(
          "Original_Pseudotimes",
          "Reest.Pseudotimes_median_original",
          "Reest.Pseudotimes_median_null"
        )
      Results[[i]][[1]] <- data[[i]]
      Results[[i]][[2]] <- Results[[i]][[3]] <- data[[i]][1:32]
      
      
      # error if the altFUN is not kmeans
      if (data[[i]]$clusterFUN != "kmeans") {
        stop("Only altFUN = kmeans is currently supported!")
      }
      
      CVresults[[i]] <- as.list(rep(0, B))
      names(CVresults[[i]]) <- paste("CV", 1:B, sep = "")
      progressions[[i]] <-
        matrix(data[[i]]$Progression[, 1], ncol = 1)
      
      for (b in 1:B) {
        mes <-
          paste(
            "| Running the CV analysis for reference run ",
            data.step2$single.batch.analysis,
            " and CV ",
            b,
            "... |",
            sep = ""
          )
        message(rep("-", length(unlist(
          strsplit(mes, "")
        ))))
        message(mes)
        message(rep("-", length(unlist(
          strsplit(mes, "")
        ))))
        
        if (B == 1) {
          a <- FluoSelection_byRun(data = data.step2, batch = batch)
        } else {
          a <- FluoSelection_byRun(data = data.step2, other = ss[b,])
        }
        data.step3 <-
          Fluo_inspection(
            data = a,
            altFUN = "kmeans",
            fixClusters = data[[i]]$fixClusters,
            SAM.sigma = data[[i]]$normal.sigma,
            k.max = data[[i]]$Kmax,
            B.kmeans = data[[i]]$B.kmeans,
            seed = seed,
            savePlot = "OFF"
          )
        
        #predict and fix
        mm <- which(mm <-
                      match(data[[i]]$samples, a$samples, nomatch = 0) == 0)
        new.signals <- data[[i]]$corrected.transformed.exprs[mm,]
        new.clusters <-
          apply(matrix(new.signals, ncol = 2),
                1,
                predict.kmeans,
                centroid = data.step3$centroids[, 2:3])
        
        ### changed here Feb 2016
        data.step3new <-
          addKpredictions(
            whole = data[[i]][1:25],
            cv = data.step3,
            new.clusters = new.clusters,
            indices = list(ss[b,], mm)
          )
        ### changed here Feb 2016
        
        path.start.full <-
          path.initiator(data = data.step3new, where = data$init.path)
        data.step31 <-
          pathEstimator(
            data.step3new,
            path.start = as.numeric(path.start.full[1]),
            path.type = data[[i]]$Path.type,
            joinedGroups = NULL
          )
        
        # keep this!
        CVresults[[b]] <-
          Fluo_modeling(
            data = data.step31,
            init.path = data.step31$Path,
            VSmethod = data[[i]]$VSmethod,
            CPmethod = data[[i]]$CPmethod,
            CPgroups = data[[i]]$CPgroups,
            CPpvalue = data[[i]]$CPsig,
            CPmingroup = data[[i]]$CPmingroup,
            seed = seed
          )
        
        # keep this!
        progressions[[i]] <-
          matrix(cbind(progressions[[i]], matrix(CVresults[[b]]$Progression[, 1], ncol =
                                                   1)), nrow = nrow(progressions[[i]]))
      }
      
      compProgressions <-
        apply(
          matrix(progressions[[i]], ncol = ncol(progressions[[i]])),
          1,
          aveDiff,
          path.type = data[[i]]$Path.type[1],
          maxPseudo = max(data[[i]]$index)
        )
      medcomp <- median(compProgressions)
      
      sorter <- sorter.leg <- c()
      if (B == 1) {
        sorter <- as.numeric(data[[i]]$batch[, 2])
        sorter.leg <- "run index"
        out <-
          setdiff(unique(as.numeric(data.step2$batch[, 2])), batch)
        main.leg <- paste("re-estimate run ", out, sep = "")
      } else {
        sorter <- data[[i]]$Updated.groups
        main.leg <- "re-estimate random leave-outs"
        sorter.leg <- "estimated cluster"
      }
      
      barplot(
        compProgressions[sort.list(sorter)],
        col = sorter[sort.list(sorter)],
        main = paste("CV analysis: ", main.leg, sep = ""),
        sub = paste(
          "The dashed line at ",
          medcomp,
          " indicates the median difference",
          sep = ""
        ),
        ylab = "Pseudotime(original) - median{Pseudotime(CV)}",
        xlab = paste("Samples (sorted by ", sorter.leg, ")", sep = "")
      )
      abline(h = medcomp, lty = 2, lwd = 1.5)
      
      
      # re-estimate progressions based on pseudo.est.method
      cut <-
        as.numeric(quantile(compProgressions, seq(0, 1, 0.01))[(q * 100 + 1)])
      pseu <-
        matrix(cbind(progressions[[i]][, 1], compProgressions, progressions[[i]][, 2:ncol(progressions[[i]])]),
               nrow = nrow(progressions[[i]]))
      pseu <-
        apply(
          matrix(pseu, ncol = ncol(pseu)),
          1,
          reestimate.pseudos.byCV,
          diff.quantile = cut,
          perc.cutoff = perc.cutoff,
          pseudotime.cutoff = pseudotime.cutoff
        )
      pseu <- estimate.new.pseudotimes(pseu)
      
      Results[[i]][[2]]$DataSorts[, 2] <- pseu[1,]
      Results[[i]][[3]]$DataSorts[, 2] <- pseu[2,]
      Results[[i]][[3]] <-
        FluoSelection_byRun(data = Results[[i]][[3]], other = which(pseu[2,] >
                                                                      0))
      
      #get the new CP analysis
      Results[[i]][[2]] <-
        cpoints(
          data = Results[[i]][[2]],
          thresh = data[[i]]$CPgroups,
          cmethod = data[[i]]$CPmethod,
          sig.level = data[[i]]$CPsig,
          Q = data[[i]]$CPmingroup,
          path.type = "other",
          seed = seed
        )
      Results[[i]][[3]] <-
        cpoints(
          data = Results[[i]][[3]],
          thresh = data[[i]]$CPgroups,
          cmethod = data[[i]]$CPmethod,
          sig.level = data[[i]]$CPsig,
          Q = data[[i]]$CPmingroup,
          path.type = "other",
          seed = seed
        )
      
    }
    
    if (savePlot != "OFF" & savePlot != "screen") {
      dev.off()
    }
    
    
    return(list(CV_results = Results, All.Progressions = progressions))
  }
