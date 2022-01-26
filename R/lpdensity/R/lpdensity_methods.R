################################################################################
#' Print Method for Local Polynomial Density Estimation and Inference
#'
#' @description The print method for local polynomial density objects.
#'
#' @param x Class "lpdensity" object, obtained from calling \code{\link{lpdensity}}.
#' @param ... Additional options.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Estimate density and report results
#' print(lpdensity(data = X, bwselect = "imse-dpi"))
#'
#' @export
print.lpdensity <- function(x, ...) {

  cat("Call: lpdensity\n\n")

  cat(paste("Sample size                                      ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation    (p=)    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated            (v=)    ", x$opt$v,        "\n", sep=""))
  if (x$opt$q > 0) {
    cat(paste("Polynomial order for confidence interval (q=)    ", x$opt$q,        "\n", sep=""))
  } else {
    cat(paste("Polynomial order for confidence interval (q=)    ", "p",            "\n", sep=""))
  }
  cat(paste("Kernel function                                  ", x$opt$kernel,   "\n", sep=""))
  if (x$opt$scale != 1) {
    cat(paste("Scaling factor                                   ",  x$opt$scale,   "\n", sep=""))
  }
  cat(paste("Bandwidth method                                 ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  cat("Use summary(...) to show estimates.\n")
}

################################################################################
#' Summary Method for Local Polynomial Density Estimation and Inference
#'
#' @description The summary method for local polynomial density objects.
#'
#' @param object Class "lpdensity" object, obtained from calling \code{\link{lpdensity}}.
#' @param ... Additional options, including (i) \code{grid} specifies a subset of grid points
#'   to display results; (ii) \code{gridIndex} specifies the indices of grid points
#'   to display results; (iii) \code{alpha} specifies the significance level; (iv)
#'   \code{CIuniform} specifies whether displaying pointwise confidence intervals (\code{FALSE}, default) or
#'   the uniform confidence band (\code{TRUE}); (v) \code{CIsimul} specifies the number of simulations used
#'   to construct critical values (default is \code{2000}).
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Estimate density and report results
#' est1 <- lpdensity(data = X, bwselect = "imse-dpi")
#' summary(est1)
#'
#' # Report results for a subset of grid points
#' summary(est1, grid=est1$Estimate[4:10, "grid"])
#' summary(est1, gridIndex=4:10)
#'
#' # Report the 99% uniform confidence band
#' set.seed(42) # fix the seed for simulating critical values
#' summary(est1, alpha=0.01, CIuniform=TRUE)
#'
#' @export
summary.lpdensity <- function(object, ...) {

  x <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }
  if (is.null(args[['CIuniform']]))   { CIuniform <- FALSE } else { CIuniform <- args[['CIuniform']] }
  if (is.null(args[['CIsimul']]))   { CIsimul <- 2000 } else { sep <- args[['CIsimul']] }

  if (is.null(args[['grid']]) & is.null(args[['gridIndex']])) {
    gridIndex <- 1:nrow(x$Estimate)
  } else if (is.null(args[['grid']]) & !is.null(args[['gridIndex']])) {
    gridIndex <- args[['gridIndex']]
    if (is.null(gridIndex)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!all(gridIndex %in% 1:nrow(x$Estimate))) {
      stop(paste("Option gridIndex incorrectly specified. Should be integers between 1 and ", nrow(x$Estimate), ".\n", sep=""))
    }
  } else {
    grid <- args[['grid']]
    if (is.null(grid)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!is.numeric(grid)) {
      stop("Option grid incorrectly specified.\n")
    } else {
      gridIndex <- rep(NA, length(grid))
      if (min(grid) < min(x$Estimate[, "grid"]) | max(grid) > max(x$Estimate[, "grid"])) {
        warning("The reporting range exceeds the original estimation range. Option summary(..., grid=) should be within the estimation range specified by lpdensity(..., grid=).\n")
      }
      for (j in 1:length(grid)) {
        gridIndex[j] <- which.min(abs(x$Estimate[, "grid"]-grid[j]))
        #gridIndex <- unique(gridIndex)
      }
    }
  }

  cat("Call: lpdensity\n\n")

  cat(paste("Sample size                              (n=)    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation    (p=)    ", x$opt$p,        "\n", sep=""))
  if (x$opt$v == 0) {
  cat(paste("Distribution function estimated          (v=)    ", x$opt$v,        "\n", sep=""))
  } else if (x$opt$v == 1) {
  cat(paste("Density function estimated               (v=)    ", x$opt$v,        "\n", sep=""))
  } else {
  cat(paste("Order of derivative estimated            (v=)    ", x$opt$v,        "\n", sep=""))
  }
  if (x$opt$q > x$opt$p) {
  cat(paste("Polynomial order for confidence interval (q=)    ", x$opt$q,        "\n", sep=""))
  } else {
  cat(paste("Polynomial order for confidence interval (q=)    ", "p",            "\n", sep=""))
  }
  cat(paste("Kernel function                                  ", x$opt$kernel,   "\n", sep=""))
  if (x$opt$scale != 1) {
  cat(paste("Scaling factor                                   ",  x$opt$scale,   "\n", sep=""))
  }
  cat(paste("Bandwidth selection method                       ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  ### compute CI
  if (CIuniform) {
    if (length(CIsimul) == 0) { CIsimul <- 2000 }
    if (!is.numeric(CIsimul) | is.na(CIsimul)) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform <- FALSE
      z_val <- qnorm(1 - alpha/2)
    } else if (ceiling(CIsimul)<2) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform <- FALSE
      z_val <- qnorm(1 - alpha/2)
    } else {
      CIsimul <- ceiling(CIsimul)
      if (x$opt$q == x$opt$p) {
        x$CovMat_q <- x$CovMat_p
        x$Estimate[, "se_q"] <- x$Estimate[, "se_p"]
      }
      corrMat <- sweep(sweep(x$CovMat_q, MARGIN=1, FUN="*", STATS=1/x$Estimate[, "se_q"]), MARGIN=2, FUN="*", STATS=1/x$Estimate[, "se_q"])
      normalSimu <- try(
        mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=corrMat),
        silent=TRUE)
      if (is.character(normalSimu)) {
        print(normalSimu)
        warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
        CIuniform <- FALSE
        z_val <- qnorm(1 - alpha/2)
      } else {
        z_val <- quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha)
      }
    }
  } else {
    z_val <- qnorm(1 - alpha/2)
  }

  if (x$opt$q == x$opt$p) {
    CI_l <- x$Estimate[, "f_p"] - x$Estimate[, "se_p"] * z_val
    CI_r <- x$Estimate[, "f_p"] + x$Estimate[, "se_p"] * z_val
  } else {
    CI_l <- x$Estimate[, "f_q"] - x$Estimate[, "se_q"] * z_val
    CI_r <- x$Estimate[, "f_q"] + x$Estimate[, "se_q"] * z_val
  }

  flagNotAllLocalSampleSizeSame <- !all(x$Estimate[, "nh"] == x$Estimate[, "nhu"])

  ### print output
  cat(paste(rep("=", 14 + 10 + 8 + 8*flagNotAllLocalSampleSizeSame + 10 + 10 + 25), collapse="")); cat("\n")

  cat(format(" ", width= 14 ))
  cat(format(" ", width= 10 ))
  cat(format(" ", width= 8 + 8*flagNotAllLocalSampleSizeSame  ))
  cat(format("Point", width= 10, justify="right"))
  cat(format("Std." , width= 10, justify="right"))
  if (x$opt$q > x$opt$p) {
    cat(format("Robust B.C."
               , width=25, justify="centre"))
    cat("\n")
  } else {
    cat(format("Conventional"
               , width=25, justify="centre"))
    cat("\n")
  }

  cat(format("Index     Grid"            , width=14, justify="right"))
  cat(format("B.W."              , width=10, justify="right"))
  cat(format("Eff.n"           , width=8 , justify="right"))
  if (flagNotAllLocalSampleSizeSame) {
    cat(format("Uniq.n"           , width=8 , justify="right"))
  }
  cat(format("Est."            , width=10, justify="right"))
  cat(format("Error"           , width=10, justify="right"))
  if (CIuniform) {
    cat(format(paste("[ Unif. ", floor((1-alpha)*100), "%", " C.I. ]", sep="")
               , width=25, justify="centre"))
  } else {
    cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep="")
               , width=25, justify="centre"))
  }
  cat("\n")

  if (all(x$Estimate[, "nh"] == x$Estimate[, "nhu"])) {
    cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

    jj <- 1
    for (j in gridIndex) {
      cat(format(toString(j), width=4))
      cat(format(sprintf("%6.4f", x$Estimate[j, "grid"]), width=10, justify="right"))
      cat(format(sprintf("%6.4f", x$Estimate[j, "bw"])  , width=10, justify="right"))
      cat(format(sprintf("%8.0f", x$Estimate[j, "nh"])  , width=8 , justify="right"))
      cat(format(sprintf("%6.4f", x$Estimate[j, "f_p"]) , width=10, justify="right"))
      cat(format(
        paste(sprintf("%6.4f", x$Estimate[j, "se_p"]), sep=""), width=10, justify="right"))
      cat(format(
        paste(sprintf("%7.4f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
      cat(format(
        paste(sprintf("%7.4f", CI_r[j]), sep=""), width=11, justify="left"))
      cat("\n")
      if (is.numeric(sep)) if (sep > 0) if (jj %% sep == 0) {
        cat(paste(rep("-", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
      }
      jj <- jj + 1
    }

    cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 14 + 10 + 8 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

    jj <- 1
    for (j in gridIndex) {
      cat(format(toString(j), width=4))
      cat(format(sprintf("%6.4f", x$Estimate[j, "grid"]), width=10, justify="right"))
      cat(format(sprintf("%6.4f", x$Estimate[j, "bw"])  , width=10, justify="right"))
      cat(format(sprintf("%8.0f", x$Estimate[j, "nh"])  , width=8 , justify="right"))
      cat(format(sprintf("%8.0f", x$Estimate[j, "nhu"]) , width=8 , justify="right"))
      cat(format(sprintf("%6.4f", x$Estimate[j, "f_p"]) , width=10, justify="right"))
      cat(format(
        paste(sprintf("%6.4f", x$Estimate[j, "se_p"]), sep=""), width=10, justify="right"))
      cat(format(
        paste(sprintf("%7.4f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
      cat(format(
        paste(sprintf("%7.4f", CI_r[j]), sep=""), width=11, justify="left"))
      cat("\n")
      if (is.numeric(sep)) if (sep > 0) if (jj %% sep == 0) {
        cat(paste(rep("-", 14 + 10 + 8 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
      }
      jj <- jj + 1
    }

    cat(paste(rep("=", 14 + 10 + 8 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
  }
}

################################################################################
#' @title Plot Method for Local Polynomial Density Estimation and Inference
#'
#' @description The plot method for local polynomial density objects.
#'
#' @param ... Class "lpdensity" object, obtained from calling \code{\link{lpdensity}}.
#' @param alpha Numeric scalar between 0 and 1, specifies the significance level for plotting
#'   confidence intervals/bands. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param type String, one of \code{"line"} (default), \code{"points"} and \code{"both"}, specifies how
#'   the point estimates are plotted. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided, they will be applied
#'   to each data series accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param grid Numeric vector, specifies a subset of grid points
#'   to plot point estimates. This option is effective only if \code{type} is \code{"points"} or
#'   \code{"both"}; or if \code{CItype} is \code{"ebar"} or
#'   \code{"all"}.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIuniform \code{TRUE} or \code{FALSE} (default), plotting either pointwise confidence intervals (\code{FALSE}) or
#'   uniform confidence bands (\code{TRUE}).
#' @param CIsimul Positive integer, specifies the number of simulations used to construct critical values (default is \code{2000}). This
#'   option is ignored if \code{CIuniform=FALSE}.
#' @param CIshade Numeric, specifies the opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIcol Color of the confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param hist \code{TRUE} or \code{FALSE} (default), specifies whether a histogram should be added to the background.
#' @param histData Numeric vector, specifies the data used to construct the histogram plot.
#' @param histBreaks Numeric vector, specifies the breakpoints between histogram cells.
#' @param histFillCol Color of the histogram cells.
#' @param histFillShade Opaqueness of the histogram cells, should be between 0 (transparent) and
#'   1. Default is 0.2.
#' @param histLineCol Color of the histogram lines.
#' @param title,xlabel,ylabel Strings, specifies the title of the plot and labels for the x- and y-axis.
#' @param legendTitle String, specifies the legend title.
#' @param legendGroups String vector, specifies the group names used in legend.
#'
#' @return
#' \item{}{A stadnard \code{ggplot} object is returned, hence can be used for further customization.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Generate a density discontinuity at 0
#' X <- X - 0.5; X[X>0] <- X[X>0] * 2
#'
#' # Density estimation, left of 0 (scaled by the relative sample size)
#' est1 <- lpdensity(data = X[X<=0], grid = seq(-2.5, 0, 0.05), bwselect = "imse-dpi",
#'   scale = sum(X<=0)/length(X))
#' # Density estimation, right of 0 (scaled by the relative sample size)
#' est2 <- lpdensity(data = X[X>0],  grid = seq(0, 2, 0.05), bwselect = "imse-dpi",
#'   scale = sum(X>0)/length(X))
#'
#' # Plot
#' plot(est1, est2, legendTitle="My Plot", legendGroups=c("Left", "Right"))
#'
#' # Plot uniform confidence bands
#' set.seed(42) # fix the seed for simulating critical values
#' plot(est1, est2, legendTitle="My Plot", legendGroups=c("Left", "Right"), CIuniform=TRUE)
#'
#' # Adding a histogram to the background
#' plot(est1, est2, legendTitle="My Plot", legendGroups=c("Left", "Right"),
#'   hist=TRUE, histBreaks=seq(-2.4, 2, 0.2), histData=X)
#'
#' # Plot point estimates for a subset of evaluation points
#' plot(est1, est2, legendTitle="My Plot", legendGroups=c("Left", "Right"),
#'   type="both", CItype="all", grid=seq(-2, 2, 0.5))
#' @export
plot.lpdensity <- function(..., alpha=NULL,
                           type=NULL, lty=NULL, lwd=NULL, lcol=NULL, pty=NULL, pwd=NULL, pcol=NULL, grid=NULL,
                           CItype=NULL, CIuniform=FALSE, CIsimul=2000, CIshade=NULL, CIcol=NULL,
                           hist=FALSE, histData=NULL, histBreaks=NULL, histFillCol=3, histFillShade=0.2, histLineCol="white",
                           title=NULL, xlabel=NULL, ylabel=NULL, legendTitle=NULL, legendGroups=NULL) {
  ########################################
  # check how many series are passed in
  ########################################

  x <- list(...)
  nfig <- length(x)
  if (nfig == 0) stop("Nothing to plot.\n")

  flagToPlot <- rep(TRUE, nfig)
  for (i in 1:length(x)) {
    # check if is a lpbwdensity object
    if (class(x[[i]])[1] == "lpbwdensity") {
      flagToPlot[i] <- FALSE
      warning(paste("Input ", i, " is an \"lpbwdensity\" object, which is not supported by the plot method.\n", sep=""))
      next
    }

    # check if there is only one grid point
    if (nrow(x[[i]]$Estimate) < 2) {
      flagToPlot[i] <- FALSE
      warning(paste("At least two grid points are needed to plot input ", i, ".\n", sep=""))
      next
    }
  }
  x <- x[flagToPlot]
  nfig <- length(x)
  if (nfig == 0) stop("Nothing to plot.\n")

  isLpbwdensity <- FALSE
  for (i in 1:nfig) {
    if (class(x[[i]])[1] == "lpbwdensity") isLpbwdensity <- TRUE
  }

  if (isLpbwdensity) stop("The plot method does not support \"lpbwdensity\" objects.\n")

  ########################################
  # error handling
  ########################################
  # alpha
  if (length(alpha) == 0) {
    alpha <- rep(0.05, nfig)
  } else if (!all(alpha>0 & alpha<1)) {
    stop("Significance level incorrectly specified.\n")
  } else {
    alpha <- rep(alpha, length.out=nfig)
  }

  # plot type
  if (length(type) == 0) {
    type <- rep("line", nfig)
  } else {
    if (!all(type%in%c("line", "points", "both"))) {
      stop("Plotting type incorrectly specified.\n")
    }
    type <- rep(type, length.out=nfig)
  }

  # CI type
  if (length(CItype) == 0) {
    CItype <- rep("region", nfig)
  } else {
    if (!all(CItype%in%c("region", "line", "ebar", "all", "none"))) {
      stop("Confidence interval type incorrectly specified.\n")
    }
    CItype <- rep(CItype, length.out=nfig)
  }

  # line style, line width, line color
  if (length(lty) == 0) {
    lty <- rep(1, nfig)
  } else {
    lty <- rep(lty, length.out=nfig)
  }
  if (length(lwd) == 0) {
    lwd <- rep(0.5, nfig)
  } else {
    lwd <- rep(lwd, length.out=nfig)
  }
  if (length(lcol) == 0) {
    lcol <- 1:nfig
  } else {
    lcol <- rep(lcol, length.out=nfig)
  }

  # point style, point width, point color
  if (length(pty) == 0) {
    pty <- rep(1, nfig)
  } else {
    pty <- rep(pty, length.out=nfig)
  }
  if (length(pwd) == 0) {
    pwd <- rep(1, nfig)
  } else {
    pwd <- rep(pwd, length.out=nfig)
  }
  if (length(pcol) == 0) {
    pcol <- lcol
  } else {
    pcol <- rep(pcol, length.out=nfig)
  }

  # CI shade, CI color
  if (length(CIshade) == 0) {
    CIshade <- rep(0.2, nfig)
  } else {
    CIshade <- rep(CIshade, length.out=nfig)
  }
  if (length(CIcol) == 0) {
    CIcol <- lcol
  } else {
    CIcol <- rep(CIcol, length.out=nfig)
  }

  # legend
  # New in v0.2.1 to handle legend
  if (length(legendTitle) == 0) {
    legendTitle <- ""
  } else {
    legendTitle <- legendTitle[1]
  }
  if (length(legendGroups) > 0) {
    legendGroups <- rep(legendGroups, length.out=nfig)
    legend_default <- FALSE
  } else {
    legend_default <- TRUE
  }

  # grid
  if (!is.null(grid)) {
    if (!is.numeric(grid)) {
      stop("Option grid incorrectly specified.\n")
    }
  }

  ########################################
  # initializing plot
  ########################################
  if (hist & !is.null(histData)) {
    histData <- as.data.frame(histData)
    colnames(histData) <- c("v1")
    if(is.null(histBreaks)) { histBreaks <- seq(from=min(histData[, 1]), to=max(histData[, 1]), length.out=21) }
    histScale <- mean(histData[, 1] >= min(histBreaks) & histData[, 1] <= max(histBreaks))
    temp_plot <- ggplot() +
      geom_histogram(data=histData, aes(x=v1, y=..density..*histScale), breaks=histBreaks, fill=histFillCol, col=histLineCol, alpha=histFillShade) +
      theme_bw() #+ theme(legend.position="none")
  } else {
    temp_plot <- ggplot() + theme_bw() #+ theme(legend.position="none")
  }


  CI_l <- CI_r <- f_p <- Sname <- v1 <- ..density.. <- NULL

  ########################################
  # looping over input models
  ########################################
  ### all colors

  col_all <- lty_all <- pty_all <- v_all <- c()
  estRangeL <- estRangeR <- c() # estimation range
  for (i in 1:nfig) {
    # get derivative order
    v_all <- c(v_all, x[[i]]$opt$v)
    # get ploting indices
    if (is.null(grid)) {
      plotIndex <- 1:nrow(x[[i]]$Estimate)
    } else {
      gridTemp <- grid[grid >= min(x[[i]]$Estimate[, "grid"]) & grid <= max(x[[i]]$Estimate[, "grid"])]
      if (length(gridTemp) == 0) {
        plotIndex <- NULL
      } else {
        plotIndex <- rep(NA, length(gridTemp))
        for (ii in 1:length(gridTemp)) {
          plotIndex[ii] <- which.min(abs(gridTemp[ii]-x[[i]]$Estimate[, "grid"]))
        }
        plotIndex <- unique(plotIndex)
      }
    }

    estRangeL <- min(estRangeL, min(x[[i]]$Estimate[, "grid"]))
    estRangeR <- max(estRangeR, max(x[[i]]$Estimate[, "grid"]))

    data_x <- data.frame(x[[i]]$Estimate[, c("grid", "f_p", "f_q", "se_p", "se_q"), drop=FALSE])

    if (x[[i]]$opt$q == x[[i]]$opt$p) {
      data_x$f_q <- data_x$f_p; data_x$se_q <- data_x$se_p
    }

    # critical value
    if (CIuniform) {
      if (length(CIsimul) == 0) { CIsimul <- 2000 }
      if (!is.numeric(CIsimul) | is.na(CIsimul)) {
        warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
        z_val <- qnorm(1 - alpha[i]/2)
      } else if (ceiling(CIsimul)<2) {
        warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
        z_val <- qnorm(1 - alpha[i]/2)
      } else {
        CIsimul <- ceiling(CIsimul)
        if (x[[i]]$opt$q == x[[i]]$opt$p) {
          x[[i]]$CovMat_q <- x[[i]]$CovMat_p
        }
        corrMat <- sweep(sweep(x[[i]]$CovMat_q, MARGIN=1, FUN="*", STATS=1/data_x$se_q), MARGIN=2, FUN="*", STATS=1/data_x$se_q)
        normalSimu <- try(
          mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=corrMat),
          silent=TRUE)
        if (is.character(normalSimu)) {
          print(normalSimu)
          warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
          z_val <- qnorm(1 - alpha[i]/2)
        } else {
          z_val <- quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha[i])
        }
      }
    } else {
      z_val <- qnorm(1 - alpha[i]/2)
    }


    data_x$CI_l <- data_x$f_q - z_val * data_x$se_q
    data_x$CI_r <- data_x$f_q + z_val * data_x$se_q

    # New in v0.2.1 to handle legend
    if (legend_default) {
      data_x$Sname <- paste("Series", i, sep=" ")
      legendGroups <- c(legendGroups, data_x$Sname)
    } else {
      data_x$Sname <- legendGroups[i]
    }

    ########################################
    # add CI regions to the plot
    if (CItype[i]%in%c("region", "all"))
      temp_plot <- temp_plot + geom_ribbon(data=data_x, aes(x=grid, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], fill=CIcol[i])

    ########################################
    # add CI lines to the plot
    if (CItype[i]%in%c("line", "all"))
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=grid, y=CI_l), linetype=2, alpha=1, col=CIcol[i]) +
      geom_line(data=data_x, aes(x=grid, y=CI_r), linetype=2, alpha=1, col=CIcol[i])
    #if (CItype[i]%in%c("line", "all"))
    #  temp_plot <- temp_plot + geom_line(data=data_x, aes(x=grid, y=CI_l), linetype=2, alpha=CIshade[i], col=CIcol[i]) +
    #  geom_line(data=data_x, aes(x=grid, y=CI_r), linetype=2, alpha=CIshade[i], col=CIcol[i])

    ########################################
    # add error bars to the plot
    if (CItype[i]%in%c("ebar", "all") & !is.null(plotIndex))
      temp_plot <- temp_plot + geom_errorbar(data=data_x[plotIndex, ], aes(x=grid, ymin=CI_l, ymax=CI_r), alpha=1, col=CIcol[i], linetype=1)
    #if (CItype[i]%in%c("ebar", "all") & !is.null(plotIndex))
    #  temp_plot <- temp_plot + geom_errorbar(data=data_x[plotIndex, ], aes(x=grid, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], col=CIcol[i], linetype=1)

    ########################################
    # add lines to the plot
    if (type[i]%in%c("line", "both")) {
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=grid, y=f_p, colour=Sname, linetype=Sname), size=lwd[i])
    }

    ########################################
    # add points to the plot
    if (type[i]%in%c("points", "both") & !is.null(plotIndex)) {
      temp_plot <- temp_plot + geom_point(data=data_x[plotIndex, ], aes(x=grid, y=f_p, colour=Sname, shape=Sname), size=pwd[i])
    }

    if (type[i] == "line") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, NA)
    } else if (type[i] == "both") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, pty[i])
    } else {
      col_all <- c(col_all, pcol[i])
      lty_all <- c(lty_all, NA)
      pty_all <- c(pty_all, pty[i])
    }
  }

  ########################################
  # change color, line type and point shape back, and customize legend
  ########################################
  # New in v0.2.1 to handle legend
  index <- sort.int(legendGroups, index.return=TRUE)$ix
  temp_plot <- temp_plot + scale_color_manual(values = col_all[index]) +
    scale_linetype_manual(values = lty_all[index]) +
    scale_shape_manual(values = pty_all[index]) +
    guides(colour=guide_legend(title=legendTitle)) +
    guides(linetype=guide_legend(title=legendTitle)) +
    guides(shape=guide_legend(title=legendTitle))

  ########################################
  # add title, x and y labs
  ########################################
  if (is.null(ylabel)) {
    if (all(v_all == v_all[1])) {
      if (v_all[1] == 0) {
        ylabel <- "Distribution function"
      } else if (v_all[1] == 1) {
        ylabel <- "Density"
      } else {
        ylabel <- paste("Density derivative (v=", v_all[1], ")", sep="")
      }
    } else {
      ylabel <- ""
    }
  }

  if (is.null(xlabel)) {
    xlabel <- ""
  }

  if (is.null(title)) {
    title <- ""
  }
  temp_plot <- temp_plot + labs(x=xlabel, y=ylabel) + ggtitle(title)

  # check plotting range vs estimation range
  if (!is.null(grid)) {
    if (min(grid) < estRangeL | max(grid) > estRangeR) {
      warning("The plotting range exceeds the original estimation range. Option plot(..., grid=) should be within the estimation range specified by lpdensity(..., grid=).\n")
    }
  }

  ########################################
  # return the plot
  ########################################
  return (temp_plot)
}

################################################################################
#' @title Plot Method for Local Polynomial Density Estimation and Inference
#'
#' @description This has been replaced by \code{\link{plot.lpdensity}}.
#'
#' @param ... Class "lpdensity" object, obtained from calling \code{\link{lpdensity}}.
#' @param alpha Numeric scalar between 0 and 1, specifies the significance level for plotting
#'   confidence intervals/bands. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param type String, one of \code{"line"} (default), \code{"points"} and \code{"both"}, specifies how
#'   the point estimates are plotted. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided, they will be applied
#'   to each data series accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param grid Numeric vector, specifies a subset of grid points
#'   to plot point estimates. This option is effective only if \code{type} is \code{"points"} or
#'   \code{"both"}; or if \code{CItype} is \code{"ebar"} or
#'   \code{"all"}.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIuniform \code{TRUE} or \code{FALSE} (default), plotting either pointwise confidence intervals (\code{FALSE}) or
#'   uniform confidence bands (\code{TRUE}).
#' @param CIsimul Positive integer, specifies the number of simulations used to construct critical values (default is \code{2000}). This
#'   option is ignored if \code{CIuniform=FALSE}.
#' @param CIshade Numeric, specifies the opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIcol Color of the confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param hist \code{TRUE} or \code{FALSE} (default), specifies whether a histogram should be added to the background.
#' @param histData Numeric vector, specifies the data used to construct the histogram plot.
#' @param histBreaks Numeric vector, specifies the breakpoints between histogram cells.
#' @param histFillCol Color of the histogram cells.
#' @param histFillShade Opaqueness of the histogram cells, should be between 0 (transparent) and
#'   1. Default is 0.2.
#' @param histLineCol Color of the histogram lines.
#' @param title,xlabel,ylabel Strings, specifies the title of the plot and labels for the x- and y-axis.
#' @param legendTitle String, specifies the legend title.
#' @param legendGroups String vector, specifies the group names used in legend.
#'
#' @return
#' \item{}{A stadnard \code{ggplot} object is returned, hence can be used for further customization.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @export
lpdensity.plot <- plot.lpdensity

################################################################################
#' Coef Method for Local Polynomial Density Estimation and Inference
#'
#' @description The coef method for local polynomial density objects.
#'
#' @param object Class "lpdensity" object, obtained by calling \code{\link{lpdensity}}.
#' @param ... Additional options.
#'
#' @return
#' \item{}{A matrix containing grid points and density estimates using p- and q-th order local polynomials.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Estimate density and report results
#' coef(lpdensity(data = X, bwselect = "imse-dpi"))
#'
#' @export
coef.lpdensity <- function(object, ...) {
  object$Estimate[, c("grid", "f_p", "f_q")]
}

################################################################################
#' Vcov Method for Local Polynomial Density Estimation and Inference
#'
#' @description The vcov method for local polynomial density objects.
#'
#' @param object Class "lpdensity" object, obtained by calling \code{\link{lpdensity}}.
#' @param ... Additional options.
#'
#' @return
#' \item{stdErr}{A matrix containing grid points and standard errors using p- and q-th order local polynomials.}
#' \item{CovMat_p}{The variance-covariance matrix corresponding to \code{f_p}.}
#' \item{CovMat_q}{The variance-covariance matrix corresponding to \code{f_q}.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Estimate density and report results
#' vcov(lpdensity(data = X, bwselect = "imse-dpi"))
#'
#' @export
vcov.lpdensity <- function(object, ...) {

  if (class(object)[1] == "lpbwdensity") stop("The vcov method does not support \"lpbwdensity\" objects.\n")

  return(list(stdErr=object$Estimate[, c("grid", "se_p", "se_q")], CovMat_p=object$CovMat_p, CovMat_q=object$CovMat_q))
}

################################################################################
#' Confint Method for Local Polynomial Density Estimation and Inference
#'
#' @description The confint method for local polynomial density objects.
#'
#' @param object Class "lpdensity" object, obtained by calling \code{\link{lpdensity}}.
#' @param parm Integer, indicating which parameters are to be given confidence intervals.
#' @param level Numeric scalar between 0 and 1, the significance level for computing confidence intervals
#' @param ... Additional options, including (i) \code{grid} specifies a subset of grid points
#'   to display the bandwidth; (ii) \code{gridIndex} specifies the indices of grid points
#'   to display the bandwidth (this is the same as \code{parm}); (iii) \code{alpha} specifies the significance level
#'   (this is \code{1-level}); (iv)
#'   \code{CIuniform} specifies whether displaying pointwise confidence intervals (\code{FALSE}, default) or
#'   the uniform confidence band (\code{TRUE}); (v) \code{CIsimul} specifies the number of simulations used
#'   to construct critical values (default is 2000).
#'
#' @return
#' \item{}{A matrix containing grid points and confidence interval end points using p- and q-th order local polynomials.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}},
#'   \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}},
#'   \code{\link{vcov.lpdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Estimate density and report 95% confidence intervals
#' est1 <- lpdensity(data = X, bwselect = "imse-dpi")
#' confint(est1)
#'
#' # Report results for a subset of grid points
#' confint(est1, parm=est1$Estimate[4:10, "grid"])
#' confint(est1, grid=est1$Estimate[4:10, "grid"])
#' confint(est1, gridIndex=4:10)
#'
#' # Report the 99% uniform confidence band
#' # Fix the seed for simulating critical values
#' set.seed(42); confint(est1, level=0.99, CIuniform=TRUE)
#' set.seed(42); confint(est1, alpha=0.01, CIuniform=TRUE)
#'
#' @export
confint.lpdensity <- function(object, parm = NULL, level = NULL, ...) {

  x <- object
  if (class(x)[1] == "lpbwdensity") stop("The confint method does not support \"lpbwdensity\" objects.\n")
  args <- list(...)

  if (!is.null(parm)) { args[['grid']] <- parm }
  if (!is.null(level)) { args[['alpha']] <- 1 - level }

  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['CIuniform']]))   { CIuniform <- FALSE } else { CIuniform <- args[['CIuniform']] }
  if (is.null(args[['CIsimul']]))   { CIsimul <- 2000 } else { sep <- args[['CIsimul']] }

  if (is.null(args[['grid']]) & is.null(args[['gridIndex']])) {
    gridIndex <- 1:nrow(x$Estimate)
  } else if (is.null(args[['grid']]) & !is.null(args[['gridIndex']])) {
    gridIndex <- args[['gridIndex']]
    if (is.null(gridIndex)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!all(gridIndex %in% 1:nrow(x$Estimate))) {
      stop("Option gridIndex incorrectly specified.\n")
    }
  } else {
    grid <- args[['grid']]
    if (is.null(grid)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!is.numeric(grid)) {
      stop("Option param/grid incorrectly specified.\n")
    } else {
      gridIndex <- rep(NA, length(grid))
      for (j in 1:length(grid)) {
        gridIndex[j] <- which.min(abs(x$Estimate[, "grid"]-grid[j]))
        #gridIndex <- unique(gridIndex)
      }
    }
  }


  Estimate <- matrix(NA, nrow=length(gridIndex), ncol=5)
  Estimate[, 1] <- x$Estimate[gridIndex, "grid"]

  if (CIuniform) {
    if (length(CIsimul) == 0) { CIsimul <- 2000 }
    if (!is.numeric(CIsimul) | is.na(CIsimul)) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform <- FALSE
      z_val <- qnorm(1 - alpha/2)
    } else if (ceiling(CIsimul)<2) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform <- FALSE
      z_val <- qnorm(1 - alpha/2)
    } else {
      CIsimul <- ceiling(CIsimul)
      if (x$opt$q == x$opt$p) {
        x$CovMat_q <- x$CovMat_p
        x$Estimate[, "se_q"] <- x$Estimate[, "se_p"]
      }
      corrMat <- sweep(sweep(x$CovMat_q, MARGIN=1, FUN="*", STATS=1/x$Estimate[, "se_q"]), MARGIN=2, FUN="*", STATS=1/x$Estimate[, "se_q"])
      normalSimu <- try(
        mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=corrMat),
        silent=TRUE)
      if (is.character(normalSimu)) {
        print(normalSimu)
        warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
        CIuniform <- FALSE
        z_val <- qnorm(1 - alpha/2)
      } else {
        z_val <- quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha)
      }
    }
  } else {
    z_val <- qnorm(1 - alpha/2)
  }

  Estimate[, 2] <- x$Estimate[gridIndex, "f_p"] - z_val * x$Estimate[gridIndex, "se_p"]
  Estimate[, 3] <- x$Estimate[gridIndex, "f_p"] + z_val * x$Estimate[gridIndex, "se_p"]

  Estimate[, 4] <- x$Estimate[gridIndex, "f_q"] - z_val * x$Estimate[gridIndex, "se_q"]
  Estimate[, 5] <- x$Estimate[gridIndex, "f_q"] + z_val * x$Estimate[gridIndex, "se_q"]

  colnames(Estimate) <- c("grid", "CI_l_p", "CI_r_p", "CI_l_q", "CI_r_q")

  #if (x$opt$p == x$opt$q) Estimate <- Estimate[, 1:3]

  return(Estimate)
}
