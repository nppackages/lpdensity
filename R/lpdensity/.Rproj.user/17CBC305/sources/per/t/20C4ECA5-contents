################################################################################
#' @title Local Polynomial Density Plot with Robust Confidence Intervals
#'
#' @description \code{lpdensity.plot} plots estimated density/derivatives. This command
#'   can also be used to plot smoothed distribution function. See
#'   Cattaneo, Jansson and Ma (2019b) for more implementation details and illustrations.
#'
#'   Companion command: \code{\link{lpdensity}} for local polynomial based density
#'   and derivatives estimation.
#'
#'   For more details, and related Stata and R packages useful for empirical analysis,
#'   visit \url{https://sites.google.com/site/nppackages/}.
#'
#' @param ... Objects returned by \code{\link{lpdensity}}.
#' @param alpha Numeric scalar between 0 and 1, the significance level for plotting
#'   confidence regions. If more than one is provided, they will be applied to data series
#'   accordingly.
#' @param type String, one of \code{"line"} (default), \code{"points"} or \code{"both"}, how
#'   the point estimates are plotted. If more than one is provided, they will be applied to data series
#'   accordingly.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to data series
#'   accordingly.
#' @param title,xlabel,ylabel Strings, title of the plot and labels for x- and y-axis.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to data series accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided, they will be applied
#'   to data series accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to data series
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to data series
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to data series
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to data series
#'   accordingly.
#' @param CIshade Numeric, opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to data series
#'   accordingly.
#' @param CIcol color for confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to data series
#'   accordingly.
#' @param legendTitle String, title of legend.
#' @param legendGroups String Vector, group names used in legend.
#'
#' @return
#' \item{}{A stadnard \code{ggplot} object is returned, hence can be used for further customization.}
#'
#' @references
#' M.D. Cattaneo, M. Jansson and X. Ma. (2019a). \href{https://arxiv.org/abs/1811.11512}{Simple Local Polynomial Density Estimators}. \emph{Journal of the American Statistical Association}, forthcoming.
#'
#' M.D. Cattaneo, M. Jansson and X. Ma. (2019b). \href{https://arxiv.org/abs/1906.06529}{\code{lpdensity}: Local Polynomial Density Estimation and Inference}. Working paper.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpdensity}}
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
#' lpdensity.plot(est1, est2, legendTitle="My Plot", legendGroups=c("Left", "Right"))
#'
#' @export
lpdensity.plot <- function(..., alpha=NULL, type=NULL, CItype=NULL,
                           title="", xlabel="", ylabel="",
                           lty=NULL, lwd=NULL, lcol=NULL, pty=NULL, pwd=NULL, pcol=NULL,
                           CIshade=NULL, CIcol=NULL, legendTitle=NULL, legendGroups=NULL) {

  ########################################
  # check how many series are passed in
  ########################################

  x <- list(...)
  nfig <- length(x)
  if (nfig == 0) stop("Nothing to plot.\n")

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

  ########################################
  # initializing plot
  ########################################
  temp_plot <- ggplot() + theme_bw() #+ theme(legend.position="none")

  CI_l <- CI_r <- f_p <- grid <- Sname <- NULL

  ########################################
  # looping over input models
  ########################################
  ### all colors
  col_all <- lty_all <- pty_all <- c()
  for (i in 1:nfig) {
    data_x <- data.frame(x[[i]]$Estimate[, c("grid", "f_p", "f_q", "se_p", "se_q")])
    z_val <- qnorm(1 - alpha[i]/2)
    if (x[[i]]$opt$q == x[[i]]$opt$p) {
      data_x$f_q <- data_x$f_p; data_x$se_q <- data_x$se_p
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
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=grid, y=CI_l), linetype=2, alpha=CIshade[i], col=CIcol[i]) +
      geom_line(data=data_x, aes(x=grid, y=CI_r), linetype=2, alpha=CIshade[i], col=CIcol[i])

    ########################################
    # add error bars to the plot
    if (CItype[i]%in%c("ebar", "all"))
      temp_plot <- temp_plot + geom_errorbar(data=data_x, aes(x=grid, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], col=CIcol[i], linetype=1)

    ########################################
    # add lines to the plot
    if (type[i]%in%c("line", "both")) {
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=grid, y=f_p, colour=Sname, linetype=Sname), size=lwd[i])
    }

    ########################################
    # add points to the plot
    if (type[i]%in%c("points", "both")) {
      temp_plot <- temp_plot + geom_point(data=data_x, aes(x=grid, y=f_p, colour=Sname, shape=Sname), size=pwd[i])
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
  temp_plot <- temp_plot + labs(x=xlabel, y=ylabel) + ggtitle(title)

  ########################################
  # return the plot
  ########################################
  return (temp_plot)
}
