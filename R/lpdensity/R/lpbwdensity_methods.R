################################################################################
#' Print Method for Local Polynomial Density Bandwidth Selection
#'
#' @description The print method for local polynomial density bandwidth selection objects.
#'
#' @param x Class "lpbwdensity" object, obtained by calling \code{\link{lpbwdensity}}.
#' @param ... Other arguments.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwdensity}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwdensity}}, \code{\link{print.lpbwdensity}}, \code{\link{summary.lpbwdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Construct bandwidth
#' print(lpbwdensity(X))
#'
#' @export
print.lpbwdensity <- function(x, ...) {

  cat("Call: lpbwdensity\n\n")

  cat(paste("Sample size                                   ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p=)    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated         (v=)    ", x$opt$v,        "\n", sep=""))
  cat(paste("Kernel function                               ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                              ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  cat("Use summary(...) to show bandwidths.\n")
}

################################################################################
#' Summary Method for Local Polynomial Density Bandwidth Selection
#'
#' @description The summary method for local polynomial density bandwidth selection objects.
#'
#' @param object Class "lpbwdensity" object, obtained by calling \code{\link{lpbwdensity}}.
#' @param ... Additional options, including (i) \code{grid} specifies a subset of grid points
#'   to display the bandwidth; (ii) \code{gridIndex} specifies the indices of grid points
#'   to display the bandwidth.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwdensity}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwdensity}}, \code{\link{print.lpbwdensity}}, \code{\link{summary.lpbwdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Construct bandwidth
#' bw1 <- lpbwdensity(X)
#' summary(bw1)
#'
#' # Display bandwidths for a subset of grid points
#' summary(bw1, grid=bw1$BW[4:10, "grid"])
#' summary(bw1, gridIndex=4:10)
#'
#' @export
summary.lpbwdensity <- function(object, ...) {

  x <- object
  args <- list(...)

  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }

  if (is.null(args[['grid']]) & is.null(args[['gridIndex']])) {
    gridIndex <- 1:nrow(x$BW)
  } else if (is.null(args[['grid']]) & !is.null(args[['gridIndex']])) {
    gridIndex <- args[['gridIndex']]
    if (is.null(gridIndex)) {
      gridIndex <- 1:nrow(x$BW)
    } else if (!all(gridIndex %in% 1:nrow(x$BW))) {
      stop(paste("Option gridIndex incorrectly specified. Should be integers between 1 and ", nrow(x$BW), ".\n", sep=""))
    }
  } else {
    grid <- args[['grid']]
    if (is.null(grid)) {
      gridIndex <- 1:nrow(x$BW)
    } else if (!is.numeric(grid)) {
      stop("Option grid incorrectly specified.\n")
    } else {
      gridIndex <- rep(NA, length(grid))
      if (min(grid) < min(x$BW[, "grid"]) | max(grid) > max(x$BW[, "grid"])) {
        warning("The reporting range exceeds the original estimation range. Option summary(..., grid=) should be within the estimation range specified by lpbwdensity(..., grid=).\n")
      }
      for (j in 1:length(grid)) {
        gridIndex[j] <- which.min(abs(x$BW[, "grid"]-grid[j]))
        #gridIndex <- unique(gridIndex)
      }
    }
  }

  cat("Call: lpbwdensity\n\n")

  cat(paste("Sample size                           (n=)    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p=)    ", x$opt$p,        "\n", sep=""))
  if (x$opt$v == 0) {
    cat(paste("Distribution function estimated       (v=)    ", x$opt$v,        "\n", sep=""))
  } else if (x$opt$v == 1) {
    cat(paste("Density function estimated            (v=)    ", x$opt$v,        "\n", sep=""))
  } else {
    cat(paste("Order of derivative estimated         (v=)    ", x$opt$v,        "\n", sep=""))
  }
  cat(paste("Kernel function                               ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth selection method                    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  ### print output
  if (all(x$BW[, "nh"] == x$BW[, "nhu"])) {
    cat(paste(rep("=", 14 + 10 + 8), collapse="")); cat("\n")

    cat(format("Index     Grid"            , width=14, justify="right"))
    cat(format("B.W."              , width=10, justify="right"))
    cat(format("Eff.n"           , width=8 , justify="right"))
    cat("\n")

    cat(paste(rep("=", 14 + 10 + 8), collapse="")); cat("\n")

    jj <- 1
    for (j in gridIndex) {
      cat(format(toString(j), width=4))
      cat(format(sprintf("%6.4f", x$BW[j, "grid"]), width=10, justify="right"))
      cat(format(sprintf("%6.4f", x$BW[j, "bw"])  , width=10, justify="right"))
      cat(format(sprintf("%8.0f", x$BW[j, "nh"])  , width=8 , justify="right"))
      cat("\n")
      if (is.numeric(sep)) if (sep > 0) if (jj %% sep == 0) {
        cat(paste(rep("-", 14 + 10 + 8), collapse="")); cat("\n")
      }
      jj <- jj + 1
    }

    cat(paste(rep("=", 14 + 10 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 14 + 10 + 8 + 8), collapse="")); cat("\n")

    cat(format("Index     Grid"            , width=14, justify="right"))
    cat(format("B.W."              , width=10, justify="right"))
    cat(format("Eff.n"           , width=8 , justify="right"))
    cat(format("Uniq.n"           , width=8 , justify="right"))
    cat("\n")

    cat(paste(rep("=", 14 + 10 + 8 + 8), collapse="")); cat("\n")

    jj <- 1
    for (j in gridIndex) {
      cat(format(toString(j), width=4))
      cat(format(sprintf("%6.4f", x$BW[j, "grid"]), width=10, justify="right"))
      cat(format(sprintf("%6.4f", x$BW[j, "bw"])  , width=10, justify="right"))
      cat(format(sprintf("%8.0f", x$BW[j, "nh"])  , width=8 , justify="right"))
      cat(format(sprintf("%8.0f", x$BW[j, "nhu"])  , width=8 , justify="right"))
      cat("\n")
      if (is.numeric(sep)) if (sep > 0) if (jj %% sep == 0) {
        cat(paste(rep("-", 14 + 10 + 8 + 8), collapse="")); cat("\n")
      }
      jj <- jj + 1
    }

    cat(paste(rep("=", 14 + 10 + 8 + 8), collapse="")); cat("\n")
  }
}

################################################################################
#' Coef Method for Local Polynomial Density Bandwidth Selection
#'
#' @description The coef method for local polynomial density bandwidth selection objects.
#'
#' @param object Class "lpbwdensity" object, obtained by calling \code{\link{lpbwdensity}}.
#' @param ... Other arguments.
#'
#' @return
#' A matrix containing grid points and selected bandwidths.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwdensity}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwdensity}}, \code{\link{print.lpbwdensity}}, \code{\link{summary.lpbwdensity}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42); X <- rnorm(2000)
#'
#' # Construct bandwidth
#' coef(lpbwdensity(X))
#'
#' @export
coef.lpbwdensity <- function(object, ...) {
  object$BW[, c("grid", "bw")]
}

