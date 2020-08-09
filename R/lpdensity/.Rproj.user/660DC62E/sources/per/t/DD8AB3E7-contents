################################################################################
#' @title Data-driven Bandwidth Selection for Local Polynomial Density Estimators
#'
#' @description  \code{\link{lpbwdensity}} implements the bandwidth selection methods for local
#'   polynomial based density (and derivatives) estimation proposed and studied
#'   in Cattaneo, Jansson and Ma (2020a) and Cattaneo, Jansson and Ma (2020b).
#'   See Cattaneo, Jansson and Ma (2020c) for more implementation details and illustrations.
#'
#'   Companion command: \code{\link{lpdensity}} for estimation and robust bias-corrected inference.
#'
#'   Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
#'   available at \url{https://sites.google.com/site/nppackages/}.
#'
#' @param data Numeric vector or one dimensional matrix/data frame, the raw data.
#' @param grid Numeric, specifies the grid of evaluation points. When set to default, grid points
#'   will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.
#' @param p Nonnegative integer, specifies the order of the local polynomial used to construct point
#'   estimates. (Default is \code{2}.)
#' @param v Nonnegative integer, specifies the derivative of the distribution function to be estimated. \code{0} for
#'   the distribution function, \code{1} (default) for the density funtion, etc.
#' @param kernel String, specifies the kernel function, should be one of \code{"triangular"}, \code{"uniform"} or
#'   \code{"epanechnikov"}.
#' @param bwselect String, specifies the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Can be (1) \code{"mse-dpi"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point); or (2) \code{"imse-dpi"} (integrated MSE-optimal bandwidth,
#'   common for all grid points); (3) \code{"mse-rot"} (rule-of-thumb bandwidth with Gaussian
#'   reference model); and (4) \code{"imse-rot"} (integrated rule-of-thumb bandwidth with Gaussian
#'   reference model).
#' @param massPoints \code{TRUE} (default) or \code{FALSE}, specifies whether point estimates and standard errors
#'   should be adjusted if there are mass points in the data.
#' @param stdVar \code{TRUE} (default) or \code{FALSE}, specifies whether the data should be standardized for
#'   bandwidth selection.
#' @param regularize \code{TRUE} (default) or \code{FALSE}, specifies whether the bandwidth should be
#'   regularized. When set to \code{TRUE}, the bandwidth is chosen such that the local region includes
#'   at least \code{nLocalMin} observations and at least \code{nUniqueMin} unique observations.
#' @param nLocalMin Nonnegative integer, specifies the minimum number of observations in each local neighborhood. This option
#'   will be ignored if \code{regularize=FALSE}. Default is \code{20+p+1}.
#' @param nUniqueMin Nonnegative integer, specifies the minimum number of unique observations in each local neighborhood. This option
#'   will be ignored if \code{regularize=FALSE}. Default is \code{20+p+1}.
#' @param Cweights Numeric vector, specifies the weights used
#'   for counterfactual distribution construction. Should have the same length as the data.
#'   This option will be ignored if \code{bwselect} is \code{"mse-rot"} or \code{"imse-rot"}.
#' @param Pweights Numeric vector, specifies the weights used
#'   in sampling. Should have the same length as the data.
#'   This option will be ignored if \code{bwselect} is \code{"mse-rot"} or \code{"imse-rot"}.
#'
#' @return
#' \item{BW}{A matrix containing (1) \code{grid} (grid point), (2) \code{bw} (bandwidth),
#'   (3) \code{nh} (number of observations in each local neighborhood), and
#'   (4) \code{nhu} (number of unique observations in each local neighborhood).}
#' \item{opt}{A list containing options passed to the function.}
#'
#' @references
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2020a. \href{https://sites.google.com/site/nppackages/lpdensity/Cattaneo-Jansson-Ma_2020_JASA.pdf}{Simple Local Polynomial Density Estimators}. \emph{Journal of the American Statistical Association}, forthcoming.
#'
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2020b. \href{https://sites.google.com/site/nppackages/lpdensity/Cattaneo-Jansson-Ma_2020_JoE.pdf}{Local Regression Distribution Estimators}. Working paper.
#'
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2020c. \href{https://sites.google.com/site/nppackages/lpdensity/Cattaneo-Jansson-Ma_2020_JSS.pdf}{lpdensity: Local Polynomial Density Estimation and Inference}. Working paper.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso Supported methods: \code{\link{coef.lpbwdensity}}, \code{\link{print.lpbwdensity}}, \code{\link{summary.lpbwdensity}}.
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
lpbwdensity <- function(data, grid=NULL,
                        p=NULL, v=NULL, kernel=c("triangular", "uniform", "epanechnikov"),
                        bwselect=c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot"),
                        massPoints=TRUE, stdVar=TRUE,
                        regularize=TRUE, nLocalMin=NULL, nUniqueMin=NULL,
                        Cweights=NULL, Pweights=NULL) {
  ################################################################################
  # Input Error Handling
  ################################################################################
  # data
  data <- as.vector(data)
  if (any(is.na(data))) {
    warning(paste(sum(is.na(data)), " missing ", switch((sum(is.na(data))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    data <- data[!is.na(data)]
  }
  n <- length(data)
  if (!is.numeric(data) | length(data)==0) {
    stop("Data has to be numeric, and cannot be empty.\n")
  }

  # grid
  if (length(grid) == 0) {
    flag_no_grid <- TRUE
    grid <- quantile(data, seq(from=0.05, to=0.95, by=0.05))
    ng <- length(grid)
  } else {
    flag_no_grid <- FALSE
    grid <- as.vector(grid)
    ng <- length(grid)
    if(!is.numeric(grid)) {
      stop("Grid points have to be numeric.\n")
    }
  }

  # bwselect
  if (length(bwselect) == 0) {
    bwselect <- "mse-dpi"
  } else {
    bwselect <- tolower(bwselect[1])
  }
  ### BACKWARD COMPATIBILITY FOR VERSION 0.1
  if (bwselect == "mse" ) { bwselect <- "mse-dpi" }
  if (bwselect == "imse") { bwselect <- "imse-dpi" }
  if (bwselect == "rot" ) { bwselect <- "mse-rot" }
  if (bwselect == "irot") { bwselect <- "imse-rot" }
  ### END BACKWARD COMPATIBILITY FOR VERSION 0.1
  if (!bwselect%in%c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot")) stop("Incorrect bandwidth selection method specified.\n")

  # p
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 2
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }

  # v
  if (length(v) == 0) {
    flag_no_v <- TRUE
    v <- min(p, 1)
  } else if ((length(v) > 1) | !(v[1]%in%c(0:20)) | (v[1]>p)) {
    stop("Derivative order v incorrectly specified.\n")
  } else {
    flag_no_v <- FALSE
  }

  # kernel
  if (length(kernel) == 0) {
    flag_no_kernel <- TRUE
    kernel <- "triangular"
  } else {
    kernel <- tolower(kernel)
    kernel <- kernel[1]
    if (!kernel%in%c("triangular", "uniform", "epanechnikov")) {
      stop("Kernel function incorrectly specified.\n")
    } else {
      flag_no_kernel <- FALSE
    }
  }

  # Cweights
  if (length(Cweights) == 0) {
    flag_no_Cweights <- TRUE
    Cweights <- rep(1, n)
  } else if (!is.numeric(Cweights)) {
    stop("Counterfactual weights incorrectly specified.\n")
  } else if (length(Cweights) != n) {
    stop("Counterfactual weights should have the same length as sample.\n")
  } else {
    flag_no_Cweights <- FALSE
  }

  # Pweights
  if (length(Pweights) == 0) {
    flag_no_Pweights <- TRUE
    Pweights <- rep(1, n)
  } else if (!is.numeric(Pweights)) {
    stop("Probability weights incorrectly specified.\n")
  } else if (length(Pweights) != n) {
    stop("Probability weights should have the same length as sample.\n")
  } else if (any(Pweights < 0)) {
    stop("Probability weights should be nonnegative.\n")
  } else{
    flag_no_Pweights <- FALSE
  }

  # massPoints
  if (length(massPoints) == 0) {
    massPoints <- TRUE
  } else if (length(massPoints) > 1 | !massPoints[1]%in%c(TRUE, FALSE)) {
    stop("Option massPoints incorrectly specified.\n")
  }

  # stdVar
  if (length(stdVar) == 0) {
    stdVar <- TRUE
  } else if (length(stdVar) > 1 | !stdVar[1]%in%c(TRUE, FALSE)) {
    stop("Option stdVar incorrectly specified.\n")
  }

  # regularize
  if (length(regularize) == 0) {
    regularize <- TRUE
  } else if (length(regularize) > 1 | !regularize[1]%in%c(TRUE, FALSE)) {
    stop("Regularization parameter incorrectly specified.\n")
  }

  # nLocalMin
  if (length(nLocalMin) == 0) { nLocalMin <- 20 + p + 1 }
  if (!is.numeric(nLocalMin) | is.na(nLocalMin)) {
    stop("Option nLocalMin incorrectly specified.\n")
  } else if (ceiling(nLocalMin) < 0) {
    stop("Option nLocalMin incorrectly specified.\n")
  } else {
    nLocalMin <- ceiling(nLocalMin)
  }

  # nUniqueMin
  if (length(nUniqueMin) == 0) { nUniqueMin <- 20 + p + 1 }
  if (!is.numeric(nUniqueMin) | is.na(nUniqueMin)) {
    stop("Option nUniqueMin incorrectly specified.\n")
  } else if (ceiling(nUniqueMin) < 0) {
    stop("Option nUniqueMin incorrectly specified.\n")
  } else {
    nUniqueMin <- ceiling(nUniqueMin)
  }

  ################################################################################
  # Sample trimming
  ################################################################################
  trim_index <- (Pweights == 0)
  if (all(trim_index)) {
    stop("All weights are zero.\n")
  } else {
    data <- data[!trim_index]
    Cweights <- Cweights[!trim_index]
    Pweights <- Pweights[!trim_index]
  }
  if (abs(sum(Cweights * Pweights)) <= .Machine$double.eps * 10) {
    stop("Composited weights (Cweights * Pweights) are numerically zero.\n")
  }

  ################################################################################
  # Bandwidth Estimation
  ################################################################################

  if (bwselect == "mse-dpi") {
    bw <- bw_MSE( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  } else if (bwselect == "imse-dpi") {
    bw <- bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  } else if (bwselect == "mse-rot") {
    bw <- bw_ROT( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  } else {
    bw <- bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  }

  BW <- matrix(NA, ncol=4, nrow=ng)
  BW[, 1] <- grid
  BW[, 2] <- bw
  rownames(BW) <- 1:ng
  colnames(BW) <- c("grid", "bw", "nh", "nhu")

  dataUnique  <- lpdensityUnique(sort(data, decreasing=FALSE))$unique

  for (i in 1:ng) {
    BW[i, 3] <- sum(abs(data       - BW[i, 1]) <= BW[i, 2])
    BW[i, 4] <- sum(abs(dataUnique - BW[i, 1]) <= BW[i, 2])
  }

  Result <- list(BW=BW,
                 opt=list(p=p, v=v, kernel=kernel, n=n, ng=ng,
                          bwselect=bwselect,
                          massPoints=massPoints, stdVar=stdVar,
                          regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                          data_min=min(data), data_max=max(data),
                          grid_min=min(grid), grid_max=max(grid)))

  class(Result) <- c("lpbwdensity", "lpdensity")

  return (Result)
}
