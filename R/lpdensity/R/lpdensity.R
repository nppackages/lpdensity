################################################################################
#' @title  Local Polynomial Density Estimation and Inference
#'
#' @description \code{\link{lpdensity}} implements the local polynomial regression based density (and derivatives)
#'   estimator proposed in Cattaneo, Jansson and Ma (2020a).  Robust bias-corrected inference methods,
#'   both pointwise (confidence intervals) and uniform (confidence bands), are also implemented
#'   following the results in Cattaneo, Jansson and Ma (2020a) and Cattaneo, Jansson and Ma (2020b).
#'   See Cattaneo, Jansson and Ma (2020c) for more implementation details and illustrations.
#'
#'   Companion command: \code{\link{lpbwdensity}} for bandwidth selection.
#'
#'   Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
#'   available at \url{https://sites.google.com/site/nppackages/}.
#'
#' @param data Numeric vector or one dimensional matrix/data frame, the raw data.
#' @param grid Numeric, specifies the grid of evaluation points. When set to default, grid points
#'   will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.
#' @param bw Numeric, specifies the bandwidth
#'   used for estimation. Can be (1) a positive scalar (common
#'   bandwidth for all grid points); or (2) a positive numeric vector specifying bandwidths for
#'   each grid point (should be the same length as \code{grid}).
#' @param p Nonnegative integer, specifies the order of the local polynomial used to construct point
#'   estimates. (Default is \code{2}.)
#' @param q Nonnegative integer, specifies the order of the local polynomial used to construct
#'   confidence intervals/bands (a.k.a. the bias correction order). Default is \code{p+1}. When set to be
#'   the same as \code{p}, no bias correction will be performed. Otherwise it should be
#'   strictly larger than \code{p}.
#' @param v Nonnegative integer, specifies the derivative of the distribution function to be estimated. \code{0} for
#'   the distribution function, \code{1} (default) for the density funtion, etc.
#' @param kernel String, specifies the kernel function, should be one of \code{"triangular"}, \code{"uniform"}, and
#'   \code{"epanechnikov"}.
#' @param scale Numeric, specifies how
#'   estimates are scaled. For example, setting this parameter to 0.5 will scale down both the
#'   point estimates and standard errors by half. Default is \code{1}. This parameter is useful if only
#'   part of the sample is employed for estimation, and should not be confused with \code{Cweights}
#'   or \code{Pweights}.
#' @param massPoints \code{TRUE} (default) or \code{FALSE}, specifies whether point estimates and standard errors
#'   should be adjusted if there are mass points in the data.
#' @param bwselect String, specifies the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Options are (1) \code{"mse-dpi"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point); (2) \code{"imse-dpi"} (integrated MSE-optimal bandwidth,
#'   common for all grid points); (3) \code{"mse-rot"} (rule-of-thumb bandwidth with Gaussian
#'   reference model); and (4) \code{"imse-rot"} (integrated rule-of-thumb bandwidth with Gaussian
#'   reference model).
#' @param stdVar \code{TRUE} (default) or \code{FALSE}, specifies whether the data should be standardized for
#'   bandwidth selection.
#' @param regularize \code{TRUE} (default) or \code{FALSE}, specifies whether the bandwidth should be
#'   regularized. When set to \code{TRUE}, the bandwidth is chosen such that the local region includes
#'   at least \code{nLocalMin} observations and at least \code{nUniqueMin} unique observations.
#' @param nLocalMin Nonnegative integer, specifies the minimum number of observations in each local neighborhood. This option
#'   will be ignored if \code{regularize=FALSE}. Default is \code{20+p+1}.
#' @param nUniqueMin Nonnegative integer, specifies the minimum number of unique observations in each local neighborhood. This option
#'   will be ignored if \code{regularize=FALSE}. Default is \code{20+p+1}.
#' @param Cweights, Numeric, specifies the weights used
#'   for counterfactual distribution construction. Should have the same length as the data.
#' @param Pweights Numeric, specifies the weights used
#'   in sampling. Should have the same length as the data.
#'
#' @return
#' \item{Estimate}{A matrix containing (1) \code{grid} (grid points), (2) \code{bw} (bandwidths),
#'   (3) \code{nh} (number of observations in each local neighborhood),
#'   (4) \code{nhu} (number of unique observations in each local neighborhood),
#'   (5) \code{f_p} (point estimates with p-th order local polynomial),
#'   (6) \code{f_q} (point estimates with q-th order local polynomial, only if option \code{q} is nonzero),
#'   (7) \code{se_p} (standard error corresponding to \code{f_p}), and (8) \code{se_q} (standard error
#'   corresponding to \code{f_q}).}
#' \item{CovMat_p}{The variance-covariance matrix corresponding to \code{f_p}.}
#' \item{CovMat_q}{The variance-covariance matrix corresponding to \code{f_q}.}
#' \item{opt}{A list containing options passed to the function.}
#'
#' @details Bias correction is only used for the construction of confidence intervals/bands, but not for point
#'   estimation. The point estimates, denoted by \code{f_p}, are constructed using local polynomial estimates
#'   of order \code{p}, while the centering of the confidence intervals/bands, denoted by \code{f_q}, are constructed
#'   using local polynomial estimates of order \code{q}. The confidence intervals/bands take the form:
#'   \code{[f_q - cv * SE(f_q) , f_q + cv * SE(f_q)]}, where \code{cv} denotes the appropriate critical value and \code{SE(f_q)}
#'   denotes an standard error estimate for the centering of the confidence interval/band. As a result,
#'   the confidence intervals/bands may not be centered at the point estimates because they have been bias-corrected.
#'   Setting \code{q} and \code{p} to be equal results on centered at the point estimate confidence intervals/bands,
#'   but requires undersmoothing for valid inference (i.e., (I)MSE-optimal bandwdith for the density point estimator
#'   cannot be used). Hence the bandwidth would need to be specified manually when \code{q=p}, and the
#'   point estimates will not be (I)MSE optimal. See Cattaneo, Jansson and Ma (2020a, 2020b) for details, and also
#'   Calonico, Cattaneo, and Farrell (2018, 2020) for robust bias correction methods.
#'
#'   Sometimes the density point estimates may lie outside of the confidence intervals/bands, which can happen
#'   if the underlying distribution exhibits high curvature at some evaluation point(s). One possible solution
#'   in this case is to increase the polynomial order \code{p} or to employ a smaller bandwidth.
#'
#' @references
#'   Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}. \emph{Journal of the American Statistical Association}, 113(522): 767-779.
#'
#'   Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020. \href{http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2020_CEopt.pdf}{Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}. Working paper.
#'
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
#' @seealso Supported methods: \code{\link{coef.lpdensity}}, \code{\link{confint.lpdensity}}, \code{\link{plot.lpdensity}}, \code{\link{print.lpdensity}}, \code{\link{summary.lpdensity}}, \code{\link{vcov.lpdensity}}.
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
#' # Plot the estimates and confidence intervals
#' plot(est1, legendTitle="My Plot", legendGroups=c("X"))
#'
#' # Plot the estimates and the 99% uniform confidence band
#' set.seed(42) # fix the seed for simulating critical values
#' plot(est1, alpha=0.01, CIuniform=TRUE, legendTitle="My Plot", legendGroups=c("X"))
#'
#' # Adding a histogram to the background
#' plot(est1, legendTitle="My Plot", legendGroups=c("X"),
#'   hist=TRUE, histData=X, histBreaks=seq(-1.5, 1.5, 0.25))
#'
#' @export
lpdensity <- function(data, grid=NULL, bw=NULL, p=NULL, q=NULL, v=NULL,
                      kernel=c("triangular", "uniform", "epanechnikov"),
                      scale=NULL, massPoints=TRUE,
                      bwselect=c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot"),
                      stdVar=TRUE,
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
    stop("Data should be numeric, and cannot be empty.\n")
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
      stop("Grid points should be numeric.\n")
    }
  }

  # bw
  if (length(bw) == 0) {
    flag_no_bw <- TRUE
    if (length(bwselect) == 0) {
      bwselect <- "mse-dpi"
    } else {
      bwselect <- tolower(bwselect[1])
      # for
    }
    ### BACKWARD COMPATIBILITY FOR VERSION 0.1
    if (bwselect == "mse" ) { bwselect <- "mse-dpi" }
    if (bwselect == "imse") { bwselect <- "imse-dpi" }
    if (bwselect == "rot" ) { bwselect <- "mse-rot" }
    if (bwselect == "irot") { bwselect <- "imse-rot" }
    ### END BACKWARD COMPATIBILITY FOR VERSION 0.1

    if (!bwselect%in%c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot")) stop("Incorrect bandwidth selection method specified.\n")
  } else if (length(bw) == 1) {
    if (!is.numeric(bw) | bw <= 0) {
      stop("Bandwidth incorrectly specified.\n")
    } else {
      flag_no_bw <- FALSE
      bw <- rep(bw, ng)
      bwselect <- "user provided"
    }
  } else {
    bw <- as.vector(bw)
    if (!is.numeric(bw)) {
      stop("Bandwidth incorrectly specified.\n")
    } else if (length(bw) != ng) {
      stop("Bandwidth has to be the same length as grid.\n")
    } else {
      bwselect <- "user provided"
      flag_no_bw <- FALSE
    }
  }

  # p
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 2
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }

  # q
  if (length(q) == 0) {
    flag_no_q <- TRUE
    q <- p + 1
  } else if ((length(q) != 1) | !(q[1]%in%c(0:20)) | (q[1]<p)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  } else {
    flag_no_q <- FALSE
  }

  # v
  if (length(v) == 0) {
    flag_no_v <- TRUE
    v <- min(1, p)
  } else if ((length(v) != 1) | !(v[1]%in%c(0:20)) | (v[1]>p)) {
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

  # scale
  if (length(scale) == 0) {
    flag_no_scale <- TRUE
    scale <- 1
  } else if (length(scale) != 1 | !is.numeric(scale) | scale <= 0) {
    stop("Scale incorrectly specified.\n")
  } else {
    flag_no_scale <- FALSE
  }

  # regularize
  if (length(regularize) == 0) {
    regularize <- TRUE
  } else if (length(regularize) > 1 | !regularize[1]%in%c(TRUE, FALSE)) {
    stop("Regularization parameter incorrectly specified.\n")
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
  # Bandwidth selection
  ################################################################################

  if (flag_no_bw & bwselect[1]=="mse-rot") {
    bw <- bw_ROT( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  }
  if (flag_no_bw & bwselect[1]=="imse-rot") {
    bw <- bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
    bw <- rep(bw, length(grid))
  }
  if (flag_no_bw & bwselect[1]=="mse-dpi") {
    bw <- bw_MSE( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
  }
  if (flag_no_bw & bwselect[1]=="imse-dpi") {
    bw <- bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
    bw <- rep(bw, length(grid))
  }

  ################################################################################
  # Point Estimation
  ################################################################################
  Temp_Result <- lpdensity_fn(data=data, grid=grid, bw=bw, p=p, q=q, v=v, kernel=kernel,
                       Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, showSE=TRUE)

  Temp_Result$Estimate[, c("f_p", "f_q", "se_p", "se_q")] <- Temp_Result$Estimate[, c("f_p", "f_q", "se_p", "se_q")] * scale
  rownames(Temp_Result$Estimate) <- 1:ng

  Temp_Result$CovMat_p <- Temp_Result$CovMat_p * scale^2
  Temp_Result$CovMat_q <- Temp_Result$CovMat_q * scale^2
  rownames(Temp_Result$CovMat_p) <- colnames(Temp_Result$CovMat_p) <- 1:ng
  rownames(Temp_Result$CovMat_q) <- colnames(Temp_Result$CovMat_q) <- 1:ng

  Result <- list(Estimate=Temp_Result$Estimate,
                 CovMat_p=Temp_Result$CovMat_p,
                 CovMat_q=Temp_Result$CovMat_q,
                 opt=list(
                   p=p, q=q, v=v, kernel=kernel, scale=scale, massPoints=massPoints, n=n, ng=ng,
                   bwselect=bwselect, stdVar=stdVar,
                   regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                   data_min=min(data), data_max=max(data),
                   grid_min=min(grid), grid_max=max(grid)
                 ))
  class(Result) <- c("lpdensity")

  ################################################################################
  # Return
  ################################################################################
  return(Result)
}
