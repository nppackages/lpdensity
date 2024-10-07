################################################################################
#' @title lpdensity: Local Polynomial Density Estimation and Inference
#'
#' @description Without imposing stringent distributional assumptions or shape restrictions,
#'   nonparametric estimation has been popular in economics and other social
#'   sciences for counterfactual analysis, program evaluation, and policy recommendations.
#'   This package implements a novel density (and derivatives) estimator based on local polynomial
#'   regressions, documented in Cattaneo, Jansson and Ma (2020, 2023).
#'
#'   \code{\link{lpdensity}} implements the local polynomial regression based density (and derivatives)
#'   estimator. Robust bias-corrected inference methods, both pointwise (confidence intervals) and
#'   uniform (confidence bands), are also implemented. \code{\link{lpbwdensity}} implements the bandwidth
#'   selection methods. See Cattaneo, Jansson and Ma (2022) for more implementation details and illustrations.
#'
#'   Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
#'   available at \url{https://nppackages.github.io/}.
#'
#' @references
#'   Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018.
#'     On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference.
#'     \emph{Journal of the American Statistical Association}, 113(522): 767-779.
#'     \doi{10.1080/01621459.2017.1285776}
#'
#'   Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2022.
#'     Coverage Error Optimal Confidence Intervals for Local Polynomial Regression.
#'     \emph{Bernoulli}, 28(4): 2998-3022.
#'     \doi{10.3150/21-BEJ1445}
#'
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2020.
#'     Simple Local Polynomial Density Estimators.
#'     \emph{Journal of the American Statistical Association}, 115(531): 1449-1455.
#'     \doi{10.1080/01621459.2019.1635480}
#'
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2022.
#'     lpdensity: Local Polynomial Density Estimation and Inference.
#'     \emph{Journal of Statistical Software}, 101(2): 1–25.
#'     \doi{10.18637/jss.v101.i02}
#'
#'   Cattaneo, M. D., M. Jansson, and X. Ma. 2023.
#'     Local Regression Distribution Estimators.
#'     \emph{Journal of Econometrics}, 240(2): 105074.
#'     \doi{10.1016/j.jeconom.2021.01.006}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats D
#' @importFrom stats integrate
#' @importFrom stats optimize
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @importFrom stats density
#' @importFrom MASS mvrnorm
#' @import ggplot2
#'
#' @aliases lpdensity-package
"_PACKAGE"
