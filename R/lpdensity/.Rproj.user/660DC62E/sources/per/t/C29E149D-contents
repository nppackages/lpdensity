################################################################################
#' @title lpdensity: Local Polynomial Density Estimation and Inference
#'
#' @description Without imposing stringent distributional assumptions or shape restrictions,
#'   nonparametric estimation has been popular in economics and other social
#'   sciences for counterfactual analysis, program evaluation, and policy recommendations.
#'   This package implements a novel density (and derivatives) estimator based on local polynomial
#'   regressions, documented in Cattaneo, Jansson and Ma (2020a) and Cattaneo, Jansson and Ma (2020b).
#'
#'   \code{\link{lpdensity}} implements the local polynomial regression based density (and derivatives)
#'   estimator. Robust bias-corrected inference methods, both pointwise (confidence intervals) and
#'   uniform (confidence bands), are also implemented. \code{\link{lpbwdensity}} implements the bandwidth
#'   selection methods. See Cattaneo, Jansson and Ma (2020c) for more implementation details and illustrations.
#'
#'   Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
#'   available at \url{https://sites.google.com/site/nppackages/}.
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
#' @importFrom MASS mvrnorm
#' @import ggplot2
#'
#' @aliases lpdensity-package
"_PACKAGE"
