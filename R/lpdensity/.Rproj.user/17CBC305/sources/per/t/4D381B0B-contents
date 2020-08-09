################################################################################
#' Supporting Function for \code{\link{lpdensity}}
#'
#' \code{lpdensity_fn} implements the local polynomial density estimator. This
#'   function is for internal use, and there is no error handling or robustness check.
#'
#' Recommend: use \code{\link{lpdensity}}.
#'
#' @param data Numeric vector or one dimensional matrix / data frame, the raw data.
#' @param grid Numeric vector or one dimensional matrix / data frame, the grid on which
#'   density is estimated.
#' @param bw Numeric vector or one dimensional matrix / data frame, the bandwidth
#'   used for estimation. Should be strictly positive, and have the same length as
#'   \code{grid}.
#' @param p Integer, the order of the local-polynomial used to construct point
#'   estimates. Should be greater than 1.
#' @param q Integer, the order of the local-polynomial used to construct point
#'   estimates. Should be greater than 1. When set to \code{0}, corresponding estimates
#'   will not be constructed.
#' @param v Integer, the derivative to be estimated. Should be nonnegative.
#' @param kernel, String, the kernel function, should be one of \code{"triangular"},
#'   \code{"uniform"} or \code{"epanechnikov"}.
#' @param Cweights Numeric vector or one dimensional matrix / data frame, the weights used
#'   for counterfactual distribution construction. Should have the same length as sample size.
#' @param Pweights Numeric vector or one dimensional matrix / data frame, the weights used
#'   in sampling. Should have the same length as sample size, and nonnegative.
#' @param showSE \code{TRUE} (default) or \code{FALSE}, whether standard errors should be computed.
#'
#' @return
#' \item{grid}{grid points.}
#' \item{bw}{bandwidth for each grid point.}
#' \item{nh}{Effective sample size for each grid point.}
#' \item{f_p}{Density estimates on the grid with local polynomial of order \code{p},
#'   with the same length as \code{grid}.}
#' \item{f_q}{Density estimates on the grid with local polynomial of order \code{q},
#'   with the same length as \code{grid}. This is reported only if \code{q} is greater than
#'   0.}
#' \item{se_p}{Standard errors corresponding to \code{hat_p}.}
#' \item{se_q}{Standard errors corresponding to \code{hat_q}. This is reported only
#'   if \code{q} is greater than 0.}
#'
#' @keywords internal
lpdensity_fn <- function(data, grid, bw, p, q, v, kernel, Cweights, Pweights, showSE=TRUE) {

  # preparation
  ii <- order(data)
  data <- data[ii]
  Cweights <- Cweights[ii]
  Pweights <- Pweights[ii]
  n    <- length(data)
  ng   <- length(grid)
  Fn   <- cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
  nh   <- rep(NA, ng)
  weights_normal <- Cweights * Pweights / sum(Cweights * Pweights) * n
  Cweights <- Cweights / sum(Cweights) * n
  Pweights <- Pweights / sum(Pweights) * n

  hat_p <- rep(NA, ng); se_p  <- rep(NA, ng)
  hat_q <- rep(NA, ng); se_q <- rep(NA, ng)

  for (j in 1:ng) {
    index_temp <- abs(data - grid[j]) <= bw[j]
    nh[j] <- sum(index_temp)

    # LHS and RHS variables
    Y_temp    <- matrix(Fn[index_temp], ncol=1)
    Xh_temp   <- matrix((data - grid[j])[index_temp], ncol=1) / bw[j]
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
    if (p == 0) Xh_p_temp <- t(Xh_p_temp)

    # weights
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / bw[j]
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / bw[j]
    } else {
      Kh_temp <- 0.75 * (1 - Xh_temp^2) / bw[j]
    }
    Kh_temp   <- Pweights[index_temp] * Kh_temp

    XhKh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
    XhKhXh_inv <- try(
      solve(t(Xh_p_temp) %*% XhKh_temp / n)
      , silent=TRUE)
    if (is.character(XhKhXh_inv)) { next }

    # point estimate
    hat_p[j] <- factorial(v) * (XhKhXh_inv %*% t(XhKh_temp) %*% Y_temp)[v+1] / bw[j]^v / n

    if (showSE) {
    # standard error estimate
    F_XhKh_temp <- matrix(Fn[index_temp], nrow=1) %*% XhKh_temp / n
    G <- XhKh_temp[nh[j]:1, ]
    for (jj in 1:ncol(G)) {
      G[, jj] <- cumsum(G[, jj]) / n - F_XhKh_temp[1, jj]
    }

    G <- sweep(G, MARGIN=1, FUN="*", STATS=weights_normal[index_temp])
    G <- t(G) %*% G / n

    index_temp_1 <- data - grid[j] < -1 * bw[j]
    index_temp_2 <- data - grid[j] >      bw[j]

    G1 <- matrix(1 - Fn[index_temp], nrow=1)
    G1 <- G1 %*% XhKh_temp / n
    G1 <- t(G1) %*% G1 * sum(weights_normal[index_temp_1]^2) / n
    G2 <- matrix(0 - Fn[index_temp], nrow=1)
    G2 <- G2 %*% XhKh_temp / n
    G2 <- t(G2) %*% G2 * sum(weights_normal[index_temp_2]^2) / n

    V <- XhKhXh_inv %*% (G+G1+G2) %*% XhKhXh_inv

    se_p[j] <- factorial(v) * sqrt( V[v+1,v+1] / (n * bw[j]^(2*v)) )
    }

    if (q > p) {
      Xh_q_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:q)))

      XhKh_temp <- sweep(Xh_q_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
      XhKhXh_inv <- try(
        solve(t(Xh_q_temp) %*% XhKh_temp / n)
        , silent=TRUE)
      if (is.character(XhKhXh_inv)) { next }

      # point estimate
      hat_q[j] <- factorial(v) * (XhKhXh_inv %*% t(XhKh_temp) %*% Y_temp)[v+1] / bw[j]^v / n

      if (showSE) {
      # standard error estimate
      F_XhKh_temp <- matrix(Fn[index_temp], nrow=1) %*% XhKh_temp / n
      G <- XhKh_temp[nh[j]:1, ]
      for (jj in 1:ncol(G)) {
        G[, jj] <- cumsum(G[, jj]) / n - F_XhKh_temp[1, jj]
      }

      G <- sweep(G, MARGIN=1, FUN="*", STATS=weights_normal[index_temp])
      G <- t(G) %*% G / n

      index_temp_1 <- data - grid[j] < -1 * bw[j]
      index_temp_2 <- data - grid[j] >      bw[j]

      G1 <- matrix(1 - Fn[index_temp], nrow=1)
      G1 <- G1 %*% XhKh_temp / n
      G1 <- t(G1) %*% G1 * sum(weights_normal[index_temp_1]^2) / n
      G2 <- matrix(0 - Fn[index_temp], nrow=1)
      G2 <- G2 %*% XhKh_temp / n
      G2 <- t(G2) %*% G2 * sum(weights_normal[index_temp_2]^2) / n

      V <- XhKhXh_inv %*% (G+G1+G2) %*% XhKhXh_inv

      se_q[j] <- factorial(v) * sqrt( V[v+1,v+1] / (n*bw[j]^(2*v)) )
      }
    }
  }

  Estimate <- cbind(grid, bw, nh, hat_p, hat_q, se_p, se_q)
  colnames(Estimate) <- c("grid", "bw", "nh", "f_p", "f_q", "se_p", "se_q")
  rownames(Estimate) <- c()

  return (Estimate)
}


