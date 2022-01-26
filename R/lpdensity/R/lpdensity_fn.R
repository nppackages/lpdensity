################################################################################
#' Supporting Function for \code{\link{lpdensity}}
#'
#' \code{lpdensity_fn} implements the local polynomial density estimator. This
#'   function is for internal use, and there is no error handling or robustness check.
#'
#' Recommend: use \code{\link{lpdensity}}.
#'
#' @param data Numeric vector or one dimensional matrix/data frame, the raw data.
#' @param grid Numeric vector or one dimensional matrix/data frame, the grid on which
#'   density is estimated.
#' @param bw Numeric vector or one dimensional matrix/data frame, the bandwidth
#'   used for estimation. Should be strictly positive, and have the same length as
#'   \code{grid}.
#' @param p Integer, nonnegative, the order of the local-polynomial used to construct point
#'   estimates.
#' @param q Integer, nonnegative, the order of the local-polynomial used to construct
#'   confidence interval (a.k.a. the bias correction order).
#' @param v Integer, nonnegative, the derivative of distribution function to be estimated. \code{0} for
#'   the distribution function, \code{1} (default) for the density funtion, etc.
#' @param kernel, String, the kernel function, should be one of \code{"triangular"},
#'   \code{"uniform"} or \code{"epanechnikov"}.
#' @param massPoints Boolean, whether whether point estimates and standard errors
#'   should be corrected if there are mass points in the data.
#' @param Cweights Numeric vector or one dimensional matrix/data frame, the weights used
#'   for counterfactual distribution construction. Should have the same length as sample size.
#' @param Pweights Numeric vector or one dimensional matrix/data frame, the weights used
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
lpdensity_fn <- function(data, grid, bw, p, q, v, kernel, Cweights, Pweights, massPoints, showSE=TRUE) {

  # preparation
  ii <- order(data, decreasing=FALSE)
  data <- data[ii]
  Cweights <- Cweights[ii]
  Pweights <- Pweights[ii]
  n    <- length(data)
  ng   <- length(grid)

  dataUnique  <- lpdensityUnique(data)
  freqUnique  <- dataUnique$freq
  indexUnique <- dataUnique$index
  dataUnique  <- dataUnique$unique
  nUnique     <- length(dataUnique)

  # whether considering mass points when constructing the empirical distribution function
  if (massPoints) {
    Fn <- rep((cumsum(Cweights * Pweights) / sum(Cweights * Pweights))[indexUnique], times=freqUnique)
  } else {
    Fn <- cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
  }

  weights_normal <- Cweights * Pweights / sum(Cweights * Pweights) * n
  Cweights <- Cweights / sum(Cweights) * n
  Pweights <- Pweights / sum(Pweights) * n

  weights_normalUnique <- cumsum(weights_normal)[indexUnique]
  if (nUnique > 1) { weights_normalUnique <- weights_normalUnique - c(0, weights_normalUnique[1:(nUnique-1)]) }

  CweightsUnique <- cumsum(Cweights)[indexUnique]
  if (nUnique > 1) { CweightsUnique <- CweightsUnique - c(0, CweightsUnique[1:(nUnique-1)]) }

  PweightsUnique <- cumsum(Pweights)[indexUnique]
  if (nUnique > 1) { PweightsUnique <- PweightsUnique - c(0, PweightsUnique[1:(nUnique-1)]) }

  hat_p <- rep(NA, ng); se_p <- rep(NA, ng); iff_p <- matrix(NA, nrow=n, ncol=ng)
  hat_q <- rep(NA, ng); se_q <- rep(NA, ng); iff_q <- matrix(NA, nrow=n, ncol=ng)
  nh <- rep(NA, ng)
  nhu <- rep(NA, ng)

  for (j in 1:ng) {
    index_temp <- abs(data - grid[j]) <= bw[j]
    nh[j]  <- sum(index_temp)
    nhu[j] <- sum(index_temp[indexUnique])

    # LHS and RHS variables
    if (massPoints) {
      Y_temp    <- matrix(Fn[indexUnique], ncol=1)
      Xh_temp   <- matrix((data[indexUnique] - grid[j]), ncol=1) / bw[j]
      Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
      if (p == 0) {
        Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
      }

      # weights
      if (kernel == "triangular") {
        Kh_temp <- ((1 - abs(Xh_temp)) / bw[j]) * index_temp[indexUnique]
      } else if (kernel == "uniform") {
        Kh_temp <- (0.5 / bw[j]) * index_temp[indexUnique]
      } else {
        Kh_temp <- (0.75 * (1 - Xh_temp^2) / bw[j]) * index_temp[indexUnique]
      }

      Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
      Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=PweightsUnique)

      XhKhXh_inv <- try(
        solve(t(Xh_p_temp[index_temp[indexUnique], , drop=FALSE]) %*% Xh_p_Kh_Pweights_temp[index_temp[indexUnique], , drop=FALSE] / n)
        , silent=TRUE)
      if (is.character(XhKhXh_inv)) { next }

      # point estimate
      hat_p[j] <- factorial(v) * (XhKhXh_inv %*% t(Xh_p_Kh_Pweights_temp[index_temp[indexUnique], , drop=FALSE]) %*% Y_temp[index_temp[indexUnique], , drop=FALSE])[v+1] / bw[j]^v / n
    } else {
      Y_temp    <- matrix(Fn, ncol=1)
      Xh_temp   <- matrix((data - grid[j]), ncol=1) / bw[j]
      Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
      if (p == 0) {
        Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
      }

      # weights
      if (kernel == "triangular") {
        Kh_temp <- ((1 - abs(Xh_temp)) / bw[j]) * index_temp
      } else if (kernel == "uniform") {
        Kh_temp <- (0.5 / bw[j]) * index_temp
      } else {
        Kh_temp <- (0.75 * (1 - Xh_temp^2) / bw[j]) * index_temp
      }

      Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
      Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=Pweights)

      XhKhXh_inv <- try(
        solve(t(Xh_p_temp[index_temp, , drop=FALSE]) %*% Xh_p_Kh_Pweights_temp[index_temp, , drop=FALSE] / n)
        , silent=TRUE)
      if (is.character(XhKhXh_inv)) { next }

      # point estimate
      hat_p[j] <- factorial(v) * (XhKhXh_inv %*% t(Xh_p_Kh_Pweights_temp[index_temp, , drop=FALSE]) %*% Y_temp[index_temp, , drop=FALSE])[v+1] / bw[j]^v / n
    }
    if (showSE) {
      # standard error estimate
      if (massPoints) {
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n

        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)
        for (jj in 1:ncol(G)) {
          G[, jj] <- (rep((cumsum(Xh_p_Kh_Pweights_temp[nUnique:1, jj]) / n)[nUnique:1], times=freqUnique) - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }

        iff_p[, j] <- (XhKhXh_inv %*% t(G))[v+1, ] * factorial(v) / sqrt(n*bw[j]^(2*v))
      } else {
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n

        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)
        for (jj in 1:ncol(G)) {
          G[, jj] <- ((cumsum(Xh_p_Kh_Pweights_temp[n:1, jj]) / n)[n:1] - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }

        iff_p[, j] <- (XhKhXh_inv %*% t(G))[v+1, ] * factorial(v) / sqrt(n*bw[j]^(2*v))
      }

    }

    if (q > p) {
      if (massPoints) {
        Xh_q_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:q)))
        Xh_q_Kh_temp <- sweep(Xh_q_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_q_Kh_Pweights_temp <- sweep(Xh_q_Kh_temp, MARGIN=1, FUN="*", STATS=PweightsUnique)

        XhKhXh_inv <- try(
          solve(t(Xh_q_temp[index_temp[indexUnique], , drop=FALSE]) %*% Xh_q_Kh_Pweights_temp[index_temp[indexUnique], , drop=FALSE] / n)
          , silent=TRUE)
        if (is.character(XhKhXh_inv)) { next }

        # point estimate
        hat_q[j] <- factorial(v) * (XhKhXh_inv %*% t(Xh_q_Kh_Pweights_temp[index_temp[indexUnique], , drop=FALSE]) %*% Y_temp[index_temp[indexUnique], , drop=FALSE])[v+1] / bw[j]^v / n
      } else {
        Xh_q_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:q)))
        Xh_q_Kh_temp <- sweep(Xh_q_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_q_Kh_Pweights_temp <- sweep(Xh_q_Kh_temp, MARGIN=1, FUN="*", STATS=Pweights)

        XhKhXh_inv <- try(
          solve(t(Xh_q_temp[index_temp, , drop=FALSE]) %*% Xh_q_Kh_Pweights_temp[index_temp, , drop=FALSE] / n)
          , silent=TRUE)
        if (is.character(XhKhXh_inv)) { next }

        # point estimate
        hat_q[j] <- factorial(v) * (XhKhXh_inv %*% t(Xh_q_Kh_Pweights_temp[index_temp, , drop=FALSE]) %*% Y_temp[index_temp, , drop=FALSE])[v+1] / bw[j]^v / n
      }


      if (showSE) {
        # standard error estimate
        if (massPoints) {
          F_Xh_q_Kh_temp <- t(Y_temp) %*% Xh_q_Kh_Pweights_temp / n

          G <- matrix(NA, ncol=ncol(Xh_q_Kh_Pweights_temp), nrow=n)
          for (jj in 1:ncol(G)) {
            G[, jj] <- (rep((cumsum(Xh_q_Kh_Pweights_temp[nUnique:1, jj]) / n)[nUnique:1], times=freqUnique) - F_Xh_q_Kh_temp[1, jj]) * weights_normal
          }
          iff_q[, j] <- (XhKhXh_inv %*% t(G))[v+1, ] * factorial(v) / sqrt(n*bw[j]^(2*v))
        } else {
          F_Xh_q_Kh_temp <- t(Y_temp) %*% Xh_q_Kh_Pweights_temp / n

          G <- matrix(NA, ncol=ncol(Xh_q_Kh_Pweights_temp), nrow=n)
          for (jj in 1:ncol(G)) {
            G[, jj] <- ((cumsum(Xh_q_Kh_Pweights_temp[n:1, jj]) / n)[n:1] - F_Xh_q_Kh_temp[1, jj]) * weights_normal
          }
          iff_q[, j] <- (XhKhXh_inv %*% t(G))[v+1, ] * factorial(v) / sqrt(n*bw[j]^(2*v))
        }
      }
    }
  }

  CovMat_p <- t(iff_p) %*% iff_p / n
  CovMat_q <- t(iff_q) %*% iff_q / n

  se_p <- sqrt(abs(diag(CovMat_p)))
  se_q <- sqrt(abs(diag(CovMat_q)))

  Estimate <- cbind(grid, bw, nh, nhu, hat_p, hat_q, se_p, se_q)
  colnames(Estimate) <- c("grid", "bw", "nh", "nhu", "f_p", "f_q", "se_p", "se_q")
  rownames(Estimate) <- c()

  return(list(Estimate=Estimate, CovMat_p=CovMat_p, CovMat_q=CovMat_q))
}


