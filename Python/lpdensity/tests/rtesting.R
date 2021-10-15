################################################################################
#' Internal function.
#'
#' Find unique elements and their frequencies in a numeric vector. This function
#'   has a similar performance as the built-in R function \code{unique}.
#'
#' @param x Numeric vector, already sorted in ascending order.
#'
#' @return
#' \item{unique}{A vector containing unique elements in \code{x}.}
#' \item{freq}{The frequency of each element in \code{x}.}
#' \item{index}{The last occurrence of each element in \code{x}.}
#'
#' @keywords internal
lpdensityUnique <- function(x) {
  n <- length(x)
  # if x has one or no element
  if (n == 0) return(list(unique=NULL, freq=c(), index=c()))
  if (n == 1) return(list(unique=x, freq=1, index=1))
  
  # else
  uniqueIndex <- c(x[2:n] != x[1:(n-1)], TRUE)
  unique <- x[uniqueIndex]
  nUnique <- length(unique)
  
  # all are distinct
  if (nUnique == n) return(list(unique=unique, freq=rep(1,length(x)), index=1:n))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=n, index=n))
  
  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1
  
  return(list(unique=unique, freq=freq, index=(1:n)[uniqueIndex]))
}

################################################################################
#' Internal function.
#'
#' Calculates density and higher order derivatives for Gaussian models.
#'
#' @param x Scalar, point of evaluation.
#' @param v Nonnegative integer, the derivative order (0 indicates cdf, 1 indicates pdf, etc.).
#' @param mean Scalar, the mean of the normal distribution.
#' @param sd Strictly positive scalar, the standard deviation of the normal distribution.
#'
#' @return
#' \item{}{Scalar.}
#'
#' @keywords internal
normal_dgps <- function(x, v, mean, sd) {
  if (v == 0) {
    return(pnorm(x, mean=mean, sd=sd))
  } else {
    temp <- expression(exp(-(x-mean)^2/(2*sd^2))/sqrt(2*pi*sd^2))
    while(v > 1) {
      temp <- D(temp, "x")
      v <- v - 1
    }
    return(eval(temp, list(x=x, mean=mean, sd=sd)))
  }
}

################################################################################
#' Internal function.
#'
#' Generate matrix.
#'
#' @param p Nonnegative integer, polynomial order.
#' @param low,up Scalar, between -1 and 1, the region of integration.
#' @param kernel String, the kernel function.
#'
#' @return
#' \item{}{A (p+1)-by-(p+1) matrix.}
#'
#' @keywords internal
Sgenerate <- function(p, low=-1, up=1, kernel="triangular") {
  S <- matrix(rep(0, (p+1)^2), ncol=(p+1))
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      if (kernel == "uniform") {
        integrand <- function(x) { x^(i+j-2)*0.5 }
      } else if (kernel == "epanechnikov") {
        integrand <- function(x) { x^(i+j-2)*0.75*(1-x^2) }
      } else {
        integrand <- function(x) { x^(i+j-2)*(1-abs(x)) }
      }
      S[i,j] <- (integrate(integrand, lower=low, upper=up)$value)
    }
  }
  return(S)
}

################################################################################
#' Internal function.
#'
#' Generate matrix.
#'
#' @param p Nonnegative integer, polynomial order.
#' @param low,up Scalar, between -1 and 1, the region of integration.
#' @param kernel String, the kernel function.
#'
#' @return
#' \item{}{A (p+1)-by-(p+1) matrix.}
#'
#' @keywords internal
Tgenerate <- function(p, low=-1, up=1, kernel="triangular") {
  S <- matrix(rep(0, (p+1)^2), ncol=(p+1))
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      if (kernel == "uniform") {
        integrand <- function(x) { x^(i+j-2) * 0.5^2 }
      } else if (kernel == "epanechnikov") {
        integrand <- function(x) { x^(i+j-2) * (0.75*(1-x^2))^2 }
      } else {
        integrand <- function(x) { x^(i+j-2) * (1-abs(x))^2 }
      }
      S[i,j] <- (integrate(integrand, lower=low, upper=up)$value)
    }
  }
  return(S)
}

################################################################################
#' Internal function.
#'
#' Generate matrix.
#'
#' @param k Nonnegative integer, extra order (usually p+1).
#' @param p Nonnegative integer, the polynomial order.
#' @param low,up Scalar, between -1 and 1, region of integration.
#' @param kernel String, the kernel function.
#'
#' @return
#' \item{}{A (p+1)-by-1 matrix.}
#'
#' @keywords internal
Cgenerate <- function(k, p, low=-1, up=1, kernel="triangular") {
  C <- matrix(rep(0, (p+1)), ncol=1)
  for (i in 1:(p+1)) {
    if (kernel == "uniform") {
      integrand <- function(x) { x^(i+k-1)*0.5 }
    } else if (kernel == "epanechnikov") {
      integrand <- function(x) { x^(i+k-1)*0.75*(1-x^2) }
    }
    else {
      integrand <- function(x) { x^(i+k-1)*(1-abs(x)) }
    }
    C[i,1] <- (integrate(integrand, lower=low, upper=up)$value)
  }
  return(C)
}

################################################################################
#' Internal function.
#'
#' Generate matrix.
#'
#' @param p Nonnegative integer, polynomial order.
#' @param low,up Scalar, between -1 and 1, the region of integration.
#' @param kernel String, the kernel function.
#'
#' @return
#' \item{}{A (p+1)-by-(p+1) matrix.}
#'
#' @keywords internal
Ggenerate <- function(p, low=-1, up=1, kernel="triangular") {
  G <- matrix(rep(0, (p+1)^2), ncol=(p+1))
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      if (kernel == "uniform") {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1)*0.25, low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j*0.25, y, up)$value
            })
          }, low, up)$value
      } else if (kernel == "epanechnikov") {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1) * 0.75^2 *
                        (1-x^2) * (1-y^2), low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j * 0.75^2 *
                          (1-x^2) * (1-y^2), y, up)$value
            })
          }, low, up)$value
      } else {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1) *
                        (1-abs(x)) * (1-abs(y)), low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j *
                          (1-abs(x)) * (1-abs(y)), y, up)$value
            })
          }, low, up)$value
      }
    }
  }
  return(G)
}

################################################################################
#' Internal function.
#'
#' Calculates rule-of-thumb bandwidth
#'
#' @param data Numeric vector, the data.
#' @param grid Numeric vector, the evaluation points.
#' @param p Integer, polynomial order.
#' @param v Integer, order of derivative.
#' @param kernel String, the kernel.
#' @param Cweights Numeric vector, the counterfactual weights.
#' @param Pweights Numeric vector, the survey sampling weights.
#' @param massPoints Boolean, whether whether point estimates and standard errors
#'   should be corrected if there are mass points in the data.
#' @param stdVar Boolean, whether the data should be standardized for bandwidth selection.
#' @param regularize Boolean, whether the bandwidth should be regularized.
#' @param nLocalMin Nonnegative integer, minimum number of observations in each local
#'   neighborhood.
#' @param nUniqueMin Nonnegative integer, minimum number of unique observations in
#'   each local neighborhood.
#'
#' @return
#' \item{}{Bandwidth sequence.}
#'
#' @keywords internal
bw_ROT  <- function(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin) {
  
  n  <- length(data)
  ng <- length(grid)
  
  dataUnique <- unique(data)
  nUnique <- length(dataUnique)
  
  if (stdVar) {
    center_temp <- mean(data)
    scale_temp <- sd(data)
    data <- (data - center_temp) / scale_temp
    dataUnique <- (dataUnique - center_temp) / scale_temp
    grid <- (grid - center_temp) / scale_temp
  }
  
  # estimate a normal reference model
  mean_hat <- sum(Cweights * Pweights * data) / sum(Cweights * Pweights)
  sd_hat   <- sqrt(sum(Cweights * Pweights * (data - mean_hat)^2) / sum(Cweights * Pweights))
  
  # normal quantities
  temp_1 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  temp_2 <- D(expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2)), "x")
  temp_3 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  temp_4 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  
  j <- p + 1
  while(j > 1) {
    temp_3 <- D(temp_3, "x")
    j <- j - 1
  }
  
  j <- p + 2
  while(j > 1) {
    temp_4 <- D(temp_4, "x")
    j <- j - 1
  }
  
  # bias estimate, no rate added, DGP constant
  bias_dgp <- matrix(NA, ncol=2, nrow=ng)
  for (j in 1:ng) {
    # this comes from a higher-order bias expansion. See Lemma 5 in the Appendix of Cattaneo, Jansson and Ma (2019a)
    bias_dgp[j, 1] <- eval(temp_3, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+1) * factorial(v)
    bias_dgp[j, 2] <- eval(temp_4, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+2) * factorial(v) +
      bias_dgp[j, 1] * eval(temp_2, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat))
  }
  
  # bias estimate, no rate added, kernel constant
  S  <- Sgenerate(       p=p, low=-1, up=1, kernel=kernel)
  C1 <- Cgenerate(k=p+1, p=p, low=-1, up=1, kernel=kernel)
  C2 <- Cgenerate(k=p+2, p=p, low=-1, up=1, kernel=kernel)
  G  <- Ggenerate(       p=p, low=-1, up=1, kernel=kernel)
  S2 <- Tgenerate(       p=p, low=-1, up=1, kernel=kernel)
  bias_dgp[, 1] <- bias_dgp[, 1] * (solve(S) %*% C1)[v+1, ]
  bias_dgp[, 2] <- bias_dgp[, 2] * (solve(S) %*% C2)[v+1, ]
  
  # variance estimate, sample size added
  sd_dgp <- matrix(NA, ncol=1, nrow=ng)
  if (v > 0) {
    for (j in 1:ng) {
      sd_dgp[j, 1] <- factorial(v) * sqrt(eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / n)
    }
    sd_dgp <- sd_dgp * sqrt(abs((solve(S) %*% G %*% solve(S))[v+1, v+1]))
  } else {
    for (j in 1:ng) {
      # this comes from a higher-order variance expansion. See Lemma 4 in the Appendix of Cattaneo, Jansson and Ma (2019a)
      sd_dgp[j, 1] <- sqrt(
        pnorm(grid[j], mean=mean_hat, sd=sd_hat) * pnorm(grid[j], mean=mean_hat, sd=sd_hat, lower.tail=FALSE) /
          dnorm(grid[j], mean=mean_hat, sd=sd_hat) / (0.5*n^2))
    }
    sd_dgp <- sd_dgp * sqrt(abs((solve(S) %*% S2 %*% solve(S))[v+1, v+1]))
  }
  
  # bandwidth
  h <- rep(NA, ng)
  for (j in 1:ng) {
    if (v > 0) {
      opt.f <- function(a) {
        a^(2*p+2-2*v) * (bias_dgp[j, 1] + a * bias_dgp[j, 2])^2 + sd_dgp[j, 1]^2 / a^(2*v - 1)
      }
    } else {
      opt.f <- function(a) {
        a^(2*p+2-2*v) * (bias_dgp[j, 1] + a * bias_dgp[j, 2])^2 + sd_dgp[j, 1]^2 / a
      }
    }
    h[j] <- optimize(opt.f, interval=c(.Machine$double.eps, max(data) - min(data)), maximum=FALSE)$minimum
  }
  
  for (j in 1:ng) {
    if (is.na(h[j])) {
      h[j] <- sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))]
    }
    if (regularize) {
      if (nLocalMin > 0) { h[j] <- max(h[j], sort(abs(data-grid[j]))[min(n, nLocalMin)]) }
      if (nUniqueMin > 0) { h[j] <- max(h[j], sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)]) }
    }
    h[j] <- min(h[j], max(abs(dataUnique-grid[j])))
  }
  
  if (stdVar) {
    h <- h * scale_temp
  }
  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates integrated rule-of-thumb bandwidth
#'
#' @param data Numeric vector, the data.
#' @param grid Numeric vector, the evaluation points.
#' @param p Integer, polynomial order.
#' @param v Integer, order of derivative.
#' @param kernel String, the kernel.
#' @param Cweights Numeric vector, the counterfactual weights.
#' @param Pweights Numeric vector, the survey sampling weights.
#' @param massPoints Boolean, whether point estimates and standard errors
#'   should be corrected if there are mass points in the data.
#' @param stdVar Boolean, whether the data should be standardized for bandwidth selection.
#' @param regularize Boolean, Whether the bandwidth should be regularized.
#' @param nLocalMin Nonnegative integer, minimum number of observations in each local neighborhood.
#' @param nUniqueMin Nonnegative integer, minimum number of unique observations in each local neighborhood.
#'
#' @return
#' \item{}{A single bandwidth.}
#'
#' @keywords internal
bw_IROT <- function(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin) {
  
  n  <- length(data)
  ng <- length(grid)
  
  dataUnique <- unique(data)
  nUnique <- length(dataUnique)
  
  if (stdVar) {
    center_temp <- mean(data)
    scale_temp <- sd(data)
    data <- (data - center_temp) / scale_temp
    dataUnique <- (dataUnique - center_temp) / scale_temp
    grid <- (grid - center_temp) / scale_temp
  }
  
  # estimate a normal reference model
  mean_hat <- sum(Cweights * Pweights * data) / sum(Cweights * Pweights)
  sd_hat   <- sqrt(sum(Cweights * Pweights * (data - mean_hat)^2) / sum(Cweights * Pweights))
  
  # normal quantities
  temp_1 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  temp_2 <- D(expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2)), "x")
  temp_3 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  temp_4 <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi*sd^2))
  
  j <- p + 1
  while(j > 1) {
    temp_3 <- D(temp_3, "x")
    j <- j - 1
  }
  
  j <- p + 2
  while(j > 1) {
    temp_4 <- D(temp_4, "x")
    j <- j - 1
  }
  
  # bias estimate, no rate added, DGP constant
  bias_dgp <- matrix(NA, ncol=2, nrow=ng)
  for (j in 1:ng) {
    # this comes from a higher-order bias expansion. See Lemma 5 in the Appendix of Cattaneo, Jansson and Ma (2019a)
    bias_dgp[j, 1] <- eval(temp_3, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+1) * factorial(v)
    bias_dgp[j, 2] <- eval(temp_4, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+2) * factorial(v) +
      bias_dgp[j, 1] * eval(temp_2, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat))
  }
  
  # bias estimate, no rate added, kernel constant
  S  <- Sgenerate(       p=p, low=-1, up=1, kernel=kernel)
  C1 <- Cgenerate(k=p+1, p=p, low=-1, up=1, kernel=kernel)
  C2 <- Cgenerate(k=p+2, p=p, low=-1, up=1, kernel=kernel)
  S2 <- Tgenerate(       p=p, low=-1, up=1, kernel=kernel)
  G  <- Ggenerate(       p=p, low=-1, up=1, kernel=kernel)
  bias_dgp[, 1] <- bias_dgp[, 1] * (solve(S) %*% C1)[v+1, ]
  bias_dgp[, 2] <- bias_dgp[, 2] * (solve(S) %*% C2)[v+1, ]
  
  # variance estimate, sample size added
  sd_dgp <- matrix(NA, ncol=1, nrow=ng)
  if (v > 0) {
    for (j in 1:ng) {
      sd_dgp[j, 1] <- factorial(v) * sqrt(eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / n)
    }
    sd_dgp <- sd_dgp * sqrt(abs((solve(S) %*% G %*% solve(S))[v+1, v+1]))
  } else {
    for (j in 1:ng) {
      # this comes from a higher-order variance expansion. See Lemma 4 in the Appendix of Cattaneo, Jansson and Ma (2019a)
      sd_dgp[j, 1] <- sqrt(
        pnorm(grid[j], mean=mean_hat, sd=sd_hat) * pnorm(grid[j], mean=mean_hat, sd=sd_hat, lower.tail=FALSE) /
          dnorm(grid[j], mean=mean_hat, sd=sd_hat) / (0.5*n^2))
    }
    sd_dgp <- sd_dgp * sqrt(abs((solve(S) %*% S2 %*% solve(S))[v+1, v+1]))
  }

  # bandwidth
  if (v > 0) {
    opt.f <- function(a) {
      a^(2*p+2-2*v) * sum((bias_dgp[, 1] + a * bias_dgp[, 2])^2) + sum(sd_dgp[, 1]^2) / a^(2*v - 1)
    }
  } else {
    opt.f <- function(a) {
      a^(2*p+2-2*v) * sum((bias_dgp[, 1] + a * bias_dgp[, 2])^2) + sum(sd_dgp[, 1]^2) / a
    }
  }
  
  h <- optimize(opt.f, interval=c(.Machine$double.eps, max(data) - min(data)), maximum=FALSE)$minimum
  if (is.na(h)) {
    h <- 0
    for (j in 1:ng) {
      h <- max(h, sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))])
    }
  }
  if (regularize) {
    for (j in 1:ng) {
      if (nLocalMin > 0) { 
        h <- max(h, sort(abs(data-grid[j]))[min(n, nLocalMin)]) 
      }
      if (nUniqueMin > 0) { 
        h <- max(h, sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)]) 
        }
    }
  }
  h <- min(h, max(abs(max(dataUnique)-min(grid)), abs(min(dataUnique)-max(grid))))
  
  if (stdVar) {
    h <- h * scale_temp
  }
  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates MSE-optimal bandwidths.
#'
#' @param data Numeric vector, the data.
#' @param grid Numeric vector, the evaluation points.
#' @param p Integer, polynomial order.
#' @param v Integer, order of derivative.
#' @param kernel String, the kernel.
#' @param Cweights Numeric vector, the counterfactual weights.
#' @param Pweights Numeric vector, the survey sampling weights.
#' @param massPoints Boolean, whether whether point estimates and standard errors
#'   should be corrected if there are mass points in the data.
#' @param stdVar Boolean, whether the data should be standardized for bandwidth selection.
#' @param regularize Boolean, whether the bandwidth should be regularized.
#' @param nLocalMin Nonnegative integer, minimum number of observations in each local neighborhood.
#' @param nUniqueMin Nonnegative integer, minimum number of unique observations in each local neighborhood.
#'
#' @return
#' \item{}{Bandwidth sequence.}
#'
#' @keywords internal
bw_MSE  <- function(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin) {
  
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
  
  if (stdVar) {
    center_temp <- mean(data)
    scale_temp <- sd(data)
    data <- (data - center_temp) / scale_temp
    dataUnique <- (dataUnique - center_temp) / scale_temp
    grid <- (grid - center_temp) / scale_temp
  }
  
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
  
  # obtain preliminary bandwidth for estimating densities
  # this is used for constructing preasymptotic matrices
  h1  <- bw_IROT(data=data, grid=grid, p=2,   v=1,   kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+2+1,  nUniqueMin=20+2+1)
  
  # obtain preliminary bandwidth for estimating F_p+1
  # this is used for constructing F_p+1
  hp1 <- bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+p+2+1, nUniqueMin=20+p+2+1)
  
  # obtain preliminary bandwidth for estimating F_p+2
  # this is used for constructing F_p+2
  hp2 <- bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+p+3+1, nUniqueMin=20+p+3+1)

  dgp_hat <- matrix(NA, ncol=2, nrow=ng) # Fp+1 and Fp+2 with normalization constants
  const_hat <- matrix(NA, ncol=3, nrow=ng)
  h <- rep(NA, ng)
  
  for (j in 1:ng) {
    # estimate F_p+2
    index_temp <- abs(data-grid[j]) <= hp2
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / hp2
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:(p+3))))
    
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / hp2
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / hp2
    } else {
      Kh_temp <- 0.75 * (1 - (Xh_temp)^2) / hp2
    }
    Kh_temp   <- Pweights[index_temp] * Kh_temp
    Y_temp   <- matrix(Fn[index_temp], ncol=1)

    
    temp <- try(
      (solve(t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)) %*% t(Xh_p_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=Kh_temp))[p+3, 1] / (hp2^(p+2))
      , silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 2] <- temp
    
    # estimate F_p+1
    index_temp <- abs(data-grid[j]) <= hp1
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / hp1
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:(p+2))))
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / hp1
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / hp1
    } else {
      Kh_temp <- 0.75 * (1 - (Xh_temp)^2) / hp1
    }
    Kh_temp   <- Pweights[index_temp] * Kh_temp
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    
    temp <- try(
      (solve(t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)) %*% t(Xh_p_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=Kh_temp))[p+2, 1]  / (hp1^(p+1)),
      silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 1] <- temp
    
    # prepare for estimating matrices
    index_temp <- abs(data-grid[j]) <= h1
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / h1
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / h1
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / h1
    } else {
      Kh_temp <- 0.75*(1 - (Xh_temp)^2) / h1
    }
    #Kh_temp   <- Pweights[index_temp] * Kh_temp
    
    # estimate Cp matrix
    if (p > 0) {
      C_p_hat <- matrix(apply(
        sweep(
          apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p_hat <- matrix(apply(
        sweep(
          matrix(apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))), nrow=1),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }
    
    # estimate Cp+1 matrix
    if (p > 0) {
      C_p1_hat <- matrix(apply(
        sweep(
          apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p1_hat <- matrix(apply(
        sweep(
          matrix(apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))), nrow=1),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }
    
    # estimate S matirx
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
    if (p == 0) {
      Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
    }
    
    S_hat <- t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Pweights[index_temp] * Kh_temp) / n
    S_hat_inv <- try(solve(S_hat), silent=TRUE)
    if (is.character(S_hat_inv)) next
    
    # estimate G matrix
    if (v == 0) {
      G_hat <- t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Pweights[index_temp] * Kh_temp^2) / n
    } else {
      # for coding clarity, here we use the full sample and the influence function approach
      if (massPoints) {
        Y_temp    <- matrix(Fn[indexUnique], ncol=1)
        Xh_temp   <- matrix((dataUnique - grid[j]), ncol=1) / h1
        Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
        if (p == 0) {
          Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
        }
        
        if (kernel == "triangular") {
          Kh_temp <- ((1 - abs(Xh_temp)) / h1) * index_temp[indexUnique]
        } else if (kernel == "uniform") {
          Kh_temp <- (0.5 / h1) * index_temp[indexUnique]
        } else {
          Kh_temp <- (0.75 * (1 - Xh_temp^2) / h1) * index_temp[indexUnique]
        }
        
        Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=PweightsUnique)
        
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n
        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)

        for (jj in 1:ncol(G)) {
          G[, jj] <- (rep((cumsum(Xh_p_Kh_Pweights_temp[nUnique:1, jj]) / n)[nUnique:1], times=freqUnique) - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }
        G_hat <- t(G) %*% G / n
      } else {
        Y_temp    <- matrix(Fn, ncol=1)
        Xh_temp   <- matrix((data - grid[j]), ncol=1) / h1
        Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
        if (p == 0) {
          Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
        }
        
        if (kernel == "triangular") {
          Kh_temp <- ((1 - abs(Xh_temp)) / h1) * index_temp
        } else if (kernel == "uniform") {
          Kh_temp <- (0.5 / h1) * index_temp
        } else {
          Kh_temp <- (0.75 * (1 - Xh_temp^2) / h1) * index_temp
        }
        
        Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=Pweights)
        
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n
        
        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)
        for (jj in 1:ncol(G)) {
          G[, jj] <- ((cumsum(Xh_p_Kh_Pweights_temp[n:1, jj]) / n)[n:1] - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }
        G_hat <- t(G) %*% G / n
      }
    }
    
    # now get all constants
    const_hat[j, 1] <- factorial(v) * (S_hat_inv%*%C_p_hat)[v+1]
    const_hat[j, 2] <- factorial(v) * (S_hat_inv%*%C_p1_hat)[v+1]
    
    if (v > 0) {
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (n*h1))
    } else {
      temp_ii <- min( max(mean(data <= grid[j]), 1/n), 1 - 1/n )
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1] / (0.5*n^2) *
                                                   h1 * temp_ii * (1 - temp_ii)))
    }

    # now optimal bandwidth
    if (v > 0) {
      opt.f <- function(a) {
        a^(2*p+2-2*v) * (dgp_hat[j, 1]*const_hat[j, 1] + a * dgp_hat[j, 2]*const_hat[j, 2])^2 + const_hat[j, 3]^2 / a^(2*v - 1)
      }
    } else {
      opt.f <- function(a) {
        a^(2*p+2-2*v) * (dgp_hat[j, 1]*const_hat[j, 1] + a * dgp_hat[j, 2]*const_hat[j, 2])^2 + const_hat[j, 3]^2 / a
      }
    }
    
    h[j] <- optimize(opt.f, interval=c(.Machine$double.eps, max(data) - min(data)), maximum=FALSE)$minimum
  }

  for (j in 1:ng) {
    if (is.na(h[j])) {
      h[j] <- sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))]
    }
    if (regularize) {
      if (nLocalMin > 0) { h[j] <- max(h[j], sort(abs(data-grid[j]))[min(n, nLocalMin)]) }
      if (nUniqueMin > 0) { h[j] <- max(h[j], sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)]) }
    }
    h[j] <- min(h[j], max(abs(dataUnique-grid[j])))
  }
  
  if (stdVar) {
    h <- h * scale_temp
  }
  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates integrated MSE-optimal bandwidth.
#'
#' @param data Numeric vector, the data.
#' @param grid Numeric vector, the evaluation points.
#' @param p Integer, polynomial order.
#' @param v Integer, order of derivative.
#' @param kernel String, the kernel.
#' @param Cweights Numeric vector, the counterfactual weights.
#' @param Pweights Numeric vector, the survey sampling weights.
#' @param massPoints Boolean, whether whether point estimates and standard errors
#'   should be corrected if there are mass points in the data.
#' @param stdVar Boolean, whether the data should be standardized for bandwidth selection.
#' @param regularize Boolean, Whether the bandwidth should be regularized.
#' @param nLocalMin Nonnegative integer, minimum number of observations in each local neighborhood.
#' @param nUniqueMin Nonnegative integer, minimum number of unique observations in each local neighborhood.
#'
#' @return
#' \item{}{A single bandwidth.}
#'
#' @keywords internal
bw_IMSE  <- function(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin) {
  
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
  
  if (stdVar) {
    center_temp <- mean(data)
    scale_temp <- sd(data)
    data <- (data - center_temp) / scale_temp
    dataUnique <- (dataUnique - center_temp) / scale_temp
    grid <- (grid - center_temp) / scale_temp
  }
  
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
  
  # obtain preliminary bandwidth for estimating densities
  # this is used for constructing preasymptotic matrices
  h1  <- bw_IROT(data=data, grid=grid, p=2,   v=1,   kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+2+1, nUniqueMin=20+2+1)
  
  # obtain preliminary bandwidth for estimating F_p+1
  # this is used for constructing F_p+1
  hp1 <- bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+p+2+1, nUniqueMin=20+p+2+1)
  
  # obtain preliminary bandwidth for estimating F_p+2
  # this is used for constructing F_p+2
  hp2 <- bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=TRUE, stdVar=TRUE, regularize=TRUE, nLocalMin=20+p+3+1, nUniqueMin=20+p+3+1)
  
  dgp_hat <- matrix(NA, ncol=2, nrow=ng) # Fp+1 and Fp+2 with normalization constants
  const_hat <- matrix(NA, ncol=3, nrow=ng)
  h <- rep(NA, ng)
  
  for (j in 1:ng) {
    # estimate F_p+2
    index_temp <- abs(data-grid[j]) <= hp2
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / hp2
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:(p+3))))
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / hp2
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / hp2
    } else {
      Kh_temp <- 0.75 * (1 - (Xh_temp)^2) / hp2
    }
    Kh_temp   <- Pweights[index_temp] * Kh_temp
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    
    temp <- try(
      (solve(t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)) %*% t(Xh_p_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=Kh_temp))[p+3, 1] / (hp2^(p+2))
      , silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 2] <- temp
    
    # estimate F_p+1
    index_temp <- abs(data-grid[j]) <= hp1
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / hp1
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:(p+2))))
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / hp1
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / hp1
    } else {
      Kh_temp <- 0.75 * (1 - (Xh_temp)^2) / hp1
    }
    Kh_temp   <- Pweights[index_temp] * Kh_temp
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    
    temp <- try(
      (solve(t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)) %*% t(Xh_p_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=Kh_temp))[p+2, 1]  / (hp1^(p+1)),
      silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 1] <- temp
    
    # prepare for estimating matrices
    index_temp <- abs(data-grid[j]) <= h1
    Xh_temp <- matrix(data[index_temp] - grid[j], ncol=1) / h1
    if (kernel == "triangular") {
      Kh_temp <- (1 - abs(Xh_temp)) / h1
    } else if (kernel == "uniform") {
      Kh_temp <- 0.5 / h1
    } else {
      Kh_temp <- 0.75*(1 - (Xh_temp)^2) / h1
    }
    #Kh_temp   <- Pweights[index_temp] * Kh_temp
    
    # estimate Cp matrix
    if (p > 0) {
      C_p_hat <- matrix(apply(
        sweep(
          apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p_hat <- matrix(apply(
        sweep(
          matrix(apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))), nrow=1),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }
    
    # estimate Cp+1 matrix
    if (p > 0) {
      C_p1_hat <- matrix(apply(
        sweep(
          apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p1_hat <- matrix(apply(
        sweep(
          matrix(apply(Xh_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))), nrow=1),
          MARGIN=2, FUN="*", STATS=Pweights[index_temp] * Kh_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }
    
    # estimate S matirx
    Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
    if (p == 0) {
      Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
    }
    
    S_hat <- t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Pweights[index_temp] * Kh_temp) / n
    S_hat_inv <- try(solve(S_hat), silent=TRUE)
    if (is.character(S_hat_inv)) next
    
    # estimate G matrix
    if (v == 0) {
      G_hat <- t(Xh_p_temp) %*% sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Pweights[index_temp] * Kh_temp^2) / n
    } else {
      # for coding clarity, here we use the full sample and the influence function approach
      if (massPoints) {
        Y_temp    <- matrix(Fn[indexUnique], ncol=1)
        Xh_temp   <- matrix((data[indexUnique] - grid[j]), ncol=1) / h1
        Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
        if (p == 0) {
          Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
        }
        
        if (kernel == "triangular") {
          Kh_temp <- ((1 - abs(Xh_temp)) / h1) * index_temp[indexUnique]
        } else if (kernel == "uniform") {
          Kh_temp <- (0.5 / h1) * index_temp[indexUnique]
        } else {
          Kh_temp <- (0.75 * (1 - Xh_temp^2) / h1) * index_temp[indexUnique]
        }
        
        Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=PweightsUnique)
        
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n
        
        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)
        for (jj in 1:ncol(G)) {
          G[, jj] <- (rep((cumsum(Xh_p_Kh_Pweights_temp[nUnique:1, jj]) / n)[nUnique:1], times=freqUnique) - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }
        G_hat <- t(G) %*% G / n
      } else {
        Y_temp    <- matrix(Fn, ncol=1)
        Xh_temp   <- matrix((data - grid[j]), ncol=1) / h1
        Xh_p_temp <- t(apply(Xh_temp, MARGIN=1, FUN=function(x) x^(0:p)))
        if (p == 0) {
          Xh_p_temp <- matrix(Xh_p_temp, ncol=1)
        }
        
        if (kernel == "triangular") {
          Kh_temp <- ((1 - abs(Xh_temp)) / h1) * index_temp
        } else if (kernel == "uniform") {
          Kh_temp <- (0.5 / h1) * index_temp
        } else {
          Kh_temp <- (0.75 * (1 - Xh_temp^2) / h1) * index_temp
        }
        
        Xh_p_Kh_temp <- sweep(Xh_p_temp, MARGIN=1, FUN="*", STATS=Kh_temp)
        Xh_p_Kh_Pweights_temp <- sweep(Xh_p_Kh_temp, MARGIN=1, FUN="*", STATS=Pweights)
        
        F_Xh_p_Kh_temp <- t(Y_temp) %*% Xh_p_Kh_Pweights_temp / n
        
        G <- matrix(NA, ncol=ncol(Xh_p_Kh_Pweights_temp), nrow=n)
        for (jj in 1:ncol(G)) {
          G[, jj] <- ((cumsum(Xh_p_Kh_Pweights_temp[n:1, jj]) / n)[n:1] - F_Xh_p_Kh_temp[1, jj]) * weights_normal
        }
        G_hat <- t(G) %*% G / n
      }
    }
    
    # now get all constants
    const_hat[j, 1] <- factorial(v) * (S_hat_inv%*%C_p_hat)[v+1]
    const_hat[j, 2] <- factorial(v) * (S_hat_inv%*%C_p1_hat)[v+1]
    
    if (v > 0) {
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (n*h1))
    } else {
      temp_ii <- min(max(mean(data <= grid[j]), 1/n), 1 - 1/n)
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1] / (0.5*n^2) *
                                                   h1 * temp_ii * (1 - temp_ii)))
    }
  }
  
  # now optimal bandwidth
  na.index <- apply(cbind(dgp_hat, const_hat), MARGIN=1, FUN=function(x) any(is.na(x)))
  dgp_hat <- dgp_hat[!na.index, , drop=FALSE]
  const_hat <- const_hat[!na.index, , drop=FALSE]
  
  if (v > 0) {
    opt.f <- function(a) {
      a^(2*p+2-2*v) * sum((dgp_hat[, 1]*const_hat[, 1] + a * dgp_hat[, 2]*const_hat[, 2])^2) + sum(const_hat[, 3]^2) / a^(2*v - 1)
    }
  } else {
    opt.f <- function(a) {
      a^(2*p+2-2*v) * sum((dgp_hat[, 1]*const_hat[, 1] + a * dgp_hat[, 2]*const_hat[, 2])^2) + sum(const_hat[, 3]^2) / a
    }
  }
  
  h <- optimize(opt.f, interval=c(.Machine$double.eps, max(data) - min(data)), maximum=FALSE)$minimum
  if (is.na(h)) {
    h <- 0
    for (j in 1:ng) {
      h <- max(h, sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))])
    }
  }
  if (regularize) {
    for (j in 1:ng) {
      if (nLocalMin > 0) { h <- max(h, sort(abs(data-grid[j]))[min(n, nLocalMin)]) }
      if (nUniqueMin > 0) { h <- max(h, sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)]) }
    }
  }
  h <- min(h, max(abs(max(dataUnique)-min(grid)), abs(min(dataUnique)-max(grid))))
  
  if (stdVar) {
    h <- h * scale_temp
  }
  
  return(h)
}




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
