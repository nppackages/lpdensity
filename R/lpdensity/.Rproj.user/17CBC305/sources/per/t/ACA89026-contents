################################################################################
#' Internal function.
#'
#' Calculates density and higher order derivatives for Gaussian models.
#'
#' @param x, Numeric, point of evaluation.
#' @param v, Numeric and nonnegative. The derivative, with 0 being cdf and 1 being pdf.'
#' @param mean,sd, the mean and standard deviation.
#'
#' @return
#' \item{}{A scalar.}
#'
#' @keywords internal
normal_dgps <- function(x, v, mean, sd) {
  if (v == 0) {
    return(pnorm(x, mean=mean, sd=sd))
  } else {
    temp <- expression(exp(-(x-mu)^2/(2*sd^2))/sqrt(2*pi))
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
#' @param p, Integer, polynomial order.
#' @param low,up, Numeric, between -1 and 1, region of integration.
#' @param kernel, String, the kernel.
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
#' @param p, Integer, polynomial order.
#' @param low,up, Numeric, between -1 and 1, region of integration.
#' @param kernel, String, the kernel.
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
#' @param k, Integer, extra order (usually p+1).
#' @param p, Integer, polynomial order.
#' @param low,up, Numeric, between -1 and 1, region of integration.
#' @param kernel, String, the kernel.
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
#' @param p, Integer, polynomial order.
#' @param low,up, Numeric, between -1 and 1, region of integration.
#' @param kernel, String, the kernel.
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
#' @param data, Numeric vector, the data.
#' @param grid, Numeric vector, the evaluation points.
#' @param p, Integer, polynomial order.
#' @param v, Integer, order of derivative.
#' @param kernel, String, the kernel.
#' @param Regularize, whether the bandwidth should be regularized.
#'
#' @return
#' \item{}{Bandwidth sequence.}
#'
#' @keywords internal
bw_ROT  <- function(data, grid, p, v, kernel, regularize) {

  n  <- length(data)
  ng <- length(grid)

  # estimate a normal reference model
  mean_hat <- mean(data)
  sd_hat   <- sd(data)

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

  # bias estimate, no rate added
  bias_dgp <- matrix(NA, ncol=2, nrow=ng)
  for (j in 1:ng) {
    bias_dgp[j, 1] <- eval(temp_3, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+1) * factorial(v)
    bias_dgp[j, 2] <- eval(temp_4, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+2) * factorial(v) +
      bias_dgp[j, 1] * eval(temp_2, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat))
  }

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
      sd_dgp[j, 1] <- factorial(v) * sqrt(
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
    if (regularize) {
      h[j] <- max(h[j], sort(abs(data-grid[j]))[20 + p + 1])
    }
  }
  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates integrated rule-of-thumb bandwidth
#'
#' @param data, Numeric vector, the data.
#' @param grid, Numeric vector, the evaluation points.
#' @param p, Integer, polynomial order.
#' @param v, Integer, order of derivative.
#' @param kernel, String, the kernel.
#' @param Regularize, whether the bandwidth should be regularized.
#'
#' @return
#' \item{}{A single bandwidth.}
#'
#' @keywords internal
bw_IROT <- function(data, grid, p, v, kernel, regularize) {

  n  <- length(data)
  ng <- length(grid)

  # estimate a normal reference model
  mean_hat <- mean(data)
  sd_hat   <- sd(data)

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

  # bias estimate, no rate added
  bias_dgp <- matrix(NA, ncol=2, nrow=ng)
  for (j in 1:ng) {
    bias_dgp[j, 1] <- eval(temp_3, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+1) * factorial(v)
    bias_dgp[j, 2] <- eval(temp_4, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / factorial(p+2) * factorial(v) +
      bias_dgp[j, 1] * eval(temp_2, list(x=grid[j], mu=mean_hat, sd=sd_hat)) / eval(temp_1, list(x=grid[j], mu=mean_hat, sd=sd_hat))
  }

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
      sd_dgp[j, 1] <- factorial(v) * sqrt(
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

  if (regularize) {
    for (j in 1:ng) {
      h <- max(h, sort(abs(data-grid[j]))[20 + p + 1])
    }
  }
  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates MSE-optimal bandwidths.
#'
#' @param data, Numeric vector, the data.
#' @param grid, Numeric vector, the evaluation points.
#' @param p, Integer, polynomial order.
#' @param v, Integer, order of derivative.
#' @param kernel, String, the kernel.
#' @param Cweights, Numeric vector, the counterfactual weights.
#' @param Pweights, Numeric vector, the survey sampling weights.
#' @param Regularize, whether the bandwidth should be regularized.
#'
#' @return
#' \item{}{Bandwidth sequence.}
#'
#' @keywords internal
bw_MSE  <- function(data, grid, p, v, kernel, Cweights, Pweights, regularize) {

  ii <- order(data)
  data <- data[ii]
  Cweights <- Cweights[ii]
  Pweights <- Pweights[ii]
  n    <- length(data)
  ng   <- length(grid)
  Fn   <- cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
  weights_normal <- Cweights * Pweights / sum(Cweights * Pweights) * n
  Cweights <- Cweights / sum(Cweights) * n
  Pweights <- Pweights / sum(Pweights) * n

  # obtain preliminary bandwidth for estimating densities
  # this is used for constructing preasymptotic matrices
  h1  <- bw_IROT(data=data, grid=grid, p=2, v=1, kernel=kernel, regularize=TRUE)

  # obtain preliminary bandwidth for estimating F_p+1
  # this is used for constructing F_p+1
  hp1 <- bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, regularize=TRUE)

  # obtain preliminary bandwidth for estimating F_p+2
  # this is used for constructing F_p+2
  hp2 <- bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, regularize=TRUE)

  dgp_hat <- matrix(NA, ncol=2, nrow=ng) # Fp+1 and Fp+2 with normalization constants
  const_hat <- matrix(NA, ncol=3, nrow=ng)
  h <- rep(NA, ng)

  for (j in 1:ng) {
    # estimate F_p+2
    index_temp <- abs(data-grid[j]) <= hp2
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) # centered
    if (kernel == "triangular") {
      K_temp <- 1 - abs(X_temp / hp2)
    } else if (kernel == "uniform") {
      K_temp <- 1
    } else {
      K_temp <- 1 - (X_temp / hp2)^2
    }
    K_temp   <- Pweights[index_temp] * K_temp
    X_temp   <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:(p+3))))
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    temp <- try(
      (solve(t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)) %*% t(X_temp) %*%
        sweep(Y_temp, MARGIN=1, FUN="*", STATS=K_temp))[p+3, 1]
      , silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 2] <- temp

    # estimate F_p+1
    index_temp <- abs(data-grid[j]) <= hp1
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) # centered
    if (kernel == "triangular") {
      K_temp <- 1 - abs(X_temp / hp1)
    } else if (kernel == "uniform") {
      K_temp <- 1
    } else {
      K_temp <- 1 - (X_temp / hp1)^2
    }
    K_temp   <- Pweights[index_temp] * K_temp
    X_temp   <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:(p+2))))
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    temp <- try(
      (solve(t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)) %*% t(X_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=K_temp))[p+2, 1],
      silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 1] <- temp

    # prepare for estimating matrices
    index_temp <- abs(data-grid[j]) <= h1
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) / h1 # centered and scaled
    if (kernel == "triangular") {
      K_temp <- (1 - abs(X_temp)) / h1
    } else if (kernel == "uniform") {
      K_temp <- 0.5 / h1
    } else {
      K_temp <- 0.75*(1 - (X_temp)^2) / h1
    }
    K_temp   <- Pweights[index_temp] * K_temp

    # estimate Cp matrix
    if (p > 0) {
      C_p_hat <- matrix(apply(
        sweep(
          apply(X_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p_hat <- matrix(apply(
        sweep(
          matrix(apply(X_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))), nrow=1),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }

    # estimate Cp+1 matrix
    if (p > 0) {
      C_p1_hat <- matrix(apply(
        sweep(
          apply(X_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p1_hat <- matrix(apply(
        sweep(
          matrix(apply(X_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))), nrow=1),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }

    # estimate S matirx
    X_temp <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:p)))
    if (p == 0) {
      X_temp <- t(X_temp)
    }
    XhKh_temp <- sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)
    S_hat <- t(X_temp) %*% XhKh_temp / n
    S_hat_inv <- try(solve(S_hat), silent=TRUE)
    if (is.character(S_hat_inv)) next

    if (v > 0) {
      # estimate G matrix
      F_XhKh_temp <- matrix(Fn[index_temp], nrow=1) %*% XhKh_temp / n
      G <- XhKh_temp[sum(index_temp):1, ]

      for (jj in 1:ncol(G)) {
        G[, jj] <- cumsum(G[, jj]) / n - F_XhKh_temp[1, jj]
      }

      G <- sweep(G, MARGIN=1, FUN="*", STATS=weights_normal[index_temp])
      G <- t(G) %*% G / n

      index_temp_1 <- data - grid[j] < -1 * h1
      index_temp_2 <- data - grid[j] >      h1

      G1 <- matrix(1 - Fn[index_temp], nrow=1)
      G1 <- G1 %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp) / n
      G1 <- t(G1) %*% G1 * sum(weights_normal[index_temp_1]^2) / n
      G2 <- matrix(0 - Fn[index_temp], nrow=1)
      G2 <- G2 %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp) / n
      G2 <- t(G2) %*% G2 * sum(weights_normal[index_temp_2]^2) / n
      G_hat <- G + G1 + G2

    } else {
      # estimate T matrix
      G_hat <- t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp^2) / n
    }

    # now get all constants
    const_hat[j, 1] <- factorial(v) * (S_hat_inv%*%C_p_hat)[v+1]
    const_hat[j, 2] <- factorial(v) * (S_hat_inv%*%C_p1_hat)[v+1]

    if (v > 0) {
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (n*h1))
    } else {
      temp_ii <- which.min(abs(data-grid[j]))
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (0.5*n^2) *
                                               h1 * Fn[temp_ii] * (1 - Fn[temp_ii]))
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
    if (regularize) {
      h[j] <- max(h[j], sort(abs(data-grid[j]))[20 + p + 1])
    }
  }

  return(h)
}

################################################################################
#' Internal function.
#'
#' Calculates integrated MSE-optimal bandwidth.
#'
#' @param data, Numeric vector, the data.
#' @param grid, Numeric vector, the evaluation points.
#' @param p, Integer, polynomial order.
#' @param v, Integer, order of derivative.
#' @param kernel, String, the kernel.
#' @param Cweights, Numeric vector, the counterfactual weights.
#' @param Pweights, Numeric vector, the survey sampling weights.
#' @param Regularize, whether the bandwidth should be regularized.
#'
#' @return
#' \item{}{A single bandwidth.}
#'
#' @keywords internal
bw_IMSE  <- function(data, grid, p, v, kernel, Cweights, Pweights, regularize) {

  ii <- order(data)
  data <- data[ii]
  Cweights <- Cweights[ii]
  Pweights <- Pweights[ii]
  n    <- length(data)
  ng   <- length(grid)
  Fn   <- cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
  weights_normal <- Cweights * Pweights / sum(Cweights * Pweights) * n
  Cweights <- Cweights / sum(Cweights) * n
  Pweights <- Pweights / sum(Pweights) * n

  # obtain preliminary bandwidth for estimating densities
  # this is used for constructing preasymptotic matrices
  h1  <- bw_IROT(data=data, grid=grid, p=2, v=1, kernel=kernel, regularize=TRUE)

  # obtain preliminary bandwidth for estimating F_p+1
  # this is used for constructing F_p+1
  hp1 <- bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, regularize=TRUE)

  # obtain preliminary bandwidth for estimating F_p+2
  # this is used for constructing F_p+2
  hp2 <- bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, regularize=TRUE)

  dgp_hat <- matrix(NA, ncol=2, nrow=ng) # Fp+1 and Fp+2 with normalization constants
  const_hat <- matrix(NA, ncol=3, nrow=ng)
  h <- rep(NA, ng)

  for (j in 1:ng) {
    # estimate F_p+2
    index_temp <- abs(data-grid[j]) <= hp2
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) # centered
    if (kernel == "triangular") {
      K_temp <- 1 - abs(X_temp / hp2)
    } else if (kernel == "uniform") {
      K_temp <- 1
    } else {
      K_temp <- 1 - (X_temp / hp2)^2
    }
    K_temp   <- Pweights[index_temp] * K_temp
    X_temp   <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:(p+3))))
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    temp <- try(
      (solve(t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)) %*% t(X_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=K_temp))[p+3, 1]
      , silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 2] <- temp

    # estimate F_p+1
    index_temp <- abs(data-grid[j]) <= hp1
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) # centered
    if (kernel == "triangular") {
      K_temp <- 1 - abs(X_temp / hp1)
    } else if (kernel == "uniform") {
      K_temp <- 1
    } else {
      K_temp <- 1 - (X_temp / hp1)^2
    }
    K_temp   <- Pweights[index_temp] * K_temp
    X_temp   <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:(p+2))))
    Y_temp   <- matrix(Fn[index_temp], ncol=1)
    temp <- try(
      (solve(t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)) %*% t(X_temp) %*%
         sweep(Y_temp, MARGIN=1, FUN="*", STATS=K_temp))[p+2, 1],
      silent=TRUE)
    if(is.character(temp)) next
    dgp_hat[j, 1] <- temp

    # prepare for estimating matrices
    index_temp <- abs(data-grid[j]) <= h1
    X_temp <- matrix(data[index_temp] - grid[j], ncol=1) / h1 # centered and scaled
    if (kernel == "triangular") {
      K_temp <- (1 - abs(X_temp)) / h1
    } else if (kernel == "uniform") {
      K_temp <- 0.5 / h1
    } else {
      K_temp <- 0.75*(1 - (X_temp)^2) / h1
    }
    K_temp   <- Pweights[index_temp] * K_temp

    # estimate Cp matrix
    if (p > 0) {
      C_p_hat <- matrix(apply(
        sweep(
          apply(X_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p_hat <- matrix(apply(
        sweep(
          matrix(apply(X_temp, MARGIN=1, FUN=function(x) x^((p+1):(2*p+1))), nrow=1),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }

    # estimate Cp+1 matrix
    if (p > 0) {
      C_p1_hat <- matrix(apply(
        sweep(
          apply(X_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    } else {
      C_p1_hat <- matrix(apply(
        sweep(
          matrix(apply(X_temp, MARGIN=1, FUN=function(x) x^((p+2):(2*p+2))), nrow=1),
          MARGIN=2, FUN="*", STATS=K_temp),
        MARGIN=1, FUN=sum) / n, ncol=1)
    }

    # estimate S matirx
    X_temp <- t(apply(X_temp, MARGIN=1, FUN=function(x) x^(0:p)))
    if (p == 0) {
      X_temp <- t(X_temp)
    }
    XhKh_temp <- sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp)
    S_hat <- t(X_temp) %*% XhKh_temp / n
    S_hat_inv <- try(solve(S_hat), silent=TRUE)
    if (is.character(S_hat_inv)) next

    if (v > 0) {
      # estimate G matrix
      F_XhKh_temp <- matrix(Fn[index_temp], nrow=1) %*% XhKh_temp / n
      G <- XhKh_temp[sum(index_temp):1, ]
      for (jj in 1:ncol(G)) {
        G[, jj] <- cumsum(G[, jj]) / n - F_XhKh_temp[1, jj]
      }

      #G <- diag(1, sum(index_temp)); G[col(G) > row(G)] <- 1
      #G <- G - matrix(Fn[index_temp], ncol=sum(index_temp), nrow=sum(index_temp), byrow=TRUE)
      #G <- G %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp) / n
      G <- sweep(G, MARGIN=1, FUN="*", STATS=weights_normal[index_temp])
      G <- t(G) %*% G / n

      index_temp_1 <- data - grid[j] < -1 * h1
      index_temp_2 <- data - grid[j] >      h1

      G1 <- matrix(1 - Fn[index_temp], nrow=1)
      G1 <- G1 %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp) / n
      G1 <- t(G1) %*% G1 * sum(weights_normal[index_temp_1]^2) / n
      G2 <- matrix(0 - Fn[index_temp], nrow=1)
      G2 <- G2 %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp) / n
      G2 <- t(G2) %*% G2 * sum(weights_normal[index_temp_2]^2) / n
      G_hat <- G + G1 + G2

    } else {
      # estimate T matrix
      G_hat <- t(X_temp) %*% sweep(X_temp, MARGIN=1, FUN="*", STATS=K_temp^2) / n
    }

    # now get all constants
    const_hat[j, 1] <- factorial(v) * (S_hat_inv%*%C_p_hat)[v+1]
    const_hat[j, 2] <- factorial(v) * (S_hat_inv%*%C_p1_hat)[v+1]
    if (v > 0) {
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (n*h1))
    } else {
      temp_ii <- which.min(abs(data-grid[j]))
      const_hat[j, 3] <- factorial(v) * sqrt(abs((S_hat_inv%*%G_hat%*%S_hat_inv)[v+1,v+1]) / (0.5*n^2) *
                                               h1 * Fn[temp_ii] * (1 - Fn[temp_ii]))
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

  if (regularize) {
    for (j in 1:ng) {
      h <- max(h, sort(abs(data-grid[j]))[20 + p + 1])
    }
  }

  return(h)
}



