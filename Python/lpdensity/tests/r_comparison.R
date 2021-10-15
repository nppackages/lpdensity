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
      if (nLocalMin > 0) { 
        h <- max(h, sort(abs(data-grid[j]))[min(n, nLocalMin)]) 
        }
      if (nUniqueMin > 0) { 
        print(sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)])
        h <- max(h, sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)]) 
        }
    }
  }
  h <- min(h, max(abs(max(dataUnique)-min(grid)), abs(min(dataUnique)-max(grid))))
  print(h)
  if (stdVar) {
    h <- h * scale_temp
  }

  return(h)
}
