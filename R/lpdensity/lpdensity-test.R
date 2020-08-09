#----------------------------------------------------------------------------------------------------
# lpdensity test file
# Cattaneo, Jansson and Ma
# June 19, 2020
#----------------------------------------------------------------------------------------------------

set.seed(42)
n <- 1000
data1 <- rnorm(n, mean=1)
data2 <- data1 * 1000
data3 <- round(data1)
data4 <- round(data2)
data5 <- lpdensityUnique(sort(data4))$unique
data5freq <- lpdensityUnique(sort(data4))$freq

lpdensityUnique <- function(x) {
  # x has one or no element
  if (length(x) == 0) return(list(unique=NULL, freq=c(), index=c()))
  if (length(x) == 1) return(list(unique=x, freq=1, index=1))

  n <- length(x)
  uniqueIndex <- c(x[2:n] != x[1:(n-1)], TRUE)
  unique <- x[uniqueIndex]
  nUnique <- sum(uniqueIndex)

  # all are distinct
  if (nUnique == n) return(list(unique=unique, freq=rep(1,length(x)), index=1:n))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=n, index=n))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, index=(1:n)[uniqueIndex]))
}

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, cdf
#----------------------------------------------------------------------------------------------------

### one grid point, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=0, v=0); summary(bw4)

### one grid point, cdf, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=1, v=0); summary(bw4)

### one grid point, cdf, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=2, v=0); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, pdf
#----------------------------------------------------------------------------------------------------

### one grid point, v=1, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=1, v=1); summary(bw4)

### one grid point, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=2, v=1); summary(bw4)

### one grid point, v=1, p=3
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-rot", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-rot", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-rot", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-rot", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-rot", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-rot", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-rot", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-rot", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="mse-dpi", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="mse-dpi", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="mse-dpi", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="mse-dpi", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, grid=0, bwselect="imse-dpi", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, grid=0, bwselect="imse-dpi", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, grid=0, bwselect="imse-dpi", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, grid=0, bwselect="imse-dpi", p=3, v=1); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, cdf
#----------------------------------------------------------------------------------------------------

### quantile grid points, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=0, v=0); summary(bw4)

### quantile grid points, cdf, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=1, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=1, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=1, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=1, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=1, v=0); summary(bw4)

### quantile grid points, cdf, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=2, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=2, v=0); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=2, v=0); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=2, v=0); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=2, v=0); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, pdf
#----------------------------------------------------------------------------------------------------

### ### quantile grid points, v=1, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=1, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=1, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=1, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=1, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=1, v=1); summary(bw4)

### quantile grid points, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=2, v=1); summary(bw4)

### quantile grid points, v=1, p=3
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-rot", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-rot", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-rot", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-rot", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-rot", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-rot", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-rot", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-rot", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="mse-dpi", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="mse-dpi", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="mse-dpi", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="mse-dpi", p=3, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(data1, bwselect="imse-dpi", p=3, v=1); summary(bw1)
bw2 <- lpbwdensity(data2, bwselect="imse-dpi", p=3, v=1); summary(bw2)
bw3 <- lpbwdensity(data3, bwselect="imse-dpi", p=3, v=1); summary(bw3)
bw4 <- lpbwdensity(data4, bwselect="imse-dpi", p=3, v=1); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, cdf, mass points
#----------------------------------------------------------------------------------------------------

### bandwidth selection with mass points, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="mse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="mse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="mse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="mse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="imse-rot", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="imse-rot", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="imse-rot", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="imse-rot", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="mse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="mse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="mse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="mse-dpi", p=0, v=0); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="imse-dpi", p=0, v=0); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="imse-dpi", p=0, v=0); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="imse-dpi", p=0, v=0); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="imse-dpi", p=0, v=0); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, pdf, mass points
#----------------------------------------------------------------------------------------------------

### bandwidth selection with mass points, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="mse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="mse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="mse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="mse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="imse-rot", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="imse-rot", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="imse-rot", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="imse-rot", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="mse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="mse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="mse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="mse-dpi", p=2, v=1); summary(bw4)

cat("\n\n\n\n\n\n\n\n\n\n")
bw1 <- lpbwdensity(round(data1    ), bwselect="imse-dpi", p=2, v=1); summary(bw1)
bw2 <- lpbwdensity(round(data1*10 ), bwselect="imse-dpi", p=2, v=1); summary(bw2)
bw3 <- lpbwdensity(round(data1*50 ), bwselect="imse-dpi", p=2, v=1); summary(bw3)
bw4 <- lpbwdensity(round(data1*500), bwselect="imse-dpi", p=2, v=1); summary(bw4)

#----------------------------------------------------------------------------------------------------
# lpbwdensity test, reporting
#----------------------------------------------------------------------------------------------------

### reporting
cat("\n\n\n\n\n\n\n\n\n\n")
coef(bw4)
print(bw4)
summary(bw4)
summary(bw4, grid=400)
summary(bw4, gridIndex=10)


#----------------------------------------------------------------------------------------------------
# lpdensity test, cdf
#----------------------------------------------------------------------------------------------------

### one grid point, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=0, v=0); summary(est4)

### one grid point, cdf, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=1, v=0); summary(est4)

### one grid point, cdf, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=2, v=0); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, pdf
#----------------------------------------------------------------------------------------------------

### one grid point, v=1, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=1, v=1); summary(est4)

### one grid point, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=2, v=1); summary(est4)

### one grid point, v=1, p=3
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-rot", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-rot", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-rot", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-rot", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-rot", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-rot", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-rot", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-rot", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="mse-dpi", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="mse-dpi", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="mse-dpi", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="mse-dpi", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, grid=0, bwselect="imse-dpi", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, grid=0, bwselect="imse-dpi", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, grid=0, bwselect="imse-dpi", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, grid=0, bwselect="imse-dpi", p=3, v=1); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, cdf
#----------------------------------------------------------------------------------------------------

### quantile grid points, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=0, v=0); summary(est4)

### quantile grid points, cdf, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=1, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=1, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=1, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=1, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=1, v=0); summary(est4)

### quantile grid points, cdf, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=2, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=0); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=2, v=0); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=2, v=0); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=2, v=0); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, pdf
#----------------------------------------------------------------------------------------------------

### ### quantile grid points, v=1, p=1
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=1, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=1, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=1, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=1, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=1, v=1); summary(est4)

### quantile grid points, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=2, v=1); summary(est4)

### quantile grid points, v=1, p=3
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-rot", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-rot", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-rot", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-rot", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-rot", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-rot", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-rot", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-rot", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="mse-dpi", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="mse-dpi", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="mse-dpi", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="mse-dpi", p=3, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(data1, bwselect="imse-dpi", p=3, v=1); summary(est1)
est2 <- lpdensity(data2, bwselect="imse-dpi", p=3, v=1); summary(est2)
est3 <- lpdensity(data3, bwselect="imse-dpi", p=3, v=1); summary(est3)
est4 <- lpdensity(data4, bwselect="imse-dpi", p=3, v=1); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, cdf, mass points
#----------------------------------------------------------------------------------------------------

### bandwidth selection with mass points, cdf, p=0
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="mse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="mse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="mse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="mse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="imse-rot", p=0, v=0); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="imse-rot", p=0, v=0); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="imse-rot", p=0, v=0); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="imse-rot", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="mse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="mse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="mse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="mse-dpi", p=0, v=0); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="imse-dpi", p=0, v=0); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="imse-dpi", p=0, v=0); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="imse-dpi", p=0, v=0); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="imse-dpi", p=0, v=0); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, pdf, mass points
#----------------------------------------------------------------------------------------------------

### bandwidth selection with mass points, v=1, p=2
cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="mse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="mse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="mse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="mse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="imse-rot", p=2, v=1); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="imse-rot", p=2, v=1); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="imse-rot", p=2, v=1); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="imse-rot", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="mse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="mse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="mse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="mse-dpi", p=2, v=1); summary(est4)

cat("\n\n\n\n\n\n\n\n\n\n")
est1 <- lpdensity(round(data1    ), bwselect="imse-dpi", p=2, v=1); summary(est1)
est2 <- lpdensity(round(data1*10 ), bwselect="imse-dpi", p=2, v=1); summary(est2)
est3 <- lpdensity(round(data1*50 ), bwselect="imse-dpi", p=2, v=1); summary(est3)
est4 <- lpdensity(round(data1*500), bwselect="imse-dpi", p=2, v=1); summary(est4)

#----------------------------------------------------------------------------------------------------
# lpdensity test, reporting
#----------------------------------------------------------------------------------------------------

### reporting
cat("\n\n\n\n\n\n\n\n\n\n")
est2 <- lpdensity(round(data1*10 ), bwselect="imse-dpi", p=2, v=1)
coef(est2)
print(est2)
vcov(est2)

summary(est2)
summary(est2, alpha=0.01)
summary(est2, alpha=0.01, CIuniform=TRUE)
summary(est2, alpha=0.01, CIuniform=TRUE, grid=c(20, 340))
summary(est2, alpha=0.01, CIuniform=TRUE, gridIndex=c(1, 5))
summary(est2, alpha=0.01, CIuniform=TRUE, grid=c(20, 340), gridIndex=c(1, 5))
summary(est2, grid=400)
summary(est2, gridIndex=10)

confint(est2, alpha=0.01, CIuniform=TRUE, grid=c(20))
confint(est2, alpha=0.01, CIuniform=TRUE, parm=c(20))

est2 <- lpdensity(round(data1*10 ), bwselect="imse-dpi", p=2, v=1, q=2)
coef(est2)
print(est2)
vcov(est2)
summary(est2)
summary(est2, alpha=0.01)
summary(est2, alpha=0.01, CIuniform=TRUE)
summary(est2, alpha=0.01, CIuniform=TRUE, grid=c(20, 340))
summary(est2, alpha=0.01, CIuniform=TRUE, gridIndex=c(1, 5))
summary(est2, alpha=0.01, CIuniform=TRUE, grid=c(20, 340), gridIndex=c(1, 5))
summary(est2, grid=400)
summary(est2, gridIndex=10)

confint(est2, alpha=0.01, CIuniform=TRUE, grid=c(20))
confint(est2, alpha=0.01, CIuniform=TRUE, parm=c(20))

#----------------------------------------------------------------------------------------------------
# lpdensity test, plotting
#----------------------------------------------------------------------------------------------------
est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=1)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=1)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

est1 <- lpdensity(data1*10000, bwselect="imse-dpi", p=2, v=1)
est2 <- lpdensity(data1*10000 + 10000, bwselect="imse-dpi", p=2, v=1)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1*10000)
plot(est1, est2, hist=TRUE, histData=data1*10000, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1*10000, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1*10000 + 10000)

est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=1, q=2)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=1, q=2)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)


est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=0, q=2)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=0, q=2)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=0, q=2)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=0, q=2)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

est1 <- lpdensity(data1, bwselect="imse-dpi", p=0, v=0, q=0)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=0, v=0, q=0)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=2, q=3)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=2, q=3)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

est1 <- lpdensity(data1, bwselect="imse-dpi", p=2, v=2, q=2)
est2 <- lpdensity(data1 + 1, bwselect="imse-dpi", p=2, v=2, q=2)
plot(est1, est2)
plot(est1, est2, hist=TRUE, histData=data1)
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue")
plot(est1, est2, hist=TRUE, histData=data1, histBreaks=seq(-1, 1, 0.2), histFillCol="red", histLineCol="blue", CIuniform=TRUE, grid=c(0, 1, 2, 2.6, 3.1), type="points", CItype="all")
plot(est1, est2, hist=TRUE, histData=data1 + 1)

#----------------------------------------------------------------------------------------------------
# lpdensity test, running speed
#----------------------------------------------------------------------------------------------------

n <- c(1000, 10000, 100000, 1000000, 10000000)
timeUsed1 <- rep(NA, length(n))
timeUsed2 <- rep(NA, length(n))

set.seed(42)
for (i in 1:length(n)) {
  data <- rnorm(n[i], sd=10)

  ptm <- proc.time()
  est <- lpdensity(data, bwselect="imse-dpi", p=2, v=1)
  timeUsed1[i] <- proc.time() - ptm

  ptm <- proc.time()
  est <- lpdensity(round(data), bwselect="imse-dpi", p=2, v=1)
  timeUsed2[i] <- proc.time() - ptm
}




