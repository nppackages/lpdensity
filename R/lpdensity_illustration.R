#--------------------------------------------------------------------------------
# lpdensity: Local Polynomial Density Estimation and Inference
# Matias D. Cattaneo, Michael Jansson, and Xinwei Ma
# Replication Code
#--------------------------------------------------------------------------------

rm(list=ls())

#--------------------------------------------------------------------------------
# Install and load the "lpdensity" and "ggplot2" packages
#--------------------------------------------------------------------------------

# install.packages("ggplot2")
# install.packages("lpdensity")

library("ggplot2")
library("lpdensity")

#----------------------------------------
# Generate data
#----------------------------------------

set.seed(42)
data <- as.data.frame(rnorm(4000, mean = -1))
data <- data[data < 0, 1, drop=FALSE]
data <- -1 * data[1:2000, 1, drop=FALSE]
colnames(data) <- c("v1")

#----------------------------------------
# Figure 1, panel (a)
#----------------------------------------

data$pdf <- dnorm(data$v1, mean = 1, sd = 1) / pnorm(0, mean = 1, sd = 1, lower.tail = FALSE)
ggplot() + geom_histogram(data=data, aes(x=v1, y=..density..), breaks=seq(0, 4, 0.2), fill=2, col="white", alpha=0.6) +
  theme_bw() + labs(x = "") +
  geom_line(data=data, aes(x=v1, y=pdf))

#----------------------------------------
# Figure 1, panel (b)
#----------------------------------------

model2 <- lpdensity(data$v1, bw = 0.5, grid = seq(0, 4, 0.05))
plot(model2, ylabel = "density") + theme(legend.position = "none")

#----------------------------------------
# lpdensity(): Estimation with bandwidth 0.5 on provided grid points
#----------------------------------------

model1 <- lpdensity(data$v1, bw = 0.5, grid = seq(0, 4, 0.5))
summary(model1)

#----------------------------------------
# lpdensity(): extracting estimation
#   results
#----------------------------------------

model1$Estimate

#----------------------------------------
# lpdensity(): conventional inference
#----------------------------------------

summary(lpdensity(data$v1, bw = 0.5, p = 2, q = 2))

#----------------------------------------
# lpdensity(): customizing screen output
#----------------------------------------

set.seed(123) # fix the random seed for critical value simulation
summary(model1, alpha = 0.01, sep = 3, grid = c(0, 0.5, 1, 2), CIuniform = TRUE)

#----------------------------------------
# lpdensity(): inconsistent density
#   estimation using partial sample
#----------------------------------------

lpdensity(data$v1[data$v1 < 1.5], bw = 0.5, grid = 1.5)$Estimate[, "f_p"]
lpdensity(data$v1[data$v1 > 1.5], bw = 0.5, grid = 1.5)$Estimate[, "f_p"]
dnorm(1.5, mean = 1, sd = 1) / pnorm(0, mean = 1, sd = 1, lower.tail = FALSE) # true density at 1.5

#----------------------------------------
# lpdensity(): consistent density
#   estimation using partial sample and
#   option "scale"
#----------------------------------------

lpdensity(data$v1[data$v1 < 1.5], bw = 0.5, grid = 1.5,
          scale = sum(data$v1 < 1.5)/2000)$Estimate[, "f_p"]
lpdensity(data$v1[data$v1 > 1.5], bw = 0.5, grid = 1.5,
          scale = sum(data$v1 > 1.5)/2000)$Estimate[, "f_p"]

#----------------------------------------
# plot(): customization
#----------------------------------------

plot(model2, CItype="line", ylabel = "density") +
  theme(legend.position = "none")
plot(model2, type="points", CItype="ebar", grid = seq(0, 4, 0.5), ylabel = "density") +
  theme(legend.position = "none")
plot(model2, hist = TRUE, histData = data$v1, histBreaks = seq(0, 4, 0.2), ylabel = "density") +
  theme(legend.position = "none")
set.seed(123) # fix the random seed for critical value simulation
plot(model2, alpha=0.1, CIuniform = TRUE, ylabel = "density") +
  theme(legend.position = "none")

#----------------------------------------
# lpbwdensity(): illustration
#----------------------------------------

model1bw <- lpbwdensity(data$v1, grid = seq(0, 4, 0.5))
summary(model1bw)

#----------------------------------------
# lpdensity(): automatic bandwidth
#   selection
#----------------------------------------

model5 <- lpdensity(data$v1, grid = seq(0, 4, 0.5), bwselect = "imse-dpi")
summary(model5)

#----------------------------------------
# lpdensity(): undersmoothing
#----------------------------------------

# Estimation and plot using IMSE bandwidth
model6bwIMSE <- lpbwdensity(data$v1, grid = seq(0, 4, 0.05),
                            bwselect = "imse-dpi")
model6 <- lpdensity(data$v1, grid = seq(0, 4, 0.05),
                    bw = model6bwIMSE$BW[, "bw"])
plot(model6, ylabel = "density") +
  theme(legend.position = "none")

# Estimation and plot using half IMSE bandwidth
model7 <- lpdensity(data$v1, grid = seq(0, 4, 0.05),
                    bw = model6bwIMSE$BW[, "bw"] / 2)
plot(model7, ylabel = "density") +
  theme(legend.position = "none")
