#----------------------------------------------------------------------------------------------------
# lpdensity test file
# Across R and Stata
# Cattaneo, Jansson and Ma
# June 19, 2020
#----------------------------------------------------------------------------------------------------

# generate data to ensure consistency of results across R and Stata
data <- matrix(NA, ncol=8, nrow=2000)
set.seed(123456)
n <- 2000
data[, 1] <- rnorm(n, mean=1)
data[, 2] <- data[, 1] * 20
data[, 3] <- round(data[, 1])
data[, 4] <- round(data[, 2])
data[1:19, 5] <- quantile(data[, 1], seq(0.05, 0.95, 0.05))
data[1:19, 6] <- quantile(data[, 2], seq(0.05, 0.95, 0.05))
data[1:6, 7] <- data[13:18, 5] + 0.01
data[1:7, 8] <- 7:1
colnames(data) <- c("v1", "v2", "v3", "v4", "g1", "g2", "rgrid", "rindex")
write.csv(data, "/Users/xinweima/Dropbox/0_Research/01_2017_LocPolDensity/Package_Stata_lpdensity/lpdensity_data.csv", row.names=FALSE, na="")

#--------------------------------------------------------------------------------
# Density estimation check
# ROT bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot"))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot"))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot"))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot"))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", p=3))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", p=3))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", p=3))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", p=3))

# IROT bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot"))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot"))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot"))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot"))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", p=3))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", p=3))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", p=3))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", p=3))

# MSE bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi"))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi"))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi"))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi"))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", p=3))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", p=3))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", p=3))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", p=3))


# IMSE bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi"))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi"))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi"))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi"))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", p=3))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", p=3))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", p=3))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", p=3))

#--------------------------------------------------------------------------------
# DF estimation check
# ROT bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=0))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=0))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=0))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=0))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=2))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=2))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=2))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=2))

# IROT bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=0))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=0))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=0))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=0))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=2))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=2))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=2))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=2))

# MSE bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=0))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=2))


# IMSE bw check
summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=1))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=1))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=0))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=0))

summary(lpbwdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=2))
summary(lpbwdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=2))

#--------------------------------------------------------------------------------
# Density estimation check
# ROT bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot"))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot"))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot"))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot"))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", p=3))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", p=3))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", p=3))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", p=3))

# IROT bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot"))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot"))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot"))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot"))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", p=3))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", p=3))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", p=3))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", p=3))

# MSE bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi"))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi"))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi"))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi"))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", p=3))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", p=3))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", p=3))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", p=3))


# IMSE bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi"))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi"))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi"))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi"))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", p=3))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", p=3))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", p=3))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", p=3))

#--------------------------------------------------------------------------------
# DF estimation check
# ROT bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=0))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=0))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=0))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=0))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=2))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=2))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-rot", v=0, p=2))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-rot", v=0, p=2))

# IROT bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=0))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=0))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=0))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=0))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=2))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=2))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-rot", v=0, p=2))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-rot", v=0, p=2))

# MSE bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=0))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=0))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=0))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=0))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=2))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=2))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="mse-dpi", v=0, p=2))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="mse-dpi", v=0, p=2))


# IMSE bw check
summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=1))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=1))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=1))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=1))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=0))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=0))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=0))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=0))

summary(lpdensity(data[, 1], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=2))
summary(lpdensity(data[, 2], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=2))
summary(lpdensity(data[, 3], grid=data[1:19, 5], bwselect="imse-dpi", v=0, p=2))
summary(lpdensity(data[, 4], grid=data[1:19, 6], bwselect="imse-dpi", v=0, p=2))
