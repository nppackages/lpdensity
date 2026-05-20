baseline_data <- function() {
  index <- seq(0, 199)
  seq(0.05, 3.95, length.out = 200) + 0.05 * sin(index)
}

baseline_weights <- function() {
  index <- seq(1, 200)
  list(
    Cweights = 1 + 0.1 * (index %% 5),
    Pweights = 1 + 0.05 * (index %% 7)
  )
}

expect_baseline_case <- function(observed, expected, case_name, tolerance = 1e-10) {
  case_expected <- expected[expected$case == case_name, colnames(observed)]
  row.names(case_expected) <- NULL
  expect_equal(
    as.data.frame(observed),
    case_expected,
    tolerance = tolerance,
    ignore_attr = TRUE,
    info = case_name
  )
}

test_that("fixed bandwidth lpdensity estimates match the numerical baseline", {
  data <- baseline_data()
  observed <- lpdensity(data, bw = 0.5, grid = seq(0, 4, 0.5))$Estimate
  expected <- read.csv(test_path("fixtures", "numerical-baseline.csv"))

  expect_equal(as.data.frame(observed), expected, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("broad lpdensity cases match numerical baselines", {
  data <- baseline_data()
  weights <- baseline_weights()
  expected <- read.csv(test_path("fixtures", "lpdensity-baseline-broad.csv"), check.names = FALSE)

  cases <- list(
    fixed_triangular_density = lpdensity(data, bw = 0.5, grid = seq(0, 4, 0.5))$Estimate,
    kernel_uniform_density = lpdensity(data, bw = 0.5, grid = seq(0, 4, 1), kernel = "uniform")$Estimate,
    kernel_epanechnikov_density = lpdensity(data, bw = 0.5, grid = seq(0, 4, 1), kernel = "epanechnikov")$Estimate,
    cdf_triangular = lpdensity(data, bw = 0.75, grid = seq(0, 4, 1), p = 1, q = 2, v = 0)$Estimate,
    weighted_triangular_density = lpdensity(
      data,
      bw = 0.65,
      grid = seq(0, 4, 1),
      Cweights = weights$Cweights,
      Pweights = weights$Pweights
    )$Estimate
  )

  for (case_name in names(cases)) {
    expect_baseline_case(cases[[case_name]], expected, case_name)
  }
})

test_that("broad lpbwdensity cases match numerical baselines", {
  data <- baseline_data()
  weights <- baseline_weights()
  expected <- read.csv(test_path("fixtures", "lpbwdensity-baseline-broad.csv"), check.names = FALSE)

  for (selector in c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot")) {
    case_name <- paste0("bw_", selector)
    observed <- as.data.frame(lpbwdensity(data, grid = seq(0, 4, 1), bwselect = selector)$BW)
    expect_baseline_case(observed, expected, case_name, tolerance = 1e-6)

    case_name <- paste0("weighted_bw_", selector)
    observed <- as.data.frame(lpbwdensity(
      data,
      grid = seq(0, 4, 1),
      bwselect = selector,
      Cweights = weights$Cweights,
      Pweights = weights$Pweights
    )$BW)
    expect_baseline_case(observed, expected, case_name, tolerance = 1e-6)
  }
})

test_that("uniform confidence band methods accept explicit simulation counts", {
  data <- baseline_data()
  fit <- lpdensity(data, bw = 0.5, grid = seq(0, 4, 0.5))

  set.seed(42)
  intervals <- confint(fit, CIuniform = TRUE, CIsimul = 10)
  expect_equal(dim(intervals), c(9, 5))
  expect_true(all(is.finite(intervals)))

  set.seed(42)
  expect_output(summary(fit, CIuniform = TRUE, CIsimul = 10), "Call: lpdensity")
})
