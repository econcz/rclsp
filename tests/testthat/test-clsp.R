test_that("CMLS model (RP-type) runs and converges", {
  skip_if_not_installed("CVXR")
  set.seed(123456789)
  
  k <- 500   # number of observations in D
  p <- 6     # number of regressors
  c <- 1.0   # sum of coefficients
  
  # design matrix D
  D        <- matrix(NA_real_, nrow = k, ncol = p)
  D[, 1]   <- 1.0
  D[, 2:p] <- matrix(rnorm(k * (p - 1)), nrow = k)
  
  # true coefficients and dependent variable
  b_true <- rnorm(p)
  b_true <- (b_true / sum(b_true)) * c
  e      <- rnorm(k)
  y      <- as.numeric(D %*% b_true + e)
  
  # constraint matrices
  b <- rbind(
    matrix(c, ncol = 1),
    matrix(0, ncol = 1, nrow = k - 2),
    matrix(0, ncol = 1, nrow = k - 1),
    matrix(y, ncol = 1)
  )
  
  C <- rbind(
    matrix(1, nrow = 1, ncol = p),
    apply(D, 2, diff, differences = 2),
    apply(D, 2, diff, differences = 1)
  )
  
  S <- rbind(
    matrix(0, nrow = 1, ncol = k - 2),
    diag(sign(diff(y, differences = 2))),
    matrix(0, nrow = k - 1, ncol = k - 2)
  )
  
  # model estimation
  model <- clsp(problem = "cmls", b = b, C = C, S = S, M = D,
                r = 1L, alpha = 1.0)
  
  # diagnostics
  expect_true(is.matrix(model$x))
  expect_true(all(is.finite(model$x)))
  expect_true(is.numeric(model$nrmse))
  
  # print summary for developer feedback
  summary(model)
  
  # test bootstrap t-test (NRMSE)
  ttest_res <- ttest(model, sample_size = 30L,
                     seed = 123456789L,
                     distribution = rnorm,
                     partial = TRUE)
  
  expect_true(is.list(ttest_res))
  expect_true("p_value" %in% names(ttest_res) || length(ttest_res) > 0)
})
