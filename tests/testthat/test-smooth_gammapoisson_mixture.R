test_that("smooth_gammapoisson_mixture", {
  set.seed(1234)

  n <- 5e3
  x <- 1 + rpois(n, 1e2)

  test_smooth_gammapoisson_mixture_1 <- function(alpha, beta) {
    y <- rnbinom(n,
                 size = alpha,
                 mu = alpha / beta * x)

    fitted <- smooth_gammapoisson_mixture(y = y,
                                          x = x,
                                          k = 1,
                                          control = control_smooth(verbose = TRUE))
    parameters <- fitted$parameters[[1]]
    print(parameters)
    expect_equal(parameters["alpha", ] / alpha, 1,
                 tolerance = 1e-1)
    expect_equal(parameters["beta", ] / beta, 1,
                 tolerance = 1e-1)
  }
  test_smooth_gammapoisson_mixture_1(3, 5)
  test_smooth_gammapoisson_mixture_1(1, 10)
  test_smooth_gammapoisson_mixture_1(10, 1)

  test_smooth_gammapoisson_mixture <- function(prior_inlier,
                                               alpha_inlier, beta_inlier,
                                               alpha_outlier, beta_outlier) {
    stopifnot(
      var_gamma(alpha_inlier, beta_inlier) < var_gamma(alpha_outlier, beta_outlier)
    )
    n_inlier <- round(n * prior_inlier)
    prob <- c(rgamma(n_inlier, alpha_inlier, beta_inlier),
              rgamma(n - n_inlier, alpha_outlier, beta_outlier))
    y <- rpois(n,
               lambda = prob * x)

    fitted <- smooth_gammapoisson_mixture(y = y,
                                          x = x,
                                          k = 2,
                                          control = control_smooth(verbose = TRUE))

    parameters <- fitted$parameters[[1]]
    print(parameters)
    expect_equal(parameters["prior", "inlier"] / prior_inlier, 1,
                 tolerance = 1e-1)
    expect_equal(parameters["alpha", "inlier"] / alpha_inlier, 1,
                 tolerance = 1e-1)
    expect_equal(parameters["beta", "inlier"] / beta_inlier, 1,
                 tolerance = 1e-1)
    expect_equal(parameters["alpha", "outlier"] / alpha_outlier, 1,
                 tolerance = 1e-1)
    expect_equal(parameters["beta", "outlier"] / beta_outlier, 1,
                 tolerance = 1e-1)
  }
  test_smooth_gammapoisson_mixture(prior_inlier = 0.25,
                                   alpha_inlier = 3,
                                   beta_inlier = 5,
                                   alpha_outlier = 10,
                                   beta_outlier = 1)
  test_smooth_gammapoisson_mixture(prior_inlier = 0.75,
                                   alpha_inlier = 3,
                                   beta_inlier = 10,
                                   alpha_outlier = 10,
                                   beta_outlier = 3)
})
