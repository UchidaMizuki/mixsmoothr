test_that("fit_gammapoisson_mixture", {
  set.seed(1234)

  n <- 1e3
  x <- rpois(n, n)

  test_fit_gammapoisson_mixture_1 <- function(alpha, beta) {
    y <- rnbinom(n,
                 size = alpha,
                 mu = alpha / beta * n)

    fitted <- fit_gammapoisson_mixture(x, y, 1,
                                       verbose = 1)
    expect_equal(fitted$parameters[[1]]["alpha", ], alpha,
                 tolerance = 1e-1)
    expect_equal(fitted$parameters[[1]]["beta", ], beta,
                 tolerance = 1e-1)
  }
  test_fit_gammapoisson_mixture_1(3, 5)
  test_fit_gammapoisson_mixture_1(1, 10)
  test_fit_gammapoisson_mixture_1(10, 1)

  test_fit_gammapoisson_mixture <- function(prior_inlier,
                                            alpha_inlier, beta_inlier,
                                            alpha_outlier, beta_outlier) {
    stopifnot(
      variance_gamma(alpha_inlier, beta_inlier) < variance_gamma(alpha_outlier, beta_outlier)
    )
    n_inlier <- round(n * prior_inlier)
    prob <- c(rgamma(n_inlier, alpha_inlier, beta_inlier),
              rgamma(n - n_inlier, alpha_outlier, beta_outlier))
    y <- rpois(n,
               lambda = prob * x)

    fitted <- fit_gammapoisson_mixture(x, y, 2,
                                       verbose = 1)

    parameters <- fitted$parameters[[1]]
    expect_equal(parameters["prior", "inlier"] / prior_inlier, 1,
                 tolerance = 5e-1)
    expect_equal(parameters["alpha", "inlier"] / alpha_inlier, 1,
                 tolerance = 5e-1)
    expect_equal(parameters["beta", "inlier"] / beta_inlier, 1,
                 tolerance = 5e-1)
    expect_equal(parameters["alpha", "outlier"] / alpha_outlier, 1,
                 tolerance = 5e-1)
    expect_equal(parameters["beta", "outlier"] / beta_outlier, 1,
                 tolerance = 5e-1)
  }
  test_fit_gammapoisson_mixture(prior_inlier = 0.25,
                                alpha_inlier = 10,
                                beta_inlier = 1,
                                alpha_outlier = 3,
                                beta_outlier = 5)
})
