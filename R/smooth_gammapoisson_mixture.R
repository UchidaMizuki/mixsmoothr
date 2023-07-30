FLXMR_gammapoisson_mixture <- function(formula = . ~ .,
                                       offset = NULL) {
  out <- methods::new("FLXMR",
                      weighted = TRUE,
                      name = "gammapoisson_mixture",
                      formula = formula,
                      offset = offset)
  out@defineComponent <- function(para) {
    predict <- function(x) {
      para$alpha / para$beta * offset
    }
    logLik <- function(x, y) {
      stats::dnbinom(y,
                     size = para$alpha,
                     mu = predict(x),
                     log = TRUE)
    }
    methods::new("FLXcomponent",
                 parameters = list(alpha = para$alpha,
                                   beta = para$beta),
                 logLik = logLik,
                 predict = predict,
                 df = para$df)
  }
  out@fit <- function(x, y, w) {
    fitted <- MASS::glm.nb(y ~ offset(log(offset)),
                           weights = w)
    out@defineComponent(para = list(alpha = fitted$theta,
                                    beta = fitted$theta / exp(stats::coef(fitted))[["(Intercept)"]],
                                    df = 2))
  }
  out
}

variance_gamma <- function(alpha, beta) {
  alpha * beta ^ 2
}

flexmix_gammapoisson_mixture <- function(k, x, y, ...) {
  fitted <- flexmix::flexmix(y ~ 1,
                             data = data_frame(y = y),
                             k = k,
                             model = FLXMR_gammapoisson_mixture(offset = x),
                             control = dots_list(minprior = 0, ...,
                                                 .named = TRUE,
                                                 .homonyms = "error"))
  parameters <- rbind(flexmix::parameters(fitted),
                      prior = flexmix::prior(fitted))
  loc_inlier <- which.min(variance_gamma(alpha = parameters["alpha", ],
                                         beta = parameters["beta", ]))
  colnames(parameters) <- vec_rep("outlier", ncol(parameters))
  colnames(parameters)[loc_inlier] <- "inlier"

  tibble::tibble(k = k,
                 parameters = list(parameters),
                 BIC = stats::BIC(fitted))
}

fit_gammapoisson_mixture <- function(x, y, k, ...) {
  fitted <- purrr::map(k,
                       flexmix_gammapoisson_mixture,
                       x = x,
                       y = y, ...)
  purrr::list_rbind(fitted)
}

smooth_gammapoisson_mixture <- function(fitted, x, y) {
  fitted <- vec_slice(fitted, which.min(fitted$BIC))
  parameters <- dplyr::first(fitted$parameters)
  alpha <- parameters["alpha", ,
                      drop = FALSE]
  beta <- parameters["beta", ,
                     drop = FALSE]
  prior <- parameters["prior", ,
                      drop = FALSE]

  responsibility <- purrr::map(seq_len(ncol(parameters)),
                               function(i) {
                                 prior[[i]] *
                                   stats::dnbinom(y,
                                                  size = alpha[[i]],
                                                  mu = alpha[[i]] / beta[[i]] * x)
                               })
  responsibility <- exec(cbind, !!!responsibility)
  type <- vec_slice(colnames(parameters),
                    apply(responsibility, 1, which.max))

  alpha_inlier <- alpha[, "inlier"]
  beta_inlier <- beta[, "inlier"]
  tibble::tibble(rate = dplyr::case_match(type,
                                          "inlier" ~ (y + alpha_inlier) / (x + beta_inlier),
                                          "outlier" ~ alpha_inlier / beta_inlier),
                 type = factor(type, c("inlier", "outlier")))
}
