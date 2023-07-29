FLXMR_gamma_poisson_mixture <- function(formula = . ~ .,
                                        offset = NULL) {
  out <- methods::new("FLXMR",
                      weighted = TRUE,
                      name = "gamma_poisson_mixture",
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

flexmix_gamma_poisson_mixture <- function(x, y, k, ...) {
  flexmix::flexmix(y ~ 1,
                   data = data_frame(y = y),
                   k = k,
                   model = FLXMR_gamma_poisson_mixture(offset = x),
                   control = dots_list(minprior = 0, ...,
                                       .named = TRUE,
                                       .homonyms = "error"))
}
