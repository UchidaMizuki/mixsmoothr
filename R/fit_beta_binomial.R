FLXMR_beta_binomial_mixture <- function(formula = . ~ .,
                                        offset = NULL) {
  out <- methods::new("FLXMR",
                      weighted = TRUE,
                      name = "beta_binomial_mixture",
                      formula = formula,
                      offset = offset)
  out@defineComponent <- function(para) {
    predict <- function(x) {
      para$alpha / (para$alpha + para$beta) * offset
    }
    logLik <- function(x, y) {
      extraDistr::dbbinom(y,
                          size = offset,
                          alpha = para$alpha,
                          beta = para$beta,
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
    minuslogl <- function(alpha = 1,
                          beta = 1) {
      -sum(extraDistr::dbbinom(y,
                               size = offset,
                               alpha = alpha,
                               beta = beta,
                               log = TRUE))
    }
    fitted <- stats4::mle(minuslogl,
                          lower = vec_rep(0, 2))
    out@defineComponent(para = list(alpha = fitted@coef[["alpha"]],
                                    beta = fitted@coef[["beta"]],
                                    df = 2))
  }
  out
}

flexmix_beta_binomial_mixture <- function(x, y, k, ...) {
  flexmix::flexmix(y ~ 1,
                   data = data_frame(y = y),
                   k = k,
                   model = FLXMR_beta_binomial_mixture(offset = x),
                   control = dots_list(minprior = 0, ...,
                                       .named = TRUE,
                                       .homonyms = "error"))
}
