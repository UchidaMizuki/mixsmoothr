FLXMR_betabinomial_mixture <- function(formula = . ~ .,
                                       offset = NULL) {
  out <- methods::new("FLXMR",
                      weighted = TRUE,
                      name = "betabinomial_mixture",
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
      -sum(w *
             extraDistr::dbbinom(y,
                                 size = offset,
                                 alpha = alpha,
                                 beta = beta,
                                 log = TRUE))
    }
    fitted <- stats4::mle(minuslogl,
                          lower = vec_rep(.Machine$double.eps, 2))
    out@defineComponent(para = list(alpha = fitted@coef[["alpha"]],
                                    beta = fitted@coef[["beta"]],
                                    df = 2))
  }
  out
}

variance_beta <- function(alpha, beta) {
  alpha * beta / (alpha + beta) ^ 2 / (alpha + beta + 1)
}

flexmix_betabinomial_mixture <- function(k, x, y, ...) {
  fitted <- flexmix::flexmix(y ~ 1,
                             data = data_frame(y = y),
                             k = k,
                             model = FLXMR_betabinomial_mixture(offset = x),
                             control = dots_list(minprior = 0, ...,
                                                 .named = TRUE,
                                                 .homonyms = "error"))
  parameters <- rbind(flexmix::parameters(fitted),
                      prior = flexmix::prior(fitted))
  loc_inlier <- which.min(variance_beta(alpha = parameters["alpha", ],
                                        beta = parameters["beta", ]))
  colnames(parameters) <- vec_rep("outlier", ncol(parameters))
  colnames(parameters)[loc_inlier] <- "inlier"

  tibble::tibble(k = k,
                 parameters = list(parameters),
                 BIC = stats::BIC(fitted))
}

#' @export
fit_betabinomial_mixture <- function(x, y, k, ...) {
  fitted <- purrr::map(k,
                       flexmix_betabinomial_mixture,
                       x = x,
                       y = y, ...)
  purrr::list_rbind(fitted)
}

#' @export
smooth_betabinomial_mixture <- function(fitted, x, y) {
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
                                   extraDistr::dbbinom(y,
                                                       size = x,
                                                       alpha = alpha[[i]],
                                                       beta = beta[[i]])
                               })
  responsibility <- exec(cbind, !!!responsibility)
  type <- vec_slice(colnames(parameters),
                    apply(responsibility, 1, which.max))

  alpha_inlier <- alpha[, "inlier"]
  beta_inlier <- beta[, "inlier"]
  tibble::tibble(rate = dplyr::case_match(type,
                                          "inlier" ~ (y + alpha_inlier) / (x + alpha_inlier + beta_inlier),
                                          "outlier" ~ alpha_inlier / (alpha_inlier + beta_inlier)),
                 type = factor(type, c("inlier", "outlier")))
}
