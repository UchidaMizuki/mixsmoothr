FLXMR_gammapoisson_mixture <- function(formula = . ~ .,
                                       offset, control) {
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
  out@fit <- function(x, y, w, component) {
    alpha <- component$alpha %||% 1
    beta <- component$beta %||% 1
    minuslogl <- function(alpha, beta) {
      -sum(w *
             stats::dnbinom(y,
                            size = alpha,
                            mu = alpha / beta * offset,
                            log = TRUE))
    }
    fitted <- stats4::mle(minuslogl,
                          start = list(alpha = alpha,
                                       beta = beta),
                          lower = vec_rep(.Machine$double.eps, 2),
                          control = purrr::compact(list(trace = control[["verbose"]],
                                                        maxit = control[["max_iter"]],

                                                        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
                                                        factr = control[["tolerance"]] / .Machine$double.eps)))
    out@defineComponent(para = list(alpha = fitted@coef[["alpha"]],
                                    beta = fitted@coef[["beta"]],
                                    df = 2))
  }
  out
}

mean_gamma <- function(alpha, beta) {
  alpha / beta
}

var_gamma <- function(alpha, beta) {
  alpha / beta ^ 2
}

flexmix_gammapoisson_mixture <- function(k, y, x, control) {
  fitted <- flexmix::flexmix(y ~ 1,
                             data = data_frame(y = y),
                             k = k,
                             model = FLXMR_gammapoisson_mixture(offset = x,
                                                                control = control[["local_control"]]),
                             control = purrr::compact(list2(iter.max = control[["max_iter"]],
                                                            minprior = .Machine$double.eps,
                                                            verbose = as.numeric(control[["verbose"]]),
                                                            !!!control[!names(control) %in% c("max_iter", "verbose", "local_control")])))
  parameters <- rbind(flexmix::parameters(fitted),
                      prior = flexmix::prior(fitted))
  loc_inlier <- which.min(var_gamma(alpha = parameters["alpha", ],
                                    beta = parameters["beta", ]))
  colnames(parameters) <- vec_rep("outlier", ncol(parameters))
  colnames(parameters)[loc_inlier] <- "inlier"

  tibble::tibble(k = k,
                 parameters = list(parameters),
                 BIC = stats::BIC(fitted))
}

#' @export
smooth_gammapoisson_mixture <- function(y, x, k,
                                        control = control_smooth()) {
  fitted <- purrr::map(k,
                       flexmix_gammapoisson_mixture,
                       y = y,
                       x = x,
                       control = control)
  structure(purrr::list_rbind(fitted),
            class = "smooth_gammapoisson_mixture")
}

#' @export
predict.smooth_gammapoisson_mixture <- function(object, y, x) {
  parameters <- object$parameters[[which.min(object$BIC)]]
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
  tibble::tibble(type = factor(type, c("inlier", "outlier")),
                 alpha = dplyr::case_match(type,
                                           "inlier" ~ y + alpha_inlier,
                                           "outlier" ~ alpha_inlier),
                 beta = dplyr::case_match(type,
                                          "inlier" ~ x + beta_inlier,
                                          "outlier" ~ beta_inlier),
                 mean = dplyr::case_match(type,
                                          "inlier" ~ mean_gamma(alpha = alpha,
                                                                beta = beta),
                                          "outlier" ~ mean_gamma(alpha = alpha,
                                                                 beta = beta)),
                 var = dplyr::case_match(type,
                                         "inlier" ~ var_gamma(alpha = alpha,
                                                              beta = beta),
                                         "outlier" ~ var_gamma(alpha = alpha,
                                                               beta = beta)))
}
