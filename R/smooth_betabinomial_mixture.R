FLXMR_betabinomial_mixture <- function(formula = . ~ .,
                                       offset, control) {
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
  out@fit <- function(x, y, w, component) {
    alpha <- component$alpha %||% 1
    beta <- component$beta %||% 1

    get_par <- function(par) {
      list(alpha = par[[1]],
           beta = par[[2]])
    }
    minuslogl <- function(par) {
      par <- get_par(par)
      -sum(w *
             extraDistr::dbbinom(y,
                                 size = offset,
                                 alpha = par$alpha,
                                 beta = par$beta,
                                 log = TRUE))
    }
    fitted <- stats::optim(c(alpha, beta), minuslogl,
                           method = "L-BFGS-B",
                           lower = vec_rep(sqrt(.Machine$double.eps), 2),
                           control = purrr::compact(list(trace = control[["verbose"]],
                                                         maxit = control[["max_iter"]],

                                                         # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
                                                         factr = control[["tolerance"]] / .Machine$double.eps)))
    par <- get_par(fitted$par)
    out@defineComponent(para = list(alpha = par[["alpha"]],
                                    beta = par[["beta"]],
                                    df = 2))
  }
  out
}

mean_beta <- function(alpha, beta) {
  alpha / (alpha + beta)
}

var_beta <- function(alpha, beta) {
  alpha * beta / (alpha + beta) ^ 2 / (alpha + beta + 1)
}

flexmix_betabinomial_mixture <- function(k, y, x, control) {
  fitted <- flexmix::flexmix(y ~ 1,
                             data = data_frame(y = y),
                             k = k,
                             model = FLXMR_betabinomial_mixture(offset = x,
                                                                control = control[["local_control"]]),
                             control = purrr::compact(list2(iter.max = control[["max_iter"]],
                                                            minprior = .Machine$double.eps,
                                                            verbose = as.numeric(control[["verbose"]]),
                                                            !!!control[!names(control) %in% c("max_iter", "verbose", "local_control")])))
  parameters <- rbind(flexmix::parameters(fitted),
                      prior = flexmix::prior(fitted))
  loc_inlier <- which.min(var_beta(alpha = parameters["alpha", ],
                                   beta = parameters["beta", ]))
  colnames(parameters) <- vec_rep("outlier", ncol(parameters))
  colnames(parameters)[loc_inlier] <- "inlier"

  tibble::tibble(k = k,
                 parameters = list(parameters),
                 BIC = stats::BIC(fitted))
}

#' @export
smooth_betabinomial_mixture <- function(y, x, k,
                                        control = control_smooth()) {
  fitted <- purrr::map(k,
                       flexmix_betabinomial_mixture,
                       y = y,
                       x = x,
                       control = control)
  structure(purrr::list_rbind(fitted),
            class = "smooth_betabinomial_mixture")
}

#' @export
predict.smooth_betabinomial_mixture <- function(object, y, x) {
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
  tibble::tibble(type = factor(type, c("inlier", "outlier")),
                 alpha = dplyr::case_match(type,
                                           "inlier" ~ y + alpha_inlier,
                                           "outlier" ~ alpha_inlier),
                 beta = dplyr::case_match(type,
                                          "inlier" ~ x - y + beta_inlier,
                                          "outlier" ~ beta_inlier),
                 mean = dplyr::case_match(type,
                                          "inlier" ~ mean_beta(alpha = alpha,
                                                               beta = beta),
                                          "outlier" ~ mean_beta(alpha = alpha,
                                                                beta = beta)),
                 var = dplyr::case_match(type,
                                         "inlier" ~ var_beta(alpha = alpha,
                                                             beta = beta),
                                         "outlier" ~ var_beta(alpha = alpha,
                                                              beta = beta)))
}
