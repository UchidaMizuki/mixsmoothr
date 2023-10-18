#' @export
control_smooth <- function(verbose = FALSE,
                           tolerance = NULL,
                           max_iter = NULL,
                           initial = NULL,
                           init_control = NULL,
                           local_control = NULL,
                           ...) {
  list2(verbose = verbose,
        tolerance = tolerance,
        max_iter = max_iter,
        initial = initial,
        init_control = init_control,
        local_control = local_control,
        others = list2(...))
}
