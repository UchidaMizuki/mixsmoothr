#' @export
control_smooth <- function(verbose = FALSE,
                           tolerance = NULL,
                           max_iter = NULL,
                           local_control = NULL,
                           ...) {
  list2(verbose = verbose,
        tolerance = tolerance,
        max_iter = max_iter,
        local_control = local_control,
        ...)
}
