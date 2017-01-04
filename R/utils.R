#' @param verbose \code{logical}, verbose output?
#' @param \ldots further arguments passed to \code{message}.
#' @noRd
.msg <- function(verbose, ...) {
  if (verbose) {
    message(...)
  }
}
