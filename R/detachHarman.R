#' @title Detach the Harman package and its shared C/C++ library code
#' @description A helper function that can be called if \code{\link{harman}}
#' had to be aborted.
#' @return    None
#' @useDynLib Harman
#' @importFrom Rcpp sourceCpp
detachHarman <- function() {
  if (is.element("Harman", .packages())){
    detach(package:Harman, unload="TRUE")
  }
}
# @example
# library(Harman)
# detachHarman()
# @export
# @details Reinitialises the OpenMP library to use available threads for then
# next time the \code{\link{harman}} function is called.
