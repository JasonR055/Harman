#' @title     Plot method for harman
#' @description Plot method for instances of \code{\link{harmanresults}}.
#' @param     x An instance of \code{harmanresults}.
#' @param     ... further plotting parameters.
#' @return    None
#' @seealso  \code{\link{harmanresults}} \code{\link{pcaPlot}}
#' @examples 
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' plot(olf.harman)
#' @export
plot.harmanresults <- function(x, ...) {
  # set xlim and ylim  to be the same
  this_pc_x <- 1
  this_pc_y <- 2
  
  params <- list(...)
  
  if('pc_x' %in% names(params)) {
    this_pc_x <- params[['pc_x']]
  } else {
    if(length(params) >= 1 && is.numeric(params[[1]])) {
      this_pc_x <- params[[1]]
    }
  }
  if('pc_y' %in% names(params)) {
    this_pc_y <- params[['pc_y']]
  } else {
    if(length(params) >= 2 && is.numeric(params[[2]])) {
      this_pc_y <- params[[2]]
    }
  }
  
  xrange <- range(c(x$original[, this_pc_x], x$corrected[, this_pc_x]))
  yrange <- range(c(x$original[, this_pc_y], x$corrected[, this_pc_y]))
  old_mfrow <- par()$mfrow
  par(mfrow=c(1, 2))
  pcaPlot(x, this='original', main='Original', xlim=xrange, ylim=yrange,
          legend=TRUE, ...)
  pcaPlot(x, this='corrected', main='Corrected', xlim=xrange, ylim=yrange,
          legend=FALSE, ...)
  par(mfrow=old_mfrow)
}
