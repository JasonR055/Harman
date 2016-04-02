#' @title     Summarizing harman results.
#' @description Summary method for class \code{\link{harmanresults}}.
#' @param object An object of class \code{harmanresults}.
#' @param ... further parameters.
#' @return    Returns an object of class \code{summary.harmanresults}.
#' @seealso  \code{\link{harmanresults}}
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' summary(olf.harman)
#' @export
summary.harmanresults <- function(object, ...) {
  # 1) % of the variance removed.
  # 2) Sequence of corrections from the 1st to last PC
  # 3) Confidence threshold
  # 4) Whether is was balanced or not
  # 5) Parameters
  # 6) Batch structure
  
  ans <- list()
  ans$totals <- c(length(levels(object$factors$expt)),
                  length(levels(object$factors$batch))
                  )
  names(ans$totals) <- c('expt', 'batch')
  ans$factor_table <- table(expt=object$factors$expt,
                            batch=object$factors$batch)
  ans$parameters <- object$parameters
  ans$correction <- object$stats$correction
  names(ans$correction) <- object$stats$dimension
  class(ans) <- c("summary.harmanresults")
  ans
}
#' @title Printing Harmanresults summaries.
#' @param x an object of class \code{summary.harmanresults}, usually, a result
#' of a call to \code{summary.harmanresults}.
#' @param ... further parameters.
#' @description Print method for \code{summary.harmanresults}.
#' @return    Prints summary information from an object of class
#' \code{summary.harmanresults}.
#' @export
print.summary.harmanresults <- function(x, ...) {
  
  cat('Total factor levels:\n\n')
  print(x$totals)
  cat('\nExperiment x Batch Design:\n\n')
  print(x$factor_table)
  cat('\nSimulation parameters:\n')
  cat(x$parameters$numrepeats, ' simulations (with seed of ',
      x$parameters$randseed, '). ForceRand is ', x$parameters$forceRand, '.\n',
      sep='')
  cat('\nHarman results with confidence limit of ', x$parameters$limit, ':\n',
      sep='')
  print(x$correction)
  cat('\nTop batch-effected PCs:\n')
  top <- sort(x$correction)
  top <- top[top < 1]
  ntop <- min(5, length(top))
  if(ntop > 0) {
    print(top[1:ntop])
  } else {
    cat('None, all have no correction with this confidence limit\n')
  }
}