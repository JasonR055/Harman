#' @title     PCA plot for harman results
#' @description Generates a Principle Component plot for an instance of
#' \code{\link{harmanresults}}.
#' @param     harmanresults An instance of \code{harmanresults}.
#' @param     pc_x integer, principle component for the plot x dimension.
#' @param     pc_y integer, principle component for the plot y dimension.
#' @param     this string, legal values are \code{original} or \code{corrected}.
#' @param     colBy string, colour the points by the experimental or batch
#' variable; legal values
#' are \code{expt} and \code{batch}. The palette function specified in
#' \code{palette} is used.
#' This parameter is overridden by \code{col}.
#' @param     pchBy string, point-type by the experimental or batch variable;
#' legal values are \code{expt} and \code{batch}. This parameter is overridden
#' by \code{pch}.
#' @param     palette string, the function to call to create a vector of
#' contiguous colours with the levels of factor in \code{colBy} steps.
#' @param     legend logical, whether to display a legend on the plot.
#' @param     col, colour vector for the points. This parameter overrides
#' \code{colBy} and \code{palette}.
#' @param     pch, integer vector giving the point type.  This parameter
#' overrides \code{pchBy}.
#' @param     ... further arguments passed to or from other methods.
#' @return    None
#' @details If a vector of colours is supplied via the \code{col} argument,
#' then a legend will not be drawn.
#' @seealso \code{\link{harmanresults}} \code{\link{plot.harmanresults}}
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(as.matrix(olf.data), expt, batch)
#' pcaPlot(olf.harman)
#' pcaPlot(olf.harman, colBy='expt')
#' pcaPlot(olf.harman, pc_x=2, pc_y=3, this='original', pch=17)
#' @export
pcaPlot <- function(harmanresults, pc_x=1, pc_y=2, this='corrected',
                    colBy='batch', pchBy='expt', palette="rainbow",
                    legend=TRUE, col, pch, ...) {

  # Sanity checking
  if(class(harmanresults) != 'harmanresults') {
    stop("Require an object of class 'harmanresults'.")
  }
  if(!(this %in% c('original', 'corrected'))) {
    stop("Require 'this' to be either the values 'original' or 'corrected'.")
  }
  if(!(colBy %in% c('batch', 'expt'))) {
    stop("Require 'colBy' to be either the values 'expt' or 'batch'.")
  }
  if(!(pchBy %in% c('batch', 'expt'))) {
    stop("Require 'pchBy' to be either the values 'expt' or 'batch'.")
  }

  scores <- harmanresults[[this]]
  draw_legend <- FALSE
  
  if(missing(col)) {
    num_levels <- length(levels(harmanresults$factors[, colBy]))
    palette <- match.fun(palette)(num_levels)
    col <- palette[harmanresults$factors[, colBy]]
    if(legend == TRUE) { draw_legend <- TRUE }
  }
  
  if(missing(pch)) {
    num_levels <- length(levels(harmanresults$factors[, pchBy]))
    #num_expt <- 1:length(levels(harmanresults$factors$expt))
    level_vector <- 1:num_levels
    #pch <- level_vector[harmanresults$factors$expt]
    pch <- level_vector[harmanresults$factors[[pchBy]]]
  }
  
  plot(x=scores[, pc_x],
       y=scores[, pc_y],
       col=col,
       pch=pch,
       xlab=colnames(scores)[pc_x],
       ylab=colnames(scores)[pc_y],
       ...)
  
  if(draw_legend == TRUE) {
    legend(x=min(scores[, pc_x]), y=max(scores[, pc_y]),
           legend=levels(harmanresults$factors[, colBy]),
           fill=palette, cex=0.7,
           bg="transparent")
  }
}
