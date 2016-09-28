#' @title     PCA before and after arrow plot for harman results
#' @description Generates an arrow plot for an instance of
#' \code{\link{harmanresults}}. The tail of the arrow is the starting point
#' (original) in principle coordinates, while the arrow head is the new point
#' (corrected) in principle coordinates. It can be observed that on principle
#' components that have undergone correction
#' (code{harmanresults$stats$correction < 1.0}), the samples within a batch will
#' be coordinately moved towards \code{0} on that priciple component.
#' @param     harmanresults an instance of \code{harmanresults}.
#' @param     pc_x integer, principle component for the plot x dimension.
#' @param     pc_y integer, principle component for the plot y dimension.
#' @param     colBy string, colour the points by the experimental or batch
#' variable; legal values are \code{expt} and \code{batch}. The palette function
#' specified in \code{palette} is used. This parameter is overridden by
#' \code{col}.
#' @param     palette string, the function to call to create a vector of
#' contiguous colours with the levels of factor in \code{colBy} steps.
#' @param     col, colour vector for the points. This parameter overrides
#' \code{palette}.
#' @param     length length of the \code{\link{arrow}} heads, default is 0.1.
#' @param     legend logical, whether to display a legend on the plot
#' @param     ... further arguments passed to or from other methods.
#' @return None
#' @details   Generates a Principle Component plot for an instance of
#' \code{harmanresults}. If a vector of colours is supplied via the \code{col}
#' argument, then a legend will not be drawn.
#' @seealso \code{\link{harmanresults}} \code{\link{plot.harmanresults}}
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' arrowPlot(olf.harman, pc_x=2, pc_y=3, length=0.2)
#' @importFrom graphics arrows points legend plot
#' @export
arrowPlot <- function(harmanresults, pc_x=1, pc_y=2, colBy='batch',
                      palette="rainbow", col, length=0.1, legend=TRUE, ...) {
    
    # Sanity checking
    if(class(harmanresults) != 'harmanresults') {
      stop("Require an object of class 'harmanresults'.")
    }
    if(!(colBy %in% c('batch', 'expt'))) {
      stop("Require 'colBy' to be either the values 'expt' or 'batch'.")
    }
    
    xrange <- range(c(harmanresults$original[, pc_x],
                      harmanresults$corrected[, pc_x]))
    yrange <- range(c(harmanresults$original[, pc_y],
                      harmanresults$corrected[, pc_y]))
    
    graphics::plot(x=NA,
                   y=NA,
                   xlab=colnames(harmanresults$original)[pc_x],
                   ylab=colnames(harmanresults$original)[pc_y],
                   xlim=xrange,
                   ylim=yrange,
                   type='n')
    
    if(missing(col)) {
      mylegend <- harmanresults$factors[, colBy]
      palette <- match.fun(palette)(length(levels(mylegend)))
      col <- palette[mylegend]
      
      if(legend == TRUE) {
        graphics::legend(x="topleft", # x=min(xrange), y=max(yrange)
                         legend=levels(mylegend),
                         fill=palette, cex=0.7)
        }
    }
    
    if(harmanresults$stats$correction[pc_x] != 1 ||
       harmanresults$stats$correction[pc_y] != 1) {
      graphics::arrows(x0=harmanresults$original[, pc_x],
                       y0=harmanresults$original[, pc_y],
                       x1=harmanresults$corrected[, pc_x],
                       y1=harmanresults$corrected[, pc_y],
                       col=col,
                       length=length, ...)
      } else {     # Case for no correction
        graphics::points(harmanresults$original[, pc_x],
                         harmanresults$original[, pc_y],
                         col=col,
                         pch=4,
                         ...)
        }
}

