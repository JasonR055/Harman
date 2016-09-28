#' @title     PCA plot
#' @description Generates a Principle Component plot for data.frames, matrices,
#' or a pre-made \code{\link{prcomp}} object.
#' @param     object data.frame, matrix or \code{prcomp} object.
#' @param     pc_x integer, principle component for the plot x dimension.
#' @param     pc_y integer, principle component for the plot y dimension.
#' @param     scale logical, whether to scale to unit variance before PCA.
#' @param     colFactor factor or vector, colour the points by this factor,
#' default is \code{NULL}. 
#' @param     pchFactor factor or vector, point-type by this factor,
#' default is \code{NULL}.
#' @param     palette string, the function to call to create a vector of
#' contiguous colours with \code{levels(colFactor)} steps.
#' @param     legend logical, whether to display a legend on the plot.
#' @param     ... further arguments passed to or from other methods.
#' @return    None
#' @details   A data.frame object will be coerced internally to a matrix.
#' Matrices must be of type \code{double} or \code{integer}. The
#' \code{prcompPlot} function will then perform a principle component analysis
#' on the data prior to plotting. The function is call
#' is \code{prcomp(t(object), retx=TRUE, center=TRUE, scale.=scale)}.
#' Instead of specifying a data.frame or matrix, a pre-made \code{prcomp} object
#' can be given to \code{prcompPlot}. In this case, care should be taken in
#' setting the appropriate value of \code{scale.}. If a vector is given to
#' \code{colFactor} or \code{pchFactor}, they will be coerced internally to
#' factors.
#' 
#' For the default \code{NULL} values of \code{colFactor} and \code{pchFactor},
#' all colours will be black and circles the point type, respectively.
#' @seealso \code{\link{prcomp}} \code{\link{rainbow}}
#' @importFrom graphics legend plot
#' @importFrom stats prcomp
#' @export
#' @examples
#' library(HarmanData)
#' data(IMR90)
#' expt <- imr90.info$Treatment
#' batch <- imr90.info$Batch
#' prcompPlot(imr90.data, colFactor=expt)
#' pca <- prcomp(t(imr90.data), scale.=TRUE)
#' prcompPlot(pca, 1, 3, colFactor=batch, pchFactor=expt, palette='topo.colors',
#' main='IMR90 PCA plot of Dim 1 and 3')

prcompPlot <- function(object, pc_x=1, pc_y=2, scale=FALSE, colFactor=NULL,
                       pchFactor=NULL, palette="rainbow", legend=TRUE, ...) {

  # Sanity check object
  if(class(object) == "data.frame") {
    object <- as.matrix(object)
  }
  if(class(object) == "matrix" & typeof(object) %in% c('double', 'integer')) {
    object <- stats::prcomp(t(object), retx=TRUE, center=TRUE, scale.=scale)
  }
  if(class(object) != "prcomp") {
    stop("Require an instance of 'prcomp', a matrix of type 'double' or
         'integer', or a data.frame coercible to such a matrix.")
  }
  
  # Sanity check col and pch factors
  if(is.null(colFactor)) {
    legend <- FALSE
    colFactor <- rep('', ncol(object$x))
  }
  if(is.null(pchFactor)) {
    pchFactor <- rep(1, ncol(object$x))
  }
  if(!is.factor(colFactor)) {
    colFactor <- factor(colFactor)
  }
  if(!is.factor(pchFactor)) {
    pchFactor <- factor(pchFactor)
  }
  
  # Sanity check col/pch lengths
  if(length(colFactor) != ncol(object$x)) {
    stop('The length of colFactor and object do not match.')
  }
  if(length(pchFactor) != ncol(object$x)) {
    stop('The length of pchFactor and object do not match.')
  }
  
  mypchs <- (1:length(levels(pchFactor)))[pchFactor]  
  
  factor_names <- levels(colFactor)
  num_levels <- length(factor_names)
  mypalette <- match.fun(palette)(num_levels)
  mycols <- mypalette[colFactor]
  
  graphics::plot(object$x[, pc_x], object$x[, pc_y],
                 xlab=paste('PC', pc_x, sep=''),
                 ylab=paste('PC', pc_y, sep=''),
                 col=mycols,
                 pch=mypchs,
                 ...)
  
  if(legend == TRUE) {
    #legend(x=min(object$x[, pc_x]), y=max(object$x[, pc_y]),
    graphics::legend(x="topleft",
                     legend=factor_names,
                     fill=mypalette,
                     cex=0.7,
                     bg="transparent")
    }
}

