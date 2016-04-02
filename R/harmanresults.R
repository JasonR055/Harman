#' @name harmanresults
#' @title Harman results object
#' @description The S3 object returned after running \code{\link{harman}}.
#' @slot factors \code{A data.frame} of the \code{expt} and \code{batch} vectors.
#' @slot parameters The harman runtime parameters. See \code{\link{harman}}
#' for details.
#' @slot stats Confidence intervals and the degree of correction for each
#' pricipal component.
#' @slot center The centering vector returned by \code{\link{prcomp}} with
#' \code{center=TRUE}.
#' @slot rotation The matrix of eigenvectors (by column) returned from
#' \code{\link{prcomp}}.
#' @slot original The original PC scores returned by \code{\link{prcomp}}.
#' @slot corrected The harman corrected PC scores.
#' @details \code{harmanresults} is the S3 object used to store the results from
#' \code{\link{harman}}.
#' This object may be presented to summary and data exploration functions such
#' as \code{\link{plot.harmanresults}}
#' and \code{\link{summary.harmanresults}} as well as the
#' \code{\link{reconstructData}} function which creates a corrected matrix of
#' data with the batch effect removed.
#' @seealso \code{\link{harman}}, \code{\link{reconstructData}},
#' \code{\link{pcaPlot}}, \code{\link{arrowPlot}}
#' @exportClass harmanresults
#' @examples
#' ## HarmanResults
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(as.matrix(olf.data), expt, batch)
#' plot(olf.harman)
#' summary(olf.harman)
#' pcaPlot(olf.harman, pc_x=2, pc_y=3)
#' pcaPlot(olf.harman, pc_x=2, pc_y=3, colBy='expt', pch=1)
#' olf.data.corrected <- reconstructData(olf.harman)
NULL