#' @title     Reconstruct corrected data from Harman results
#' @description Method which reverts the PCA factorisation for instances of
#' \code{\link{harmanresults}}. This allows the original or corrected data to be
#' returned back from the PCA domain into the original data domain.
#' @param object An instance of \code{harmanresults}.
#' @param this string, legal values are \code{original} or \code{corrected}.
#' @return \code{matrix} of data
#' @seealso \code{\link{harman}} \code{\link{harmanresults}}
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' olf.data.corrected <- reconstructData(olf.harman)
#' @export
reconstructData <- function(object, this='corrected')  {
  
  # The Matlab code from which this function is built:
  # reconstructedmicroarray = corrected_scorebatch_matrix*coeff'+\
  # ones(n,1)*means_initialprobesets;
  # corrected_scorebatch_matrix == object[[this]]
  # coeff == t(object$rotation), so coeff' == object$rotation
  # n == number of samples
  # ones(n,1) == matrix(1, n, 1)
  # means_initialprobesets == object$center
  
  if(class(object) != 'harmanresults') {
    stop(paste("Require an instance of 'harmanresults', not class \'",
               class(object), "\'.", sep=""))
  }
  
  if(!(this %in% c('original', 'corrected'))) {
    stop("Require 'this' to be either the values 'original' or 'corrected'.")
  }
  
  n <- ncol(object$rotation)
  #n <- nrow(object$rotation)
  #zeros <- matrix(0, n, 1)
  ones <- matrix(1, n, 1)
  
  # Add the extra column of zeros.
  #scores <- cbind(object[[this]], zeros)
  scores <- object[[this]]
  
  # corrected_scorebatch_matrix*coeff'
  centred_scores_matrix <- scores %*% t(object$rotation)
  # ones(n,1)*means_initialprobesets
  means_matrix <- ones %*% matrix(object$center, nrow=1)
  # Now add the means and transpose
  reconstructed <- t(centred_scores_matrix + means_matrix)
  reconstructed
}
