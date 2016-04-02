#' @name callHarman
#' @title Wrapper function to call the shared C/C++ library code
#' @description This wrapper should probably not be addressed directly except
#' for debugging. Instead use \code{\link{harman}}. Input of PCA scores and the
#' experiment structure (treatments and batches) and returns a batch corrected
#' version of the PCA scores matrix
#' @param pc_data_scores 2D NumericMatrix of PCA scores data (from the
#' \code{prcomp$x} slot), rows = samples, cols = PC scores
#' @param group The structure of the experiment, consisting of batch numbers and
#' treatment numbers forming 2 rows or columns (HarmanMain works out which).
#' Each entry for a sample describes what batch it came from and what treatment
#' it was given. Has to be integer formated data.
#' @param limit A double precsion value indicating the limit of confidence in
#' which to stop removing a batch effect
#' @param numrepeats The number of repeats in which to run the simulated batch
#' mean distribution estimator. Probably should be greater than 100,000.
#' @param randseed Random seed to pass to the random number generator (0 for use
#' default from system time)
#' @param forceRand Force algorithm
#' @param printInfo Print update information to screen 
#' @return SEXP  R list: scores.corrected  = harman_res_list["corrected_scores"]
#'                             correction        = harman_res_list["correction"]
#'                             confidence        = harman_res_list["confidence"]
#' @note A data matrix with samples in columns must be transposed before PCA
#' analysis and these scores in turn are tweaked a little before handing over
#' to .callHarman. See the example below.
#' @useDynLib Harman
#' @importFrom Rcpp sourceCpp
.callHarman <- function(pc_data_scores, group, limit, numrepeats, randseed,
                        forceRand, printInfo) {
  
  .Call("HarmanMain", pc_data_scores, as.matrix(group), as.double(limit),
        as.integer(numrepeats), as.integer(randseed),
        as.logical(forceRand), as.logical(printInfo), PACKAGE = "Harman")
}

