#' @title     Shift beta values from 0 and 1 to avoid infinite M values
#' @description A convienance function for methylation data.
#' @param     betas matrix, beta values.
#' @param     shiftBy numeric, the amount to shift values of \code{0} and
#' \code{1} by.
#' @return    None
#' @examples
#' betas <- seq(0, 1, by=0.05)
#' range(betas)
#' newBetas <- shiftBetas(betas, shiftBy=1e-4)
#' newBetas
#' range(newBetas)
#' @export
shiftBetas <- function(betas, shiftBy=1e-4) {
  
  if(sum(is.nan(betas) > 0)) { stop('Beta matrix contains NaN values')}
  if(sum(is.na(betas) > 0)) { stop('Beta matrix contains NA values')}
  
  is_zero <- betas == 0
  is_one <- betas == 1
  
  if(sum(is_zero) > 0) {
    betas[is_zero] <- betas[is_zero] + shiftBy
  }
  
  if(sum(is_one) > 0) {
    betas[is_one] <- betas[is_one] - shiftBy
  }
  
  if(sum(is_zero) == 0 & sum(is_one) == 0) {
    cat("No beta values were found to be 0 or 1. No shifts made.\n")
  }
  betas
}

