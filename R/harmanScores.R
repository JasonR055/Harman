#' @title     A Principal components prcomp function tweaked for Harman
#' @description A tweaking of \code{stats::prcomp} such that for the svd, the
#' transpose of u is used instead of v when the number of assays is less than
#' the number of samples.
#' @param     x matrix, data matrix of values to perform PCA on.
#' @return    scores, a prcomp-like object with rotation, scores and the center
#' values. The scores are corrected, but all three are needed later to
#' reconstruct the data.
harmanScores <- function(x) {

  # x <- cars
  # x <- t(cars)
  # library(HarmanData)
  # x <- olf.data
  # x <- olf.data[1:60, 1:20]
  # x <- olf.data[1:16, 1:20]
  # x <- df16_scaled
  x <- as.matrix(x)

  if(nrow(x) < ncol(x)) {
    # Special case for less assays than samples, need to use u' instead of v
    
    # For now, until the code is complete, we quit.
    stop("Cannot presently handle cases where the number of features is less
         than samples")
    
    tx <- t(x)
    x <- t(scale(tx, center = TRUE, scale = FALSE))
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0)) {
      stop("cannot rescale a constant/zero column to unit variance")
    }
    
    s <- svd(x)
    
    rotation <- s$u
    #dimnames(rotation) <- list(colnames(x),
    #                           paste0("PC", seq_len(ncol(rotation))))
    scores <- t(x) %*% rotation
  } else {
    # Standard case like in prcomp()
    
    x <- scale(t(x), center = TRUE, scale = FALSE)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0)) {
      stop("cannot rescale a constant/zero column to unit variance")
    }
    
    s <- svd(x, nu = 0)
    rotation <- s$v
    dimnames(rotation) <- list(colnames(x),
                               paste0("PC", seq_len(ncol(rotation))))
    scores <- x %*% rotation
  }

  sdev <- s$d / sqrt(max(1, nrow(x) - 1))

  r <- list(sdev = sdev, rotation = rotation,
            center = if(is.null(cen)) FALSE else cen,
            scale = if(is.null(sc)) FALSE else sc,
            scores = scores)
  class(r) <- "scores"
  r
}

