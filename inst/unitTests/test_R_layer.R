library(Harman)
library(HarmanData)
library(RUnit)

#####  Global Data  #####

data(OLF)
expt <- olf.info$Treatment
batch <- olf.info$Batch
res <- harman(as.matrix(olf.data), expt, batch)
resf <- harman(as.matrix(olf.data), expt, batch, forceRand=TRUE)


#####  Tests  #####

test.harman <- function() {

  resdf <- harman(olf.data, expt, batch)
  checkEquals(resdf, res)
  checkException(harman(as.character(as.matrix(olf.data)), expt, batch))
  checkException(harman(as.matrix(olf.data), expt, batch, limit=1.1))
  checkException(harman(as.matrix(olf.data), expt, batch, limit=-1))
  checkException(harman(as.matrix(olf.data), expt, batch, limit='limit'))
  checkException(harman(as.matrix(olf.data), expt, batch, randseed='bad'))
  checkException(harman(as.matrix(olf.data), expt[1:5], batch))
  checkException(harman(as.matrix(olf.data), expt, batch[1:5]))
  checkException(harman(as.matrix(olf.data), expt=rep(1, ncol(olf.data)), batch))
  checkException(harman(as.matrix(olf.data), expt, batch=rep(1, ncol(olf.data))))

  badb <- batch
  badb[1] <- NA
  checkException(harman(as.matrix(olf.data), expt, batch=badb))
  
  is_corrected_pc <- res$stats$correction < 1
  checkEquals(res$original[, !is_corrected_pc], res$corrected[, !is_corrected_pc])
}


test.reconstructData <- function() {
  
  remade <- reconstructData(res, this='original')
  checkEquals(remade, as.matrix(olf.data), msg='reconstructData test')
}


test.arrowPlot <- function() {
  
  checkException(arrowPlot(res, 1, 2, colBy='exp'))
  
}


test.pcaPlot <- function() {
  
  checkException(pcaPlot(res, colBy='bad'))
  checkException(pcaPlot(res, pchBy='bad'))
  checkException(pcaPlot(res, this='incorrect'))
  
}


test.forceRand <- function() {

  is_same_correction <- resf$stats$correction == res$stats$correction
  checkEquals(res$corrected[, is_same_correction], resf$corrected[, is_same_correction])
}


test.two_batches <- function() {

  two_batch <- olf.info$Batch %% 2
  res2b <- harman(as.matrix(olf.data), expt=expt, batch=two_batch)
  is_same_correction <- resf$stats$correction == res2b$stats$correction
  checkEquals(res$corrected[, is_same_correction], res2b$corrected[, is_same_correction])
}


test.three_batches <- function() {
  
  three_batch <- olf.info$Batch
  three_batch[three_batch == 4] <- 3
  
  res3b <- harman(as.matrix(olf.data), expt=expt, batch=three_batch)
  is_same_correction <- resf$stats$correction == res3b$stats$correction
  checkEquals(res$corrected[, is_same_correction], res3b$corrected[, is_same_correction])
}

