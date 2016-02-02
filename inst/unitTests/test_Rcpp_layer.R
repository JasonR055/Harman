library(Harman)
library(HarmanData)
library(RUnit)


#####  Tests  #####

test.callHarman <- function() {

  data(IMR90)
  pca <- prcomp(t(imr90.data), retx=TRUE, center=TRUE)
  pc_data_scores <- pca$x[, 1:(ncol(pca$x) - 1)]
  expt <- as.factor(imr90.info$Treatment)
  batch <- as.factor(imr90.info$Batch)
  group <- as.matrix(data.frame(expt=as.integer(expt), batch=as.integer(batch)))
  res1 <- Harman:::.callHarman(pc_data_scores, group, limit=0.95, numrepeats=1e5, randseed=64, forceRand=FALSE, printInfo=FALSE)
  res2 <- harman(imr90.data, expt=imr90.info$Treatment, batch=imr90.info$Batch)
  checkEquals(res1$corrected_scores, as.numeric(res2$corrected))
  
}


test.forceRand <- function() {
  
  data(IMR90)
  imr1 <- harman(as.matrix(imr90.data), expt=imr90.info$Treatment, batch=imr90.info$Batch, numrepeats=1e5, randseed=42, forceRand=FALSE)
  imr2 <- harman(as.matrix(imr90.data), expt=imr90.info$Treatment, batch=imr90.info$Batch, numrepeats=1e5, randseed=42, forceRand=TRUE)
  
  is_same_correction <- imr1$stats$correction == imr2$stats$correction
  # Within a tolerance of 0.01
  checkEquals(abs(round(imr1$stats$correction - imr2$stats$correction, 2)) <= 0.01, rep(TRUE, length(imr1$stats$correction)))
  #all.equal(res1$corrected, res2$corrected, tolerance = .Machine$double.eps ^ 0.5)
  checkEquals(imr1$original, imr2$original, msg='forceRand IMR90 test: original scores', tolerance = .Machine$double.eps ^ 0.5)
  checkEquals(imr1$corrected[, is_same_correction], imr2$corrected[, is_same_correction], msg='forceRand IMR90 test: corrected scores', tolerance = .Machine$double.eps ^ 0.5)

  data(NPM)
  npm1 <- harman(as.matrix(npm.data), expt=npm.info$Treatment, batch=npm.info$Batch, numrepeats=1e5, randseed=42, forceRand=FALSE)
  npm2 <- harman(as.matrix(npm.data), expt=npm.info$Treatment, batch=npm.info$Batch, numrepeats=1e5, randseed=42, forceRand=TRUE)

  is_same_correction <- npm1$stats$correction == npm2$stats$correction
  # Within a tolerance of 0.01
  checkEquals(abs(round(npm1$stats$correction - npm2$stats$correction, 2)) <= 0.01, rep(TRUE, length(npm1$stats$correction)))
  #all.equal(res1$corrected, res2$corrected, tolerance = .Machine$double.eps ^ 0.5)
  checkEquals(npm1$original, npm2$original, msg='forceRand NPM test: original scores', tolerance = .Machine$double.eps ^ 0.5)
  checkEquals(npm1$corrected[, is_same_correction], npm2$corrected[, is_same_correction], msg='forceRand NPM test: corrected scores', tolerance = .Machine$double.eps ^ 0.5)
}

