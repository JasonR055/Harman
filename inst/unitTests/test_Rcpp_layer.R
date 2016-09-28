library(Harman)
library(HarmanData)
library(RUnit)


#####  Tests  #####

test.callHarman <- function() {

  data(IMR90)
  pca <- prcomp(t(imr90.data), retx=TRUE, center=TRUE)
  #pc_data_scores <- pca$x[, 1:(ncol(pca$x) - 1)]
  pc_data_scores <- pca$x
  expt <- as.factor(imr90.info$Treatment)
  batch <- as.factor(imr90.info$Batch)
  group <- as.matrix(data.frame(expt=as.integer(expt), batch=as.integer(batch)))
  res1 <- Harman:::.callHarman(pc_data_scores, group, limit=0.95,
                               numrepeats=1e5, randseed=64, forceRand=FALSE,
                               printInfo=FALSE)
  res2 <- harman(imr90.data, expt=imr90.info$Treatment, batch=imr90.info$Batch)
  checkEquals(res1$corrected_scores, as.numeric(res2$corrected))
  
}


test.forceRand <- function() {
  
  precision_limit <- 0.02
  
  data(IMR90)
  imr <- list()
  imr[['F']] <- harman(imr90.data, expt=imr90.info$Treatment,
                       batch=imr90.info$Batch, numrepeats=1e5, randseed=42,
                       forceRand=FALSE)
  
  imr[['T']] <- harman(imr90.data, expt=imr90.info$Treatment,
                       batch=imr90.info$Batch, numrepeats=1e5, randseed=42,
                       forceRand=TRUE)
  
  data(NPM)
  npm <- list()
  npm[['F']] <- harman(npm.data, expt=npm.info$Treatment,
                       batch=npm.info$Batch, numrepeats=1e5, randseed=42,
                       forceRand=FALSE)
  
  npm[['T']] <- harman(npm.data, expt=npm.info$Treatment,
                       batch=npm.info$Batch, numrepeats=1e5, randseed=42,
                       forceRand=TRUE)
  
  data(OLF)
  olf <- list()
  olf[['F']] <- harman(olf.data, expt=olf.info$Treatment,
                       batch=olf.info$Batch, numrepeats=3e5, randseed=42,
                       forceRand=FALSE)
  
  olf[['T']] <- harman(olf.data, expt=olf.info$Treatment,
                       batch=olf.info$Batch, numrepeats=3e5, randseed=42,
                       forceRand=TRUE)
  
  results <- list(imr=imr, npm=npm, olf=olf)
  
  for(i in 1:length(results)) {

    res <- results[[i]]
    this_set <- names(results)[i]
    
    is_same_amount <- res[['F']]$stats$correction == res[['T']]$stats$correction

    # Corrected scores are within a tolerance of precision_limit
    x <- abs(res[['F']]$stats$correction - res[['T']]$stats$correction)
    checkEquals(x  <= precision_limit, rep(TRUE, length(x)),
                msg=paste('forceRand', this_set, 'test: rounding tolerance'))
    
    # Original scores should the same
    checkEquals(res[['F']]$original, res[['T']]$original,
                msg=paste('forceRand', this_set, 'test: original scores'),
                tolerance = .Machine$double.eps ^ 0.5)

    # The PCs which are corrected by the same amount should be the same
    checkEquals(res[['F']]$corrected[, is_same_amount],
                res[['T']]$corrected[, is_same_amount],
                msg=paste('forceRand', this_set, 'test: corrected scores'),
                tolerance = .Machine$double.eps ^ 0.5)
    
  }

}

