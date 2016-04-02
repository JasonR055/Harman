# library(xlsx)
# library(Harman)

# ## Load R data and compute scores

# data(NPM)
# expt <- npm.info$Treatment
# batch <- npm.info$Batch
# npm.harman <- harman(npm.data, expt, batch)

# ## Load Matlab corrected scores

# xls <- read.xlsx(file.path(paste0(path.package('Harman'),'/unitTests/corrected_scorebatch_nonpreg95.xls')), 1, header=FALSE)
# xls <- as.matrix(xls)
# dimnames(xls) <- dimnames(npm.harman$corrected) 

# ## Correlate, storing the results in a list of vectors (lv)

# lv <- list()

# for(i in 1:ncol(npm.harman$corrected)) {
  # colvec <- npm.harman$corrected[, i]
  # lv[[i]] <- apply(xls, 2, cor, colvec)
# }

# ## Munge list of vectors into a matrix

# corr_matrix <- do.call("cbind", lv)
# colnames(corr_matrix) <- rownames(corr_matrix)

# write.csv(round(corr_matrix, 3), file='PearsonCorrelationMatrix_Matlab_with_R.csv')
# write.csv(npm.harman$stats, file='Stats_NPM_R.csv')




