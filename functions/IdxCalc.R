# Function to be optimized to obtain best weights/coefficients
# for the index linear combination

# I think using CV on it would be ideal...
# And maximizing the mean...
IdxCalc <- function(coefs, prxy, target, matG, vFolds){
  
  traits <- prxy$traits
  weights <- prxy$weights # BLUE calculation weights, not
  # the coefficients we seek for the index
  
  # Vector where each element corresponds to the index for a genotype
  Idx <- as.matrix(traits) %*% coefs
  
  # Weight vector for the new index variable
  wtIdx <- 1/((1/as.matrix(weights)) %*% (coefs^2))
  
  # DF with the index values as "BLUEs", and the index weights as weights
  IdxDF <- data.frame(genotype = prxy$genotypes,
                      BLUE = Idx,
                      weight = wtIdx)
  
  # CV parameters
  nrep <- length(vFolds)
  k <- length(vFolds[[1]]) # all equal length
  
  # Vector of correlations (one for each repetition of k-fold CV)
  cors <- numeric(nrep)
  for(r in 1:nrep){
    # Data frame to save results of each repetition
    gpDF <- data.frame()
    for(f in 1:k){
      # I want to "mask" all the BLUEs in the testing rows
      Idx_test <- IdxDF
      Idx_test[Idx_test$genotype %in% vFolds[[r]][[f]], "BLUE"] <- NA
      
      Idx_test <- Idx_test[Idx_test$genotype %in% rownames(matG), ]
      matG <- matG[as.character(Idx_test$genotype), 
                   as.character(Idx_test$genotype)]
      
      GBLUP <- asreml(fixed = BLUE ~ 1,
                      # Variance structure of the genotypes
                      random = ~ vm(genotype, matG),
                      weights = weight,
                      residual = ~ idv(units),
                      data = Idx_test)
  
      # Predicted values
      predVals <- predict(GBLUP, classify = "genotype")$pvals
  
      # Filtering the predicted values for only those present in the
      # current fold
      predVals <- predVals[predVals$genotype %in% vFolds[[r]][[f]],]
  
      # Merge the predicted (GEBV) values to the original 
      # training dataset keeping only the rows relevant 
      # to the current fold
      predMerged <- merge(predVals, IdxDF[, c("genotype", "BLUE")], 
                          by = "genotype")
  
      # Naming the GEBV column accordingly
      colnames(predMerged)[2] <- "GEBV"
  
      # Data frame merging the proxy trait with the target trait
      # via the common genotypes, plus the relevant GEBVs and BLUEs
      predMerged <- merge(predMerged |> select(genotype, GEBV), 
                    target |> select(genotype, BLUE), 
                    by = "genotype")
      
      # Data frame with BLUEs and GEBVs for the current repetition
      gpDF <- rbind(gpDF, predMerged)
    }
    cors[r] <- cor(gpDF$GEBV, gpDF$BLUE)
  }
  
  # Returns the negative mean correlation (across all repetitions)
  # Negative because the optimization function minimizes the target
  return(-mean(cors))
}

