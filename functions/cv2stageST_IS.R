# This function is a slightly different version of the original cv2stage
# created to accommodate an external dataset corresponding to the 
# target trait

cv2stageST_IS <- function(dataset, tgtset, matG, vFolds){
  
  # Filter dataset for only genotypes present in the G matrix
  dataset <- dataset[dataset$genotype %in% rownames(matG), ]
  
  # Then filter G for only genotypes present in dataset
  # I wonder if this also orders the elements of G accordingly...
  Gfilt <- matG[as.character(dataset$genotype), 
                as.character(dataset$genotype)]
  
  # CV parameters
  nrep <- length(vFolds)
  k <- length(vFolds[[1]]) # all equal length
  
  # Vector of accuracies for each repetition
  accs <- numeric()
  
  for (j in 1:nrep){
    # Data frame to store the results
    gpDF <- data.frame()
    
    # Loop over folds (basically loops over validation folds)
    # It will usually be a 80/20 split, so 5 folds
    # 4 for training, 1 for testing
    for(f in 1:k){
      # Create a copy of the train_set for each iteration,
      # otherwise there will be problems (I don't want to
      # alter the original train dataset)
      trainFold <- dataset
      # Mask the BLUEs for genotypes present in the f-th validation
      # fold, f = 1, 2, ..., k, so they are absent from training the model
      trainFold[trainFold$genotype %in% vFolds[[j]][[f]], "BLUE"] <- NA
      
      GBLUPmodel <- asreml(fixed = BLUE ~ 1,
                           # Variance structure of the genotypes
                           random = ~ vm(genotype, Gfilt),
                           weights = weight,
                           residual = ~ idv(units),
                           data = trainFold)
      
      # Predicted values
      predVals <- predict(GBLUPmodel, classify = "genotype")$pvals
      
      # Filtering the predicted values for only those present in the
      # (current) validation fold
      predVals <- predVals[predVals$genotype %in% vFolds[[j]][[f]], ]
      
      # Merge the predicted (GEBV) values to the original 
      # training dataset keeping only the rows relevant 
      # to the current fold
      predMerged <- merge(predVals, dataset[, c("genotype", "BLUE")], 
                          by = "genotype")
      
      # Naming the GEBV column accordingly
      colnames(predMerged)[2] <- "GEBV"
      
      # Append the rows with the GEBVs and BLUEs to the data frame storing
      # the results of the genomic prediction
      gpDF <- rbind(gpDF, predMerged)
      
    }
    # Data frame merging the proxy trait with the target trait
    # via the common genotypes, plus the relevant GEBVs and BLUEs
    gpDF <- merge(gpDF |> select(genotype, GEBV), 
                            tgtset |> select(genotype, BLUE), 
                            by = "genotype")
    
    accs[j] <- cor(gpDF$GEBV, gpDF$BLUE)
  }
  
  return(accs)
}
