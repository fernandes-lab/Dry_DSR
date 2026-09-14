# Function to perform cross-validation on the dataset 
# with the BLUEs (adjusted means)
# The arguments are the adjusted means dataset, the G matrix, and the
# number of folds k

# This function is for standard genomic selection, other versions will be
# developed later

cv2stageST <- function(dataset, matG, vFolds){
  
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
    # Data frame to store the results for each repetition
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
    accs[j] <- cor(gpDF$GEBV, gpDF$BLUE)
  }
  
 return(accs)
}
