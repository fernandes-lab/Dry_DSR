
cv2stgBGLR_ST_IS <- function(dataset, tgtset, matG, vFolds, nIt, brnIn){
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
  
  # Storing genotype information for CV
  genotype <- dataset$genotype
  
  for (i in 1:nrep){
    
    # Data frame to save results of the CV for each rep
    results <- data.frame()
    
    for(f in 1:k){
      # The BLUEs here are for the index variable obtained previously
      
      trainData <- dataset
      trainData[trainData$genotype %in% valFolds[[i]][[f]], "BLUE"] <- NA
      
      fit <- BGLR(
        y = trainData$BLUE,
        ETA = list(
          kernel = list(K = Gfilt, model = "RKHS")
        ),
        nIter = nIt,
        burnIn = brnIn,
        saveAt = ""
      )
      
      predVals <- as.data.frame(cbind(as.data.frame(genotype), fit$yHat))
      predVals <- predVals[predVals$genotype %in% valFolds[[i]][[f]], ]
      
      predMerged <- merge(predVals |> select(genotype, pred = `fit$yHat`), 
                          dataset |> select(genotype, BLUE),
                          # Keeping the BLUE here is not really necessary
                          # As we care about the target trait BLUE only
                          by = "genotype")
      
      results <- rbind(results, predMerged)
    }
    
    # Data frame merging the proxy trait with the target trait
    # via the common genotypes, plus the relevant GEBVs and BLUEs
    gpDF <- merge(results |> select(genotype, pred), 
                  tgtset |> select(genotype, BLUE), 
                  by = "genotype")
    
    accs[i] <- cor(gpDF$pred, gpDF$BLUE)
    
  }
  return(accs)
}







