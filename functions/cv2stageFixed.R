# Alternative version of CV with "mlid0051837994" as a fixed effect: 

cv2stageFixed <- function(dataset, matG, k){
  genotypes <- dataset$genotype
  
  # Number of distinct genotypes in the dataset
  n <- length(genotypes)
  
  ## Create folds:
  
  # Fold (index) assignment for each genotype
  # Split the 1:n sequence into k folds and them randomly arranges them
  # across the range of the dataset
  folds <- sample(cut(1:n, breaks = k, labels = FALSE))
  
  # Genotype validation folds
  # Each member of the list is a subset of the genotypes column
  # pertaining to that specific fold assignment
  valFolds <- lapply(1:k, function(i) genotypes[folds == i])
  
  # Data frame to store the results
  gpDF <- data.frame()
  
  # Loop over folds (basically loops over validation folds)
  # It will usually be a 80/20 split, so 5 folds
  # 4 for training, 1 for testing
  for(f in 1:k){
    trainData <- dataset
    # Masks the BLUEs for genotypes present in the f-th validation
    # fold, f = 1, 2, ..., k, so they are absent from training the model
    trainData[trainData$genotype %in% valFolds[[f]], "BLUE"] <- NA
    
    
    # Filter trainData for only genotypes present in the G matrix
    trainData <- trainData[trainData$genotype %in% rownames(matG), ]
    
    # Then filter G for only genotypes present in trainData
    # I wonder if this also orders the elements of G accordingly...
    Gfilt <- matG[as.character(trainData$genotype), as.character(trainData$genotype)]
    
    # Simple modeling structure (for now)
    
    GBLUPmodel <- asreml(fixed = BLUE ~ snpMajor,
                         # Variance structure of the genotypes
                         random = ~ vm(genotype, Gfilt),
                         weights = weight,
                         residual = ~ idv(units),
                         data = trainData)
    
    # Predicted values
    predVals <- predict(GBLUPmodel, classify = "genotype")$pvals
    
    # Filtering the predicted values for only those present in the
    # (current) validation fold
    predVals <- predVals[predVals$genotype %in% valFolds[[f]], ]
    
    # Merge the predicted (GEBV) values to the original dataset
    # keeping only the rows relevant to the current fold
    predMerged <- merge(predVals, dataset[, c("genotype", "BLUE")], 
                        by = "genotype")
    
    # Naming the GEBV column accordingly
    colnames(predMerged)[2] <- "GEBV"
    
    # Append the rows with the GEBVs and BLUEs to the data frame storing
    # the results of the genomic prediction
    gpDF <- rbind(gpDF, predMerged)
    
  }
  
  # Return the "genomic prediction" data frame with GEBVs and BLUEs
  return(gpDF)
}
