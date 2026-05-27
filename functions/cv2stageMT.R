# Multi-trait cross validation

cv2stageMT <- function(dataset, matG, k){
  # Unique genotypes in the dataset
  genotypes <- unique(dataset$genotype)
  
  # Number of distinct genotypes in the dataset
  n <- length(genotypes)
  
  ## Create folds:
  
  # Fold assignment for each genotype
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
    # This should mask the BLUEs for both traits
    trainData[trainData$genotype %in% valFolds[[f]], "BLUE"] <- NA
    
    # Filter trainData for only genotypes present in the G matrix
    trainData <- trainData[trainData$genotype %in% rownames(matG), ]
    
    # Filter G according to the genotypes found in the training dataset
    Gfilt <- G[rownames(G) %in% trainData$genotype,
               colnames(G) %in% trainData$genotype]
    
    # Bivariate GBLUP model
    MT_GBLUPmodel <- asreml(fixed = BLUE ~ trait,
                            random = ~ corgh(trait):vm(genotype, Gfilt),
                            weights = weight,
                            residual = ~ dsum(~ units | trait),
                            data = trainData)
    
    # Predict for each trait "within" a genotype
    predVals <- predict(MT_GBLUPmodel, classify = "genotype:trait")$pvals
    
    # Filtering the predicted values for only those present in the
    # (current) validation fold
    predVals <- predVals[predVals$genotype %in% valFolds[[f]], ]
    
    # Stack one trait on top of the other in the predVals dataframe
    predVals <- predVals |>
      arrange(trait)
    
    # Merge predVals with dataset with the BLUEs
    predMerged <- merge(predVals[,c("genotype","trait", "predicted.value")], 
                        longMT_IS[,c("trait","genotype","BLUE")], 
                        by=c("genotype", "trait"))
    
    # Renaming predicted values column to GEBVs
    predMerged <- predMerged |>
      rename(GEBV = predicted.value)
    
    # Pivoting predMerged to wide format for better understading
    predMerged <- predMerged |>
      pivot_wider(
        id_cols = genotype,
        names_from = trait, 
        values_from = c(GEBV, BLUE)
      )
    
    # Append the rows with the GEBVs and BLUEs to the data frame storing
    # the results of the genomic prediction
    gpDF <- rbind(gpDF, predMerged)
    
  }
  
  # Return the "genomic prediction" data frame with GEBVs and BLUEs
  return(gpDF)
}
