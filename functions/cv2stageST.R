# Function to perform cross-validation on the dataset 
# with the BLUEs (adjusted means)
# The arguments are the adjusted means dataset, the G matrix, and the
# number of folds k

# This function is for standard genomic selection, other versions will be
# developed later

cv2stageST <- function(dataset, matG, k, nrep){
  
  # Filter dataset for only genotypes present in the G matrix
  dataset <- dataset[dataset$genotype %in% rownames(matG), ]
  
  # Then filter G for only genotypes present in dataset
  # I wonder if this also orders the elements of G accordingly...
  Gfilt <- matG[as.character(dataset$genotype), 
                as.character(dataset$genotype)]
  
  # Getting genotypes
  genotypes <- dataset$genotype
    
  # Number of distinct genotypes in the dataset
  n <- length(genotypes)
  
  # We will perform 5-fold CV on the data
  # The CV will be performed 10 times
  
  # Vector of accuracies for each repetition
  accs <- numeric()
  
  for (j in 1:nrep){
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
      # Create a copy of the train_set for each iteration,
      # otherwise there will be problems (I don't want to
      # alter the original train dataset)
      trainFold <- dataset
      # Mask the BLUEs for genotypes present in the f-th validation
      # fold, f = 1, 2, ..., k, so they are absent from training the model
      trainFold[trainFold$genotype %in% valFolds[[f]], "BLUE"] <- NA
    
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
      predVals <- predVals[predVals$genotype %in% valFolds[[f]], ]
    
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
  
  # Let's build a small interval for the CV accuracies
  # Assuming normality of the accuracies and a 95% interval
  Z <- qnorm(0.025, lower.tail = F)
  IC <- mean(accs) + c(-1, 1) * sd(accs)/sqrt(nrep)
  
 return(IC)
}
