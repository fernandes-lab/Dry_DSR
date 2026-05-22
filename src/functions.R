# Function to obtain adjusted means assuming a fixed model structure and different
# response variables/datasets -> for field/ragdoll experiments

# Note for field data: 
# Mixed model assuming block has no effect (the genotypes are not carried from
# one block to another within a single replication, so estimating block effect
# is meaningless)
# Actually, the AIC and BIC decreased with the inclusion of block 
# in the field experiment, so I will leave block as an optional argument

adjMeans <- function(dataset, resp_name, blck = NULL){
  
  response <- dataset[[resp_name]]
  
  # If we want to include blocking effect (in the field experiment case)
  if(!is.null(blck)){
    block <- dataset[[blck]]
  
    model <- asreml(fixed = response ~ replication + genoID,
                    random = ~ block,
                    residual = ~ idv(units),
                    data = dataset)
  }else{
    model <- asreml(fixed = response ~ replication + genoID,
                    residual = ~ idv(units),
                    data = dataset)
  }
  
  # Obtaining BLUEs
  adjms <- predict(model, classify = "genoID")$pvals
  
  # Squared prediction error (?) matrix:
  # varcov structure of the estimated means across genotypes
  # we expects different genotypes to be uncorrelated in the absence
  # of genetic information
  # vcov <- predict(model, classify = "genoID", vcov = TRUE)
  # mat <- vcov$vcov
  
  # Converting mat to conventional numeric matrix format
  # mat <- as.matrix(mat)
  
  # Naming the columns and rows of the matrix according to the genotypes
  # dimnames(mat) <- list(adjms$genoID, adjms$genoID)
  
  # Heatmap to visually assess the correlation structure between genotypes
  # without any genomic input
  # heatmap(mat) # May be uncommented if need be
  
  # Per Piepho's paper, since we have an almost perfectly balanced design, in
  # a single environment, calculating weights by the inverse of the squared 
  # standard error of the genotype means is reasonable
  
  # Note: commented line results might be returned in future 
  # iterations of the function
  
  # Thus, weights are calculated as 1/SE^2:
  adjms$weight <- 1/(adjms$std.error^2)
  
  # Keeping only the information needed for GBLUP
  adjms <- adjms |>
    select(genoID, predicted.value, weight) |>
    rename(genotype = genoID, BLUE = predicted.value)
  
  return(adjms)
}

#------------------------------------------------------------------------------#
# Function to perform cross-validation on the dataset 
# with the BLUEs (adjusted means)
# The arguments are the adjusted means dataset, the G matrix, and the
# number of folds k

# This function is for standard genomic selection, other versions will be
# developed later

cv2stage <- function(dataset, matG, k){
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
    
    GBLUPmodel <- asreml(fixed = BLUE ~ 1,
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

#-------------------------------------------------------------
# Alternative version with "mlid0051837994" as a fixed effect: 

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

#------------------------------------------------------------------------------#
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




