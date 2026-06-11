# Multi-trait cross validation
# This function is limited to two traits

# For multi-trait asreml, the data must be in long format to
# specify the weights correctly
# We first remove the columns associated with the field experiment
# Then pivot to longer format and finally pivot the result to wider format
# In the end, we want a column for the trait, another for the BLUEs, and another
# for the weights

# Auxiliary data frames to help organize the process of pivoting to long format
# We will pivot each individually to long format, then merge them

cv2stageMT_IS <- function(dataset1, dataset2, tgtset, matG, k, nrep){
  
  # Dataset1 represents the first proxy trait (mesocotyl)
  # Dataset2 represents the second proxy trait (coleoptile)
  # tgtset represents the target trait (emergence)
  dataset0 <- merge(dataset1 |> select(genotype, RagMeso = BLUE, 
                                          wtMeso = weight),
                 dataset2 |> select(genotype, RagColeo = BLUE,
                                           wtColeo = weight), 
                 by = "genotype") |>
    merge(tgtset |> select(genotype, FieldEmer = BLUE,
                                  wtEmerg = weight), 
          by = "genotype") |>
    droplevels()
  
  # The merge with the target trait dataset was merely to ensure
  # consistency of genotypes across the datasets
  dataset <- dataset0 |>
    select(-c(FieldEmer, wtEmerg))
  
  # For multi-trait GBLUP with weights, the data must be in long format
  # We go through a few auxiliary steps to achieve that in a clean(er)
  # way
  
  # Data frame with BLUEs only
  traitAux1 <- dataset |>
    select(genotype, RagMeso, RagColeo)
  
  traitAux1 <- traitAux1 |>
    pivot_longer(
      !genotype,
      names_to = "trait",
      # Captures only the "Meso" or "Coleo" part
      names_pattern = "Rag(.*)",
      values_to = "BLUE"
    )
  
  # Data frame with weights only
  traitAux2 <- dataset |>
    select(genotype, wtMeso, wtColeo)
  
  traitAux2 <- traitAux2 |>
    pivot_longer(
      !genotype,
      names_to = "trait",
      # Captures only the "Meso" or "Coleo" part
      names_pattern = "wt(.*)",
      values_to = "weight"
    )
  
  # Final dataset in long format
  longData <- merge(traitAux1, traitAux2, by = c("genotype", "trait"))
  
  # Converting "trait" column to factor
  longData <- longData |>
    mutate(trait = as.factor(trait))
  
  # Arrange the data so that all rows with a given trait are followed by all
  # rows with the other
  longData <- longData |>
    arrange(trait)
  
  # Unique genotypes in the dataset
  genotypes <- unique(longData$genotype)
  
  # Number of distinct genotypes in the dataset
  n <- length(genotypes)
  
  # Vector of accuracies for each repetition
  accs <- numeric()
  
  for (j in 1:nrep){
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
    trainData <- longData
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
                        longData[,c("trait","genotype","BLUE")], 
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
  
  # Data frame merging the (main) proxy trait with the target trait
  # via the common genotypes, plus the relevant GEBVs and BLUEs
  gpDF <- merge(gpDF |> select(genotype, GEBV_Meso), 
                dataset0 |> select(genotype, FieldEmer), 
                by = "genotype")
  
  accs[j] <- cor(gpDF$GEBV_Meso, gpDF$FieldEmer)
  }
  
  return(accs)
}
