cv2stageIdx <- function(dataset1, dataset2, tgtset, matG, k, nrep){
  
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
  
  # Filter dataset for only genotypes present in the G matrix
  dataset0 <- dataset0[dataset0$genotype %in% rownames(matG), ]
  
  # Then filter G for only genotypes present in dataset
  # I wonder if this also orders the elements of G accordingly...
  Gfilt <- matG[as.character(dataset0$genotype), 
                as.character(dataset0$genotype)]
  
  # Getting genotypes
  genotypes <- dataset0$genotype
  
  # Number of distinct genotypes in the dataset
  n <- length(genotypes)
  
  # Remove columns related to the field experiment
  # and build a dataset to be fed into the index building
  # algorithm
  preIdx <- dataset0 |>
    select(-c(FieldEmer, wtEmerg))
  
  # Standardize the trait columns in preIdx
  # So their combination does not unfairly favor the one
  # with larger variance solely due to scale
  preIdx <- preIdx |>
    mutate_at(c("RagMeso", "RagColeo"), function(x) scale(x))
  
  # The weight columns refer to the estimation errors when obtaining
  # BLUEs, so they will be kept the same
  
  # List to store prediction accuracies for each weight setup
  # Each element is gonna be the accuracy vector
  # for one CV repetition
  accsWt <- vector(mode = "list")
  
  # Weights to be tried for the traits
  # Just a simple max grid approach
  wIdx <- seq(0, 1, by = 0.01)

  for (w in wIdx){
    # Index variable
    IdxVar <- w * preIdx$RagMeso + (1 - w) * preIdx$RagColeo
    
    # Index variable weight for GBLUP
    wtIdx <- 1/((w^2)/preIdx$wtMeso + ((1 - w)^2)/preIdx$wtColeo)
    
    # Generating data frame to be fed to cv2stage function:
    # The data frame must be in genotype-BLUE-weight format
    IdxDF <- data.frame(genotype = preIdx$genotype,
                        BLUE = IdxVar,
                        weight = wtIdx)
    
    # We will perform 5-fold CV on the data
    # The CV will be performed 10 times
    
    # Vector of accuracies for the repetitions
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
        trainFold <- IdxDF
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
        predMerged <- merge(predVals, IdxDF[, c("genotype", "BLUE")], 
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
    
    # One accs vector for each weight configuration
    # We need to update accsWt now with the accuracies
    # for this CV cycle
    # Each element index corresponds to the weight given to
    # mesocotyl
    accsWt[[as.character(w)]] <- accs
    
  }
  
  # Each element of accsWt has a vector containing the accuracies
  # for a given weight setup
  # We will choose the index corresponding 
  # to the highest mean accuracy
  # Obtaining best index weight combination
  avgAccWt <- lapply(accsWt, mean) # avgAccWt is a list of averages
  
  bestIdxWt <- wIdx[which.max(avgAccWt)]
  
  # Chosen vector of accuracies:
  bestAccs <- accsWt[[as.character(bestIdxWt)]]
  
  # Let's build a small interval for the CV accuracies
  # Assuming normality of the accuracies and a 95% interval
  Z <- qnorm(0.025, lower.tail = F)
  IC <- mean(bestAccs) + c(-1, 1) * sd(bestAccs)/sqrt(nrep)
  
  return(IC)
}





