IdxSel2tFixed <- function(dataset1, dataset2, tgtset, snpData, matG, 
                        k){
  
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
  
  # Now we must add the column with the major SNP dosages across
  # the genotypes
  
  # Column with only the major effect SNP
  snpMajor <- snpData[, colnames(snpData) == "mlid0051837994"]
  
  # Keeping the genotype information
  snpMajor <- cbind(rownames(snpData), snpMajor)
  colnames(snpMajor)[1] <- "genotype"
  rownames(snpMajor) <- NULL
  
  # Vector to store prediction accuracies for each weight setup
  accsWt <- numeric()
  
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
    
    # Merging IdxDF to snpMajor
    IdxMajor <- merge(IdxDF, snpMajor, by = "genotype")
    
    # snpMajor column must be numeric
    IdxMajor <- IdxMajor |>
      mutate(snpMajor = as.numeric(snpMajor))
  
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
      trainFold <- IdxMajor
      # Mask the BLUEs for genotypes present in the f-th validation
      # fold, f = 1, 2, ..., k, so they are absent from training the model
      trainFold[trainFold$genotype %in% valFolds[[f]], "BLUE"] <- NA
      
      GBLUPmodel <- asreml(fixed = BLUE ~ snpMajor,
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
      predMerged <- merge(predVals, IdxMajor[, c("genotype", "BLUE")], 
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
    
    accsWt <- c(accsWt, cor(gpDF$GEBV, gpDF$BLUE))
  }
  
  # Each element of accsWt has a vector containing the accuracy
  # for a given weight setup
  # We will choose the index corresponding 
  # to the highest accuracy
  # Obtaining best index weight combination
  
  bestIdxWt <- wIdx[which.max(accsWt)] 
  bestAcc <- max(accsWt)
  
  return(list(bestWt = bestIdxWt, bestIdxAcc = bestAcc))
}





