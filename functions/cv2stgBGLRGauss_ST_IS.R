
cv2stgBGLRGauss_ST_IS <- function(dataset, tgtset, matG, vFolds, nIt, brnIn){
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
  
  # Building Gaussian kernel
  # We first need a distance matrix built from the G matrix
  # (Gfilt in this case)
  # D[i, j] =  G[i ,i] + G[j, j] -  2*G[i, j]
  
  # Vectorized calculation:
  n = nrow(Gfilt)
  D = outer(diag(Gfilt), rep(1, n)) + 
    outer(rep(1, n), diag(Gfilt)) - 2*Gfilt 
  
  #---------------------------------------
  # The Gaussian kernel is computed as:
  
  # K(i, j) = exp(-h * d(i, j)^2)
  
  # where h is the bandwidth parameter
  #---------------------------------------
  
  # We will use a grid of h values, all scaled by the median
  # distance in the D matrix (so the h values make sense within
  # the scale of the data)
  
  d_median = median(D)
  
  h_values = c(0.25, 0.5, 1, 2, 4)/d_median
  
  # Compute one kernel matrix per h value
  
  # The exponentiation is an element-wise operation!
  # List of Gaussian kernels (matrices)
  K_list <- lapply(h_values, function(h) exp(-h * D))
  
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
          kernel_1 = list(K = K_list[[1]], model = "RKHS"),
          kernel_2 = list(K = K_list[[2]], model = "RKHS"),
          kernel_3 = list(K = K_list[[3]], model = "RKHS"),
          kernel_4 = list(K = K_list[[4]], model = "RKHS"),
          kernel_5 = list(K = K_list[[5]], model = "RKHS")
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



