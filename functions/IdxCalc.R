# Function to be optimized to obtain best weights/coefficients
# for the index linear combination

# I think using CV on it would be ideal...
# And maximizing the mean...
IdxCalc <- function(coefs, prxy, target, matG, train_ind){
  
  traits <- prxy$traits
  weights <- prxy$weights # BLUE calculation weights, not
  # the coefficients we seek for the index
  
  # Vector where each element corresponds to the index for a genotype
  Idx <- as.matrix(traits) %*% coefs
  
  # Weight vector for the new index variable
  wtIdx <- 1/((1/as.matrix(weights)) %*% (coefs^2))
  
  # DF with the index values as "BLUEs", and the index weights as weights
  IdxDF <- data.frame(genotype = prxy$genotypes,
                      BLUE = Idx,
                      weight = wtIdx)
  
  # I want to "mask" all the BLUEs in the testing rows
  Idx_test <- IdxDF
  Idx_test[-train_ind, "BLUE"] <- NA
  
  GBLUP <- asreml(fixed = BLUE ~ 1,
                  # Variance structure of the genotypes
                  random = ~ vm(genotype, matG),
                  weights = weight,
                  residual = ~ idv(units),
                  data = Idx_test)
  
  # Predicted values
  predVals <- predict(GBLUP, classify = "genotype")$pvals
  
  # Filtering the predicted values for only those present in the
  # testing set
  predVals <- predVals[-train_ind, ]
  
  # Merge the predicted (GEBV) values to the original 
  # training dataset keeping only the rows relevant 
  # to the current fold
  predMerged <- merge(predVals, IdxDF[, c("genotype", "BLUE")], 
                      by = "genotype")
  
  # Naming the GEBV column accordingly
  colnames(predMerged)[2] <- "GEBV"
  
  # Data frame merging the proxy trait with the target trait
  # via the common genotypes, plus the relevant GEBVs and BLUEs
  gpDF <- merge(predMerged |> select(genotype, GEBV), 
                target |> select(genotype, BLUE), 
                by = "genotype")
  
  return(-cor(gpDF$GEBV, gpDF$BLUE))
}

