library(here)
library(tidyverse)
library(BGLR)

# In this script we will implement a Bayesian GBLUP much the same way
# as was implemented in 2-RKHS.R, however we will use a Gaussian instead
# of a linear kernel built around the G matrix.
# This approach might be able to capture non-linear effects

# Loading folds list for repeated (10 times) 5-fold CV
load(file = here("output", "valFolds.RData"))

# G matrix:
load(here("output", "G.RData"))

# BLUEs data:
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)
# Note: the weights derived from the BLUEs' prediction error will be
# used as features too

# Combined dataset with all features and the response
all_DF <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                         wtMeso = weight),
                adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                          wtColeo = weight), 
                by = "genotype") |> 
  merge(adjRagdollRoot |> select(genotype, RagRoot = BLUE, 
                                 wtRoot = weight),
        by = "genotype") |>
  merge(adjRagdollShoot |> select(genotype, RagShoot = BLUE, 
                                  wtShoot = weight),
        by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                                wtEmer = weight), 
        by = "genotype") |>
  droplevels()

rm(adjFieldEmerg, adjRagdollColeo, adjRagdollMeso, adjRagdollRoot,
   adjRagdollShoot)

# Now we make sure the genotypes in the comprehensive dataset match
# those present in the G matrix, and vice-versa

# Filter dataset for only genotypes present in the G matrix
# Safety step as it's unlikely that the dataset won't be fully
# represented in the G matrix
all_DF <- all_DF[all_DF$genotype %in% rownames(G), ]

# Then filter G for only genotypes present in dataset
# I wonder if this also orders the elements of G accordingly...
Gfilt <- G[as.character(all_DF$genotype), 
           as.character(all_DF$genotype)]
rm(G)

# -------------------------------------------------------

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

# The kernel used is a weighted average of all kernels relative
# to each h value, weighting each based on their respective
# variances

# Accuracy vector with the accuracies for each k-fold CV rep
accs <- numeric()

# Number of CV repetitions
nrep <- 10

# Number of CV folds per repetition
k <- 5

# Storing genotype information for CV
genotype <- all_DF$genotype

for (i in 1:nrep){
  
  # Data frame to save results of the CV for each rep
  results <- data.frame()
  
  for(f in 1:k){
    # We omit the response only (BGLR predicts NA responses by default)
    # In a real scenario, we would have all the proxy traits and the
    # full genomic information, needing to predict only the target trait
    # Reminder that we are emulating a genomic prediction indirect selection
    # scenario
    
    trainData <- all_DF
    trainData[trainData$genotype %in% valFolds[[i]][[f]], "FieldEmer"] <- NA
    
    fit <- BGLR(
      y = trainData$FieldEmer,
      ETA = list(
        proxies = list(X = 
                         trainData |> select(RagMeso, RagColeo, RagRoot,
                                             RagShoot), model = "BRR"),
        kernel_1 = list(K = K_list[[1]], model = "RKHS"),
        kernel_2 = list(K = K_list[[2]], model = "RKHS"),
        kernel_3 = list(K = K_list[[3]], model = "RKHS"),
        kernel_4 = list(K = K_list[[4]], model = "RKHS"),
        kernel_5 = list(K = K_list[[5]], model = "RKHS")
      ),
      nIter = 10000,
      burnIn = 2000,
      saveAt = ""
    )
    
    predVals <- as.data.frame(cbind(as.data.frame(genotype), fit$yHat))
    predVals <- predVals[predVals$genotype %in% valFolds[[i]][[f]], ]
    
    predMerged <- merge(predVals |> select(genotype, pred = `fit$yHat`), 
                        all_DF |> select(genotype, FieldEmer),
                        by = "genotype")
    
    results <- rbind(results, predMerged)
  }
  
  accs[i] <- cor(results$pred, results$FieldEmer)
}

save(accs, file = here("output", "accs_Gaussian.RData"))



