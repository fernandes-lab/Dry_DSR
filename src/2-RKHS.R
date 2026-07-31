library(here)
library(tidyverse)
library(BGLR)

# The main idea here is to capture non-linear relationships between the
# variables

# I will use the BLUEs data as input, but this can be further discussed

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

#----------------------------------------------------------------

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
      genomic = list(K = Gfilt, model = "RKHS")
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


# With the proxy traits only, this approach performed similarly to
# the index GBLUP approach, whereas with the Gaussian G kernel alone it
# performed almost as well as with both the kernel and the proxies, notably
# even better than the baseline model (should be further scrutinized...)


load(file = here("output", "accs_List.RData"))
accs_List[["accRKHS"]] <- accs
save(accs_List, file = here("output", "accs_List.RData"))










