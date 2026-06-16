library(here)
library(tidyverse)
library(BGLR)

# The main idea here is to capture non-linear relationships between the
# variables

# I will use the BLUEs data as input, but this can be further discussed

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

n <- nrow(all_DF)

# Accuracy vector with the accuracies for each k-fold CV rep
accs <- numeric()

# Number of cross-validation folds
k <- 5

# Storing genotype information for CV
genotype <- all_DF$genotype

for (i in 1:10){
  # Mix and divide
  folds <- cut(seq(1, n), breaks = k, labels = FALSE)
  folds <- sample(folds)
  
  # Validation groups
  valFolds <- lapply(1:k, function(l) genotype[folds == l])
  
  # Data frame to save results of the CV for each rep
  results <- data.frame()
  
  for(f in 1:k){
  # We omit the response only (BGLR predicts NA responses by default)
  # In a real scenario, we would have all the proxy traits and the
  # full genomic information, needing to predict only the target trait
  # Reminder that we are emulating a genomic prediction indirect selection
  # scenario
    
  trainData <- all_DF
  trainData[trainData$genotype %in% valFolds[[f]], "FieldEmer"] <- NA
  
  fit <- BGLR(
    y = trainData$FieldEmer,
    ETA = list(
      proxies = list(X = 
                       trainData |> select(RagMeso, RagColeo, RagRoot,
                                         RagShoot), model = "BRR"),
      genomic = list(K = Gfilt, model = "RKHS")
    ),
    nIter = 1000,
    burnIn = 200,
    saveAt = ""
  )
  
  predVals <- as.data.frame(cbind(as.data.frame(genotype), fit$yHat))
  predVals <- predVals[predVals$genotype %in% valFolds[[f]], ]
  
  predMerged <- merge(predVals |> select(genotype, pred = `fit$yHat`), 
                      all_DF |> select(genotype, FieldEmer),
                      by = "genotype")
  
  results <- rbind(results, predMerged)
  }
  
  accs[i] <- cor(results$pred, results$FieldEmer)
}

load(file = here("output", "accs_List.RData"))
accs_List[["accsRKHS"]] <- accs
save(accs_List, file = here("output", "accs_List.RData"))










