library(ranger)
library(tidyverse)

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
try_DF <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
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
try_DF <- try_DF[RF_DF$genotype %in% rownames(G), ]

# Then filter G for only genotypes present in dataset
# I wonder if this also orders the elements of G accordingly...
Gfilt <- G[as.character(RF_DF$genotype), 
              as.character(RF_DF$genotype)]
rm(G)

#----------------------------------------------------------------
# Trying some things
# I think there is some overfitting going on right now
n <- nrow(try_set)

accs <- numeric()

# An attempt at repeated 5-fold CV
for (i in 1:10){
# Train/test split
# A single 80/20 split for now (repeated 5-fold CV to come later)

train_Ind <- sample(1:n, floor(0.80*n))

# We omit the response only (BGLR predicts NA responses by default)
# In a real scenario, we would have all the proxy traits and the
# full genomic information, needing to predict only the target trait
# Reminder that we are emulating a genomic prediction indirect selection
# scenario
yNA <- try_DF$FieldEmer
yNA[-train_Ind] <- NA

library(BGLR)

fit <- BGLR(
  y = yNA,
  ETA = list(
    proxies = list(X = 
                     try_DF |> select(RagMeso, RagColeo, RagRoot,
                                       RagShoot), model = "BRR"),
    genomic = list(K = Gfilt, model = "RKHS")
  ),
  nIter = 10000,
  burnIn = 2000,
  saveAt = ""
)

yHat_test <- fit$yHat[-train_Ind]
accs <- c(accs, cor(yHat_test, try_DF$FieldEmer[-train_Ind]))
}

save(accs, file = here("output", "accRKHS.RData"))
