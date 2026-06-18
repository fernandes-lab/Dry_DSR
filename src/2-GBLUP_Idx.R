library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(nloptr)

# Setting a seed
set.seed(199927)

# Sourcing index selection custom function
# and single-trait IS function
source(file = here("functions", "IdxCalc.R"))
source(file = here("functions", "cv2stageST_IS.R"))

# G matrix:
load(here("output", "G.RData"))

# Experimental data (BLUEs)
# Loads lab proxy traits and field emergence
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)

# Single dataset with all four ragdoll experiment proxy traits
# and the target field emergence trait. The merging is basically
# to make sure the genotypes conform throughout all datasets
# Note: IS stands for indirect selection

IS_DF <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
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
                                wtEmerg = weight), 
        by = "genotype") |>
  droplevels()

# Matching dataset to G matrix' genotypes
IS_DF <- IS_DF[IS_DF$genotype %in% rownames(G), ]
Gfilt <- G[as.character(IS_DF$genotype), as.character(IS_DF$genotype)]

# Separate the target trait BLUE and weight column from the IS_DF dataset
targetDF <- IS_DF |> select(genotype, BLUE = FieldEmer, 
                            weight = wtEmerg)
IS_DF <- IS_DF |> select(-c(FieldEmer, wtEmerg))

# Note: merging into the unified dataset was done merely to ensure
# consistency across the different data sources

# Four lab proxy traits
proxyTraits <- IS_DF |> select(RagMeso, RagColeo, RagRoot, RagShoot)

# Proxy traits' corresponding BLUE weights (inverse of pred error)
proxyWeights <- IS_DF |> select(wtMeso, wtColeo, wtRoot, wtShoot)

proxy <- list(genotypes = IS_DF$genotype,
              traits = proxyTraits, 
              weights = proxyWeights)

# Splitting the dataset into validation sets for 5-fold CV
k <- 5 # number of folds
nrep <- 5 # 5-fold CV reps
n <- nrow(IS_DF)

# Each element is assigned a number between 1 and k
folds <- cut(seq(1, n), breaks = k, labels = FALSE)

# List of folds for each repetition
# Each element of the list represents a single repetition
# of k-fold CV, and each element is itself a list of folds
valFolds <- vector(mode = "list")

genotypes <- IS_DF$genotype

for(r in 1:nrep){
  aux <- sample(folds) # different sample each time
  
  # Validation groups (each element in the list is a 5-fold split)
  valFolds[[r]] <- lapply(1:k, function(i) genotypes[aux == i])
}

# Based on previous runs
# The coefficients/weights represent the relative contribution
# of each proxy trait to the index
starting_coefs <- c(0.5, 0.25, 0.125, 0.125)

# SLSQP function from nloptr library
coefOptim <- slsqp(
  x0 = starting_coefs, # starting at equal weights for each trait
  fn = IdxCalc, # sourced IdxCalc function
  prxy = proxy,
  target = targetDF,
  matG = Gfilt, 
  vFolds = valFolds,
  lower = rep(0, 4), # lower bound for the coefficients
  upper = rep(1, 4),
  heq = function(w){sum(w)-1} # coefficients should add up to 1
)

# Obtaining the best weights and their corresponding accuracy
bestWt <- coefOptim$par
accBestWt <- -coefOptim$value

# We can now build a data frame with the index variable and properly
# access the predictive ability of indirect selection (IS) using 
# repeated k-fold cross-validation (CV)
# The index can be treated as a single proxy variable and thus fed into
# the cv2stageST_IS function

# Vector where each element corresponds to the index for a genotype
Idx <- as.matrix(proxyTraits) %*% bestWt

# Weight vector for the new index variable
wtIdx <- 1/((1/as.matrix(proxyWeights)) %*% (bestWt^2))

# DF with the index values as "BLUEs", and the index weights as weights
IdxDF <- data.frame(genotype = IS_DF$genotype,
                    BLUE = Idx,
                    weight = wtIdx)

# Performing repeated 5-fold CV to assess the prediction accuracy
# with the selected index
accIdx <- cv2stageST_IS(IdxDF, adjFieldEmerg, G, k = 5, nrep = 10)

load(file = here("output", "accs_List.RData"))

accs_List[["accIdx"]] <- accIdx

save(accs_List, file = here("output", "accs_List.RData"))
