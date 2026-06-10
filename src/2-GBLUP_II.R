library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(nloptr)

# Setting a seed
set.seed(199927)

# Sourcing index selection custom function
source(file = here("functions", "IdxCalc.R"))

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
G <- G[as.character(IS_DF$genotype), as.character(IS_DF$genotype)]

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

# Splitting the dataset into training/test sets
# The accuracy will be measured on the test set
# The test set will remain the same throughout the coefs optimization
# To ensure it's a consistent "hold-out"
# I will do a 70/30 split
# Sampling row indices:
n <- nrow(IS_DF)
trnInd <- sample(1:n, floor(0.70*n))
# Note: the above is done outside the function because the split is
# only done once

# Based on previous runs
# The coefficients/weights represent the relative contribution
# of each proxy trait to the index
starting_coefs <- c(0.5, 0.25, 0.125, 0.125)

# SLSQP function from nloptr library
coefOptim <- slsqp(
  x0 = starting_coefs, # starting at equal weights for each trait
  fn = IdxCalc,
  prxy = proxy,
  target = targetDF,
  matG = G, 
  train_ind = trnInd,
  lower = rep(0, 4), # lower bound for the coefficients
  upper = rep(1, 4),
  heq = function(w){sum(w)-1} # coefficients should add up to 1
)

# Obtaining the best weights and their corresponding accuracy
bestWt <- coefOptim$par
accBestWt <- -coefOptim$value


