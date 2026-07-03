library(here)
library(dplyr)
library(tidyr)

# Setting a seed
set.seed(199927)

# G matrix:
load(here("output", "G.RData"))

# Experimental data (BLUEs)
# Loads lab proxy traits and field emergence
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)

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

# Matching dataset to G matrix's genotypes
IS_DF <- IS_DF[IS_DF$genotype %in% rownames(G), ]
Gfilt <- G[as.character(IS_DF$genotype), as.character(IS_DF$genotype)]

# Splitting the dataset into validation sets for 5-fold CV
k <- 5 # number of folds
nrep <- 10 # 5-fold CV reps
n <- nrow(IS_DF)

# Each element is assigned a number between 1 and k
folds <- cut(seq(1, n), breaks = k, labels = FALSE)

# List of folds for each repetition
# Each element of the list represents a single repetition
# of k-fold CV, and each element is itself a list of folds
valFolds <- vector(mode = "list")

genotypes <- IS_DF$genotype

# Each rep is a different k-fold CV split
for(r in 1:nrep){
  aux <- sample(folds) # different sample each time
  
  # Validation groups (each element in the list is a 5-fold split)
  valFolds[[r]] <- lapply(1:k, function(i) genotypes[aux == i])
}

# To validate all models the same way
# Not the same folds used to obtain the optimal index for
# index selection
save(valFolds, file = here("output", "valFolds.RData"))


