library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(nloptr)
library(ggplot2)


# Loading functions from the functions folder
sapply(list.files(path = here("functions"), 
                  pattern = "\\.R$", full.names = T), source)

# Setting a seed
set.seed(199927)

###############################################################
##                Loading experimental data                  ##
###############################################################

# G matrix:
load(here("output", "G.RData"))

# Ragdoll experiment
load(here("data", "expRagdoll.RData"))

# Field experiment
load(here("data", "expField.RData"))

# Filtering field data for only "Deep" (> 8 cm) treatment
expFieldDeep <- expField |> filter(depth == "Deep")

###############################################################
##                Adjusted means (first stage)               ##
###############################################################

# Here genotype is a fixed effect since we are interested in each genotype's
# individual effect
# Note: Sandeep claimed spacial correction didn't influence the results 
# significantly, hence we can skip that

#################### Field adjusted means ######################

# Our only trait of interest in the field experiment is % emergence

# Mixed model assuming block has no effect (the genotypes are not carried from
# one block to another within a single replication, so estimating block effect
# is meaningless)

adjFieldEmerg <- adjMeans(expFieldDeep, "emergence", blck = "block")
# save(adjFieldEmerg, file = here("output", "adjFieldEmerg.RData"))

##################### Ragdoll adjusted means ##################

# Our two traits of interest are mesocotyl and coleoptile length
# Same procedure as for % emergence in the field experiment

adjRagdollMeso <- adjMeans(expRagdoll, "mesocotyl")
# save(adjRagdollMeso, file = here("output", "adjRagdollMeso.RData"))

adjRagdollColeo <- adjMeans(expRagdoll, "coleoptile")
# save(adjRagdollColeo, file = here("output", "adjRagdollColeo.RData"))

# For the sake of building an index combining all responses
# collected in the ragdoll experiment, we will also obtain
# adjusted means for root length and shoot length

adjRagdollRoot <- adjMeans(expRagdoll, "rootlength")
# save(adjRagdollRoot, file = here("output", "adjRagdollRoot.RData"))

adjRagdollShoot <- adjMeans(expRagdoll, "shootlength")
# save(adjRagdollShoot, file = here("output", "adjRagdollShoot.RData"))

###############################################################
##                     GBLUP (second stage)                  ##
###############################################################

### Checking genetic variance structure:
# sum(diag(G))/nrow(G)

# A value well above 1 points to a high level of inbreeding in the population
# which tracks because rice is predominantly self-pollinating

### Heatmap illustrating the genetic covariance structure:
# heatmap(G)
# We can see big correlation regions
# PCA would probably be interesting for this dataset
# especially when conducting GWAS

# List to store prediction accuracies for each modeling approach
# Each element of the list is itself a list of accuracies,
# one for each repetition of k-fold CV
accs_List <- vector(mode = "list")

# Experimental data (BLUEs)
# Loads lab proxy traits and field emergence
# Note: in case the above parts of the code were run in a different
# Instance
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)

######################### Single-trait GP ##########################

#------------ Field emergence (standard selection) ----------#

# -> This is our baseline model <-

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

accEmerField <- cv2stageST(adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accField"]] <- accEmerField

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

accMesoIS <- cv2stageST_IS(adjRagdollMeso, adjFieldEmerg, G, k = 5,
                        nrep = 10)

# To evaluate the prediction accuracy for the indirect selection approach,
# we will assess the correlation between the lab mesocotyl GEBVs and the field 
# emergence BLUEs. For that, we have to filter the genotypes so that only those
# common to both datasets are left
# Remember: our target trait is emergence, so we will divide by its heritability
# (in the field)

accs_List[["accMesoIS"]] <- accMesoIS

#------------ Ragdoll coleoptile (indirect selection - IS) -----#

accColeoIS <- cv2stageST_IS(adjRagdollColeo, adjFieldEmerg, G, k = 5,
                           nrep = 10)

accs_List[["accColeoIS"]] <- accColeoIS

######################### Multi-trait GP ############################

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

accIS_ML_CL <- cv2stageMT_IS(adjRagdollMeso, adjRagdollColeo,
                          adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accMT_IS"]] <- accIS_ML_CL

######################### Index GP ##############################

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
  fn = IdxCalc, # sourced IdxCalc function
  prxy = proxy,
  target = targetDF,
  matG = Gfilt, 
  train_ind = trnInd,
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

accs_List[["accIdx"]] <- accIdx

###############################################################
##                Assessing model accuracies                 ##
###############################################################

save(accs_List, file = here("output", "accs_List.RData"))

# Plotting accuracies calculated so far
accs <- data.frame(
  Model = names(accs_List),
  Accuracy = unlist(lapply(accs_List, mean))
)

# Bar chart
ggplot(accs, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Prediction Accuracies (Baseline + IS)", y = "Accuracy") +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = round(Accuracy, 2), vjust = -0.08))

# Wilcox test to compare vectors of accuracies (alternative hypothesis:
# first vector is <greater or less> than second vector):

# Comparing non-index approaches to index approach
with(accs_List, wilcox.test(accMesoIS, accIdx, alternative = "less"))
with(accs_List, wilcox.test(accColeoIS, accIdx, alternative = "less"))
with(accs_List, wilcox.test(accMT_IS, accIdx, alternative = "less"))

