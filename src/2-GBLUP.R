library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(pls)
library(ggplot2)
library(glmnet)

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
accs_List <- vector(mode = "list")

######################### Single-trait GP #####################

#------------ Field emergence (standard selection) ----------#

# -> This is our baseline model <-

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

accEmerField <- cv2stageST(adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accField"]] <- accEmerField

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

# Coleoptile will be used for two-trait GP later

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

######################### Multi-trait GP ########################

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

accIS_ML_CL <- cv2stageMT_IS(adjRagdollMeso, adjRagdollColeo,
                          adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accMT_IS"]] <- accIS_ML_CL








#----------------------------------------------------------------
# SECTION SUBJECT TO BIG CHANGES

###################### Index variable GP ########################

# The criteria for choosing the best linear combination will be the
# correlation between the GEBVs of the trait and field emergence

# Returns the weight of mesocotyl in the index linear combination
# as well as the accuracy calculated via simple CV
Idx2t <- IdxSel2t(adjRagdollMeso, adjRagdollColeo,
                  adjFieldEmerg, G, k = 5)

# Now, with the weights, we can use single-trait GBLUP for indirect
# selection with the index variable

# Building data frame to be fed to single-trait indirect selection
# function
preIdx <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                     wtMeso = weight),
                  adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                     wtColeo = weight), 
                  by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                         wtEmerg = weight), 
        by = "genotype") |>
  droplevels()

# Remove columns related to the field experiment
# and build a dataset to be fed into the index building
# algorithm
preIdx <- preIdx |>
  select(-c(FieldEmer, wtEmerg))

# Standardize the trait columns in preIdx
# So their combination does not unfairly favor the one
# with larger variance solely due to scale
preIdx <- preIdx |>
  mutate_at(c("RagMeso", "RagColeo"), function(x) scale(x))

#-------------------------------------------------------------

w <- Idx2t$bestWt

# Index variable
IdxVar <- w * preIdx$RagMeso + (1 - w) * preIdx$RagColeo

# Index variable weight for GBLUP
wtIdx <- 1/((w^2)/preIdx$wtMeso + ((1 - w)^2)/preIdx$wtColeo)

# Generating data frame to be fed to cv2stage function:
# The data frame must be in genotype-BLUE-weight format
IdxDF <- data.frame(genotype = preIdx$genotype,
                    BLUE = IdxVar,
                    weight = wtIdx)
rm(IdxVar, wtIdx)

accIdx <- cv2stageST_IS(IdxDF, adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accIdx"]] <- accIdx

#------------------ Major QTL as fixed effect -------------------#

# Note: this section makes use of the index variable built/selected
# previously

# The following code relies on the "3-GWAS.R" script
# "mlid0051837994" was identified as a SNP with roughly
# 10% as the percentage of variance explained (PVE)
# Rex Bernardo's paper points out that QTLs with >=
# 10% PVE and traits with around 80% heritability
# provide minor improvements to the model


# The construction of the new G matrix without the 
# fixed effect SNP is performed in the "1-GenoAnalysis.R"
# script

# Index selection will be performed again because the best index
# without the major QTL as a fixed effectmight be different from 
# the one with it

# Loading the new G matrix
load(file = here("output", "G_NoMajor.RData"))

# Now we must add the column with the major SNP dosages across
# the genotypes

# Loading SNP dosage per genotype after pruning
load(file = here("output", "snpPruned.RData"))

Idx2tFixed <- IdxSel2tFixed(adjRagdollMeso, adjRagdollColeo,
                  adjFieldEmerg, snpPruned, G_NoMajor, k = 5)

wF <- Idx2tFixed$bestWt # Kinda funky, need to double-check

# Using the same preIdx dataset from before:

# Index variable
IdxVar <- wF * preIdx$RagMeso + (1 - wF) * preIdx$RagColeo

# Index variable weight for GBLUP
wtIdx <- 1/((wF^2)/preIdx$wtMeso + ((1 - wF)^2)/preIdx$wtColeo)

# Generating data frame to be fed to cv2stage function:
# The data frame must be in genotype-BLUE-weight format
IdxDF <- data.frame(genotype = preIdx$genotype,
                    BLUE = IdxVar,
                    weight = wtIdx)
rm(IdxVar, wtIdx)

accIdxMajor <- cv2stageST_IS_Fixed(IdxDF, adjFieldEmerg, snpPruned, 
                              G_NoMajor, k = 5, nrep = 10)

accs_List[["accIdxMajor"]] <- accIdxMajor
  
###############################################################
##                    Saving model accuracies                ##
###############################################################

save(accs_List, file = here("output", "accs_List.RData"))

# Plotting accuracies calculated so far
accs <- data.frame(
  Model = names(accs_List),
  Accuracy = unlist(accs_List)
)

# Bar chart
ggplot(accs, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Prediction Accuracies (Baseline + IS)", y = "Accuracy") +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = round(Accuracy, 2), vjust = -0.08))

# Clean plot later
# The best index so far was optimized in a grid search, not through
# a model...

