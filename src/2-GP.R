library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(ggplot2)
library(asremlPlus)


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
# # TBD # #

# Experimental data (BLUEs)
# Loads lab proxy traits and field emergence
# Note: in case the above parts of the code were run in a different
# instance
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)

# Loading folds list for repeated (10 times) 5-fold CV
load(file = here("output", "valFolds.RData"))

# List to store accuracy values for Reml
accGP_Reml <- list()

# Same for linear kernel
accGP_LK <- list()

# Gaussian kernel list
accGP_GK <- list()

######################### Single-trait GP ##########################
# All models will have their respective Reml GBLUP + Linear kernel
# + Gaussian kernel approaches

#------------ Field emergence (standard selection) ----------#

# -> This is our baseline model <-

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

accField <- cv2stageST(adjFieldEmerg, G, valFolds)

accGP_Reml["accField"] <- accField

accFieldLK <- cv2stgBGLR_ST(adjFieldEmerg, G, valFolds, nIt = 10000, 
                            brnIn = 2000)

accGP_LK["accFieldLK"] <- accFieldLK

accFieldGK <- cv2stgBGLRGauss_ST(adjFieldEmerg, G, valFolds, 
                                        nIt = 10000, brnIn = 2000)

accGP_GK["accFieldGK"] <- accFieldGK

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

accMesoIS <- cv2stageST_IS(adjRagdollMeso, adjFieldEmerg, G, 
                           valFolds)

accGP_Reml["accMesoIS"] <- accMesoIS

accMesoLK_IS <- cv2stgBGLR_ST_IS(adjRagdollMeso,
                                      adjFieldEmerg, G, valFolds, nIt = 10000,
                                      brnIn = 2000)

accGP_LK["accMesoLK"] <- accMesoLK_IS


accMesoGK_IS <- cv2stgBGLRGauss_ST_IS(adjRagdollMeso,
                                                adjFieldEmerg, G, valFolds, nIt = 10000,
                                                brnIn = 2000)

accGP_GK["accMesoGK"] <- accMesoGK_IS


# To evaluate the prediction accuracy for the indirect selection approach,
# we will assess the correlation between the lab mesocotyl GEBVs and the field 
# emergence BLUEs. For that, we have to filter the genotypes so that only those
# common to both datasets are left
# Remember: our target trait is emergence

#------------ Ragdoll coleoptile (indirect selection - IS) -----#

accColeoIS <- cv2stageST_IS(adjRagdollColeo, adjFieldEmerg, G, 
                            valFolds)

accGP_Reml["accColeoIS"] <- accColeoIS

accColeoLK_IS <- cv2stgBGLR_ST_IS(adjRagdollColeo,
                                       adjFieldEmerg, G, valFolds, nIt = 10000,
                                       brnIn = 2000)

accGP_LK["accColeoLK"] <- accColeoLK_IS

accColeoGK_IS <- cv2stgBGLRGauss_ST_IS(adjRagdollColeo,
                                                 adjFieldEmerg, G, valFolds, nIt = 10000,
                                                 brnIn = 2000)

accGP_GK["accColeoGK"] <- accColeoGK_IS

######################### Multi-trait GP ############################

# Same as before, all models will have their respective Reml GBLUP 
# + Linear kernel + Gaussian kernel approaches

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

accIS_ML_CL <- cv2stageMT_IS(adjRagdollMeso, adjRagdollColeo,
                          adjFieldEmerg, G, valFolds)

accGP_Reml["accMT_IS"] <- accIS_ML_CL

# A bit of pre-processing for the kernel models
MT_DF <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE),
               adjRagdollColeo |> select(genotype, RagColeo = BLUE), 
               by = "genotype") |> 
  droplevels()

accLK_MT_IS <- cv2stgBGLR_MT_IS(MT_DF, adjFieldEmerg, 
                                  G, valFolds, nIt = 10000, 
                                  brnIn = 2000)

accGP_LK["accMT_LK"] <- accLK_MT_IS

accGK_MT_IS <- cv2stgBGLRGauss_MT_IS(MT_DF, adjFieldEmerg, 
                                            G, valFolds, nIt = 10000, 
                                            brnIn = 2000)

accGP_GK["accMT_GK"] <- accGK_MT_IS

# Saving accuracy lists
save(accGP_Reml, file = here("output", "accGP_Reml.RData"))
save(accGP_LK, file = here("output", "accGP_LK.RData"))
save(accGP_GK, file = here("output", "accGP_GK.RData"))



