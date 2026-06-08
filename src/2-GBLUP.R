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

accEmerField <- cv2stage(adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accField"]] <- accEmerField

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

# Coleoptile will be used for two-trait GP later

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

accMesoIS <- cv2stageIS(adjRagdollMeso, adjFieldEmerg, G, k = 5,
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

accIS_ML_CL <- cv2stageMT(adjRagdollMeso, adjRagdollColeo,
                          adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accMT_IS"]] <- accIS_ML_CL

###################### Index variable GP ########################

# The criteria for choosing the best linear combination will be the
# correlation between the GEBVs of the trait and field emergence
accIdx_ML_CL <- cv2stageIdx(adjRagdollMeso, adjRagdollColeo,
                           adjFieldEmerg, G, k = 5, nrep = 10)

accs_List[["accIdx"]] <- accIdx_ML_CL

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

# Accuracy remained basically the same...
# Repeating CV multiple times within the respective functions 
# would be ideal!!

accIdxMajor <- cv2stageIdxMajor(adjRagdollMeso, adjRagdollColeo,
                                adjFieldEmerg, snpPruned, G_NoMajor, 
                                k = 5, nrep = 10)

accs_List[["accIdxMajor"]] <- accIdxMajor

#--------- Adding shoot length and root length to index -------#

# We will use the preIdx dataset and merge it with adjRagdollRoot
# and adjRagdollShoot

expandPreIdx <- merge(preIdx, adjRagdollRoot |> 
            select(genotype, RagRoot = BLUE, wtRoot = weight), 
            by = "genotype") |>
            merge(adjRagdollShoot |> 
            select(genotype, RagShoot = BLUE, wtShoot = weight),
            by = "genotype") |>
            droplevels()

# I will save the above dataset because it may speed up future
# developments of this analysis
save(expandPreIdx, file = here("output", "expandPreIdx.RData"))

# Since now we have four traits, the max grid approach used earlier
# is not ideal. 
# One alternative approach is to regress emergence on the four variables
# and obtain the coefficients from that. 
# Standard LS is less computationally expensive than running GBLUP 
# in an loop, especially for four traits, so we will work with that
# for now
# For added simplicity, the BLUE weights won't be incorporated,
# which is reasonable given their relatively small ranges

# To account for the target trait prediction errors when building
# the BLUEs, we will use weighted least squares (WLS)

# Single dataset with the target field trait and the four ragdoll 
# traits

LS_Idx <- merge(expandPreIdx |>
                  select(genotype, meso = RagMeso, coleo = RagColeo,
                         root = RagRoot, shoot = RagShoot), 
                adjFieldEmerg |> 
                  select(genotype, emerg = BLUE, emrWt = weight),
                by = "genotype") |>
  droplevels()

# The proxy trait coefficients are obtaining through OLS
modCoefIdx <- lm(emerg ~ meso + coleo + root + shoot, 
                 weights = emrWt,
                 data = LS_Idx)

# plot(modCoefIdx) # Quick residual diagnostics

# Index variable vector (turned into a data frame column)
IdxVar4t <- as.matrix(LS_Idx |> select(-c(genotype, emerg, emrWt))) %*%
                      modCoefIdx$coefficients[-1] |>
                      as.data.frame()

# Calculating index weight vector from the 4 proxy traits weight
# vectors:

# Data frame with only the proxy traits' weights
wtDF <- expandPreIdx |> select(wtMeso, wtColeo, wtRoot, wtShoot)

# Converting wtDF to a matrix and multiplying it by the squared coefficient
# vector (excluding the intercept)
# The index weights will be the inverse of the above calculation
wtIdx4t <- 1 / (as.matrix(wtDF) %*%  (modCoefIdx$coefficients[-1]^2))

# Appending genotype names to index variable vector
IdxVar4t <- cbind(LS_Idx$genotype, IdxVar4t) |>
            rename(genotype = `LS_Idx$genotype`, index = V1)

# Adding index weights column to the IdxVar4t data frame
IdxVar4t <- cbind(IdxVar4t, wtIdx4t)

# Changing weight column name accordingly
IdxVar4t <- IdxVar4t |> rename(weight = wtIdx4t)

# Obtaining the GEBVs for the new index:
# Renaming the index column to "BLUE" to make it compatible
# with the cv2stage function, even though 
# that's not exactly what it means
IdxVar4t <- IdxVar4t |> rename(BLUE = index)
IdxGBLUP4t <- cv2stage(IdxVar4t, G, k = 5)

# Joining GBLUP dataset to field emergence BLUEs
IS_Idx4t <- merge(IdxGBLUP4t |> select(genotype, GEBV), 
                adjFieldEmerg |> select(genotype, BLUE), 
                by = "genotype")

# Calculating accuracy for the new, broader, index
accIS_Idx4t <- cor(IS_Idx4t$GEBV, IS_Idx4t$BLUE)/
  sqrt(h2CullisEmerField) 

# The meso + coleo index combination continues to yield the best
# model prediction accuracy (or is it predictive ability?)
accs_List[["accIdx4t"]] <- accIS_Idx4t

#----------------- Partial Least Squares Approach -------------#

# We will generate a component variable from the four lab variables
# and see if it acts as a good predictor for field emergence
# This approach is unlikely to bring any improvement to the accuracy
# because the proxy traits are not highly correlated

# Including the weights does not seem to affect the results much, 
# so I will perform PLS without accounting for weights (for now)

# LS_Idx will be reused as it contains the lab proxy traits BLUEs as
# well as the BLUEs for field emergence

# The function for PLSR natively performs cross-validation (CV)
# with 10 folds by default

modPLS4t <- plsr(emerg ~ meso + coleo + root + shoot, 
                 data = LS_Idx, scale = TRUE, validation = "CV")

summary(modPLS4t)
validationplot(modPLS4t)
# 1 component seems to be enough, although the % variance
# explained is not even 30%

# As an index variable, we will use the first component, building
# it from the coefficients returned by the model

# Index variable vector (turned into a data frame column)
Idx4tPLS <- as.matrix(LS_Idx |> select(-c(genotype, emerg, emrWt))) %*%
  modPLS4t[["coefficients"]][1:4] |>
  as.data.frame()

# Appending genotype names to index variable vector
Idx4tPLS <- cbind(LS_Idx$genotype, Idx4tPLS) |>
  rename(genotype = `LS_Idx$genotype`, index = V1)

# Renaming index column to be compatible with the GBLUP CV function
# and running the CV
Idx4tPLS <- Idx4tPLS |> rename(BLUE = index)
IdxGBLUP4tPLS <- cv2stage(Idx4tPLS, G, k = 5)

# Joining GBLUP dataset to field emergence BLUEs
IS_Idx4tPLS <- merge(IdxGBLUP4tPLS |> select(genotype, GEBV), 
                  adjFieldEmerg |> select(genotype, BLUE), 
                  by = "genotype")

# Calculating accuracy for the new, broader, index
accIS_Idx4tPLS <- cor(IS_Idx4tPLS$GEBV, IS_Idx4tPLS$BLUE)/
  sqrt(h2CullisEmerField)

accs_List[["accIdx4tPLS"]] <- accIS_Idx4tPLS

# The accuracy was really close to the index combination between
# only mesocotyl and coleoptile, the best among the IS scenarios
# so far 
# I wonder if somehow putting the weights into the model would
# improve it...
# Weighting actually caused a decrease in accuracy, so it won't be
# implemented

# ---------------- Ridge Regression Approach -----------------#

# The variables perceived as less relevant to the response will
# have their coefficients shrunk to near 0 
# We will figure out how to deal with weights later...

# Dataset with all four traits' BLUEs and weights, and the same
# for field emergence
RidgeDF <- merge(expandPreIdx |>
                  select(genotype, meso = RagMeso, wtMeso, 
                         coleo = RagColeo, wtColeo,
                         root = RagRoot, wtRoot,
                         shoot = RagShoot, wtShoot), 
                adjFieldEmerg |> 
                  select(genotype, emerg = BLUE, emrWt = weight),
                by = "genotype") |>
  droplevels()

# CV to find the best lambda value (the function natively performs
# 10-fold CV):

y <- RidgeDF$emerg
x <- data.matrix(RidgeDF |> select(meso, coleo, root, shoot))

cv_ridge <- cv.glmnet(x, y, alpha = 0)

# Finding best lambda by minimizing the mean squared error
# (MSE)
(best_lambda <- cv_ridge$lambda.min)

# To obtain the coefficients for the index, we will run a model
# with the best lambda

best_ridge <- glmnet(x, y, alpha = 0, lambda = best_lambda)
coef(best_ridge)

# Calculating the index 
IdxRidge <- as.matrix(RidgeDF|> select(-c(genotype, wtMeso,
                                          wtColeo, wtRoot,
                                          wtShoot, emerg,
                                          emrWt))) %*%
  coef(best_ridge)[-1] |>
  as.data.frame()

# Appending genotype names to index variable vector
IdxRidge <- cbind(RidgeDF$genotype, IdxRidge) |>
  rename(genotype = `RidgeDF$genotype`, index = V1)

# Renaming index column to be compatible with the GBLUP CV function
# and running the CV
IdxRidge <- IdxRidge |> rename(BLUE = index)
IdxGBLUPRidge <- cv2stage(IdxRidge, G, k = 5)

# Matching GEBVs and emergence BLUEs to the same genotypes
IdxRidgeConsol <- merge(IdxGBLUPRidge |> select(genotype, GEBV),
                        RidgeDF |> select(genotype, emerg),
                        by = "genotype")

# Calculating accuracy for the new index
accIS_IdxRidge <- cor(IdxRidgeConsol$GEBV, IdxRidgeConsol$emerg)/
  sqrt(h2CullisEmerField)

accs_List[["accIS_IdxRidge"]] <- accIS_IdxRidge

###############################################################
##                    Saving model accuracies                ##
###############################################################

save(accs_List, file = here("output", "modelAccs.RData"))

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

