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

# Loading Cullis heritability for predictive ability assessment
# Note: heritability related to field emergence only
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

################# Subpopulation Information ###################

# Keeping subpopulation information for possible future uses
# First we must map each genotype to a specific subpopulation

# We will consolidate the genotype-subpop rows from both data sources 
# into a single dataset
popAux1 <- expField |>
  select(genoID, subpop) |>
  distinct(genoID, .keep_all = TRUE)

popAux2 <- expRagdoll |>
  select(genoID, subpop) |>
  distinct(genoID, .keep_all = TRUE)

genoPopMap <- merge(popAux1 |> select(genoID), popAux2, by = "genoID")
rm(popAux1, popAux2)

genoPopMap <- genoPopMap |>
  rename(genotype = genoID) |>
  droplevels()

# Saving map to memory
save(genoPopMap, file = here("output", "genoPopMap.RData"))

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
accs_List <- vector(mode = "list", length = 6)
#names(accs_List) <- c("accField", "accMesoIS", "accMT_IS", 
                     #"accIdx", "accIdxFixed", "accIdx4t")

######################### Single-trait GP #####################

#------------ Field emergence (standard selection) ----------#

# -> This is our baseline model <-

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

GP_EmerField <- cv2stage(adjFieldEmerg, G, k = 5)

# Predictive ability as the ratio between the correlation of GEBVs and BLUEs
# and the Cullis heritability for emergence in the field
accFieldEmergence <- cor(GP_EmerField$GEBV, GP_EmerField$BLUE)/
  sqrt(h2CullisEmerField)

accs_List[["accField"]] <- accFieldEmergence

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

# Coleoptile will be used for two-trait GP later

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

GP_MesoRagdoll <- cv2stage(adjRagdollMeso, G, k = 5)

# To evaluate the prediction accuracy for the indirect selection approach,
# we will assess the correlation between the lab mesocotyl GEBVs and the field 
# emergence BLUEs. For that, we have to filter the genotypes so that only those
# common to both datasets are left
# Remember: our target trait is emergence, so we will divide by its heritability
# (in the field)

# Data frame with the common genotypes, plus the relevant GEBVs and BLUEs
GP_MesoRagdoll <- merge(GP_MesoRagdoll |> select(genotype, GEBV), 
              adjFieldEmerg |> select(genotype, BLUE), 
              by = "genotype")

# Calculating prediction accuracy:
# Note: the value indicates that Sandeep only calculated the correlation
# without dividing by Cullis heritability (he calculated predictive )
accIS_MESO <- cor(GP_MesoRagdoll$GEBV, GP_MesoRagdoll$BLUE)/
  sqrt(h2CullisEmerField) 

accs_List[["accMesoIS"]] <- accIS_MESO

######################### Multi-trait GP ########################

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

# Assessing the correlation between coleoptile and mesocotyl in the ragdoll
# experiment:
# To assess the correlation between each of the above two traits and emergence 
# in the field, let's consolidate the information in a single data frame
MT_IS <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                                wtMeso = weight),
                       adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                                 wtColeo = weight), 
                       by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                                wtEmerg = weight), 
        by = "genotype") |>
        droplevels()

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

# First of all, for multi-trait asreml, the data must be in long format to
# specify the weights correctly
# We first remove the columns associated with the field experiment
# Then pivot to longer format and finally pivot the result to wider format
# In the end, we want a column for the trait, another for the BLUEs, and another
# for the weights
longMT_IS <- MT_IS |>
            select(-c(FieldEmer, wtEmerg))

# Auxiliary data frames to help organize the process of pivoting to long format
# We will pivot each individually to long format, then merge them

# Data frame with BLUEs only
traitAux1 <- longMT_IS |>
            select(genotype, RagMeso, RagColeo)

traitAux1 <- traitAux1 |>
            pivot_longer(
              !genotype,
              names_to = "trait",
              # Captures only the "Meso" or "coleo" part
              names_pattern = "Rag(.*)",
              values_to = "BLUE"
            )

# Data frame with weights only
traitAux2 <- longMT_IS |>
  select(genotype, wtMeso, wtColeo)

traitAux2 <- traitAux2 |>
  pivot_longer(
    !genotype,
    names_to = "trait",
    # Captures only the "Meso" or "coleo" part
    names_pattern = "wt(.*)",
    values_to = "weight"
  )

# Merging the data frames back into a single one
longMT_IS <- merge(traitAux1, traitAux2, by = c("genotype", "trait"))
rm(traitAux1, traitAux2)

# Converting "trait" column to factor
longMT_IS <- longMT_IS |>
             mutate(trait = as.factor(trait))

# Arrange the data so that all rows with a given trait are followed by all
# rows with the other
longMT_IS <- longMT_IS |>
  arrange(trait)

# From this point onwards, it is useful to have a function 
# for cross-validation (CV)
# We can use longMT_IS and the G matrix as arguments for the function
GP_Meso_Coleo <- cv2stageMT(longMT_IS, G, k = 5)

# Data frame with the common genotypes, plus the relevant GEBVs and BLUEs
# FieldEmer represents the field emergence BLUEs
GP_Meso_Coleo <- merge(GP_Meso_Coleo |> select(genotype, GEBV_Meso), 
              MT_IS |> select(genotype, FieldEmer), 
              by = "genotype")

accIS_ML_CL <- cor(GP_Meso_Coleo$GEBV_Meso, GP_Meso_Coleo$FieldEmer)/
  sqrt(h2CullisEmerField) 

accs_List["accMT_IS"] <- accIS_ML_CL

###################### Index variable GP ########################

# We will use the same MT_IS dataset built for the multi-trait
# genomic prediction approach
# It contains the BLUEs for mesocotyl and coleoptile from the lab 
# experiment, and for emergence in the field experiment

# The criteria for choosing the best linear combination will be the
# correlation between the GEBVs of the trait and field emergence

# Remove columns related to the field experiment
# and build a dataset to be fed into the index building
# algorithm
preIdx <- MT_IS |>
          select(-c(FieldEmer, wtEmerg))

# Weights to be tried for the traits
wIdx <- seq(0, 1, by = 0.01)

# Standardize the trait columns in preIdx
# So their combination does not unfairly favor the one
# with larger variance solely due to scale
preIdx <- preIdx |>
          mutate_at(c("RagMeso", "RagColeo"), function(x) scale(x))

# The weight columns refer to the estimation errors when obtaining
# BLUEs, so they will be kept the same

# Vector to store prediction accuracies for each weight setup
accs <- numeric()

for (i in wIdx){
  # Index variable
  IdxVar <- i * preIdx$RagMeso + (1 - i) * preIdx$RagColeo
  
  # Index variable weight for GBLUP
  wtIdx <- 1/((i^2)/preIdx$wtMeso + ((1 - i)^2)/preIdx$wtColeo)
  
  # Generating data frame to be fed to cv2stage function:
  # The data frame must be in genotype-BLUE-weight format
  IdxDF <- data.frame(genotype = preIdx$genotype,
                      BLUE = IdxVar,
                      weight = wtIdx)
  
  # Performing GBLUP with the index variable as response
  IdxGBLUP <- cv2stage(IdxDF, G, k = 5)
  
  # Joining GBLUP dataset to field emergence BLUEs
  IS_Idx <- merge(IdxGBLUP |> select(genotype, GEBV), 
                adjFieldEmerg |> select(genotype, BLUE), 
                by = "genotype")
  
  # Calculating prediction accuracy for indirect selection
  # with index variable
  accIS_Idx <- cor(IS_Idx$GEBV, IS_Idx$BLUE)/
    sqrt(h2CullisEmerField) 
  
  accs <- c(accs, accIS_Idx)
}

# Obtaining best index weight combination
bestIdxWt <- wIdx[which.max(accs)] # 0.6 mesocotyl / 0.4 coleoptile
bestAcc <- max(accs)

accs_List[["accIdx"]] <- bestAcc

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

# Loading the new G matrix
load(file = here("output", "G_NoMajor.RData"))

# The approach will be to use the index variable selected
# above:

IdxVar <- bestIdxWt * preIdx$RagMeso + (1-bestIdxWt) * preIdx$RagColeo

wtIdx <- 1/((bestIdxWt^2)/preIdx$wtMeso + ((1-bestIdxWt)^2)/preIdx$wtColeo)

IdxDF <- data.frame(genotype = preIdx$genotype,
                    BLUE = IdxVar,
                    weight = wtIdx)

# Now we must add the column with the major SNP dosages across
# the genotypes

# Loading SNP dosage per genotype after pruning
load(file = here("output", "snpPruned.RData"))

# Column with only the major effect SNP
snpMajor <- snpPruned[, colnames(snpPruned) == "mlid0051837994"]

# Keeping the genotype information
snpMajor <- cbind(rownames(snpPruned), snpMajor)
colnames(snpMajor)[1] <- "genotype"
rownames(snpMajor) <- NULL

# Merging IdxDF to snpMajor
IdxMajor <- merge(IdxDF, snpMajor, by = "genotype")

# snpMajor column must be numeric
IdxMajor <- IdxMajor |>
            mutate(snpMajor = as.numeric(snpMajor))

# Performing GBLUP with the index variable as response
IdxGBLUP_Major <- cv2stageFixed(IdxMajor, G_NoMajor, k = 5)

# Joining GBLUP dataset to field emergence BLUEs
IS_IdxMajor <- merge(IdxGBLUP_Major |> select(genotype, GEBV), 
                adjFieldEmerg |> select(genotype, BLUE), 
                by = "genotype")

# Calculating prediction accuracy for indirect selection
# with index variable and major SNP fixed effect
accIS_IdxMajor <- cor(IS_IdxMajor$GEBV, IS_IdxMajor$BLUE)/
  sqrt(h2CullisEmerField) 

accs_List[["accIdxFixed"]] <- accIS_IdxMajor

# Accuracy remained basically the same...
# Repeating CV multiple times within the respective functions 
# would be ideal!!

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


