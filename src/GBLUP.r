library(here)
library(asreml)
library(dplyr)

# At first, we will perform GBLUP without GWAS (to select major QTL)

source(here("src", "functions.R"))

# Setting a seed
set.seed(199927)

##################### Loading G matrix #########################################

load(here("output", "G.RData"))

###################### Loading experimental data ###############################

load(here("data", "expRagdoll.RData"))

load(here("data", "expField.RData"))

# Filtering field data for only "Deep" (> 8 cm) treatment
expFieldDeep <- expField |> filter(depth == "Deep")
rm(expField)

###################### Adjusted means (first stage) ############################

# Here genotype is a fixed effect since we are interested in each genotype's
# individual effect
# Note: Sandeep claimed spacial correction didn't influence the results 
# significantly, hence we can skip that

######## Field adjusted means ########

# Our only trait of interest in the field experiment is % emergence

# Mixed model assuming block has no effect (the genotypes are not carried from
# one block to another within a single replication, so estimating block effect
# is meaningless)

mFieldFix <- asreml(fixed = emergence ~ replication + genoID,
                    #random = ~ block,
                    residual = ~ idv(units),
                    data = expFieldDeep)

# Including block as a random effect barely changes the results
summary(mFieldFix)

# Obtaining BLUEs
adjField <- predict(mFieldFix, classify = "genoID")$pvals

# Prediction errors matrix:
# varcov structure of the estimated means across genotypes
# we expects different genotypes to be uncorrelated in the absence
# of genetic information
vcovField <- predict(mFieldFix, classify = "genoID", vcov = TRUE)
matField <- vcovField$vcov

# Converting matField to conventional numeric matrix format
matField <- as.matrix(matField)

# Naming the columns and rows of the matrix according to the genotypes
dimnames(matField) <- list(adjField$genoID, adjField$genoID)

# Heatmap to visually assess the correlation structure between genotypes
# without any genomic input
heatmap(matField) # Basically a straight line, implying no correlation

# Per Piepho's paper, since we have an almost perfectly balanced design, in
# a single environment, calculating weights by the inverse of the squared 
# standard error of the genotype means is reasonable

# Thus, weights are calculated as 1/SE^2:
adjField$weight <- 1/(adjField$std.error^2)

# Keeping only the information needed for GBLUP
adjField <- adjField |>
          select(genoID, predicted.value, weight) |>
          rename(genotype = genoID, BLUE = predicted.value)

save(adjField, file = here("output", "adjField.RData"))

######## Ragdoll adjusted means ########

# Our two traits of interest are mesocotyl and coleoptile length
# Same procedure as for % emergence in the field experiment

#### Mesocotyl model:

mRagdollFixMeso <- asreml(fixed = mesocotyl ~ replication + genoID,
                          residual = ~ idv(units),
                          data = expRagdoll)
summary(mRagdollFixMeso)

# Obtaining BLUEs
adjRagdollMeso <- predict(mRagdollFixMeso, classify = "genoID")$pvals

# Error varcov matrix
vcovRagdollMeso <- predict(mRagdollFixMeso, classify = "genoID", vcov = TRUE)
matRagdollMeso <- vcovRagdollMeso$vcov
matRagdollMeso <- as.matrix(matRagdollMeso)

dimnames(matRagdollMeso) <- list(adjRagdollMeso$genoID, adjRagdollMeso$genoID)

heatmap(matRagdollMeso) # no visually significant correlation structure

adjRagdollMeso$weight <- 1/(adjRagdollMeso$std.error^2)

adjRagdollMeso <- adjRagdollMeso |>
  select(genoID, predicted.value, weight) |>
  rename(genotype = genoID, BLUE = predicted.value)

save(adjRagdollMeso, file = here("output", "adjRagdollMeso.RData"))

#### Coleoptile model:

mRagdollFixColeo <- asreml(fixed = coleoptile ~ replication + genoID,
                          residual = ~ idv(units),
                          data = expRagdoll)
summary(mRagdollFixColeo)

# Obtaining BLUEs
adjRagdollColeo <- predict(mRagdollFixColeo, classify = "genoID")$pvals

# Error varcov matrix
vcovRagdollColeo <- predict(mRagdollFixColeo, classify = "genoID", vcov = TRUE)
matRagdollColeo <- vcovRagdollColeo$vcov
matRagdollColeo<- as.matrix(matRagdollColeo)

dimnames(matRagdollColeo) <- list(adjRagdollColeo$genoID, adjRagdollColeo$genoID)

heatmap(matRagdollColeo) # no visually significant correlation structure

adjRagdollColeo$weight <- 1/(adjRagdollColeo$std.error^2)

adjRagdollColeo <- adjRagdollColeo |>
  select(genoID, predicted.value, weight) |>
  rename(genotype = genoID, BLUE = predicted.value)

save(adjRagdollColeo, file = here("output", "adjRagdollColeo.RData"))

