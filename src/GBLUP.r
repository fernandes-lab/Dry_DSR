library(here)
library(asreml)
library(dplyr)
library(tidyr)

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

adjFieldEmerg <- adjMeans(expFieldDeep, "emergence", blck = "block")
save(adjFieldEmerg, file = here("output", "adjFieldEmerg.RData"))
rm(expFieldDeep)

######## Ragdoll adjusted means ########

# Our two traits of interest are mesocotyl and coleoptile length
# Same procedure as for % emergence in the field experiment

adjRagdollMeso <- adjMeans(expRagdoll, "mesocotyl")
save(adjRagdollMeso, file = here("output", "adjRagdollMeso.RData"))

adjRagdollColeo <- adjMeans(expRagdoll, "coleoptile")
save(adjRagdollColeo, file = here("output", "adjRagdollColeo.RData"))
rm(expRagdoll)

############## Subpop Strats ############################

# First we must map each genotype to a specific subpopulation
# Loading the original experimental datasets with the subpopulation information:

load(here("output", "expField.Rdata"))
load(here("output", "expRagdoll.Rdata"))

# We will consolidate the genotype-subpop rows from both data sources 
# into a single dataset
popAux1 <- expField |>
  select(genoID, subpop) |>
  distinct(genoID, .keep_all = TRUE)
rm(expField)

popAux2 <- expRagdoll |>
  select(genoID, subpop) |>
  distinct(genoID, .keep_all = TRUE)
rm(expRagdoll)

genoPopMap <- merge(popAux1 |> select(genoID), popAux2, by = "genoID")
rm(popAux1, popAux2)

genoPopMap <- genoPopMap |>
  rename(genotype = genoID) |>
  droplevels()

# Saving map to memory
save(genoPopMap, file = here("output", "genoPopMap.RData"))

######################## GBLUP (second stage) ##################################

### Loading the G matrix
load(here("output", "G.RData"))

### Checking genetic variance structure:

sum(diag(G))/nrow(G)
# A value well above 1 points to a high level of inbreeding in the population
# which tracks because rice is predominantly self-pollinating

heatmap(G)
# We can see big correlation regions
# PCA would probably be interesting for this dataset
# especially when conducting GWAS

########## Single-trait GP ##########

#### Field emergence (standard selection)

load(file = here("output", "adjFieldEmerg.RData"))

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

rstlsEmerField <- cv2stage(adjFieldEmerg, G, k = 5)

# Loading Cullis heritability for predictive ability assessment
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

# Predictive ability as the ratio between the correlation of GEBVs and BLUEs
# and the Cullis heritability for emergence in the field
accFieldEmergence <- cor(rstlsEmerField$GEBV, rstlsEmerField$BLUE)/
  sqrt(h2CullisEmerField)                  

#### Ragdoll mesocotyl (indirect selection - IS)

# Coleoptile will be used for two-trait GP later

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

load(file = here("output", "adjRagdollMeso.RData"))

rstlsMesoRagdoll <- cv2stage(adjRagdollMeso, G, k = 5)

# To evaluate the prediction accuracy for the indirect selection approach,
# we will assess the correlation between the lab mesocotyl GEBVs and the field 
# emergence BLUEs. For that, we have to filter the genotypes so that only those
# common to both datasets are left
# Remember: our target trait is emergence, so we will divide by its heritability
# (in the field)

# Data frame with the common genotypes, plus the relevant GEBVs and BLUEs
ISdf <- merge(rstlsMesoRagdoll |> select(genotype, GEBV), 
              adjFieldEmerg |> select(genotype, BLUE), by = "genotype")

# Calculating prediction accuracy:
# Note: the value indicates that Sandeep only calculated the correlation
# without dividing by Cullis heritability (he calculated predictive )
accIS_MESO <- cor(ISdf$GEBV, ISdf$BLUE)/
  sqrt(h2CullisEmerField) 

########## Multi-trait GP ##########

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

load(file = here("output", "adjRagdollColeo.RData"))
load(file = here("output", "adjRagdollMeso.RData"))
load(file = here("output", "adjFieldEmerg.RData"))

# Assessing the correlation between coleoptile and mesocotyl in the ragdoll
# experiment:
# To assess the correlation between each of the above two traits and emergence 
# in the field, let's consolidate the information in a single data frame
MT_ISdf <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                                wtMeso = weight),
                       adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                                 wtColeo = weight), 
                       by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                                wtEmerg = weight), 
        by = "genotype") |>
        droplevels()

rm(adjFieldEmerg, adjRagdollMeso, adjRagdollColeo)

with(MT_ISdf, cor(RagMeso, FieldEmer)) # Mesocotyl w/ emergence
with(MT_ISdf, cor(RagColeo, FieldEmer)) # Coleoptile w/ emergence
with(MT_ISdf, cor(RagColeo, RagMeso)) # Coleoptile w/ mesocotyl

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

# Multi-trait GBLUP:

# First of all, for multi-trait asreml, the data must be in long format to
# specify the weights correctly
# We first remove the columns associated with the field experiment
# Then pivot to longer format and finally pivot the result to wider format
# In the end, we want a column for the trait, another for the BLUEs, and another
# for the weights
longMT_IS <- MT_ISdf |>
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

# Loading the G matrix again
load(here("output", "G.RData"))

# From this point onwards, it is useful to have a function 
# for cross-validation (CV)
# We can use longMT_IS and the G matrix as arguments for the function
CV_MTdf <- cv2stageMT(longMT_IS, G, k = 5)

# Loading Cullis heritability for predictive ability assessment
load(file = here("output", "cullisHeritField.RData"))

# Data frame with the common genotypes, plus the relevant GEBVs and BLUEs
# FieldEmer represents the field emergence BLUEs
MLCL_ISdf <- merge(CV_MTdf |> select(genotype, GEBV_Meso), 
              MT_ISdf |> select(genotype, FieldEmer), by = "genotype")

# Assessing prediction accuracy for mesocotyl + coleoptile (ragdoll) relative
# to field emergence, with mesocotyl as the primary proxy trait:
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

accIS_ML_CL <- cor(MLCL_ISdf$GEBV_Meso, MLCL_ISdf$FieldEmer)/
  sqrt(h2CullisEmerField) 

### TO DO: 
# Stratified sampling may be integrated into the original functions
# GWAS for major/minor QTL split

############# Index combining mesocotyl and coleoptile ################

# We will use the same MT_ISdf dataset built for the multi-trait
# genomic prediction approach
# It contains the BLUEs for mesocotyl and coleoptile from the lab 
# experiment, and for emergence in the field experiment

# The criteria for choosing the best linear combination will be the
# correlation between the GEBVs of the trait and field emergence

# Remove columns related to the field experiment
# and build a dataset to be fed into the index building
# algorithm
preIdx <- MT_ISdf |>
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

# Loading field emergence heritability for accuracy scaling
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

# Loading G matrix (again)
load(here("output", "G.RData"))

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

#######################################################
#                 Major QTL as fixed effect           #
#######################################################

# The following code relies on the GWAS.R script
# "mlid0051837994" was identified as a SNP with roughly
# 10% as the percentage of variance explained (PVE)
# Rex Bernardo's paper points out that QTLs with >=
# 10% PVE and traits with around 80% heritability
# provide minor improvements to the model
# However, I am unsure how well this translates to 
# an indirect selection scenario, as the GWAS was
# performed with respect to field emergence

# The construction of the new G matrix without the 
# fixed effect SNP is performed in the GenoAnalysis.R
# script

# Loading the new G matrix
load(file = here("output", "G_NoMajor.RData"))

# The approach will be to use the index variable selected
# above:

IdxVar <- 0.6 * preIdx$RagMeso + 0.4 * preIdx$RagColeo

wtIdx <- 1/((0.6^2)/preIdx$wtMeso + (0.4^2)/preIdx$wtColeo)

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
load(file = here("output", "adjFieldEmerg.RData"))
IS_IdxMajor <- merge(IdxGBLUP_Major |> select(genotype, GEBV), 
                adjFieldEmerg |> select(genotype, BLUE), 
                by = "genotype")

# Loading field emergence heritability for accuracy scaling
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

# Calculating prediction accuracy for indirect selection
# with index variable and major SNP fixed effect
accIS_IdxMajor <- cor(IS_IdxMajor$GEBV, IS_IdxMajor$BLUE)/
  sqrt(h2CullisEmerField) 

# Accuracy remained basically the same...