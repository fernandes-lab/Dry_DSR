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

############## Stratified Cross-Validation Setup ############################

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

# Stratified CV approach:
load(file = here("output", "genoPopMap.RData"))
rstlsEmerField <- str_cv2stage(adjFieldEmerg, G, k = 5, genoPopMap)

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

# Stratified CV approach:
load(file = here("output", "genoPopMap.RData"))
rstlsMesoRagdoll <- str_cv2stage(adjRagdollMeso, G, k = 5, genoPopMap)

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

# Stratified CV approach:
CV_MTdf <- str_cv2stageMT(longMT_IS, G, k = 5, genoPopMap)

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



