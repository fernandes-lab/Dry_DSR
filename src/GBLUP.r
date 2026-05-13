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

# Calculating prediction accuract:
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

# Multi-trait CV:

# First of all, for multi-trait asreml, the data must be in long format to
# specify the weights correctly
# We first remove the columns associated with the field experiment
# Then pivot to longer format and finally pivot the result to wider format
# In the end, we want a column for the trait, another for the BLUEs, and another
# for the weights
longMT_ISdf <- MT_ISdf |>
               select(-c(FieldEmer, wtEmerg)) |>
               pivot_longer(
                 !genotype, 
                 names_to = c("typeOfValue", "trait"),
                 names_pattern = "(wt|Rag)(.*)",
                 values_to = "blueOrweight") |>
                 pivot_wider(
                   id_cols = c(genotype, trait),
                   names_from = typeOfValue,
                   values_from = blueOrweight
                 )

# Renaming the columns appropriately
longMT_ISdf <- longMT_ISdf |>
               select(genotype, Rag, wt, trait) |>
               rename(BLUE = Rag, weight = wt) |>
               mutate(trait = as.factor(trait))
               
# Ordering the data according to the traits (Mesocotyl and Coleoptile)
# All coleoptile values come before all mesocotyl values
longMT_ISdf <- longMT_ISdf |>
              arrange(trait)

genotypes <- MT_ISdf$genotype

# Number of distinct genotypes in the dataset
n <- length(genotypes)

## Create folds:
k = 5

# Fold assignment for each genotype
# Split the 1:n sequence into k folds and them randomly arranges them
# across the range of the dataset
folds <- sample(cut(1:n, breaks = k, labels = FALSE))

# Genotype validation folds
# Each member of the list is a subset of the genotypes column
# pertaining to that specific fold assignment
valFolds <- lapply(1:k, function(i) genotypes[folds == i])

# Data frame to store the results
gpDF <- data.frame()

# Loop over folds (basically loops over validation folds)
# It will usually be a 80/20 split, so 5 folds
# 4 for training, 1 for testing
for(f in 1:k){
  trainData <- MT_ISdf
  # Masks the BLUEs for genotypes present in the f-th validation
  # fold, f = 1, 2, ..., k, so they are absent from training the model
  # In this case, the BLUEs for both traits (mesocotyl and coleoptile)
  # are masked
  trainData[trainData$genotype %in% valFolds[[f]], 
            c("RagMeso", "RagColeo")] <- NA
  
  # Filter trainData for only genotypes present in the G matrix
  trainData <- trainData[trainData$genotype %in% rownames(G), ]
  
  # Then filter G for only genotypes present in trainData
  # I wonder if this also orders the elements of G accordingly...
  Gfilt <- G[as.character(trainData$genotype), as.character(trainData$genotype)]
  
  # Blueprint for CV in two-trait scenarios (example of multi-trait)
  MT_GBLUPmodel <- asreml(fixed = BLUE ~ trait,
                       random = ~ corgh(trait):vm(genotype, Gfilt),
                       weights = weight,
                       residual = ~ dsum(~ units | trait),
                       data = longMT_ISdf)
  
  # Predicted values
  predVals <- predict(MT_GBLUPmodel, classify = "trait:genotype")$pvals
  
  # Filtering the predicted values for only those present in the
  # (current) validation fold
  predVals <- predVals[predVals$genotype %in% valFolds[[f]], ]
  
  # Merge the predicted (GEBV) values to the original dataset
  # keeping only the rows relevant to the current fold
  predMerged <- merge(predVals, MT_ISdf[, c("genotype", "RagMeso", "RagColeo")], 
                      by = "genotype")
  
  # Naming the GEBV column accordingly
  colnames(predMerged)[2:3] <- c("GEBVmeso", "GEBVcoleo")
  
  # Append the rows with the GEBVs and BLUEs to the data frame storing
  # the results of the genomic prediction
  gpDF <- rbind(gpDF, predMerged)
  
}






