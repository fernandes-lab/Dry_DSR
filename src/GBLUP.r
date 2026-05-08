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

#### Field emergence

load(file = here("output", "adjFieldEmerg.RData"))

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

rstlsEmerField <- cv2stage(adjFieldEmerg, G, k = 5)

# Loading Cullis heritability for predictive ability assessment
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence

# Predictive ability as the ratio between the correlation of GEBVs and BLUEs
# and the Cullis heritability for emergence in the field
PAfieldEmergence <- cor(rstlsEmerField$GEBV, rstlsEmerField$BLUE)/
  sqrt(h2CullisEmerField)                  

#### Ragdoll mesocotyl (will double as an indirect selection)

# Coleoptile will be used for two-trait GP later

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

load(file = here("output", "adjRagdollMeso.RData"))

rstlsMesoRagdoll <- cv2stage(adjRagdollMeso, G, k = 5)

# How do I go about indirect selection????

