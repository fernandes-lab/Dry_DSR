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

adjFieldEmerg <- adjMeans(expFieldDeep, "emergence")
save(adjFieldEmerg, file = here("output", "adjFieldEmerg.RData"))

######## Ragdoll adjusted means ########

# Our two traits of interest are mesocotyl and coleoptile length
# Same procedure as for % emergence in the field experiment

adjRagdollMeso <- adjMeans(expRagdoll, "mesocotyl")
save(adjRagdollMeso, file = here("output", "adjRagdollMeso.RData"))

adjRagdollColeo <- adjMeans(expRagdoll, "coleoptile")
save(adjRagdollColeo, file = here("output", "adjRagdollColeo.RData"))

