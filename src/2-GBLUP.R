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
accs_List <- vector(mode = "list")

# Experimental data (BLUEs)
# Loads lab proxy traits and field emergence
# Note: in case the above parts of the code were run in a different
# instance
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), 
       load, .GlobalEnv)

# Loading folds list for repeated (10 times) 5-fold CV
load(file = here("output", "valFolds.RData"))

######################### Single-trait GP ##########################

#------------ Field emergence (standard selection) ----------#

# -> This is our baseline model <-

# Calling function that performs CV and returns a data frame with the GEBVs
# and BLUEs

accEmerField <- cv2stageST(adjFieldEmerg, G, valFolds)

accs_List[["accField"]] <- accEmerField

#------------ Ragdoll mesocotyl (indirect selection - IS) -----#

# We will also have to eventually assess (for indirect selection)
# how the GEBVs in the ragdoll experiment correlate with the BLUEs for
# emergence in the field

accMesoIS <- cv2stageST_IS(adjRagdollMeso, adjFieldEmerg, G, 
                           valFolds)

# To evaluate the prediction accuracy for the indirect selection approach,
# we will assess the correlation between the lab mesocotyl GEBVs and the field 
# emergence BLUEs. For that, we have to filter the genotypes so that only those
# common to both datasets are left
# Remember: our target trait is emergence, so we will divide by its heritability
# (in the field)

accs_List[["accMesoIS"]] <- accMesoIS

#------------ Ragdoll coleoptile (indirect selection - IS) -----#

accColeoIS <- cv2stageST_IS(adjRagdollColeo, adjFieldEmerg, G, 
                            valFolds)

accs_List[["accColeoIS"]] <- accColeoIS

######################### Multi-trait GP ############################

## Basically a multi-trait indirect selection

# To improve indirect selection, we will do multi-trait prediction with ragdoll
# mesocotyl + coleoptile, mesocotyl being the primary trait

# Even though mesocotyl has higher heritability than coleoptile, we will use 
# mesocotyl as the primary trait (for CV) because its correlation with field
# emergence is higher than the coleoptile's correlation with field emergence

accIS_ML_CL <- cv2stageMT_IS(adjRagdollMeso, adjRagdollColeo,
                          adjFieldEmerg, G, valFolds)

accs_List[["accMT_IS"]] <- accIS_ML_CL

###############################################################
##                Assessing model accuracies                 ##
###############################################################

# save(accs_List, file = here("output", "accs_List.RData"))

# The assessment uses output from this source file, from 2-GBLUP_Idx,
# and from 2-RKHS.R

load(file = here("output", "accs_List.RData"))

# Testing the accuracy vectors for normality to check whether t-tests are
# viable
lapply(accs_List, shapiro.test)

lapply(accs_List, function(a) {
  qqnorm(a)
  qqline(a)
  })
# Normality seems like a reasonable assumption


# Model mean accuracies and standard deviations for these
accs <- data.frame(
  Model = names(accs_List),
  Accuracy = unlist(lapply(accs_List, mean)),
  StdDev = unlist(lapply(accs_List, sd))
)

# Bar chart
ggplot(accs, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Mean Accuracies", y = "Accuracy", x = "") +
  geom_errorbar(aes(ymin = Accuracy - StdDev, ymax = Accuracy + StdDev)) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = round(Accuracy, 2), vjust = -0.7)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# In order to compare the model accuracies, we will perform analysis of variance,
# and Tukey HSD for pairwise differences

accDF <- as.data.frame(accs_List)

# Adding a column to represent the 5-fold CV repetition that groups 
# the accuracy values across models
accDF$repID <- row.names(accDF)

# Converting data frame to long format to prepare for aov
accDF <- accDF |>
  pivot_longer(-repID, names_to = "model", values_to = "accuracy") |>
  mutate(model = as.factor(model), repID = as.factor(repID))

# Accuracy is modeled in terms of the model plus a fluctuation due to the
# specific 5-fold CV repetition
modAOV <- aov(accuracy ~ model + Error(repID), accDF)
summary(modAOV)

# Since "model" has a significant effect, there is evidence that at least 
# two models differ significantly

# Note: HSD corresponds to "honestly significant difference"
# For TukeyHSD, we fit a mixed model to the accuracy ~ model relationship, so
# repID can be treated as a random effect (the Error notation from aov has poor
# compatibility with TukeyHSD)


modLMM <- asreml(
  fixed = accuracy ~ model,
  random = ~ factor(repID), 
  data = accDF
)

pred_obj <- predictPlus(
  # Differences across models
  classify = "model",
  asreml.obj = modLMM,
  # Wald tests whether the differences are statistically equal to 0
  wald.tab = wald.asreml(modLMM, denDF = "numeric")$Wald
)

# Note: LSD stands for "least significant differences"
# Pairwise comparisons
pred_obj$differences    # LSD-based differences
pred_obj$LSD            # LSD threshold

# The assigned LSD is 0.006192499 (as seen in pred_obj$LSD), and all values
# in pred_obj$differences are larger than this value (in module), hence
# all the models differ significantly from each other

sum(abs(pred_obj$differences) < 0.006192499)
# There are 6 differences below the threshold, but they correspond to
# the models compared to themselves

# Now, to rank the models, we can add up how many negative differences are
# found in each row, which correspond to how many models perform better than
# the one represented in that given row (that + 1 is their rank)

rowSums(pred_obj$differences < 0) + 1 # model ranking
