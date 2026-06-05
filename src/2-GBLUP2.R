library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(nloptr)

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

# Experimental data (BLUEs)
# Loads lab traits and field emergence
lapply(list.files(path = here("output"), 
                  pattern = "adj.*.RData", full.names = T), load, 
       .GlobalEnv)

# Loading Cullis heritability for predictive ability assessment
# Note: heritability related to field emergence only
load(file = here("output", "cullisHeritField.RData"))
h2CullisEmerField <- h2CullisField$emergence
rm(h2CullisField)

#--------- Index with equal weights for meso and coleo --------------#

# Following the same pipeline used for the index construction before
MT_IS <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                        wtMeso = weight),
               adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                         wtColeo = weight), 
               by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                                wtEmerg = weight), 
        by = "genotype") |>
  droplevels()

preIdx <- MT_IS |>
  select(-c(FieldEmer, wtEmerg))

IdxVar <- 0.5 * preIdx$RagMeso + 0.5 * preIdx$RagColeo

wtIdx <- 1/((0.5^2)/preIdx$wtMeso + (0.5^2)/preIdx$wtColeo)

preIdx <- preIdx |>
  mutate_at(c("RagMeso", "RagColeo"), function(x) scale(x))

IdxDF <- data.frame(genotype = preIdx$genotype,
                    BLUE = IdxVar,
                    weight = wtIdx)

IdxGBLUP <- cv2stage(IdxDF, G, k = 5)

IS_Idx <- merge(IdxGBLUP |> select(genotype, GEBV), 
                adjFieldEmerg |> select(genotype, BLUE), 
                by = "genotype")

# Accuracy is 0.62, smaller than the 0.6/0.4 split
accIS_Idx <- cor(IS_Idx$GEBV, IS_Idx$BLUE)/
  sqrt(h2CullisEmerField) 

#----------------------------------------------------------------#


##################################################################
##                   Optimization for four traits               ##
##################################################################

# I will try to implement an algorithm called sequential least squares
# quadratic programming (SLSQP) to search over the four lab traits
# - mesocotyl, coleoptile, root length, and shoot length - with the
# predictive ability being my objective function
# Basically, here goes nothing :)

blues4t <- merge(adjRagdollMeso |> select(genotype, RagMeso = BLUE, 
                                        wtMeso = weight),
               adjRagdollColeo |> select(genotype, RagColeo = BLUE,
                                         wtColeo = weight), 
               by = "genotype") |> # what comes b4 the pipe is the 1st
                                    # argument of the next merge
          merge(adjRagdollRoot |> select(genotype, RagRoot = BLUE, 
                                 wtRoot = weight),
              by = "genotype") |>
          merge(adjRagdollShoot |> select(genotype, RagShoot = BLUE, 
                                 wtShoot = weight),
              by = "genotype") |>
  merge(adjFieldEmerg |> select(genotype, FieldEmer = BLUE,
                                wtEmerg = weight), 
        by = "genotype") |>
  droplevels()

# Standardize the trait columns in blues_nt
# So their combination does not unfairly favor the one
# with larger variance solely due to scale
blues4t <- blues4t |>
  mutate_at(c("RagMeso", "RagColeo", "RagRoot", "RagShoot"), 
            function(x) scale(x))

# For the objective function, the code used for the meso + coleo index
# will be partially reused

objAcc <- function(coefs, blues, GM){
 
  # Four lab proxy traits
  proxyTraits <- blues |> select(RagMeso, RagColeo, RagRoot, RagShoot)
  
  Idx <- as.matrix(proxyTraits) %*% coefs
  
  # Proxy traits' corresponding BLUE weights (inverse of pred error)
  proxyWeights <- blues |> select(wtMeso, wtColeo, wtRoot, wtShoot)
  
  # Weight vector for the new index variable
  wt <- as.matrix(proxyWeights) %*% (coefs^2)
  
  IdxDF <- data.frame(genotype = blues$genotype,
                      BLUE = Idx,
                      weight = wt)
  
  # Performing GBLUP with the index variable as response
  IdxGBLUP <- cv2stage(IdxDF, GM, k = 5)
  
  # Joining GBLUP dataset to field emergence BLUEs
  emerg <- blues |> select(genotype, BLUE = FieldEmer)
  IS_Idx <- merge(IdxGBLUP |> select(genotype, GEBV), 
                  emerg, by = "genotype")
  
  # Calculating prediction accuracy for indirect selection
  # with index variable
  accIS_Idx <- cor(IS_Idx$GEBV, IS_Idx$BLUE)
  
  return(-accIS_Idx) # slsqp minimizes the obj function
                     # but I want the maximum accuracy
  
}

# SLSQP function from nloptr library
coefOptim <- slsqp(
  x0 = rep(0.25, 4), # starting at equal weights for each trait
  fn = objAcc,
  blues = blues4t,
  GM = G, 
  lower = rep(0, 4), # lower bound for the coefficients
  upper = rep(1, 4),
  heq = function(w){sum(w)-1} # coefficients should add up to 1
)

# Code ran in a little over 1 hour

# Assessing the new predictive ability:
coefOptim$value # returns the negative of what we actually want

# mesocotyl coleoptile root length shoot length
coefOptim$par # weights: 0.5157933 0.2345079 0.1248494 0.1248494

# Loading list of predictive abilities:
load(file = here("output", "modelAccs.RData"))

# Adding the newly calculated accuracy to the list
accs_List[["accIdx4tSLSQP"]] <- (-coefOptim$value)/h2CullisEmerField 
