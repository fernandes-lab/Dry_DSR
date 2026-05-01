library(here)
library(dplyr)
library(stringr)
library(readxl)
library(asreml)
library(ggplot2)
library(tidyr)

#-------------------------------------------------------------------------------------------------------#
#### Basic exploratory analysis and cleaning of the data ####

# Field trial raw data
exp_field <- read.csv(here("data/sandeepOnly/rawData", "Compiled_RDP1_fieldTrial_data_v3.0.csv"))

# Checking dataset
dim(exp_field)
str(exp_field)
summary(exp_field)
# 2515 observations of 34 variables (in total)

# Checking which of the "metadata" columns have only 1 factor 
# thus being liable for exclusion
exp_field |>
        summarise(across(everything(), ~ n_distinct(.)))

# We can see that SiteYear, SiteName and Sown have only one unique value, so we can remove them
expField <- exp_field |>
            select(-c(SiteYear, SiteName, Sown))

# TrialType and TrialNumber seem to represent the same kind of information, so let's verify that
sum(with(expField, substr(TrialType, 1, 8) == TrialNumber))
# They are equivalent on all 2515 observations, 
# so I will consolidate them in a single column with different category names (deep and shallow)
expField <- expField |>
            select(-TrialNumber) |>
            rename(depth = TrialType) |>  
            mutate(depth = recode(depth, "A19DS-02_Stress" = "Deep", "A19DS-01_control" = "Shallow"))

# Removing columns that are immediately known to not be relevant to our analysis
expField <- expField |>
            select(-c(Day_5:emergence_count_days, Shoot.dry.wt..g.:TrialUnitComment))

# Renaming columns to simpler names (replicate becomes replication to stay the same as ragdoll)
# and setting all of them to lower case
# besides converting the data types to relevant ones, if necessary
expField <- expField |> 
            rename(plotnum = PlotBarcode, mesocotyl = mesocotyl_length, rootlength = root_length, shootlength = shoot_length, seedcount = Seed_count, replication  = Replicate) |>
            mutate(across(c("depth", "block", "Row", "replication", "Genotype", "Pedigree"), as.factor)) |>
            rename_with(tolower)
str(expField)
summary(expField)

# Further renaming for better mutual understanding with genoInfo:
expField <- expField |>
            rename(genonumber = genotype, genotype = pedigree)

# Saving before any major changes
save(expField, file = here("data", "expField.Rdata"))

# Checking how many times each genotype appears in each replication
with(expField, table(genotype, replication) > 0) |>
        rowSums() |>
        table()
# 503 genotypes appear in 3 replications and 1 appears in 2 (One genotype was replicated only once)
# Bear in mind that the deep sowing treatment had 3 replications while the shallow sowing one had 2

# Checking for NAs
expField |>
    is.na() |>
    colSums()
# Since there are NAs only in the response variables, we can leave them be because
# the genotypic matrix is going to share information across different genotypes,
# and we know it has no NA

# Checking number of blocks per replication
with(expField, table(block, replication, depth))
# We can see that roughly all the genotypes appear in each rep, being equally split across
# three incomplete blocks in each rep. Furthermore, the genotypes are planted "deeply" in
# all 3 reps, whereas the shallow treatment is only applied in the first 2 reps


# Checking distribution of variables
png(here("output/plots","FieldExpDists.png"), width=800, height=600)
expField |>
  pivot_longer(cols = 8:13) |>
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free")
dev.off()

# Emergence is quite left skewed, whereas rootlength
# and shootlength seem approximately normal
# soe is somewhat left skewed, and mesocotyl seems
# right skewed(?)

# We will need to link the field phenotypes to the genotypic dataset just like we did with the
# ragdoll phenotypes. For that we will also need the genoInfo dataset created in PhenoRagdoll:
load(file = here("data", "genoInfo.Rdata"))

# Matching expField to genoInfo through the genonumber column
# There are missing values in the genonumber column of genoInfo, 
# but there is no trivial way to retrieve them
tgtIdx <- match(expField$genonumber, genoInfo$genonumber)
sum(is.na(tgtIdx))

# Add the genoID column to the expRagdoll dataset (for eventual use in GBLUP)
expField <- expField |>
  mutate(genoID = genoInfo$genoID[tgtIdx], subpop = genoInfo$subpop[tgtIdx])

# Checking for NAs in expField again

expField |>
  is.na() |>
  colSums()

# Let's load the accessions to check if the experimental dataset has samples from all of them
load(file = here("data", "snpAccessions.RData"))

# Elements present in expField but not in the snp Accessions
setdiff(unique(expField$genoID), snpAccessions)

# Filter expField for only genotypes found in the snpAccessions dataset
expField <- expField |>
  filter(genoID %in% snpAccessions)

# Dropping levels no longer represented:
expField <- expField |>
  droplevels()

# Checking for NAs in the full consolidated experimental dataset
expField |>
  is.na() |>
  colSums()

# NAs only in the responses now, which is manageable with GBLUP
# We have 454 genotypes left, which is acceptable
# (The ceiling is 470 given the SNP dataset)

# Saving for posterity
save(expField, file = here("data", "expField.Rdata"))

#--------------------------------------------------------------------------------------------------#

#### Using basic mixed models to extract heritability measures for each trait ####

load(file = here("data", "expField.Rdata"))

# We filter for plots that were assigned the planting depth of 8 cm, as our interest
# lies in assessing emergence at significant depths
expFieldDeep <- expField |>
  filter(depth == 'Deep') |>
  droplevels() # in case some genotype is not represented in the "deep" treatment

# In this case, our traits of interest are soe, emergence, mesocotyl, rootlength and shootlength

# Creating list of traits to analyze
traits <- colnames(expField)[colnames(expField) %in% c("soe", "emergence", "mesocotyl", "rootlength", "shootlength")]

### Estimating heritability with genotype as random for variance component extraction #### 

## Useful to gauge whether genomic prediction (GP) is useful
## (if heritability is low, GP might not be effective)

# Heritability with Cullis method
cullisHerit <- list()

# Heritability in an alternative way
altHerit <- list()

# The classification will be by genoID since these values are the ones present in the genotypic dataset
for (trait in traits) {
  model <- asreml(as.formula(paste(trait, "~ replication")),
                  random = ~ genoID, 
                  residual = ~ idv(units), 
                  data = expFieldDeep)
  pred <- predict(model, classify = "genoID")
  
  ###  Cullis heritability calculation
  
  # Average standard error of difference between the predicted means
  avsed <- pred$avsed
  
  # Average variance
  aved <- avsed^2
  
  # Variance component associated with the genoID variable
  vg <- summary(model)$varcomp["genoID", "component"]
  
  cullisHerit[[trait]] <- 1 - (aved / (2*vg))
  
  ### Alternative way to calculate heritability - Based on error plot variance
  # Divide by 3 due to the number of reps -> average residual error for a given genotype
  altHerit[[trait]] <- vpredict(model, h2 ~ V1 / (V1 + (V2/3))) 
}

# Cullis seems to yield smaller heritabilities























