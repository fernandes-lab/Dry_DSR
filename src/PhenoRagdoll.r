library(here)
library(asreml)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(readxl)

# According to the slides provided by Sandeep, there were 417 (RDP1) + 91 (3K RG) = 508 total accessions, 
# meaning we have about 508 genotypes before any filtering

#-------------------------------------------------------------------------------------------------------#
#### Basic exploratory analysis and cleaning of the data ####

# Experimental information dataset
exp_ragdoll <- read.csv2(file = here("data/sandeepOnly/rawData", "Compiled_TestTube_data_means_V3.csv"), 
                            header = TRUE, sep = ",")

dim(exp_ragdoll)
str(exp_ragdoll)
summary(exp_ragdoll)
# 1971 observations of 15 variables (in total)

# Some variables need to be converted to meaningful types
# Also changing to CamelBack notation for better organization
# Furthermore, let's reorder the columns and rename some of them
# Also keep everything lowercase
# ML_means1 and RL_means1 are after removing outliers from the plots
# I do not have access to the specific procedures so I will just remove both
# SiteName and TrialType have only one level each, so it is pointless to keep them
# The columns siteyear and replication represent the same information as well
expRagdoll <- exp_ragdoll |>
                # 'Genotype' becomes the first column
                relocate(Genotype) |>
                select(-c(SiteName, TrialType, ML_means1, RL_means1)) |>
                mutate(across(Genotype:Row, as.factor)) |>
                mutate(across(seed_count:SL_means, as.numeric)) |> 
                rename(coleoptile = CL_means, mesocotyl = ML_means, rootlength = RL_means, shootlength = SL_means) |>
                # Turn column names to lower case
                rename_with(tolower) |>
                rename(seedcount = seed_count, plotnum  = plotbarcode) |>
                select(-siteyear)

# Each plot corresponds to five seeds of a given genotype
str(expRagdoll)
summary(expRagdoll)

# Checking how many times each genotype appears in each replication
# Each genotype was supposed to be replicated 4 times
with(expRagdoll, table(genotype, replication) > 0) |>
        rowSums() |>
        table()
# We can see that 472 genotypes were replicated 4 times
# The experiment is almost a perfect randomized complete block design (if each rep is to be treated as a block)
# Also, accession is not necessarily the same as genotype, though it often is
# Checking for NAs
expRagdoll |>
    is.na() |>
    colSums()
# There are 178 NAs in total, of which 133 come from the mesocotyl column and 24 from the coleoptile column, besides others
# Reminder that there are 1971 observations of 10 variables in total (133 is less than 10% of 1971, for instance)

# Checking distribution of variables
png(here("output/plots","RagDollExpDists.png"), width=800, height=600)
expRagdoll |>
    pivot_longer(cols = 7:10) |>
    ggplot(aes(x = value)) +
    geom_density() +
    facet_wrap(~name, scales = "free")
dev.off()
# Coleoptile and shootlength seem to have approximately normal distributions
# Mesocotyl is right-skewed, and rootlength is bimodal

# LS_means_ragdoll has the genotype IDs
# Let's get the relevant columns from it
# We need the IDs to link the ragdoll phenotype data to the genotypic data
# Also renaming columns for easier understanding and manipulation
genoInfo <- read_excel(here("data/sandeepOnly/rawData", "LS_means_ragdoll.xlsx")) |>
                select(SEQ, Genotype, Sub_pop, 'HDRA genotype assay ID') |>
                rename(genotype = Genotype, genonumber = SEQ, subpop = Sub_pop, genoID = 'HDRA genotype assay ID') |>
                mutate(across(everything(), as.factor))

# Checking for NAs in genoInfo:
genoInfo |>
    is.na() |>
    colSums()

# Build indices to match the experimental dataset to the genotype codes dataset using
# shared genotype names
# "This guy in expRagdoll$genotype is in position target_idx in genoInfo$genotype" -> tgtIdx is a vector of such positions 
tgtIdx <- match(expRagdoll$genotype, genoInfo$genotype)
sum(is.na(tgtIdx))
# There are 4 NAs, meaning 4 genotypes in expRagdoll could not be found among the 504 available in genoInfo


# Add the genoID column to the expRagdoll dataset (for eventual use in GBLUP)
# Then clean genoID of NA values (there is a "NA" string as well)
# Finally, drop levels that were lost after removing NAs
expRagdoll <- expRagdoll |>
                mutate(genoID = genoInfo$genoID[tgtIdx], subpop = genoInfo$subpop[tgtIdx]) |>
                filter(!genoID %in% c("NA", NA)) |>
                droplevels()

# Correct some genotype names so they match the 
# genotypic dataset (e.g. IRIS 313-11949 becomes IRIS_313-11949)
expRagdoll <- expRagdoll |> 
    mutate(genoID = as.factor(str_replace_all(genoID, "IRIS ", "IRIS_")))

# Checking for NAs in the full consolidated experimental dataset
expRagdoll |>
    is.na() |>
    colSums()
# NAs only in the responses. The G matrix will hopefully bridge this gap
# by leveraging the shared information across genotypes

# Saving for posterity
save(expRagdoll, file = here("data", "expRagdoll.Rdata"))

#--------------------------------------------------------------------------------------------------#

#### Using basic mixed models to extract heritability measures for each trait ####

# In this case, our traits of interest are coleoptile, mesocotyl, rootlength and shootlength

# Creating list of traits to analyze (CL_means through SL_means)
traits <- colnames(expRagdoll)[colnames(expRagdoll) %in% c("coleoptile", "mesocotyl", "rootlength", "shootlength")]
          
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
                    data = expRagdoll)
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
    # Divide by 4 due to the number of reps -> average residual error for a given genotype
    altHerit[[trait]] <- vpredict(model, h2 ~ V1 / (V1 + (V2/4))) 
}

# Both heritability estimation methods yielded similar resuls to Sandeep's


