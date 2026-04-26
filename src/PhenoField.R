library(here)
library(dplyr)
library(stringr)
library(readxl)
library(asreml)

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

# genotype and pedigree seem to provide the same information in different formats
# pedigree should be removed, since genotype can be linked to the ragdoll and thus the
# genotypic dataset, and we can obtain the genotype names from expRagdoll
expField <- expField |>
            select(-pedigree)

# Checking how many times each genotype appears in each replication
with(expField, table(genotype, replication) > 0) |>
        rowSums() |>
        table()
# 503 genotypes appear in 3 replications and 1 appears in 2
# Bear in mind that the deep sowing treatment had 3 replications while the shallow sowing one had 2

# Checking for NAs
expField |>
    is.na() |>
    colSums()











###############################################################################

# Checking the number of genotypes in the pedigree column in expField
unique(compiled$Pedigree)
# We can see there are 502 unique genotypes

# Let's check for NAs
sum(is.na(compiled$Pedigree))
# No NAs

# compiled$Replicate <- as.numeric(gsub("Rep_", "", compiled$Replicate))
# List of traits to analyze
colnames(compiled)

# Our traits of interest are Mesocotyl length, Root length, Shoot length, Speed of Emergence, % Emergence, Shoot dry weight,
# Root dry weight, SOE index and Emergence index
traits <- colnames(compiled)[c(21, 23, 24, 25, 26, 27, 28, 32, 33)]

# Taking a look at the column data types and entries
str(compiled)

# Converting to adequate data types (mainly design columns to factors)
compiled <- compiled |>
            mutate(across(c(SiteYear:TrialNumber, block, Row, Replicate, Genotype, Pedigree), as.factor))

unique(compiled$SiteYear)
unique(compiled$SiteName)
unique(compiled$TrialType)
unique(compiled$TrialNumber)

## Matching genotype IDs to pedigree/genotype information:
# Genotype in the Compiled_RDP1 dataset will be matched to SEQ in the LS_means dataset so the former can have the HDRA assay ID column
# The ID column is important for matching phenotypic data with genotypic data

ls_means <- read_excel(here("data/sandeepOnly/rawData", "LS_means_ragdoll.xlsx"))

# Same approach as the one used in ProtoPhenotypesRagdoll.r, but now with the field trial data
# Here, instead of replacing the Genotype column, we will create a new one for the assay IDs

# The match function gives us the index of the SEQ in ls_means that corresponds to each Genotype in compiled
# Reminder that the first row does not count as it is the header, so the first SEQ is in row 2
target_idx <- match(compiled$Genotype, ls_means$SEQ)

# The GenoID column will have the assay IDs that correspond to each genotype in the field trial dataset, which will be useful for matching with genotypic data later on
compiled$GenoID <- ls_means$'HDRA genotype assay ID'[target_idx]

# Besides the conversions already done above, let's go ahead and do some more manipulations
# SiteYear and SiteName only have one type of entry, and TrialType and TrialNumber give the exact same information

compiled <- compiled |> 
            mutate(GenoID = as.factor(GenoID)) |>
            select(-c(SiteYear, SiteName, PlotBarcode, TrialUnitComment)) |>
            filter(!is.na(GenoID) & GenoID != "NA")


# Let's also rename some variables for better readability, and drop redundant/unused columns:

compiled <- compiled |>
            select(-c(TrialType, Sown, Day_5:Day_14, Seed_count)) |>
            rename(Depth = TrialNumber, ShootWeight = Shoot.dry.wt..g., RootWeight = root.dry.wt..g.)

# save(compiled, file = here("data", "exp_field_filtered.RData"))

# For running the code again:
# This object is read as "compiled"
load(file = here("data", "exp_field_filtered.RData"))

# Note: this is actually a split-plot design (genotypes within depth treatments)
# Depth treatments are whole plots and genotypes (pedigrees) are subplots
# Note2: there were 2 reps for shallow sowing and 3 reps for deep sowing

### Estimating heritability with genotype as random for variance component extraction #### 

## Note: this is basically a copy-paste from the ragdoll experiment script

## Useful to gauge whether genomic prediction (GP) is useful
## (if heritability is low, GP might not be effective)

# Heritability with Cullis method
cullisHerit <- list()

# Heritability in an alternative way
altHerit <- list()

for (trait in traits) {
    model <- asreml(as.formula(paste(trait, "~ Depth + Replicate + Replicate:TrialNumber")),
                    random = ~ Pedigree, 
                    residual = ~ idv(units), 
                    data = compiled)
    pred <- predict(model, classify = "Pedigree")

    # Cullis heritability calculation
    av_sed <- pred$avsed
    aved <- av_sed^2
    vg <- summary(model)$varcomp["Pedigree", "component"]
    cullisHerit[[trait]] <- 1 - (aved / (2*vg))

    # Alternative way to calculate heritability - Based on error plot variance
    altHerit[[trait]] <- vpredict(model, h2 ~ V1 / (V1 + (V2/4))) 
}

cullisHerit
altHerit

#-----------------------------------------------------------------------------------------------------#
# The results imply the model is a bit off, hence the need to further discuss the experimental design
# Blocks are nested within replicates
# Most of the genotypes are not repeated across blocks, appearing only 3 times each in an experiment with 9 blocks
# Most blocks do appear in both depth (whole plot) treatments, but not all of them -> I will not consider this nested
# Pedigrees/genotypes treated as fixed effects for BLUEs estimation
# Even though they were treated as random for heritability estimation

#------------------------------------------------------------------------#
# # Trying to understand the study design
# compiled <- compiled |> mutate(block = as.numeric(block))
# require(desplot)
# png(filename = here("output", "TestDesplot.png"), width = 800, height = 600, units = "px")
# desplot(compiled, Replicate ~ block + Row,
#         col=Row, text=Replicate, cex=1, aspect=511/176,
#         out1=TrialNumber, out2=block, 
#         out2.gpar=list(col = "gray50", lwd = 1, lty = 1))
# dev.off()
#------------------------------------------------------------------------#

# Blocks were not used for heritability estimation as they were deemed not relevant, but I am inserting them here again
# Each block
traits <- colnames(compiled |> 
                select(SOE, Emergence, mesocotyl_length, root_length, shoot_length, SOE_index, Emergence_index, 
                ShootWeight, RootWeight))

BLUEsField <- list()

# Auxilliary list to store the names of the genotypes
auxGenoField <- c()

# Another auxilliary list to store the prediction errors
# (to be used to weight the values in GBLUP)
auxErrorField <- list()

for (trait in traits) {
    # Genotype effect is fixed now as we want to get the adjusted means for each genotype (BLUEs)
    model <- asreml(as.formula(paste(trait, "~ Replicate + Depth + GenoID + GenoID:Depth")),
                    random = ~ Replicate:block,
                    residual = ~ idv(units),
                    workspace = "2gb",
                    maxit = 100,
                    data = compiled)
    pred <- predict(model, classify = "GenoID")
    if(trait == traits[1]){
        auxGenoField <- pred$pvals$GenoID
    }
    BLUEsField[[trait]] <- pred$pvals$predicted.value

    auxErrorField[[trait]] <- pred$pvals$std.error
}

# The BLUEs calculated do not account for spatial corrections, which will be done later (Sandeep claimed there was no significant spatial variation))
names(auxErrorField) <- paste0(names(BLUEsField), "_Error")

adjMeansField <- cbind(as.data.frame(auxGenoField), as.data.frame(BLUEsField), as.data.frame(auxErrorField))
save(adjMeansField, file = here("output", "adjMeansField.RData"))




