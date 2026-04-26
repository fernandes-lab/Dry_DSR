library(here)
library(asreml)
library(dplyr)

# The commented part below contains some data manipulation that was done to prepare the data for the GBLUP model, 
# but it has been commented out because it has already been done and saved to disk.

# # Loading field data, ragdoll data, and genomic data
# Preprocessed in allTraitsPhenoField.R, allTraitsPhenoRagdoll.R, and SNP_data.R, respectively
load(file = here("output", "adjMeansField.RData"))
load(file = here("output", "adjMeansRagdoll.Rdata"))
load(file = here("output", "snpDosesPruned.RData"))

dim(adjMeansField)
dim(adjMeansRagdoll)
dim(snpDosesPruned)

# Rename aux columns to Genotype for easier understanding
# as well as standard error columns
adjMeansField <- adjMeansField |> 
                    rename(Genotype = `auxGenoField`)

adjMeansRagdoll <- adjMeansRagdoll |>
                    rename(Genotype = `auxGenoRag`)

# And saving that to disk

save(adjMeansField, file = here("output", "adjMeansField.RData"))
save(adjMeansRagdoll, file = here("output", "adjMeansRagdoll.Rdata"))

str(adjMeansField)
str(adjMeansRagdoll)
str(snpDosesPruned)

# Converting snpDosesPruned to a data frame
# The rownames of the data frame are the genotype IDs
snpDosesPruned <- as.data.frame(snpDosesPruned)

save(snpDosesPruned, file = here("output", "snpDosesPruned.RData"))

# Loading the G matrix:
# Reminder: the G matrix contains covariance information between genotypes based on their SNP data
load(file = here("output", "G_pruned.RData"))
str(G)

# Two separate genomic predictions models will be built and compared
# One for the field and another for the ragdoll (lab) experiment

############################# FIELD ###############################
# Matching the field to the genomic (snpDosesPruned) dataset:
# Or rather, filtering the G matrix according to the genotypes in the field dataset
G_Field <- G[as.factor(rownames(G)) %in% adjMeansField$Genotype, as.factor(colnames(G)) %in% adjMeansField$Genotype]
dim(G_Field)

# 364 matches only because Genotype in adjMeansField is formatted a bit differently
# Let's fix the column in adjMeansField with the same format as the rownames of the G matrix
# The problem is in the genotypes that start with IRIS, as they have an empty space instead of an
# underscore
adjMeansField$Genotype1 <- as.factor(gsub(" ", "_", adjMeansField$Genotype))
adjMeansField$Genotype1

# Having both Genotype and Genotype1 is redundant, so let's drop Genotype and rename Genotype1 to Genotype
adjMeansField <- adjMeansField |>
                    select(-Genotype) |>
                    rename(Genotype = Genotype1)

save(adjMeansField, file = here("output", "adjMeansField.RData"))

#---------------------------------------------------------------------#
# Loading field data, ragdoll data, and genomic data after changes
# Loading field data, ragdoll data, and genomic data
# Preprocessed in allTraitsPhenoField.R, allTraitsPhenoRagdoll.R, and SNP_data.R, respectively
load(file = here("output", "adjMeansField.RData"))
load(file = here("output", "adjMeansRagdoll.Rdata"))
load(file = here("output", "snpDosesPruned.RData"))

# Loading G matrix (reminder that this matrix is post-pruning)
load(file = here("output", "G_pruned.RData"))

#---------------------------------- Univariate models for FIELD traits -----------------------------------#

# We must filter the adjMeansField dataset to keep only the genotypes that are in the G matrix
adjMeansField_G <- adjMeansField[adjMeansField$Genotype %in% rownames(G), ]

# After filtering out observations, unused factor levels (or genotypes) must be dropped from the Genotype column
adjMeansField_G <- droplevels(adjMeansField_G)

# G has 470 genotypes, of which 454 are in the fields dataset, already filtered for those in common with G
# However G must be subset as well:
G_Field <- G[as.factor(rownames(G)) %in% adjMeansField_G$Genotype, as.factor(colnames(G)) %in% adjMeansField_G$Genotype]

# Running the GBLUP model for the field data:
# The model will be run for each trait separately, and the predictions will be stored in a
# list of data frames, one for each trait
# SOE_index and Emergence_index will be dropped as their adjusted means yielded NA values for all genotypes, so they cannot be used in the model
# Genotype is also evidently not a response trait we are interested in

# We must make 'Genotype' the first column of adjMeansField_G
# And match its order with that of the rownames (each representing a genotype) in the G matrix
adjMeansField_G <- adjMeansField_G |>
                    relocate(Genotype) |>
                    arrange(match(Genotype, rownames(G_Field)))

# Let's save these core data to disk (adjusted means data for field and filtered genomic information matrix for field experiment)
save(adjMeansField_G, file = here("output", "adjMeansField_G.RData"))
save(G_Field, file = here("output", "G_Field.RData"))

# For calculating the numerical estimation of breeding values, we run the GBLUP
# For now, let's run the GBLUP without any cross-validation just to rank the genotypes according to their predicted breeding values
# Soon we will go through the proper pipeline of assessing the model's predictive accuracy

# Reloading necessary files to run GBLUP
load(file = here("output", "adjMeansField_G.RData"))
load(file = here("output", "G_Field.RData"))

# We must weight the adjusted values using their standard error to use in the GBLUP
weightsField <- adjMeansField_G |> 
                select(SOE_Error:RootWeight_Error) |> 
                # weights are basically the inverse of the variances
                mutate(across(everything(), ~ 1/(.x^2))) |>
                rename_with(~ gsub("Error", "Weight", .x), ends_with("Error"))


# Selecting traits to run GBLUP
# We are interested only in emergence related traits (SOE and Emergence)
traitsField <- colnames(adjMeansField_G)[(colnames(adjMeansField_G) %in% c("SOE", "Emergence"))]

# We do not need the standard errors in the adjMeansField_G dataset, but rather the weights, so let's replace that
adjMeansField_G <- adjMeansField_G |>
                select(-SOE_Error:-RootWeight_Error) |>
                cbind(weightsField)

# Subsetting adjMeansField_G for only the traits we are interested in and their corresponding weights
# We will keep adjMeansField_G around for possible future inclusion of other traits
adjMeansField_GBLUP <- adjMeansField_G |>
                        select("Genotype", "SOE", "SOE_Weight", "Emergence", "Emergence_Weight")

# Let's drop the NAs from adjMeansField_GBLUP
# 448 genotypes remain
adjMeansField_GBLUP <- adjMeansField_GBLUP |>
                        na.omit() |>
                        droplevels()

# Further filtering the G matrix to match the dataset without NAs
G_Field_NA <- G_Field[as.factor(rownames(G_Field)) %in% adjMeansField_GBLUP$Genotype, 
as.factor(rownames(G_Field)) %in% adjMeansField_GBLUP$Genotype]

# List to store predicted breeding values for each trait in the field dataset
# As well as their related information
GBLUP_Field <- list()

# Running the GBLUP model for each trait in the field dataset (SOE and Emergence)
for (trait in traitsField) {
    
    if (trait == "SOE"){
    # Different weight columns for different traits
    model <- asreml(as.formula(paste(trait, "~ 1")),
                    random = ~ vm(Genotype, G_Field_NA),
                    residual = ~ idv(units),
                    weights = SOE_Weight,
                    data = adjMeansField_GBLUP)
    }else{
        model <- asreml(as.formula(paste(trait, "~ 1")),
                    random = ~ vm(Genotype, G_Field_NA),
                    residual = ~ idv(units),
                    weights = Emergence_Weight,
                    data = adjMeansField_GBLUP)
    }

    pred <- predict(model, classify = "Genotype")
    
    GBLUP_Field[[trait]] <- pred$pvals
}

# GBLUP predictions saved without cross-validation and predictive accuracy measured
# Once all of that is done, its rankings on the full dataset will be compared to ragdoll rankings
save(GBLUP_Field, file = here("output", "GBLUP_Field.RData"))

#---------------------------------- Univariate models for RAGDOLL traits -----------------------------------#

# Proceeding in an analogous way to the field data, but with the ragdoll dataset instead
G_Ragdoll <- G[as.factor(rownames(G)) %in% adjMeansRagdoll$Genotype, as.factor(colnames(G)) %in% adjMeansRagdoll$Genotype]

adjMeansRagdoll_G <- adjMeansRagdoll[adjMeansRagdoll$Genotype %in% rownames(G_Ragdoll), ]
adjMeansRagdoll_G <- droplevels(adjMeansRagdoll_G)

# 462 genotypes matched, a little over the 454 for the field experiment


# Matching the order of genotypes in adjMeansRagdoll_G with that of the rownames in G_Ragdoll
adjMeansRagdoll_G <- adjMeansRagdoll_G |>
                    arrange(match(Genotype, rownames(G_Ragdoll)))

# Let's save these core data to disk (adjusted means data for ragdoll and filtered genomic information matrix for ragdoll experiment)
save(adjMeansRagdoll_G, file = here("output", "adjMeansRagdoll_G.RData"))
save(G_Ragdoll, file = here("output", "G_Ragdoll.RData"))

# Reloading necessary files to run GBLUP
load(file = here("output", "adjMeansRagdoll_G.RData"))
load(file = here("output", "G_Ragdoll.RData"))

# We must weight the adjusted values using their standard error to use in the GBLUP
weightsRagdoll <- adjMeansRagdoll_G |> 
                select(CL_means_Error:SL_means_Error) |> 
                # weights are basically the inverse of the variances
                mutate(across(everything(), ~ 1/(.x^2))) |>
                rename_with(~ gsub("Error", "Weight", .x), ends_with("Error"))

# Selecting traits to run GBLUP (excluding those with only NAs)
# We are interested only in Mesocotyl and Coleoptile length
traitsRagdoll <- colnames(adjMeansRagdoll_G)[(colnames(adjMeansRagdoll_G) %in% c("CL_means", "ML_means"))]

# We do not need the standard errors in the adjMeansRagdoll_G dataset, but rather the weights, so let's replace that
adjMeansRagdoll_G <- adjMeansRagdoll_G |>
                select(-CL_means_Error:-SL_means_Error) |>
                cbind(weightsRagdoll)

# Subsetting adjMeansRagdoll_G for only the traits we are interested in and their corresponding weights
# We will keep adjMeansRagdoll_G around for possible future inclusion of other traits
adjMeansRagdoll_GBLUP <- adjMeansRagdoll_G |>
                        select("Genotype", "CL_means", "CL_means_Weight", "ML_means", "ML_means_Weight")

# Let's drop the NAs from adjMeansRagdoll_GBLUP
# 462 genotypes remain (all of them) -> no need to filter G_Ragdoll matrix
adjMeansRagdoll_GBLUP <- adjMeansRagdoll_GBLUP |>
                        na.omit() |>
                        droplevels()     

# List for GBLUPs of traits in the field dataset
GBLUP_Ragdoll <- list()

# Running the GBLUP model for each trait in the ragdoll dataset:
for (trait in traitsRagdoll) {
    if (trait == "CL_means"){
    model <- asreml(as.formula(paste(trait, "~ 1")),
                    random = ~ vm(Genotype, G_Ragdoll),
                    residual = ~ idv(units),
                    weights = CL_means_Weight,
                    data = adjMeansRagdoll_GBLUP)
    }else{
        model <- asreml(as.formula(paste(trait, "~ 1")),
                    random = ~ vm(Genotype, G_Ragdoll),
                    residual = ~ idv(units),
                    weights = ML_means_Weight,
                    data = adjMeansRagdoll_GBLUP)
    }

    pred <- predict(model, classify = "Genotype")
    
    GBLUP_Ragdoll[[trait]] <- pred$pvals
}

# GBLUP predictions saved without cross-validation and predictive accuracy measured
# Once all of that is done, its rankings on the full dataset will be compared to field rankings
save(GBLUP_Ragdoll, file = here("output", "GBLUP_Ragdoll.RData"))


#---------------------------------- Multivariate Analysis for Ragdoll ------------------------------------------------#

# Let's do that following Salvador Gezan's example
# mesocotyl + coleoptile as indirect selection for seedling emergence in field
# Maybe isolate the depth (deep sowing) to compare to the lab experiment output
# I wonder if I did adjusted means wrong

# Loading adjusted means for ragdoll data and G matrix with matched genotypes
# Saved previously in this same script
load(here("output", "adjMeansRagdoll_G.RData"))
load(here("output", "G_Ragdoll.RData"))

summary(adjMeansRagdoll_G)
# The CL_means (coleoptile) and ML_means (mesocotyl) columns will be used jointly in a bivariate model
# They are already on similar scales, as recommended by Salvador Gezan
# It seems we can use variance components from the univariate models as starting values


# GBLUP_Ragdoll includes a list of predicted values for each trait in isolation
# As one univariate model was run for each
load(here("output", "GBLUP_Ragdoll.RData"))
summary(GBLUP_Ragdoll)
# GBLUP_Ragdoll$ML_means$predicted.value
# GBLUP_Ragdoll$CL_means$predicted.value

# I wonder if I have to re-run adjusted means with two variables...

# Bivariate model with no starting values:

# Bear in mind that the traits are correlated via the same genotype
# So that must be taken into account in the residuals
# Trait is a new factor, native to asreml, indicating that each trait has its own mean
# vm defines the correlation structure between the genotypes, as opposed to id (no correlation, identity matrix), 
# and other forms of correlation
# corgh structure for traits since two traits in the same genotype are not independent
# multiplying a covariance structure between the genotypes as well
# For the residuals, the correlation between traits make it so that the residuals are correlated across traits
# Note: corgh and us are synonymous, but the former is modeled after correlations and the latter after covariances
model_BI_1 <- asreml(fixed = cbind(ML_means, CL_means) ~ trait,
                    random = ~ corgh(trait):vm(Genotype, G_Ragdoll),
                    residual = ~ id(units):corgh(trait),
                    data = adjMeansRagdoll_G)
summary(model_BI_1)$varcomp
# The model did not crash without starting values...it did converge!

# How do I incorporate weights into the bivariate model???

# Calculating heritabilities for both traits:
vpredict(model_BI_1, heritML1 ~ V2/(V2+V6))
vpredict(model_BI_1, heritML1 ~ V3/(V3+V7))
# Do these heritabilities represent how well the model takes in the genetic information?
# Wait... did the inclusion of CL actually decrease the heritability?

# Predictions for bivariate model
predsBI1 <- predict.asreml(model_BI_1, classify = "trait:Genotype")$pvals
head(predsBI1)
# The prediction errors are quite high still...

# Comparing field genotype ranking to ML+CL genotype ranking
# I assume it is assuming along these following lines ...?
# load(here("output", "GBLUP_Field.RData"))
# sort_by(GBLUP_Field$Emergence, GBLUP_Field$Emergence$predicted.value, decreasing = TRUE)

# The next steps will be to implement cross-validation for both the field and ragdoll GBLUP models, and
# possibly introduce GWAS for separating SNPs with large effects from those with small effects
# Check out GWAS in the GWAS.R script