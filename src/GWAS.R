library(here)
library(dplyr)
library(tibble)

if (!require(GAPIT)){
  devtools::install_github("jiabowang/GAPIT", force = T)
}
library(GAPIT)

#########################################################
#                    Preparing data                     #
#########################################################

# Phenotypic data (454 genotypes)
# Already corrected for experimental conditions
myY <- load(file = here("output", "adjFieldEmerg.RData"))
rm(adjFieldEmerg)
myY <- myY |>
     select(-weight) |>
     rename(taxa = genotype, emerg = BLUE)

# Genomic numeric data (post-pruning)
myGD <- load(file = here("output", "snpPruned.RData"))
rm(snpPruned)
myGD <- myGD |>
      as.data.frame() |>
      rownames_to_column(var = "taxa")

# Genetic map file
myGM <- load(file = here("data", "metadata.RData"))
rm(metadata)
myGM <- myGM |>
  select(-c(allele, cm)) |>
  rownames_to_column(var = "snpID")

# Metadata is not yet filtered post-pruning
# Filtering GM according to GD
myGM <- myGM[myGM$snpID %in% colnames(myGD), ]


#########################################################
#                 Running GAPIT                         #
#########################################################

gwas <- GAPIT(
        Y = myY,
        GD = myGD,
        GM = myGM,
        PCA.total = 5,
        model = "BLINK"
)

# Results stored in "./output/GWAS" folder
