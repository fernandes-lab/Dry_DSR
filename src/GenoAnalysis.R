library(here)
library(data.table)
library(ASRgenomics)
library(SNPRelate)
library(AGHmatrix)
library(dplyr)
library(janitor)
library(stringr)


#### Building genomic relationship matrix (based on shared genetic material - SNPs) ####
# Read SNP data
# We also want each SNP as a column and each accession as a row
# This is the full SNP dataset with 2.9 million SNPs,
# per https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2023.1172816/full
snpData <- t(fread(here("data/sandeepOnly", "sandeep1_numeric.txt")))
snpData[1:10, 1:10]
dim(snpData)
# The first five rows are metadata, while the remainder are genotypes

# The very first row is the SNP names
# The fifth row has only NAs, so it may be excluded
unique(snpData[5, ])
snpData <- snpData[-5, ]

# Promoting first row to column names
# Now the column names correspond to SNP names
snpData <- row_to_names(snpData, row_num = 1)

# Storing metadata in a separate data frame for later use
# allele, chromosome and position information
metadata <- snpData[1:3, ]
# Remove whitespaces from metadata rows
metadata <- apply(metadata, 2, str_trim)

# Saving metadata file
save(metadata, file = here("data", "metadata.RData"))

# Removing metadata rows from main SNP dataset
# Removing spaces from SNP doses and turning them to numbers
snpData <- snpData[-1:-3, ]

# Remove leading and trailing whitespaces
snpData <- apply(snpData, 2, str_trim)

# Convert SNP doses to numeric values
# We want to keep row names intact, hence the list syntax
snpData[] <- lapply(snpData, as.numeric)

# Save the processed SNP data for later use
save(snpData, file = here("data", "snpData.RData"))

# Extracting rownames from snpData to do some checking with the phenotypic data
# Saving to memory
snpAccessions <- rownames(snpData)
save(snpAccessions, file = here("Thesis/Dry_DSR/CarlosReDo/data", "snpAccessions.RData"))

#------------------------------------------------------------------------------------------------------------------------------#
### LD pruning
# We want to remove SNPs that are in high linkage disequilibrium (LD) with each other, as they provide redundant information
# Load original SNP dataset, as we will need metadata to perform LD pruning

## First of all, let's convert the SNP dosages to (0, 1, 2) format:
load(file = here("data", "snpData.RData"))

snpData <- snpData + 1

# Saving snpData with updated dose values 
save(snpData, file = here("data", "snpData.RData"))

## Onto LD pruning itself:

# Gds file path for processing data
outGds <- here("output", "geno.gds")

# Creating gds file
snpgdsCreateGeno(outGds, genmat = snpData,
                 sample.id = rownames(snpData),
                 snp.id = colnames(snpData),
                 snp.chromosome = metadata$chr,
                 snp.position = metadata$pos,
                 snpfirstdim = F)

# File ready for pruning
genofile <- snpgdsOpen(here("output", "outGds.gds"))

# Defining an LD threshold of 0.8
# When two SNPs have correlation coefficient greater than 0.8, one of them gets removed
thr <- 0.8
sel_pruned <- snpgdsLDpruning(genofile, ld.threshold = thr, start.pos = "first", method = "corr", verbose = T)

# List of selected SNPs after pruning
# This list has 12 sublists, one for each chromosome
# save(sel_pruned, file = here("output", "sel_pruned.RData"))
load(file = here("output", "sel_pruned.RData"))

# Get the list of selected SNPs after pruning
# in a format that has all SNP names in a single vector/list of one dimension, instead of a list of vectors (one for each chromosome)
sel_pruned <- unname(unlist(sel_pruned))

# Now we filter snpData to keep only the pruned SNPs
load(file = here("data", "snpData.RData"))
snpPruned <- snpData[, sel_pruned]
# snpPruned has the doses of all the SNPs that remained after pruning

# Verifying the dimension of the pruned SNP dataset
dim(snpPruned)

# Saving it to memory
save(snpPruned, file = here("output", "snpPruned.RData"))

# The above SNP dataset is what we will use to build the genomic information matrix
# Building the genomic information matrix with the pruned SNP dataset:

load(file = here("output", "snpPruned.RData"))

# No further pruning, so the minimum allele frequency (maf) is set to 0
G <- Gmatrix(snpPruned, maf = 0)
save(G, file = here("output", "G.RData"))


#----------------------------------------------------------------------------------------#


# The result below gives us an idea of the inbreeding coefficient
# It is close to 1, which means our population is highly inbred, or homozygous
# Which is to be expected given the preferred mating system (selfing) of rice
# sum(diag(G_matrix)/nrow(G_matrix)) - 1

# In the heatmap of the G matrix, darker shades indicate more highly correlated individuals,
# with specially darker shades along the diagonal, since it's individuals with themselves
# The genetic basis of rice is more narrow, which ties back to its high inbreeding coefficient
# (All the correlation observed is genetic, as it comes from the G matrix)
# Consider doing PCA later to see the population structure!!