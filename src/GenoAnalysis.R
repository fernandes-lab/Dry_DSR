library(here)
library(data.table)
library(ASRgenomics)
library(AGHmatrix)
library(dplyr)
library(janitor)
library(stringr)
library(tibble)
library(SNPRelate)

# To create the genomic information matrix, we must use the SNPRelate package
# hosted on Bioconductor

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SNPRelate")

################################################################################
#### Building genomic relationship matrix (based on shared genetic material - SNPs)
# Read SNP data
# We also want each SNP as a column and each accession as a row
# This is the full SNP dataset with 2.9 million SNPs,
# per https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2023.1172816/full
# Local workflow to extract metadata and SNP doses while the cluster is down:

# Reading snp file
snpData <- fread(here("data/sandeepOnly", "sandeep1_numeric.txt"))
snpData[1:10, 1:10]

# Extracting metadata
metadata <- snpData[, 2:5]
snpData <- snpData[, -c(2:5)] # removing from original SNP dataset

# Converting first column into row names
snpData <- column_to_rownames(snpData, var = "snp")

# Converting data frame to matrix for better memory usage:
snpMatrix <- as.matrix(snpData)

# Transposing matrix:
snpMatrix <- t(snpMatrix)

# Each row in metadata is a SNP, and we do not transpose it
# Let's remove the whitespaces from metadata
rownames(metadata) <- colnames(snpMatrix) # re-add snp names
metadata <- apply(metadata, 1, str_trim)

# Metadata may be saved as a data frame
metadata <- as.data.frame(metadata)
save(metadata, file = here("data", "metadata.RData"))

# Before saving the matrix, let's convert the snp dosages to 0, 1, 2
snpMatrix <- snpMatrix + 1
save(snpMatrix, file = here("data", "snpMatrix.RData"))

# Proceed to LD pruning section after this (some changes will be made there)

#------------------------------------------------------------------------------------------------------------------------------#
################################ LD pruning ###################################################
# We want to remove SNPs that are in high linkage disequilibrium (LD) with each other, as they provide redundant information
# Load original SNP dataset, as we will need metadata to perform LD pruning

## First of all, let's convert the SNP dosages to (0, 1, 2) format:
# load(file = here("data", "snpData.RData"))
# snpData <- snpData + 1

# Loading metadata again:
load(file = here("data", "metadata.RData"))

# And SNP matrix:
load(file = here("data", "snpMatrix.RData"))

# Saving snpData with updated dose values (maybe not necessary)
# save(snpData, file = here("data", "snpData.RData"))

## Onto LD pruning itself:

# Gds file path for processing data
outGds <- here("output", "geno.gds")

# Creating gds file
snpgdsCreateGeno(outGds, genmat = snpMatrix,
                 sample.id = rownames(snpMatrix),
                 snp.id = colnames(snpMatrix),
                 snp.chromosome = metadata$chr,
                 snp.position = metadata$pos,
                 snpfirstdim = F)

# File ready for pruning
genofile <- snpgdsOpen(outGds)

# Defining an LD threshold of 0.8
# When two SNPs have correlation coefficient greater than 0.8, one of them gets removed
thr <- 0.8
sel_pruned <- snpgdsLDpruning(genofile, ld.threshold = thr, start.pos = "first", method = "corr", verbose = T)

# List of selected SNPs after pruning
# This list has 12 sublists, one for each chromosome
save(sel_pruned, file = here("output", "sel_pruned.RData"))
# load(file = here("output", "sel_pruned.RData"))

# Get the list of selected SNPs after pruning
# in a format that has all SNP names in a single vector/list of one dimension, instead of a list of vectors (one for each chromosome)
sel_pruned <- unname(unlist(sel_pruned))

# Now we filter snpData to keep only the pruned SNPs
# load(file = here("data", "snpData.RData"))
snpPruned <- snpMatrix[, sel_pruned]
# snpPruned has the doses of all the SNPs that remained after pruning

# Verifying the dimension of the pruned SNP dataset
# dim(snpPruned)

# Saving it to memory
save(snpPruned, file = here("output", "snpPruned.RData"))

# The above SNP dataset is what we will use to build the genomic information matrix
# Building the genomic information matrix with the pruned SNP dataset:

# load(file = here("output", "snpPruned.RData"))

# No further pruning, so the minimum allele frequency (maf) is set to 0
G <- Gmatrix(snpPruned, maf = 0)
save(G, file = here("output", "G.RData"))

#----------------------------------------------------------------------------------------#


# The result below gives us an idea of the inbreeding coefficient
sum(diag(G)/nrow(G)) - 1
# It is close to 1, which means our population is highly inbred, or homozygous
# Which is to be expected given the preferred mating system (selfing) of rice

# In the heatmap of the G matrix, darker shades indicate more highly correlated individuals,
# with specially darker shades along the diagonal, since it's individuals with themselves
# The genetic basis of rice is more narrow, which ties back to its high inbreeding coefficient
# (All the correlation observed is genetic, as it comes from the G matrix)
# Consider doing PCA later to see the population structure!!