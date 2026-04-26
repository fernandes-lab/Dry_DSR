library(here)
library(data.table)
library(ASRgenomics)
library(SNPRelate)
library(AGHmatrix)
# library(dplyr)

#### Building genomic relationship matrix (based on shared genetic material - SNPs) ####
snps <- fread(here("data/sandeepOnly", "sandeep1_numeric.txt")) # Read SNP data
snps <- t(snps) # We want each SNP as a column and each accession as a row
colnames(snps) <- snps[1, ] # Set the first row as column names
snps <- snps[-(1:5), ] # Remove first five rows (metadata), leaving only SNP dosage

save(snps, file = here("data", "snps.RData")) # Save the processed SNP data for later use

#----------------------------------------------------------------------------------------#

load(file = here("data", "snps.RData"))

snps_numeric <- apply(snps, 2, as.numeric) # SNP dosages are coded as characters, but we need numeric values

row.names(snps_numeric) <- row.names(snps) # So row names are preserved after numeric conversion

save(snps_numeric, file = here("data", "snps_numeric.Rdata"))

# Calculating genomic relationship matrix (kinship matrix or GRM)

load(file = here("data", "snps_numeric.Rdata"))

snps_numeric <- snps_numeric + 1 # SNP dosages coded to 0, 1 and 2

# Generating a list where the G matrix lies
# Calculate GRM using VanRaden method
G_list <- G.matrix(M = snps_numeric, 
                    method = "VanRaden") 

# Save the G matrix for later use
save(G_list, file = here("data", "G_list.Rdata")) 

# Extract the G matrix from the list
# The G matrix informs us on how the genotypes are genetically correlated
# Among themselves, enabling them to borrow information from one another
# Which helps when not all genotypes can be tested in all locations
G_matrix <- G_list$G 

# The result below gives us an idea of the inbreeding coefficient
# It is close to 1, which means our population is highly inbred, or homozygous
# Which is to be expected given the preferred mating system (selfing) of rice
# sum(diag(G_matrix)/nrow(G_matrix)) - 1

# In the heatmap below, darker shades indicate more highly correlated individuals,
# with specially darker shades along the diagonal, since it's individuals with themselves
# The genetic basis of rice is more narrow, which ties back to its high inbreeding coefficient
# (All the correlation observed is genetic, as it comes from the G matrix)
# Consider doing PCA later to see the population structure!!
heatmap(G_matrix)

#------------------------------------------------------------------------------------------------------------------------------#
### LD pruning
# We want to remove SNPs that are in high linkage disequilibrium (LD) with each other, as they provide redundant information
# Load original SNP dataset, as we will need metadata to perform LD pruning
snps <- fread(here("data/sandeepOnly", "sandeep1_numeric.txt"))

## First of all, let's convert the SNP dosages to (0, 1, 2) format:

# One data structure for metadata
metadata <- snps[, 1:5]

# Another data structure for SNP dosages
snpDoses <- snps[, -1:-5]

# Convert to (0, 1, 2) format
snpDoses <- snpDoses + 1 

# We will work with the snpDoses dataset, keeping the other part as metadata
# The first column of metadata will be used as column names for the trasnposed snpDoses dataset
snpDoses <- t(snpDoses)
colnames(snpDoses) <- metadata$snp

# Saving snpDoses for later use in LD pruning and genomic prediction
save(snpDoses, file = here("data", "snpDoses.RData"))

## Onto LD pruning itself:

# Gds file path for processing data
snpsPrePruned <- here("output", "snpsPrePruned.gds")

# Creating gds file
snpgdsCreateGeno(snpsPrePruned, genmat = snp_chrPos_Doses, 
                sample.id = rownames(snp_chrPos_Doses), 
                snp.id = colnames(snp_chrPos_Doses),
                snp.chromosome = metadata$chr,
                snp.position = metadata$pos, 
                snpfirstdim = F)

# File ready for pruning
genofile <- snpgdsOpen(here("output", "snpsPrePruned.gds"))

# Defining an LD threshold of 0.8
# When two SNPs have correlation coefficient greater than 0.8, one of them gets removed
thr <- 0.8
snpsPruned <- snpgdsLDpruning(genofile, ld.threshold = thr, start.pos = "first", method = "corr", verbose = T)

# List of selected SNPs after pruning
# This list has 12 sublists, one for each chromosome
# save(snpsPruned, file = here("output", "snpsPruned.RData"))
load(file = here("output", "snpsPruned.RData"))

# Get the list of selected SNPs after pruning
# in a format that has all SNP names in a single vector/list of one dimension, instead of a list of vectors (one for each chromosome)
snpsPruned <- unname(unlist(snpsPruned)) 

# Now we filter snpDoses to keep only the pruned SNPs
load(file = here("data", "snpDoses.RData"))
snpDosesPruned <- snpDoses[, snpsPruned]
# snpDosesPruned has the doses of all the SNPs that remained after pruning

# Verifying the dimension of the pruned SNP dataset
dim(snpDosesPruned)

# Saving it to memory
save(snpDosesPruned, file = here("output", "snpDosesPruned.RData"))

# The above SNP dataset is what we will use to build the genomic information matrix
# Building the genomic information matrix with the pruned SNP dataset:

load(file = here("output", "snpDosesPruned.RData"))

# No further pruning, so the minimum allele frequency (maf) is set to 0
G <- Gmatrix(snpDosesPruned, maf = 0)
save(G, file = here("output", "G_pruned.RData"))

# Phenotypic analysis starts at allTraitsPhenoRagdoll.R