# Function to obtain adjusted means assuming a fixed model structure and different
# response variables/datasets -> for field/ragdoll experiments

# Note for field data: 
# Mixed model assuming block has no effect (the genotypes are not carried from
# one block to another within a single replication, so estimating block effect
# is meaningless)
# Actually, the AIC and BIC decreased with the inclusion of block 
# in the field experiment, so I will leave block as an optional argument

adjMeans <- function(dataset, resp_name, blck = NULL){
  
  response <- dataset[[resp_name]]
  
  # If we want to include blocking effect (in the field experiment case)
  if(!is.null(blck)){
    block <- dataset[[blck]]
  
    model <- asreml(fixed = response ~ replication + genoID,
                    random = ~ block,
                    residual = ~ idv(units),
                    data = dataset)
  }else{
    model <- asreml(fixed = response ~ replication + genoID,
                    residual = ~ idv(units),
                    data = dataset)
  }
  
  # Obtaining BLUEs
  adjms <- predict(model, classify = "genoID")$pvals
  
  # Squared prediction error (?) matrix:
  # varcov structure of the estimated means across genotypes
  # we expects different genotypes to be uncorrelated in the absence
  # of genetic information
  # vcov <- predict(model, classify = "genoID", vcov = TRUE)
  # mat <- vcov$vcov
  
  # Converting mat to conventional numeric matrix format
  # mat <- as.matrix(mat)
  
  # Naming the columns and rows of the matrix according to the genotypes
  # dimnames(mat) <- list(adjms$genoID, adjms$genoID)
  
  # Heatmap to visually assess the correlation structure between genotypes
  # without any genomic input
  # heatmap(mat) # May be uncommented if need be
  
  # Per Piepho's paper, since we have an almost perfectly balanced design, in
  # a single environment, calculating weights by the inverse of the squared 
  # standard error of the genotype means is reasonable
  
  # Note: commented line results might be returned in future 
  # iterations of the function
  
  # Thus, weights are calculated as 1/SE^2:
  adjms$weight <- 1/(adjms$std.error^2)
  
  # Keeping only the information needed for GBLUP
  adjms <- adjms |>
    select(genoID, predicted.value, weight) |>
    rename(genotype = genoID, BLUE = predicted.value)
  
  return(adjms)
}




