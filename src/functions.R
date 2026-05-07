# Function to obtain adjusted means assuming a fixed model structure and different
# response variables/datasets -> for field/ragdoll experiments

adjMeans <- function(dataset, resp_name){
  
  response <- dataset[[resp_name]]
  
  model <- asreml(fixed = response ~ replication + genoID,
                      residual = ~ idv(units),
                      data = dataset)
  
  # Obtaining BLUEs
  adjms <- predict(model, classify = "genoID")$pvals
  
  # Prediction errors matrix:
  # varcov structure of the estimated means across genotypes
  # we expects different genotypes to be uncorrelated in the absence
  # of genetic information
  vcov <- predict(model, classify = "genoID", vcov = TRUE)
  mat <- vcov$vcov
  
  # Converting mat to conventional numeric matrix format
  mat <- as.matrix(mat)
  
  # Naming the columns and rows of the matrix according to the genotypes
  dimnames(mat) <- list(adjms$genoID, adjms$genoID)
  
  # Heatmap to visually assess the correlation structure between genotypes
  # without any genomic input
  heatmap(mat)
  
  # Per Piepho's paper, since we have an almost perfectly balanced design, in
  # a single environment, calculating weights by the inverse of the squared 
  # standard error of the genotype means is reasonable
  
  # Thus, weights are calculated as 1/SE^2:
  adjms$weight <- 1/(adjms$std.error^2)
  
  # Keeping only the information needed for GBLUP
  adjms <- adjms |>
    select(genoID, predicted.value, weight) |>
    rename(genotype = genoID, BLUE = predicted.value)
    
    # Extract dataset name
    data_name <- deparse(substitute(dataset))
    
    # Extracts the first 5 characters of the dataset name after "exp"
    initDataset <- str_sub(data_name, start = 4, end = 8) |> str_to_title()
    
    # Extracts the first 5 characters of the response name
    initResp <- str_sub(resp_name, start = 1, end = 5) |> str_to_title()
    
    path = paste0("adj", initDataset, initResp, ".RData")
  
  save(adjms, file = here("output", path))
}

