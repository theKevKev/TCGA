library(dplyr)

#' Clean a TCGA methylation data set
#'
#' @param data TCGA methylation data
#'
#' @return Clean data set.
#' 
clean_methy <- function(data,
                      zerosAsNA = TRUE,
                      removeZeros = TRUE,
                      zero_cutoff = 0.3,
                      knn = TRUE,
                      norm = TRUE) {
  
  ## Remove columns with var = 0, (removes 0-cols too)
  data <- data[, !apply(data, 2, function(x){var(x, na.rm = TRUE)}) == 0]
  
  if(removeZeros) {
    ## Remove columns that with more than x% of zeros
    data <- data[, colSums(data == 0, na.rm = TRUE)/nrow(data) < zero_cutoff]
  }
  
  if(zerosAsNA) {
    data[data == 0] <- NA
  }
  
  if (knn) {
    ## KNN imputation
    data <- t(data)
    data <- impute::impute.knn(data)
    data <- t(data$data)
  } else {
    stop("No Imputation Method Provided")
  }
  
  if (norm) {
    data <- scale(data)
    data <- as.data.frame(data)
  }
  
  return(data)
}

# ================== Example Usage ================== 
# pure <- read.table("/Users/home/Desktop/Cancer Datasets/aml/methy")
# data <- as.data.frame(t(pure))
# cleaned_data <- clean_methy(data)
