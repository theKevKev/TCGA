library(dplyr)

#' Clean a TCGA expression data set
#'
#' @param data TCGA expression data
#'
#' @return Clean data set.
#' 
clean_exp <- function(data,
                      removeUnknownMarkers = TRUE,
                      zerosAsNA = TRUE,
                      removeZeros = TRUE,
                      zero_cutoff = 0.3,
                      knn = TRUE,
                      norm = TRUE) {
  # Standardize Gene Formatting 
  colnames(data) <- colnames(data) %>% 
    sub("\\|", ".", .) %>% 
    sub("\\?", "X.", .)
  
  # Remove Unknown Marker Columns
  if (removeUnknownMarkers) {
    data <- data[!grepl("X..", colnames(data), fixed = TRUE)]
  }
  
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
  
  # Aggregate Gene Column Data 
  if (removeUnknownMarkers) {
    data <- remove_duplicate_genes(data)
  } else {
    # remove_duplicate_genes_with_X..(data) # Not Implemented!!
    stop("Unknown Markers not processed")
  }
  
  if (norm) {
    data <- scale(data)
    data <- as.data.frame(data)
  }
  
  return(data)
  
}

#' Assimilates duplicate gene expression data through taking median
#'
#' @param exp gene expression data in standardized GENE.1234 format
#'
#' @return
#' @export
#'
#' @examples
remove_duplicate_genes <- function(exp) {
  library(data.table)
  # Transpose the data frame
  exp_t <- as.data.frame(t(exp))
  exp_t$Gene <- rownames(exp_t)
  
  print(paste("start", length(rownames(exp_t))))
  
  # Clean gene names
  exp_t$Gene <- sub("\\..*", "", exp_t$Gene)
  
  # Identify duplicates
  duplicated_genes <- exp_t$Gene[duplicated(exp_t$Gene)]
  unique_duplicated_genes <- unique(duplicated_genes)
  
  if (length(unique_duplicated_genes) > 0) {
    print("Duplicate genes:")
    print(unique_duplicated_genes)
  } else {
    print("No duplicate genes found.")
  }
  
  # Group by gene names and calculate the median
  exp_t_median <- as.data.table(exp_t)
  exp_t_median <- exp_t_median[, lapply(.SD, median, na.rm = TRUE), by = Gene]
  exp_t_median <- as.data.frame(exp_t_median)
  
  print(paste("end", length(rownames(exp_t_median))))
  
  # Transpose back to original orientation
  exp_cleaned <- as.data.frame(t(exp_t_median))
  colnames(exp_cleaned) <- exp_cleaned[1, ]
  exp_cleaned <- exp_cleaned[-1, ]
  exp_cleaned[] <- lapply(exp_cleaned, as.numeric)
  
  return (exp_cleaned)
}

# ================== Example Usage ================== 
# pure <- read.table("/Users/home/Desktop/Cancer Datasets/aml/exp")
# data <- as.data.frame(t(pure))
# cleaned_data <- clean_exp(data)
