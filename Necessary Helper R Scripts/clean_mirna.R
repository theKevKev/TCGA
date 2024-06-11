library(dplyr)

#' Clean a TCGA expression data set
#'
#' @param data TCGA expression data
#' @param onlyHumanMiRNA whether or not to keep EBV and HCMV mirna data
#' @param aggregateDifferentLocuses Set this to TRUE if we only care about the 
#' total amount of mirna and not which locuses they came from (will take the 
#' median of the readings (do we want to sum instead?), which is preferred over
#' summing given normalization to avoid outliers)
#'
#' @return Clean data set.
#' 
clean_mirna <- function(data,
                      onlyHumanMiRNA = FALSE, 
                      zerosAsNA = TRUE,
                      removeZeros = TRUE,
                      zero_cutoff = 0.3,
                      knn = TRUE,
                      aggregateDifferentLocuses = FALSE, 
                      norm = TRUE) {
  # Standardize Gene Formatting 
  colnames(data) <- colnames(data) %>% 
    gsub("\\.", "-", .) %>% 
    sub("miR", "mir", .)
  
  # Remove Unknown Marker Columns
  if (onlyHumanMiRNA) {
    data <- data[grepl("hsa-let", colnames(data), fixed = TRUE) | grepl("hsa-mir", colnames(data), fixed = TRUE)]
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
  if (aggregateDifferentLocuses) {
    data <- combine_different_locuses(data)
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
combine_different_locuses <- function(mirna) {
  library(data.table)
  # Transpose the data frame
  mirna_t <- as.data.frame(t(mirna))
  mirna_t$microRNA <- rownames(mirna_t)
  
  print(paste("start", length(rownames(mirna_t))))
  
  # Clean mirna names
  clean_mirna_name <- function(mirna_name) {
    # If not human mirna (different naming convention) don't clean name
    if (!grepl("hsa-let", mirna_name, fixed = TRUE) & 
        !grepl("hsa-mir", mirna_name, fixed = TRUE)) {
      return (mirna_name)
    }
    
    # Check if the string contains "3p" or "5p"
    suffix <- ""
    if (grepl("3p", mirna_name)) {
      suffix <- "-3p"
    } else if (grepl("5p", mirna_name)) {
      suffix <- "-5p"
    }
    
    # Use sub to replace the string after the third hyphen with an empty string
    truncated_string <- sub("^(([^-]+-){2}[^-]+).*", "\\1", mirna_name)
    
    # Concatenate the truncated string with the suffix
    result_string <- paste0(truncated_string, suffix)
    
    return(result_string)
  }
  
  mirna_t$microRNA <- sapply(mirna_t$microRNA, clean_mirna_name)
  
  # Identify duplicates
  duplicated_mirnas <- mirna_t$microRNA[duplicated(mirna_t$microRNA)]
  unique_duplicated_mirnas <- unique(duplicated_mirnas)
  
  if (length(unique_duplicated_mirnas) > 0) {
    print("Duplicate miRNAs:")
    print(unique_duplicated_mirnas)
  } else {
    print("No duplicate miRNAs found.")
  }
  
  # Group by miRNA names and calculate the median
  mirna_t_median <- as.data.table(mirna_t)
  mirna_t_median <- mirna_t_median[, lapply(.SD, median, na.rm = TRUE), by = microRNA]
  mirna_t_median <- as.data.frame(mirna_t_median)
  
  print(paste("end", length(rownames(mirna_t_median))))
  
  # Transpose back to original orientation
  mirna_cleaned <- as.data.frame(t(mirna_t_median))
  colnames(mirna_cleaned) <- mirna_cleaned[1, ]
  mirna_cleaned <- mirna_cleaned[-1, ]
  mirna_cleaned[] <- lapply(mirna_cleaned, as.numeric)
  
  return (mirna_cleaned)
}

# ================== Example Usage ================== 
# pure <- read.table("/Users/home/Desktop/Cancer Datasets/gbm/mirna")
# data <- as.data.frame(t(pure))
# cleaned_data <- clean_mirna(data, aggregateDifferentLocuses = TRUE)
