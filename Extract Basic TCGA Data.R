library(dplyr)

# This should be the directory of the user's path where the TCGA files are stored
working.directory <- "/Users/home/Desktop/Cancer Datasets"

#' Extract TCGA Data
#'
#' @param working.directory a string of the path where the TCGA files are stored
#' @param cancer.types the set of cancers to use, as a list of strings
#' @param data.types the set of files to extract per cancer, as a list of strings
#'
#' @return a list of data frames of the raw data
tcga_data_from_files <- function(
    working.directory, 
    cancer.types = c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma"), 
    data.types = c("exp", "methy", "mirna", "survival")
) {
  data.list <- list()
  
  for (cancer in cancer.types) {
    for (data in data.types) {
      file.path <- file.path(working.directory, cancer, data)
      data.name <- paste(cancer, data, sep = ".")
      data.list[[data.name]] <- read.table(file.path)
    }
    print(paste('Finished reading: ', cancer))
  }
  
  return(data.list)
}

rename_tcga_data <- function(tcga.data) {
  for (name in names(tcga.data)) {
    raw.data <- tcga.data[[name]]
    
    if(grepl("survival", name, fixed = TRUE)) {
      # Set the column names using the second row
      colnames(raw.data) <- raw.data[1,]
      raw.data <- raw.data[-1,]
    } else { # non-survival data table
      # Transpose for consistency
      raw.data <- t(raw.data)
      raw.data <- as.data.frame(raw.data)
    }
    
    tcga.data[[name]] <- raw.data
  }
  
  return(tcga.data)
}

tcga.data <- tcga_data_from_files(working.directory)
renamed.tcga.data <- rename_tcga_data(tcga.data)
