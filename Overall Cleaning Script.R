library(dplyr)


#' Load Cleaning Functions
#'
#' Loads the necessary different cleaning functions for TCGA expression, 
#' methylation, and miRNA data to homogenize and clean all files. 
load_cleaning_functions <- function() {
  source("clean_exp.R")
  source("clean_methy.R")
  source("clean_mirna.R")
  print("Loading Functions Completed")
}

#' Extract and Clean TCGA Data
#'
#' @param cancers A list of strings naming which cancers to pull data from. 
#' Options include: aml, breast, colon, gbm, kidney, liver, lung, melanoma, 
#' ovarian, and sarcoma. 
#' @param datasets A list of strings naming which data tables to pull for each
#' cancer. Options include: exp, methy, mirna, survival (corresponds to gene 
#' expression, methylation, and micro RNA)
#' @param base_path The base path directory in which the TCGA data is stored 
#' (assumes all TCGA data was inserted into one folder), 
#' ex. /.../Desktop/Cancer Data
#'
#' @return A collection of cleaned data sets, represented as a list of list of 
#' data frames, arranged by cancer and cleaned according to data type
#' 
#' @examples 
#' cleaned_data_list <- extract_and_clean_data(cancers, datasets, base_path)
extract_and_clean_data <- function(cancers, datasets, base_path) {
  load_cleaning_functions()
  
  cleaned_data_list <- list()
  
  # Iterate over each dataset
  for (cancer_name in cancers) {
    # Create folder for cancer data
    cancer <- list()
    
    for (data_name in datasets) {
      # Construct the file path for the current dataset
      file_path <- paste(base_path, cancer_name, data_name, sep = "/")
      
      # Read the dataset
      pure_data <- read.table(file_path)
      
      # Clean the data based on the type of file
      cleaned_data <- switch(data_name, 
        "exp"      = {
          print("Reading expression data...")
          # Convert the data to a data frame and transpose it
          data <- as.data.frame(t(pure_data))
          
          # Clean the data using the clean_exp function and return it
          clean_exp(data)
        }, 
        "methy"    = {
          print("Reading methylation data...")
          # Convert the data to a data frame and transpose it
          data <- as.data.frame(t(pure_data))
          
          # Clean the data using the clean_exp function
          clean_methy(data)
        }, 
        "mirna"    = {
          print("Reading micro RNA data...")
          # Convert the data to a data frame and transpose it
          data <- as.data.frame(t(pure_data))
          
          # Clean the data using the clean_exp function
          clean_mirna(data)
        }, 
        "survival" = {
          print("Reading survival data...")
          # Set the column names using the second row
          data <- pure_data
          colnames(data) <- data[1,]
          data <- data[-1,]
          
          # Formatting of PatientID
          data$PatientID <- data$PatientID %>% 
            toupper() %>% 
            chartr(old = "-", new = ".") 
          # %>% sub("^(([^.]*\\.){2}[^.]*).*", "\\1", .) 
          # ^ if we want to remove the "month" parameter
          
          # Remove Duplicates
          cleaned_data <- data %>% 
            distinct() %>% 
            filter(Survival != "" & Death != "")
          
          cleaned_data
        }
      )
      
      cancer[[data_name]] <- cleaned_data
      print(paste0("Finished Reading: ", data_name))
    }
    
    # Store the cleaned data in the list
    cleaned_data_list[[cancer_name]] <- cancer
    print(paste0("Finished Reading: ", cancer_name))
  }
  return (cleaned_data_list)
}

## ================================================
## ============ Example Implementation ============
## ================================================
# cancers <- c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
# datasets <- c("exp", "methy", "mirna", "survival")
# base_path <- "/Users/home/Desktop/Cancer Datasets"
# 
# cleaned_data_list <- extract_and_clean_data(cancers, datasets, base_path)