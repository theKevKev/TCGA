library(dplyr)

load_cleaning_functions <- function() {
  source("clean_exp.R")
  source("clean_methy.R")
  source("clean_mirna.R")
}

cancers <- c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
datasets <- c("exp", "methy", "mirna", "survival")
cleaned_data_list <- list()
base_path <- "/Users/home/Desktop/Cancer Datasets"

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
        # Convert the data to a data frame and transpose it
        data <- as.data.frame(t(pure_data))
        
        # Clean the data using the clean_exp function and return it
        clean_exp(data)
      }, 
      "methy"    = {
        # Convert the data to a data frame and transpose it
        data <- as.data.frame(t(pure_data))
        
        # Clean the data using the clean_exp function
        clean_methy(data)
      }, 
      "mirna"    = {
        # Convert the data to a data frame and transpose it
        data <- as.data.frame(t(pure_data))
        
        # Clean the data using the clean_exp function
        clean_mirna(data)
      }, 
      "survival" = {
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
  }
  
  # Store the cleaned data in the list
  cleaned_data_list[[cancer_name]] <- cancer
}

