datasets <- c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
cleaned_data_list <- list()
base_path <- "/Users/home/Desktop/Cancer Datasets/"

# Iterate over each dataset
for (dataset in datasets) {
  # Construct the file path for the current dataset
  file_path <- paste0(base_path, dataset, "/exp")
  
  # Read the dataset
  pure <- read.table(file_path)
  
  # Convert the data to a data frame and transpose it
  data <- as.data.frame(t(pure))
  
  # Clean the data using the clean_exp function
  cleaned_data <- clean_exp(data)
  
  # Store the cleaned data in the list
  cleaned_data_list[[dataset]] <- cleaned_data
}
