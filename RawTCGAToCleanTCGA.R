library(dplyr)

# This should be the directory of the user's path where the TCGA files are stored
working.directory <- "/Users/home/Desktop/Cancer Datasets"

raw_data_from_file_paths <- function(
    cancer, 
    working.directory, 
    data.types = c("exp", "methy", "mirna", "survival")
) {
  data.list <- list()
  
  for (data in data.types) {
    file.path <- file.path(working.directory, cancer, data)
    data.name <- paste(cancer, data, sep = ".")
    data.list[[data.name]] <- read.table(file.path)
  }
  print(paste('Finished reading: ', cancer))
  
  return (data.list)
}

fix_patient_id_names <- function(data) {
  for (name in names(data)) {
    data.set <- data[[name]]
    
    if (grepl("survival", name, fixed = TRUE)) {
      # Set the column names using the second row
      colnames(data.set) <- data.set[1,]
      data.set <- data.set[-1,]

      # Formatting of PatientID
      data.set$PatientID <- data.set$PatientID %>% 
        toupper() %>% 
        chartr(old = "-", new = ".") %>% 
        sub("^(([^.]*\\.){2}[^.]*).*", "\\1", .)
      
      data.set <- data.set %>% 
        distinct() %>% 
        filter(Survival != "" & Death != "")
    } else { # non-survival data table
      # Transpose for consistency
      data.set <- t(data.set)
      data.set <- as.data.frame(data.set)
    }
    
    data[[name]] <- data.set
  }
  
  return (data)
}

formatted.data <- list()
cancers <- c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
for (cancer in cancers) {
  # Extract data from files
  raw.data <- raw_data_from_file_paths(cancer, working.directory)
  
  # Patient ID formatting
  raw.data <- fix_patient_id_names(raw.data)
  # clean data
  
  # Place in bucket or smth
  formatted.data[[cancer]] <- raw.data
}
