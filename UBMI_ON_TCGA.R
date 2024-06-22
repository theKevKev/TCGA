source("Overall Cleaning Script.R")
source("Compute UBMI.R")
source("Plot Survival.R")

base_path <- "/Users/home/Desktop/Cancer Datasets"
output_path <- "/Users/home/TCGA/Cancer Plots"

cleaned_data <- extract_and_clean_data(base_path = base_path)

ubmi_list <- compute_ubmi(cleaned_data, output_path)

sorted_ubmi_results <- ubmi_list[1:(length(ubmi_list) - 1)]
file_path <- ubmi_list[[length(ubmi_list)]] # Possibly NULL, if save_data = FALSE

plot_survival(cleaned_data, sorted_ubmi_results, file_path = file_path)
