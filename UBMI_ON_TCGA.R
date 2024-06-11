library(tools)
library(ubmi)
library(patchwork)
library(ggplot2)
source("Overall Cleaning Script.R")

cancers <- c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
datasets <- c("exp", "methy", "mirna", "survival")
base_path <- "/Users/home/Desktop/Cancer Datasets"
output_path <- "/Users/home/TCGA/Cancer Plots"

# Create the main folder if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

timely_path <- file.path(output_path, format(Sys.time(), "%Y-%m-%d %H-%M-%S"))
dir.create(timely_path)

cleaned_data <- extract_and_clean_data(cancers, datasets, base_path)

# Initialize lists to store silhouette scores and UBMI objects
silhouette_scores <- list()
ubmi_results <- list()
  
for (cancer in cancers) {
  # Extract Data
  exp <- cleaned_data[[cancer]][["exp"]]
  methy <- cleaned_data[[cancer]][["methy"]]
  mirna <- cleaned_data[[cancer]][["mirna"]]
  
  cleaned_ubmi <- ubmi(list(exp, methy, mirna), combine_omics = TRUE)
  
  silhouette_scores[[cancer]] <- cleaned_ubmi@silhouette_score
  ubmi_results[[cancer]] <- cleaned_ubmi
  
  factor_plot <- ubmi::plot_factors(cleaned_ubmi) + 
    ggtitle(toTitleCase(paste0(cancer, ": Factor Plot")))
  
  metagene_plot <- ubmi::plot_metagenes(cleaned_ubmi) + 
    ggtitle(toTitleCase(paste0(cancer, ": Metagenes")))
  
  ubmi_grid_plot <- ubmi::plot_ubmi_grid(cleaned_ubmi) + 
    plot_annotation(title = (toTitleCase(paste(cancer, ": UBMI Grid"))))
  
  subfolder_path <- file.path(timely_path, cancer)
  
  # Create the subfolder if it doesn't exist
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path)
  }
  
  # Save the plots in the subfolder
  ggsave(file.path(subfolder_path, paste0(cancer, "_factor_plot.png")), plot = factor_plot, width = 10, height = 6)
  ggsave(file.path(subfolder_path, paste0(cancer, "_metagene_plot.png")), plot = metagene_plot, width = 10, height = 6)
  ggsave(file.path(subfolder_path, paste0(cancer, "_ubmi_grid_plot.png")), plot = ubmi_grid_plot, width = 15, height = 9)
  
  # Print statements to confirm that the plots are generated and saved
  print(paste("Plots saved for:", cancer))
}

# Create a data frame with silhouette scores
silhouette_scores_df <- data.frame(
  Cancer = names(silhouette_scores),
  SilhouetteScore = unlist(silhouette_scores)
)

# Sort the data frame by silhouette score
silhouette_scores_df <- silhouette_scores_df[order(silhouette_scores_df$SilhouetteScore, decreasing = TRUE), ]

# Print the ranking
print(silhouette_scores_df)

write.csv(silhouette_scores_df, file.path(timely_path, "silhouettes"), row.names = FALSE)

# Store the UBMI objects in a named list
sorted_ubmi_results <- ubmi_results[silhouette_scores_df$Cancer]

save(sorted_ubmi_results, file = file.path(timely_path, "ubmi_results"))