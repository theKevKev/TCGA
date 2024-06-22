library(dplyr)
library(stringr)
library(rlang)

#' Lenient Left Join
#'
#' Function to perform the join with an option for partial matching, capable of
#' adding one column from the second data frame. Attempts a basic left join of 
#' the second data frame, but if no matches occur or allow_partial_match = TRUE, 
#' performs a secondary join on remaining elements looking for PatientIDs ('by') 
#' that are substrings in the second column
#' 
#' @param df1 Main data frame, should have a column with a title string 
#' designated using the 'by' parameter
#' @param df2 Data frame with column to append, should have one column with 
#' title string 'by' and another with title string 'on'. 
#' @param by String to match joining, should be in both df1 and df2. Items in 
#' df2 under column 'by' can be partial strings of items in df1$by
#' @param on String of column to append to first data frame, present in df2
#' @param allow_partial_match toggleable, if TRUE, fills in blanks of initial 
#' join with the first parameter found in df2 that is a partial substring. 
#' Warning: if an item in df2$by is a short string or single character, could 
#' potentially match with too many items in df1. 
#'
#' @return a data frame, essentially df1 joined with df2$on
survival_left_join <- function(df1, df2, 
                               by = "PatientID", 
                               on = "Survival", 
                               allow_partial_match = FALSE) {
  # Step 1: Perform the regular left join
  df2 <- df2 %>% 
    select(by, on)
  
  result <- df1 %>%
    left_join(df2, by = by)
  
  # Check if there are no matches in the initial join
  if (all(is.na(result[[on]]))) {
    allow_partial_match <- TRUE
  }
  
  if (!allow_partial_match) {
    # If partial matching is not allowed, return the result as is
    return(result)
  }
  
  find_survival <- function(patient_id) {
    matching_survivals <- df2[[on]][str_detect(patient_id, df2[[by]])]
    if (length(matching_survivals) > 0) {
      return(matching_survivals[1])  # Return the first matching survival value
    } else {
      return(NA)  # Return NA if no match is found
    }
  }
  
  # Perform left join with substring matching
  result <- df1 %>%
    mutate(!!on := sapply(!!sym(by), find_survival))
  
  return(result)
}

#' Plot Overlay Survival
#'
#' Plots the patients in a factor plot of 2 dimensions (UBMI1, UBMI2), with
#' labeling by survival rather than by cluster
#'
#' @param surv_by_cluster A list generated during the plot_survival function 
#' call, should have a column named "clust" holding cluster identification, a 
#' column named "Survival". Survival should be numeric type
#' @param ubmi_object A UBMI Object as designated by the UBMI package
#' @param cancer A string for identification when saving data
#' @param file_path path to save file to, unneeded if save_data = FALSE
#' @param save_data toggle of whether to save data, if FALSE, will plot to 
#' console
plot_overlay_survival <- function(surv_by_cluster, ubmi_object, cancer, file_path, save_data = TRUE) {
  # Extract Survival column
  surv_vector <- as.numeric(surv_by_cluster$Survival)
  max_survival <- max(surv_vector, na.rm = TRUE)
  breaks <- c(-Inf, 0, as.integer(seq(0, max_survival, length.out = 8)) + 1, Inf)
  categories <- cut(surv_vector, breaks, include.lowest = TRUE)
  
  # Create factor plot with survival labels imposed
  survival_plot <- ubmi::plot_factors(ubmi_object, ad_hoc_label = categories) +
    ggtitle(toTitleCase(paste0(cancer, ": Survival Plot")))
  
  if (save_data) {
    if (!dir.exists(file_path)) {
      dir.create(file_path)
    }
    
    subfolder_path <- file.path(file_path, cancer)
    
    # Create the subfolder if it doesn't exist
    if (!dir.exists(subfolder_path)) {
      dir.create(subfolder_path)
    }
    
    ggsave(file.path(subfolder_path, paste0(cancer, "_survival_plot.png")), plot = survival_plot, width = 10, height = 6)
  } else {
    survival_plot
  }
}


#' Plot Kaplan-Meier Curves
#' 
#' Function that creates and saves (or plots) Kaplan-Meier Survival Curves
#'
#' @param surv_by_cluster A list generated during the plot_survival function 
#' call, should have a column named "clust" holding cluster identification, a 
#' column named "Survival", and a column named "Death". Survival and Death 
#' should follow standard conventions and be numeric type
#' @param cancer A string for identification when saving data
#' @param most_different_clusters togglable, if TRUE, will only plot the two 
#' survival curves with the minimum and maximum median survival
#' @param file_path path to save file to, unneeded if save_data = FALSE
#' @param save_data toggle of whether to save data, if FALSE, will plot to 
#' console
plot_kaplan_meier <- function(surv_by_cluster, 
                              cancer, 
                              most_different_clusters = FALSE, 
                              file_path, 
                              save_data = TRUE) {
  library(survival)
  library(survminer)
  library(ggplot2)
  
  noise_removed <- surv_by_cluster %>% 
    dplyr::filter(clust != 0) %>% # 0 represents noise
    tidyr::drop_na()
  
  if (most_different_clusters) {
    extremes <- noise_removed %>%
      dplyr::group_by(clust) %>%
      dplyr::summarise(median_surv = median(Survival, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(median_surv == min(median_surv) | median_surv == max(median_surv)) %>%
      dplyr::pull(clust)
    
    noise_removed <- noise_removed %>%
      dplyr::filter(clust %in% extremes)
  }
  
  noise_removed$clust <- factor(noise_removed$clust, levels = sort(unique(noise_removed$clust)))
  
  model_event <- survfit(Surv(Survival, Death) ~ clust, data = noise_removed)
  pval <- signif(surv_pvalue(model_event, data = noise_removed)$pval, digits = 3)
  
  kaplan_meier_plot <- ggsurvplot(
    model_event,
    palette = viridis::viridis(length(levels(noise_removed[["clust"]])), end = 0.8),
    surv.median.line = ifelse(most_different_clusters, "hv", "none"),
    data = noise_removed,
    legend.title = "Cluster",
    legend.labs = levels(noise_removed$clust)
  )
  
  if (most_different_clusters) {
    kaplan_meier_plot$plot <- kaplan_meier_plot$plot +
      labs(
        title = paste0("Cluster Distribution (N = ", 
                       paste0(table(noise_removed$clust), collapse = "/"), 
                       ")"),
        subtitle = paste0("pval = ", pval)
      )
  } else {
    kaplan_meier_plot$plot <- kaplan_meier_plot$plot +
      labs(
        title = paste0("Cluster Distribution (", 
                       length(levels(noise_removed$clust)), 
                       " clusters)"),
        subtitle = paste0("pval = ", pval)
      )
  }
  
  if (save_data) {
    if (!dir.exists(file_path)) {
      dir.create(file_path)
    }
    
    subfolder_path <- file.path(file_path, cancer)
    
    # Create the subfolder if it doesn't exist
    if (!dir.exists(subfolder_path)) {
      dir.create(subfolder_path)
    }
    
    if (most_different_clusters) {
      path <- file.path(subfolder_path, paste0(cancer, "_kaplan_meier_plot_most_different.png"))
    } else {
      path <- file.path(subfolder_path, paste0(cancer, "_kaplan_meier_plot.png"))
    }
    
    ggsave(path, plot = kaplan_meier_plot$plot, width = 10, height = 6)
    
  } else {
    kaplan_meier_plot
  }
}


#' Plot Survival
#'
#' @param cleaned_data cleaned data object generated by a cleaning script (of 
#' type list of lists)
#' @param sorted_ubmi_results UBMI results generated by "Compute UBMI.R", a 
#' list of ubmi objects sorted by silhouette score
#' @param only_first_n an integer for which of the first n object you would like
#' to plot, by default equal to 0 includes all plots
#' @param file_path output file path for survival plots, can be NA if save_data 
#' is FALSE
#' @param save_data whether to save data to a file
#' @param overlay_plot whether to include the overlay survival plot
#' @param kaplan_meier_plot whether to include the Kaplan-Meier (regular and 
#' only extremes) plot, organized by cluster
plot_survival <- function(cleaned_data, 
                          sorted_ubmi_results, 
                          only_first_n = 0, 
                          overlay_plot = TRUE, 
                          kaplan_meier_plot = TRUE, 
                          file_path, 
                          save_data = TRUE) {
  if (is.null(file_path)) {
    save_data <- FALSE
    warning("File path was NULL, save_data converted to FALSE")
  }
  
  if (!overlay_plot & !kaplan_meier_plot) {
    stop("What do you even want to plot then?")
  }
  
  num_plots <- length(sorted_ubmi_results)
  if (only_first_n > 0) {
    num_plots <- min(num_plots, only_first_n)
  }
  
  library(tibble)
  for (i in 1:num_plots) {
    cancer <- names(sorted_ubmi_results)[[i]]
    ubmi_object <- sorted_ubmi_results[[i]]
    surv <- cleaned_data[[cancer]][["survival"]]
    
    surv$Survival <- as.numeric(surv$Survival)
    surv$Death <- as.numeric(surv$Death)
    
    # Create Cluster | Survival | Death table
    clusters <- ubmi_object@factors %>% 
      select(clust) %>% 
      rownames_to_column(var = "PatientID")
    
    surv_by_cluster <- survival_left_join(clusters, surv, on = "Survival") %>%
      survival_left_join(., surv, on = "Death")
    
    # Overlay Plot
    if (overlay_plot) {
      plot_overlay_survival(surv_by_cluster, ubmi_object, cancer, file_path, save_data)
      print(paste("Overlay survival plot saved for:", cancer))
    }
    
    # Kaplan-Meier Plot
    if (kaplan_meier_plot) {
      plot_kaplan_meier(surv_by_cluster, cancer, most_different_clusters = TRUE, file_path, save_data)
      plot_kaplan_meier(surv_by_cluster, cancer, most_different_clusters = FALSE, file_path, save_data)
      print(paste("Kaplan-Meier plots saved for:", cancer))
    }
  }
  
  print("All Plots Computed! ")
}

