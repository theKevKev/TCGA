library(umap)
library(plyr)
library(ubmi)
library(factoextra)

working_directory = "/Users/home/Desktop/Cancer Datasets"

# Uncleaned UMAP
data.list <- list()

# Get Uncleaned Data
for (data in c("exp", "methy", "mirna")) {
  file.path <- file.path(working_directory, "sarcoma", data)
  data.name <- paste("sarcoma", data, sep = ".")
  data.list[[data.name]] <- read.table(file.path)
}

# Remove Null Values
for (data_name in names(data.list)) {
  data <- data.list[[data_name]]
  data[is.na(data)] <- 0
  data <- as.data.frame(t(data))
  data.list[[data_name]] <- data
}

# Inner join the three data sets
uf1 <- data.list[["sarcoma.exp"]]
uf2 <- data.list[["sarcoma.methy"]]
uf3 <- data.list[["sarcoma.mirna"]]

uncleaned_ubmi <- ubmi(list(uf1, uf2, uf3))

uf1$row_name <- rownames(uf1)
uf2$row_name <- rownames(uf2)
uf3$row_name <- rownames(uf3)

uncleaned_combined <- join_all(list(uf1, uf2, uf3), by = 'row_name', type = 'inner')
rownames(uncleaned_combined) <- uncleaned_combined$row_name
uncleaned_combined$row_name <- NULL

uncleaned_umap <- umap(uncleaned_combined)
uf123 <- as.data.frame(uncleaned_umap$layout)

uncleaned_pca <- prcomp(uncleaned_combined)
fviz_eig(uncleaned_pca)
fviz_pca_ind(uncleaned_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             geom = 'point'
)

# Collect Survival Data
uf123$PatientID <- rownames(uf123)
uf4 <- read.table(file.path(working_directory, "sarcoma", "survival"))
colnames(uf4) <- uf4[1,]
uf4 <- uf4[-1,]

# Formatting of PatientID
uf4$PatientID <- uf4$PatientID %>% 
  toupper() %>% 
  chartr(old = "-", new = ".") 

uf5 <- merge(uf123, uf4, by = 'PatientID')
uf5$Survival <- as.numeric(uf5$Survival)
uf5$Survival <- as.numeric(lapply(uf5$Survival, function(x) log(x)))

uf <- data.frame(x = uncleaned_umap$layout[,1],
                 y = uncleaned_umap$layout[,2], 
                 Species = uf5$Death)

ggplot(uf, aes(x, y, colour = Species)) +
  geom_point()

## Plotting Uncleaned UBMI Results
df <- data.frame(x = uncleaned_ubmi@factors[,1],
                 y = uncleaned_ubmi@factors[,2], 
                 Species = uncleaned_ubmi@factors[,3])

ggplot(df, aes(x, y, colour = Species)) +
  geom_point()


## ============================================================================
## ================== Analysis of Cleaned and Combined Data ===================
## ============================================================================


# Get Cleaned Data for Sarcoma
source("Overall Cleaning Script.R")
sarcoma_data <- extract_and_clean_data("sarcoma", c("exp", "methy", "mirna", "survival"), working_directory)

# Inner join the three data sets
df1 <- sarcoma_data[["sarcoma"]][["exp"]]
df2 <- sarcoma_data[["sarcoma"]][["methy"]]
df3 <- sarcoma_data[["sarcoma"]][["mirna"]]

cleaned_ubmi <- ubmi(list(df1, df2, df3))

df1$row_name <- rownames(df1)
df2$row_name <- rownames(df2)
df3$row_name <- rownames(df3)

sarcoma_combined <- join_all(list(df1, df2, df3), by = 'row_name', type = 'inner')
rownames(sarcoma_combined) <- sarcoma_combined$row_name
sarcoma_combined$row_name <- NULL

# UMAP
sarcoma_umap <- umap(sarcoma_combined)
df123 <- as.data.frame(sarcoma_umap$layout)

# PCA
cleaned_pca <- prcomp(sarcoma_combined)
fviz_eig(cleaned_pca)
fviz_pca_ind(cleaned_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             geom = 'point'
)

# Collect Survival Data
df123$PatientID <- rownames(df123)
df4 <- sarcoma_data[["sarcoma"]][["survival"]]
df5 <- merge(df123, df4, by = 'PatientID')
df5$Survival <- as.numeric(df5$Survival)
df5$Survival <- as.numeric(lapply(df5$Survival, function(x) log(x)))

df <- data.frame(x = sarcoma_umap$layout[,1],
                 y = sarcoma_umap$layout[,2], 
                 Species = df5$Death)

ggplot(df, aes(x, y, colour = Species)) +
  geom_point()

## Plotting Cleaned UBMI Results
df <- data.frame(x = cleaned_ubmi@factors[,1],
                 y = cleaned_ubmi@factors[,2], 
                 Species = cleaned_ubmi@factors[,3])

ggplot(df, aes(x, y, colour = Species)) +
  geom_point()
