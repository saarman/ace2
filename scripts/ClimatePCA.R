# Load libraries
library(ggplot2)
library(ggfortify)  # For PCA visualization
library(factoextra) # For enhanced PCA plots
library(dplyr)      # Data manipulation
library(readr)

# Load each data set
hybrid_data <- read_csv("./data/summary_stats_5y_hyb.csv")
pipiens_data <- read_csv("./data/summary_stats_5y_pip.csv")
quinquefasciatus_data <- read_csv("./data/summary_stats_5y_qui.csv")

# Add species labels to each dataset
hybrid_data <- hybrid_data %>% mutate(species = "Hybrid")
pipiens_data <- pipiens_data %>% mutate(species = "Culex pipiens")
quinquefasciatus_data <- quinquefasciatus_data %>% mutate(species = "Culex quinquefasciatus")

# Merge all data into one dataframe
climate_data <- bind_rows(hybrid_data, pipiens_data, quinquefasciatus_data)

# Check structure
str(climate_data)


# Select only numeric climate variables (remove non-climate columns)
# Leaving latitude in as a proxy for autumn/spring photoperiod
climate_vars <- climate_data %>% select(-species, -lon, -year, -ID)

# Standardize climate variables (PCA works best with scaled data)
climate_vars_scaled <- scale(climate_vars)

# Run PCA
pca_result <- prcomp(climate_vars_scaled, center = TRUE, scale. = TRUE)

# Check variance explained by PCs
summary(pca_result)

# Scree plot to visualize explained variance
fviz_eig(pca_result)

# Add PCA results to dataset
pca_data <- as.data.frame(pca_result$x)
pca_data$species <- climate_data$species  # Retain species labels

# PCA Scatter Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = species)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA of Climate Data for *Culex pipiens*, *Culex quinquefasciatus*, and Hybrids",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% variance)")) +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "purple"))  # Adjust colors
