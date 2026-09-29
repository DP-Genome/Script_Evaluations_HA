#### PCA ####
library (ggplot2)
library(plotly)
library(uwot)

# STEP 1: prepare data
evec_file <- read.table("~/Desktop/merged_ref_labeled.evec",
                       header = FALSE, comment.char = "#", fill = TRUE)
pca_data = evec_file # copy data
# Assign column names
colnames(pca_data) <- c("ID", paste0("PC", 1:10), "Population", "Subpopulation")

# STEP 2: loop through PC pairs and plot
for (firstPC in seq(1, 9, by = 2)) {
  secondPC <- firstPC + 1
  # Define column names as strings
  x_col <- paste0("PC", firstPC)
  y_col <- paste0("PC", secondPC)
  
  # Create the plot
  pca_plot <- ggplot(pca_data, aes(x = .data[[x_col]], 
                               y = .data[[y_col]], 
                               color = Population, 
                               text = paste("ID:", ID, "<br>Pop:", Subpopulation))) +
    geom_point(size = 1, alpha = 0.8) +
    theme_minimal() +
    labs(title = paste("HGDP PCA -", x_col, "vs", y_col),
         x = paste("Principal Component", firstPC),
         y = paste("Principal Component", secondPC)) +
    theme(plot.title = element_text(hjust = 0.5))
  # print static plot
  print(pca_plot)
  
  # Convert to plotly
  int_plot_pca <- ggplotly(pca_plot, tooltip = "text")
  print(int_plot_pca) 
}

# df_pcs should contain only the PC columns
neighbors = 35
min_distance = 0.1
umap_result <- umap(pca_data, 
                    n_neighbors = neighbors,
                    min_dist = min_distance,
                    metric = "euclidean", 
                    verbose = TRUE)

# Convert to dataframe for plotting
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
# add subpopulation data
umap_df$subpop <- pca_data$Subpopulation
# plot
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,
                         col = pca_data$Population,
                         text = paste("<br>SubPop:", subpop))) +
  geom_point(size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = paste ("UMAP of 10 PCs - ",neighbors,
                      " neighbors and min distance of ",min_distance),
       x = "UMAP1",
       y = "UMAP2")
print (umap_plot)
# create interactive plot
int_plot_umap <- ggplotly(umap_plot, tooltip = "text")
print (int_plot_umap)

