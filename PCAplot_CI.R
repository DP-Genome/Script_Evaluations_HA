#### PCA ####
library (ggplot2)
library(plotly)

# STEP 1: prepare evec data
evec_file <- read.table("~/Desktop/merged_ref_CI_labeled.evec",
                       header = FALSE, comment.char = "#", fill = TRUE)
pca_data = evec_file # copy data
# Assign column names
colnames(pca_data) <- c("ID", paste0("PC", 1:2), "Population", "Subpopulation")


# STEP 2: read CI data
ci_file <- "~/Desktop/expoCase"
ci_lines <- readLines(ci_file)

# Sample names
sample_lines <- ci_lines[grepl("^sample:", ci_lines)]
sample_names <- trimws(sub("^sample:\\s*", "", sample_lines))

# Ellipse coordinates
ell_lines <- ci_lines[grepl("^ellcoords:", ci_lines)]
ci_data <- do.call(rbind, lapply(ell_lines, function(x) {
    x <- sub("^ellcoords:\\s*", "", x)
    x <- sub("\\s*::.*$", "", x)
    as.numeric(strsplit(trimws(x), "\\s+")[[1]])}))
ci_data <- as.data.frame(ci_data)

colnames(ci_data) <- c("PC1","PC2","semi_major","semi_minor","angle")
ci_data$ID <- sample_names

# generate function
make_ellipse <- function(x, y, a, b, angle, n = 100) {
  theta <- seq(0, 2 * pi, length.out = n)
  
  # Ellipse before rotation
  ellipse <- rbind(a * cos(theta),b * sin(theta))
  
  # Rotation matrix
  rotation <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)),
    nrow = 2,byrow = TRUE)
  
  # Rotate and translate
  coords <- rotation %*% ellipse
  data.frame(PC1 = coords[1, ] + x, PC2 = coords[2, ] + y)
}

# create per sample using function above
ellipse_data <- do.call(
  rbind,
  lapply(1:nrow(ci_data), function(i) {
    e <- make_ellipse(
      x = ci_data$PC1[i],
      y = ci_data$PC2[i],
      a = ci_data$semi_major[i],
      b = ci_data$semi_minor[i],
      angle = ci_data$angle[i]
    )
    e$ID <- ci_data$ID[i]
    e
  })
)

# STEP 3: PLOT
x_col <- paste0("PC", 1)
y_col <- paste0("PC", 2)
  
pca_plot <- ggplot(pca_data, aes(x = .data[[x_col]], 
                           y = .data[[y_col]], 
                           color = Population, 
                           text = paste("ID:", ID, "<br>Pop:", Subpopulation))) +
geom_point(size = 1, alpha = 0.8) +
theme_minimal() +
labs(title = paste("HGDP PCA -", x_col, "vs", y_col),
     x = paste("Principal Component 1"),
     y = paste("Principal Component 2")) +
theme(plot.title = element_text(hjust = 0.5)) +
# add the confidence ellipses here
geom_path(data = ellipse_data,
     aes(x = PC1, y = PC2, group = ID), inherit.aes = FALSE,
     linewidth = 0.3, alpha = 0.35
)

print(pca_plot)
  
# Convert to plotly
int_plot_pca <- ggplotly(pca_plot, tooltip = "text")
print(int_plot_pca) 



