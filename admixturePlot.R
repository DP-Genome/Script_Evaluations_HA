#### ADMIXTURE ####
library(RColorBrewer)

k_values <- 5:14

Q_list <- lapply(k_values, function(k) {
  
  Q_file_path <- paste0("~/Desktop/merged_ref_pruned.", k, ".Q.withPop")
  if (!file.exists(Q_file_path)) return(NULL)
  
  Q_data <- read.table(Q_file_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(Q_data) <- c("sampleID", "subpopulation", "population", "gender", paste0("K", 1:k))
  
  # --- ORDER DATA ---
  # Separate HGDP and target samples to ensure HGDP are grouped by population
  # while target samples stay at the end
  is_target <- grepl(":", Q_data$sampleID) | !grepl("HGDP", Q_data$sampleID)
  Q_hgdp <- Q_data[!is_target, ]
  Q_targets <- Q_data[is_target, ]
  
  # Order HGDP by population
  Q_hgdp <- Q_hgdp[order(Q_hgdp$population, Q_hgdp$subpopulation), ]
  
  # Re-combine
  Q_data <- rbind(Q_hgdp, Q_targets)
  
  # --- SPLIT DATA FOR PLOTTING ---
  total_rows <- nrow(Q_data)
  n_targets <- nrow(Q_targets)
  
  # Main group (Left): HGDP samples
  Q_left <- Q_data[1:(total_rows - n_targets), ]
  # Subset group (Right): Target samples
  Q_right <- Q_data[(total_rows - n_targets + 1):total_rows, ]
  
  # Prepare Matrices without labels
  mat_left <- t(as.matrix(Q_left[, paste0("K", 1:k)]))
  mat_right <- t(as.matrix(Q_right[, paste0("K", 1:k)]))
    
  # --- PLOTTING SETUP ---
  # Create a layout: Left panel gets 70% width, Right panel gets 30% width
  layout(matrix(c(1, 2), nrow = 1), widths = c(8, 2))
  par(mar = c(5, 4, 4, 0)) # Reduce right margin for left plot
  
  cols <- colorRampPalette(brewer.pal(min(k, 12), "Paired"))(k)
  
  # Plot 1: Main (Left)
  barplot(mat_left, col = cols, border = NA, space = 0,
          ylab = "Ancestry proportion", xaxt = "n",
          main = paste("Admixture plot K =", k))
  # Labels for Left
  pop_l <- table(Q_left$population)
  axis(1, at = cumsum(pop_l)-k-(pop_l/2), labels = names(pop_l), tick = F, cex.axis = 0.5)

  # Plot 2: Subset (Right)
  par(mar = c(5, 0, 4, 2)) # margins for right plot
  b <- barplot(mat_right, col = cols, border = "white", lwd = 0.5, space = 0.15,
          xaxt = "n", yaxt = "n", main = "Target Samples")
   # Labels for Right (Sample IDs)
  axis(1, at = b, labels = Q_right$sampleID, 
       tick = F, las=2, cex.axis = 0.6)
  
  return(Q_data)
})