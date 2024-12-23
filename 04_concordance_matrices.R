#### run freemuxlet - unsupervised, then correlate with smart-seq data 
### run comparison per array 

## also include doublet category 

# Load necessary libraries
library(Seurat) # Assumed you're using Seurat objects
library(dplyr)
library(tidyr)
library(caret) # For confusion matrix 
library(psych) # For Cohen's kappa

# Assuming 'hto.all' is your Seurat object
# First, ensure 'sample_id' is present in meta data
if("sample_id" %in% names(hto.all@meta.data)) {
  # Split the data by 'sample_id'
  data_by_sample <- split(hto.all@meta.data, hto.all@meta.data$sample_id)
  
  # Initialize an empty list to store results for each sample_id
  results_list <- list()
  
  # Iterate over each subset of data
  for(sample_id in names(data_by_sample)) {
    current_data <- data_by_sample[[sample_id]]
    
    # Follow the steps as you described, adapting them for 'current_data'
    labels_comparison <- data.frame(
      'antibody_labels' = current_data$hash.ID,
      'genetic_labels' = current_data$SNG.BEST.GUESS,
      row.names = rownames(current_data)
    )
    
    labels_comparison$antibody_labels <- ifelse(labels_comparison$antibody_labels %in% c("Negative", "unmapped"), "AMB",
                                                ifelse(labels_comparison$antibody_labels == "Doublet", "DBL", labels_comparison$antibody_labels))
    
    labels_comparison$antibody_labels <- ifelse(labels_comparison$antibody_labels == "Pt2", "Pt3",
                                                ifelse(labels_comparison$antibody_labels == "Pt3", "Pt2", labels_comparison$antibody_labels))
    
    labels_comparison <- labels_comparison %>%
      filter(!(antibody_labels %in% c("AMB", "DBL") | genetic_labels %in% c("AMB", "DBL")))
    
    total_cells <- nrow(labels_comparison)
    concordant_cells <- sum(labels_comparison$antibody_labels == labels_comparison$genetic_labels)
    concordance_rate <- concordant_cells / total_cells * 100
    
    # Create concordance matrix
    concordance_matrix <- table(labels_comparison$antibody_labels, labels_comparison$genetic_labels)
    
    # Cohen's Kappa
    kappa_stat <- cohen.kappa(as.matrix(labels_comparison[,c('antibody_labels','genetic_labels')]))
    
    # Collecting results
    results_list[[sample_id]] <- list(
      TotalCells = total_cells,
      ConcordantCells = concordant_cells,
      ConcordanceRate = concordance_rate,
      ConcordanceMatrix = concordance_matrix,
      KappaStat = kappa_stat$kappa
    )
    
    # Optionally: Print results for each sample_id
    cat(paste("Results for sample_id:", sample_id, "\n"))
    cat("Total Cells:", total_cells, "\n")
    cat("Concordant Cells:", concordant_cells, "\n")
    cat("Concordance Rate:", concordance_rate, "%", "\n")
    cat("Cohen's Kappa for Agreement:", kappa_stat$kappa, "\n")
    print(concordance_matrix)
    cat("\n\n")
  }
  
  # You can access results for each sample_id from 'results_list'
  # Example: results_list$sample_id1
  # Plotting and further analysis can proceed here
  
} else {
  cat("sample_id is not present in the metadata.")
}
