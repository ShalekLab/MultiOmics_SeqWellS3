### clear console 
ls()
rm(list=ls())

setwd("/Users/ddr13/Dropbox (MIT)/2023 Seq-Well Protocols Paper/protocol_paper/new_samples")

library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(clustree)
library(ggpubr)
library(patchwork)

metadata  <- read.csv("metadata3.csv", stringsAsFactors = FALSE)

# Function to process a sample (Part 1)
process_sample_part1 <- function(sample, umis, experiment, stim, project, id, id_factor, sample_id) {
  
  print(paste("Processing for sample", sample, "has started."))
  
  # Construct file paths
  mtx_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/matrix.mtx")
  features_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/features.tsv")
  cells_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/barcodes.tsv")
  hto_path <- paste0("outputs/HTO_", sample, "/umi_count/")
  
  # Load UMI matrix
  umis_filtered <- ReadMtx(mtx = mtx_path, features = features_path, cells = cells_path)
  
  # Find cells with non-zero expression and filter matrix
  keep_cells <- colnames(umis_filtered[, colSums(umis_filtered) > 0])
  umis_filtered <- umis_filtered[, keep_cells]
  
  # Load reads matrix
  umis_hto <- Read10X(hto_path, gene.column = 1)
  
  # Filter HTO data for cells with non-zero total counts
  keep_cells_hto <- colnames(umis_hto[, colSums(umis_hto) > 0])
  umis_hto <- umis_hto[, keep_cells_hto]
  
  # Select cell barcodes detected by both RNA and HTO
  joint_bcs <- intersect(colnames(umis_filtered), colnames(umis_hto))
  umis_filtered <- umis_filtered[, joint_bcs]
  
  # Generate metadata and perform operations
  num_cells <- length(colnames(umis_filtered))
  exp_data <- data.frame(exp = rep(experiment, num_cells),
                         stim = rep(stim, num_cells),
                         project = rep(project, num_cells),
                         id = rep(id, num_cells),
                         id_factor = rep(id_factor, num_cells),
                         sample_id = rep(sample_id, num_cells))
  
  # Setup Seurat object
  hashtag_object <- CreateSeuratObject(counts = umis_filtered, meta.data = exp_data, project = project)
  hashtag_object <- NormalizeData(hashtag_object)
  hashtag_object <- FindVariableFeatures(hashtag_object, selection.method = "vst")
  hashtag_object <- ScaleData(hashtag_object, features = VariableFeatures(hashtag_object))
  hashtag_object[["percent.mt"]] <- PercentageFeatureSet(hashtag_object, pattern = "^MT-")
  hashtag_object[["percent.ribo"]] <- PercentageFeatureSet(hashtag_object, pattern = "^RP[SL][[:digit:]]")
  hashtag_object <- subset(hashtag_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  
  
  # Save interim results if needed before Part 2
  interim_filename <- paste0(sample, "_interim_hashtag_object.rds")
  saveRDS(hashtag_object, file = interim_filename)
  print(paste("Interim Seurat object for sample", sample, "saved to file:", interim_filename))
  
  return(hashtag_object)
}

# Iterate over each row and process samples
results <- lapply(1:nrow(metadata), function(i) {
  row <- metadata[i, ]
  assign(paste0(row$sample, ".hashtag"), process_sample_part1(row$sample, row$umis, row$experiment, row$stim, row$project, row$id, row$id_factor, row$sample_id), envir = .GlobalEnv)
})

process_sample_part2 <- function(sample, umis) {
  # Load the interim Seurat object saved at the end of Part 1
  interim_filename <- paste0(sample, "_interim_hashtag_object.rds")
  hashtag_object <- readRDS(file = interim_filename)
  
  print(paste("Processing for sample", sample, "has started."))
  
  # Assuming 'umis_hto' needs to be loaded similarly to Part 1
  hto_path <- paste0("outputs/HTO_", sample, "/umi_count/")
  umis_hto <- Read10X(hto_path, gene.column = 1)
  rownames(umis_hto) <- gsub("-.*", "", rownames(umis_hto))
  
  # Filter for cells with non-zero total counts in the HTO data
  total_counts_hto <- Matrix::colSums(umis_hto)
  keep_cells_hto <- names(total_counts_hto[total_counts_hto > 0])
  umis_hto <- umis_hto[, keep_cells_hto]
  
  # Ensure matching cell barcodes between RNA and HTO assays
  joint_bcs <- intersect(colnames(hashtag_object), keep_cells_hto)
  hashtag_object <- subset(hashtag_object, cells = joint_bcs)
  umis_hto <- umis_hto[, joint_bcs]
  
  # Update the Seurat object with the filtered HTO data
  hashtag_object[["HTO"]] <- CreateAssayObject(counts = umis_hto)
  
  # Ensure HTO assay is properly set before HTODemux
  hashtag_object <- HTODemux(hashtag_object, assay = "HTO")
  
  # Print the completion of processing for the sample
  print(paste("Processing for sample", sample, "is now complete."))
  
  # Save the Seurat object to a file
  output_filename <- paste0(umis, ".hashtag.rds")
  saveRDS(hashtag_object, file = output_filename)
  print(paste("Saved Seurat object for sample", sample, "to file:", output_filename))
  
  return(hashtag_object)
}

# Script to iterate over Seurat objects from Part 1 and process each in Part 2
results <- lapply(1:nrow(metadata), function(i) {
  row <- metadata[i, ]
  sample <- row$sample
  umis <- row$umis
  # Process each sample using Part 2 function
  hashtag_object <- process_sample_part2(sample, umis)
  assign(paste0(sample, ".hashtag"), hashtag_object, envir = .GlobalEnv)
})

# Now call `process_all_samples_part2(metadata)` to process all samples through Part 2
####### RIDGE PLOTS ######
## GEX 1
library(cowplot)
# Ridge plot HTO
Idents(GEX_1.hashtag) <- 'HTO_maxID'
# Ensure that the identity variable is a factor and its levels are ordered alphabetically
GEX_1.hashtag$HTO_maxID <- factor(GEX_1.hashtag$HTO_maxID, levels = rev(sort(unique(GEX_1.hashtag$HTO_maxID))))

Idents(GEX_1.hashtag) <- "HTO_maxID"
ridge_list <- RidgePlot(GEX_1.hashtag, assay = "HTO", features = rownames(GEX_1.hashtag[["HTO"]]), 
                        ncol = 1, same.y.lims = TRUE, combine = FALSE)


# Apply the x-axis limits to each plot in the list
ridge_list <- lapply(ridge_list, function(p) {
  p + coord_cartesian(xlim = c(0, 750), expand = FALSE)
})

# Combine the plots into one plot
ridge_combined_GEX_1 <- plot_grid(plotlist = ridge_list, align = 'v', ncol = 1)

# Print the combined plot
ridge_combined_GEX_1

## GEX 2
# Ridge plot HTO
Idents(GEX_2.hashtag) <- 'HTO_maxID'
# Ensure that the identity variable is a factor and its levels are ordered alphabetically
GEX_2.hashtag$HTO_maxID <- factor(GEX_2.hashtag$HTO_maxID, levels = rev(sort(unique(GEX_2.hashtag$HTO_maxID))))

Idents(GEX_2.hashtag) <- "HTO_maxID"
ridge_list <- RidgePlot(GEX_2.hashtag, assay = "HTO", features = rownames(GEX_2.hashtag[["HTO"]]), 
                        ncol = 1, same.y.lims = TRUE, combine = FALSE)


# Apply the x-axis limits to each plot in the list
ridge_list <- lapply(ridge_list, function(p) {
  p + coord_cartesian(xlim = c(0, 750), expand = FALSE)
})

# Combine the plots into one plot
ridge_combined_GEX_2 <- plot_grid(plotlist = ridge_list, align = 'v', ncol = 1)

# Print the combined plot
ridge_combined_GEX_2
####### ADD ADT DATA #######
process_sample_adt <- function(sample, umis) {
  # Load the interim Seurat object saved at the end of Part 1
  interim_filename <- paste0(umis, ".hashtag.rds")
  hashtag_object <- readRDS(file = interim_filename)
  
  print(paste("Processing for sample", sample, "has started."))
  
  # Assuming 'umis_adt' needs to be loaded similarly to Part 1
  adt_path <- paste0("outputs/ADT_", sample, "/umi_count/")
  umis_adt <- Read10X(adt_path, gene.column = 1)
  rownames(umis_adt) <- gsub("-.*", "", rownames(umis_adt))
  
  # Filter for cells with non-zero total counts in the adt data
  total_counts_adt <- Matrix::colSums(umis_adt)
  keep_cells_adt <- names(total_counts_adt[total_counts_adt > 0])
  umis_adt <- umis_adt[, keep_cells_adt]
  
  # Ensure matching cell barcodes between RNA and adt assays
  joint_bcs <- intersect(colnames(hashtag_object), keep_cells_adt)
  hashtag_object <- subset(hashtag_object, cells = joint_bcs)
  umis_adt <- umis_adt[, joint_bcs]
  
  common_barcodes <- intersect(colnames(hashtag_object), keep_cells_adt)
  # Filter Seurat object to keep only common barcodes
  if (length(common_barcodes) > 0) {
    
    # Update the Seurat object with the filtered adt data
    hashtag_object[["ADT"]] <- CreateAssayObject(counts = umis_adt)
    
  } else {
    warning(paste("No common barcodes found for sample", sample))
  }
  
  # Print the completion of processing for the sample
  print(paste("Processing for sample", sample, "is now complete."))
  
  # Save the Seurat object to a file
  output_filename <- paste0(umis, ".hashtag.rds")
  saveRDS(hashtag_object, file = output_filename)
  print(paste("Saved Seurat object for sample", sample, "to file:", output_filename))
  
  return(hashtag_object)
}

# Script to iterate over Seurat objects from Part 2 and process ADTs
results <- lapply(1:nrow(metadata), function(i) {
  row <- metadata[i, ]
  sample <- row$sample
  umis <- row$umis
  # Process each sample using function
  hashtag_object <- process_sample_adt(sample, umis)
  assign(paste0(sample, ".hashtag"), hashtag_object, envir = .GlobalEnv)
})

# At this point, each Seurat object in your environment has ADT data added to it, where possible

# Function to extract nCount_RNA and return a data frame
extract_nCount_RNA <- function(seurat_obj, sample_name) {
  df <- data.frame(nCount_RNA = seurat_obj$nCount_RNA)
  df$sample <- sample_name
  return(df)
}

# Create an empty list to store data frames
df_list <- list()

# Iterate over each row in the metadata dataframe
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  seurat_obj_name <- paste0(sample_name, ".hashtag")
  
  # Check if the Seurat object exists in the global environment
  if (exists(seurat_obj_name, envir = .GlobalEnv)) {
    seurat_obj <- get(seurat_obj_name, envir = .GlobalEnv)
    df_list[[sample_name]] <- extract_nCount_RNA(seurat_obj, sample_name)
  } else {
    warning(paste("Seurat object not found for sample:", sample_name))
  }
}

# Combine all data frames into a single data frame
combined_df <- do.call(rbind, df_list)

# Plot reads per cell for each sample
ggplot(combined_df, aes(x = nCount_RNA, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")

pdf(paste0("reads_per_cell_array.pdf"), height = 15, width = 15)
ggplot(combined_df, aes(x = nCount_RNA, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")
dev.off()

# Function to extract nCount_HTO and return a data frame
extract_nCount_HTO <- function(seurat_obj, sample_name) {
  df <- data.frame(nCount_HTO = seurat_obj$nCount_HTO)
  df$sample <- sample_name
  return(df)
}

# Create an empty list to store data frames
df_list <- list()

# Iterate over each row in the metadata dataframe
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  seurat_obj_name <- paste0(sample_name, ".hashtag")
  
  # Check if the Seurat object exists in the global environment
  if (exists(seurat_obj_name, envir = .GlobalEnv)) {
    seurat_obj <- get(seurat_obj_name, envir = .GlobalEnv)
    df_list[[sample_name]] <- extract_nCount_HTO(seurat_obj, sample_name)
  } else {
    warning(paste("Seurat object not found for sample:", sample_name))
  }
}

# Combine all data frames into a single data frame
combined_df <- do.call(rbind, df_list)

# Plot reads per cell for each sample
ggplot(combined_df, aes(x = nCount_HTO, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")

pdf(paste0("reads_per_cell_HTO.pdf"), height = 15, width = 15)
ggplot(combined_df, aes(x = nCount_HTO, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")
dev.off()

# Function to extract nCount_ADT and return a data frame
extract_nCount_ADT <- function(seurat_obj, sample_name) {
  df <- data.frame(nCount_ADT = seurat_obj$nCount_ADT)
  df$sample <- sample_name
  return(df)
}

# Create an empty list to store data frames
df_list <- list()

# Iterate over each row in the metadata dataframe
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  seurat_obj_name <- paste0(sample_name, ".hashtag")
  
  # Check if the Seurat object exists in the global environment
  if (exists(seurat_obj_name, envir = .GlobalEnv)) {
    seurat_obj <- get(seurat_obj_name, envir = .GlobalEnv)
    df_list[[sample_name]] <- extract_nCount_ADT(seurat_obj, sample_name)
  } else {
    warning(paste("Seurat object not found for sample:", sample_name))
  }
}

# Combine all data frames into a single data frame
combined_df <- do.call(rbind, df_list)

# Plot reads per cell for each sample
ggplot(combined_df, aes(x = nCount_ADT, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")

pdf(paste0("reads_per_cell_ADT.pdf"), height = 15, width = 15)
ggplot(combined_df, aes(x = nCount_ADT, fill = sample)) +
  geom_histogram(bins = 30, alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Reads Per Cell for Each Sample",
       x = "Reads per Cell (log scale)",
       y = "Frequency",
       fill = "Sample") +
  facet_wrap(~sample, scales = "free_y")
dev.off()

# Initialize a list to store the processed Seurat objects
seurat_objects_list <- list()
# Initialize a vector to store the id_factors for each object
id_factors <- c()

# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  umis <- metadata$umis[i]
  
  # Construct the file name, load the Seurat object, process it
  file_name <- paste0(umis, ".hashtag.rds")
  singlet_object <- readRDS(file = file_name)
  
  DefaultAssay(singlet_object) <- "RNA"
  Idents(singlet_object) <- "orig.ident"
  singlet_object <- NormalizeData(singlet_object, verbose = FALSE)
  singlet_object <- FindVariableFeatures(singlet_object, selection.method = "vst", nfeatures = 2000)
  
  # Extract the 'id_factor' from the Seurat object metadata
  if("id_factor" %in% names(singlet_object@meta.data)) {
    id_factor <- unique(singlet_object@meta.data$id_factor)
    if(length(id_factor) != 1) {
      stop("Each Seurat object must have a unique 'id_factor' for all its cells.")
    }
    id_factors <- c(id_factors, id_factor)
  } else {
    stop("Column 'id_factor' not found in Seurat object metadata.")
  }
  
  # Store the processed object in the list with its sample name as the key
  seurat_objects_list[[sample_name]] <- singlet_object
}

# Prepare the objects for merging by removing the first object's name and id_factor from the lists
first_object <- seurat_objects_list[[1]]
remaining_objects <- seurat_objects_list[-1]
remaining_ids <- id_factors[-1]

# Perform the merge in one step using the prepared lists of objects and IDs
merged_seurat_object <- merge(first_object, y = remaining_objects, add.cell.ids = id_factors, project = "MergedSeurat")

# Save the merged Seurat object
saveRDS(merged_seurat_object, file = "merged_seurat_object_HTO_preprocessing.rds")

hto.all <- readRDS(file = paste0("merged_seurat_object_HTO_preprocessing.rds"))

########################## ANALYSIS ########################
library(cowplot)
# Ridge plot HTO
Idents(hto.all) <- 'HTO_maxID'
# Ensure that the identity variable is a factor and its levels are ordered alphabetically
hto.all$HTO_maxID <- factor(hto.all$HTO_maxID, levels = rev(sort(unique(hto.all$HTO_maxID))))

Idents(hto.all) <- "HTO_maxID"
ridge_list <- RidgePlot(hto.all, assay = "HTO", features = rownames(hto.all[["HTO"]]), 
                        ncol = 1, same.y.lims = TRUE, combine = FALSE)


# Apply the x-axis limits to each plot in the list
ridge_list <- lapply(ridge_list, function(p) {
  p + coord_cartesian(xlim = c(0, 750), expand = FALSE)
})

# Combine the plots into one plot
ridge_combined <- plot_grid(plotlist = ridge_list, align = 'v', ncol = 1)

# Print the combined plot
ridge_combined

# Counts of each stim condition for all singlets by patient
stim.cond.counts.HTO.all.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.HTO.all.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.HTO.all.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  theme_classic2()
print(stim.cond.counts.HTO.all.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.HTO.all.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.HTO.all.array)

pdf(paste0("bar_stim_cond_all_HTO_sample.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.HTO.all.array, stim.cond.proportion.HTO.all.array, common.legend = TRUE)
dev.off()

pdf(paste0("bar_stim_cond_count_HTO_sample.pdf"), height = 8, width = 8)
stim.cond.counts.HTO.all.array
dev.off()