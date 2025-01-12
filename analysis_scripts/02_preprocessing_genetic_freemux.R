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

metadata <- read.csv("metadata_genetic2.csv", stringsAsFactors = FALSE)

# Function to process a sample (Part 1)
process_sample_part1 <- function(sample, umis, experiment, stim, project, id, id_factor, sample_id) {
  
  print(paste("Processing for sample", sample, "has started."))
  
  # Construct file paths
  mtx_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/matrix.mtx")
  features_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/features.tsv")
  cells_path <- paste0("starsolo_out/", sample, "/Solo.out/GeneFull/filtered/barcodes.tsv")
  freemuxlet_path <- paste0("freemux_Popscle/", sample, "/demuxlet_merged.clust1.samples")
  
  # Load UMI matrix
  umis_filtered <- ReadMtx(mtx = mtx_path, features = features_path, cells = cells_path)
  
  #Load freemuxlet
  freemuxlet_delim <- read.delim(freemuxlet_path, check.names = FALSE, header = TRUE, row.names = 1)
  rownames(freemuxlet_delim) <- freemuxlet_delim$BARCODE

  # Function to replace 2, 1, 0 with Patient1, Patient2, and Patient3 respectively
  replace_patient_terms <- function(column) {
    # Replace each value based on specific matches
    column <- ifelse(column == "2", "Patient1", column)
    column <- ifelse(column == "1", "Patient2", column)
    column <- ifelse(column == "0", "Patient3", column)
    
    return(column)
  }
  
  # Apply the function only to columns 5, 7, 12, 14, and 17 of demuxlet
 columns_to_modify <- c(5, 7, 12, 14, 17)
freemuxlet_delim[columns_to_modify] <- lapply(freemuxlet_delim[columns_to_modify], replace_patient_terms)
  
  # Find cells with non-zero expression and filter matrix
  keep_cells <- colnames(umis_filtered[, colSums(umis_filtered) > 0])
  umis_filtered <- umis_filtered[, keep_cells]
  joint.bcs <- intersect(colnames(umis_filtered), freemuxlet_delim$BARCODE)
  umis_filtered <- umis_filtered[, joint.bcs]
  freemuxlet <- freemuxlet_delim[joint.bcs, ]
  
  # Generate metadata and perform operations
  num_cells <- length(colnames(umis_filtered))
  exp_data <- data.frame(sample = rep(sample, num_cells),
                         exp = rep(experiment, num_cells),
                         stim = rep(stim, num_cells),
                         project = rep(project, num_cells),
                         id = rep(id, num_cells),
                         id_factor = rep(id_factor, num_cells))
  exp_data <- cbind(exp_data, freemuxlet)
  
  for(i in 1:nrow(exp_data)) {
    if(exp_data$DROPLET.TYPE[i] %in% c("DBL", "AMB")) {
      exp_data$SNG.BEST.GUESS[i] <- exp_data$DROPLET.TYPE[i]
    }
  }
  
  # Setup Seurat object
  hashtag_object <- CreateSeuratObject(counts = umis_filtered, meta.data = exp_data, project = project)
  hashtag_object <- NormalizeData(hashtag_object)
  hashtag_object <- FindVariableFeatures(hashtag_object, selection.method = "mean.var.plot")
  hashtag_object <- ScaleData(hashtag_object, features = VariableFeatures(hashtag_object))
  hashtag_object[["percent.mt"]] <- PercentageFeatureSet(hashtag_object, pattern = "^MT-")
  hashtag_object[["percent.ribo"]] <- PercentageFeatureSet(hashtag_object, pattern = "^RP[SL][[:digit:]]")
 # hashtag_object <- subset(hashtag_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  #save freemuxlet metadata 
  #create the data.frame
  freemux_meta <- data.frame(hashtag_object@meta.data)
  freemux_meta <- freemux_meta[,-c(1:8,28:29)]
  
  # Construct a path for the CSV file
  csv_file_path <- paste0(sample, "_freemux_metadata_all.csv")
  
  # Use write.csv to save the data frame to a CSV file
  write.csv(freemux_meta, csv_file_path, row.names = TRUE)
  
  # Save interim results if needed before Part 2
  interim_filename <- paste0(sample, "_interim_hashtag_object_freemuxlet.rds")
  saveRDS(hashtag_object, file = interim_filename)
  print(paste("Interim Seurat object for sample", sample, "saved to file:", interim_filename))
  
  return(hashtag_object)
  
}

# Iterate over each row and process samples
results <- lapply(1:nrow(metadata), function(i) {
  row <- metadata[i, ]
  assign(paste0(row$sample, ".hashtag"), process_sample_part1(row$sample, row$umis, row$experiment, row$stim, row$project, row$id, row$id_factor, row$sample_id), envir = .GlobalEnv)
})

# Doublet Removal
library(DoubletFinder)
library(readr)
# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  
  sample <- metadata$sample[i]
  umis <- metadata$umis[i]
  
  # Construct the filename for the Seurat object
  filename <- paste0(sample, "_interim_hashtag_object_freemuxlet.rds")
  
  # Load the Seurat object from the RDS file
  seurat_object <- readRDS(filename)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  seurat_object <- RunUMAP(seurat_object, dims = 1:10)
  
  # Define the conversion mapping
  conversion_mapping <- c(
    "SNG" = "Singlet",
    "DBL" = "Doublet",
    "AMB" = "Doublet")
  
  df <- as.data.frame(seurat_object@meta.data)
  # Convert the 'id' column values and add them as a separate column
  df <- df %>%
    mutate(GT = as.factor(conversion_mapping[DROPLET.TYPE]))
  # Extract the 'GT' column
  new_meta <- subset(df, select = c('GT'))
  # Add metadata back to the Seurat object
  seurat_object <- AddMetaData(seurat_object, new_meta)
  
  sweep.res.list <- paramSweep(seurat_object, PCs = 1:10, sct = FALSE)
  gt.calls <- seurat_object@meta.data[rownames(sweep.res.list[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
  bcmvn <- find.pK(sweep.stats)
  plot <- ggplot(bcmvn, aes(pK, BCmetric)) +
    geom_point()
  
  ggsave(plot, filename = paste0(sample, "_pKvBCmetric.png"))
  
  ## Run DoubletFinder with varying classification stringencies 
  seurat_object <- doubletFinder(seurat_object, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])), nExp = 0.1482, reuse.pANN = FALSE, sct = FALSE)
  
  doublets <- as.data.frame(cbind(colnames(seurat_object), seurat_object@meta.data[,grepl(paste0("pANN_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat_object@meta.data))], seurat_object@meta.data[,grepl(paste0("DF.classifications_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat_object@meta.data))]))
  
  colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
  
  doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)
  
  ### Calculate number of doublets and singlets ###
  summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
  colnames(summary) <- c("Classification", "Droplet N")
  write_delim(summary, paste0(sample, "_DoubletFinder_doublet_summary.tsv"), "\t")
  
  # Optionally save the singlet object
  saveRDS(seurat_object, file = paste0(sample, ".doublet.hashtag.rds"))
  
}

# Normalize and Identify Variable Features for each dataset
# Assuming you are using the same feature set for all objects

# Initialize a list to store the processed Seurat objects
seurat_objects_list <- list()
# Initialize a vector to store the id_factors for each object
id_factors <- c()

# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  
  # Construct the file name, load the Seurat object, process it
  file_name <- paste0(sample_name, ".doublet.hashtag.rds")
  seurat_object <- readRDS(file = file_name)
  
  # Extract the 'id_factor' from the Seurat object metadata
  if("id_factor" %in% names(seurat_object@meta.data)) {
    id_factor <- unique(seurat_object@meta.data$id_factor)
    if(length(id_factor) != 1) {
      stop("Each Seurat object must have a unique 'id_factor' for all its cells.")
    }
    id_factors <- c(id_factors, id_factor)
  } else {
    stop("Column 'id_factor' not found in Seurat object metadata.")
  }
  
  # Store the processed object in the list with its sample name as the key
  seurat_objects_list[[sample_name]] <- seurat_object
}

# Prepare the objects for merging by removing the first object's name and id_factor from the lists
first_object <- seurat_objects_list[[1]]
remaining_objects <- seurat_objects_list[-1]
remaining_ids <- id_factors[-1]

# Perform the merge in one step using the prepared lists of objects and IDs
merged_seurat_object <- merge(first_object, y = remaining_objects, add.cell.ids = id_factors, project = "MergedSeurat")

merged_seurat_object2 <- subset(merged_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

all_meta <- merged_seurat_object2@meta.data

labels_comparison <- data.frame(
  'id_factor' = all_meta$id_factor,
  'barcode' = all_meta$BARCODE,
  'droplet_type' = all_meta$DROPLET.TYPE,
  'best_guess' = all_meta$SNG.BEST.GUESS,
  row.names = rownames(all_meta))

# Count occurrences of each term in text_column
term_counts <- labels_comparison %>%
  group_by(best_guess) %>%
  summarise(count = n(), .groups = 'drop')

# Print the counts of each term
print(term_counts)

# Save the merged Seurat object
saveRDS(merged_seurat_object, file = "merged_seurat_object_genetic_freemux.rds")

gen.all <- readRDS(file = paste0("merged_seurat_object_genetic_freemux.rds"))

########################## ANALYSIS ########################
# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=SNG.BEST.GUESS)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  theme_classic2()
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=SNG.BEST.GUESS)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_all_sample_all_genetic_freemux.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

pdf(paste0("bar_stim_cond_counts_only_all_genetic_freemux.pdf"), height = 8, width = 8)
stim.cond.counts.sing.array
dev.off()

#### Prep genetic data for integration with gen 
# Define a function to create a table and add CLASS column
process_data <- function(sample) {
  
  # Construct the filename for the Seurat object
  filename <- paste0(sample, "_interim_hashtag_object_freemuxlet.rds")
  
  # Load the Seurat object from the RDS file
  seurat_object <- readRDS(filename)
  freemuxlet <- seurat_object@meta.data
  
  cat(paste("Table:", sample, "\n"))
  print(table(freemuxlet$SNG.BEST.GUESS, freemuxlet$DROPLET.TYPE))
  
  # Convert the table to a data frame for easier saving
  table_df <- as.data.frame.matrix(table(freemuxlet$SNG.BEST.GUESS, freemuxlet$DROPLET.TYPE))
  # Construct a path for the CSV file
  csv_file_path <- paste0(sample, "_genetic_assignment_table.csv")
  # Use write.csv to save the data frame to a CSV file
  write.csv(table_df, csv_file_path, row.names = TRUE)
  
  freemuxlet$CLASS <- freemuxlet$SNG.BEST.GUESS
  freemuxlet$CLASS[freemuxlet$DROPLET.TYPE == "DBL"] <- "DBL"
  
  return(table(freemuxlet$SNG.BEST.GUESS, freemuxlet$DROPLET.TYPE))
}

# Initialize an empty data frame to store aggregated data
aggregated_data <- data.frame()

# Inside your loop or function, after creating each table
for(i in 1:nrow(metadata)) {
  row <- metadata[i, ]
  sample_name <- row$sample 
  
  # Process the data to get the table
  table_data <- process_data(sample_name) 
  
  # Convert the table to a data frame and add a sample prefix to row names
  table_df <- as.data.frame.matrix(table_data)
  row_names_with_prefix <- paste(sample_name, rownames(table_df), sep = "_")
  rownames(table_df) <- row_names_with_prefix
  
  # Append this sample's data to the aggregated data frame
  aggregated_data <- rbind(aggregated_data, table_df)
}

# Write the processed metadata to a CSV file
write.csv(aggregated_data, "all_genetic_assignment_table_freemux.csv", row.names = TRUE)

# Initialize a list to store singlet objects
singlet_objects <- list()

# Loop over each sample
for (i in 1:nrow(metadata)) {
  # Extract the current sample name
  sample <- metadata$sample[i]
  
  # Construct the filename for the Seurat object
  filename <- paste0(sample, "_interim_hashtag_object_freemuxlet.rds")
  
  # Load the Seurat object from the RDS file
  seurat_object <- readRDS(filename)
  seurat_object@meta.data$DROPLET.TYPE <- as.factor(seurat_object@meta.data$DROPLET.TYPE)
  #Set identity to be droplet type
  Idents(seurat_object) <- "DROPLET.TYPE"
  
  # Keep only singlet cells
  sample_hashtag.singlet <- subset(seurat_object, idents = "SNG")
  
  ## Save df with genetic demulti metadata
  #create the data.frame
  freemux_meta <- data.frame(sample_hashtag.singlet@meta.data)
  freemux_meta <- freemux_meta[,-c(1:8,28:29)]
  
  # Construct a path for the CSV file
  csv_file_path <- paste0(sample, "_freemux_metadata_sing.csv")
  
  # Use write.csv to save the data frame to a CSV file
  write.csv(freemux_meta, csv_file_path, row.names = TRUE)
  
  # Add the processed singlet object to the list
  singlet_objects[[sample]] <- sample_hashtag.singlet
  
  # Optionally save the singlet object
  saveRDS(sample_hashtag.singlet, file = paste0(sample, ".hashtag.singlet.freemuxlet.rds"))
  
}

# Normalize and Identify Variable Features for each dataset
# Assuming you are using the same feature set for all objects

# Initialize a list to store the processed Seurat objects
seurat_objects_list <- list()
# Initialize a vector to store the id_factors for each object
id_factors <- c()

# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  
  # Construct the file name, load the Seurat object, process it
  file_name <- paste0(sample_name, ".hashtag.singlet.freemuxlet.rds")
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
saveRDS(merged_seurat_object, file = "merged_seurat_object_genetic_freemux.rds")

gen.all <- readRDS(file = paste0("merged_seurat_object_genetic_freemux.rds"))

########################## ANALYSIS ########################
# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=stim)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  theme_classic2()
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=SNG.BEST.GUESS)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(gen.all@meta.data, aes(x=sample, fill=SNG.BEST.GUESS)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_all_sample_singlet_genetic_freemux.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

pdf(paste0("bar_stim_cond_counts_only_singlet_genetic_freemux.pdf"), height = 8, width = 8)
stim.cond.counts.sing.array
dev.off()

# Normalization and identification of variable features
gen.all <- NormalizeData(gen.all, verbose = FALSE)
gen.all <- FindVariableFeatures(gen.all, verbose = FALSE)

# Run the standard workflow for visualization and clustering
gen.all <- ScaleData(gen.all, verbose = FALSE)
gen.all <- RunPCA(gen.all, npcs = 30, verbose = FALSE, reduction.name = "pca.unintegrated")
gen.all <- RunUMAP(gen.all, reduction = "pca.unintegrated", dims = 1:15, reduction.name = "umap.unintegrated")
gen.all <- RunTSNE(gen.all, reduction = "pca.unintegrated", dims = 1:15, check_duplicates = FALSE, reduction.name = "tsne.unintegrated")
gen.all <- FindNeighbors(gen.all, reduction = "pca.unintegrated", dims = 1:15)
gen.all <- FindClusters(gen.all, resolution = 0.6)

Idents(gen.all) <- "hash.ID"
gen.id <- subset(gen.all, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
gen.id <- JoinLayers(gen.id)
gen.id
# Normalization and identification of variable features
gen.id <- NormalizeData(gen.id, verbose = FALSE)
gen.id <- FindVariableFeatures(gen.id, verbose = FALSE)

# Run the standard workflow for visualization and clustering
gen.id <- ScaleData(gen.id, verbose = FALSE)
gen.id <- RunPCA(gen.id, npcs = 30, verbose = FALSE, reduction.name = "pca.unintegrated")
gen.id <- RunUMAP(gen.id, reduction = "pca.unintegrated", dims = 1:15, reduction.name = "umap.unintegrated")
gen.id <- RunTSNE(gen.id, reduction = "pca.unintegrated", dims = 1:15, check_duplicates = FALSE, reduction.name = "tsne.unintegrated")
gen.id <- FindNeighbors(gen.id, reduction = "pca.unintegrated", dims = 1:15)
gen.id <- FindClusters(gen.id, resolution = 0.6)

saveRDS(gen.id, file = paste0("gen.sing.rds"))

# Visualization
p1 <- DimPlot(gen.id, reduction = "umap.unintegrated", group.by = "id")
p2 <- DimPlot(gen.id, reduction = "umap.unintegrated", group.by = "stim")
p3 <- DimPlot(gen.id, reduction = "umap.unintegrated", group.by = "exp")
p4 <- DimPlot(gen.id, reduction = "umap.unintegrated", label = TRUE, repel = TRUE)

p1 

p2 

p3

p4

pdf("umap_unintegrated_mapped.pdf", width = 10)
ggarrange(p1, p2, p3, p4, nrow =2, ncol=2)
dev.off()
