library(cowplot)
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
library(DoubletFinder)
library(readr)
# library(caret) - load later to ensure no issues with running Seurat

metadata <- read.csv("metadata3.csv", stringsAsFactors = FALSE)

## We want to see which cells we can save with genetic demultiplexing 
# Initialize a list to store singlet objects
combined_objects <- list()

# Loop over each sample
for (i in 1:nrow(metadata)) {
  # Extract the current sample name
  sample <- metadata$sample[i]
  umis <- metadata$umis[i]
  
  # Construct the filename for the Seurat object
  filename <- paste0(umis, ".hashtag.rds")
  
  # Load the Seurat object from the RDS file
  seurat_object <- readRDS(filename)
  
  # Keep only singlet cells
  sample_hashtag.nonsing <- subset(seurat_object@meta.data, hash.ID %in% c('Negative', 'unmapped'))
  
  #Load gen_freemux data 
  filename2 <- paste0(sample, "_freemux_metadata_sing.csv")
  gen_freemux_obj <- read.csv(filename2, header = T, row.names = 1)
  
  #Determine barcodes that exist in both files - these are the instances in which we see singlet genetic demultiplexing calling where it does not exist for HTO 
  joint_bcs <- intersect(rownames(sample_hashtag.nonsing), rownames(gen_freemux_obj))
  
  # Filter Seurat object to keep only common barcodes
  if (length(joint_bcs) > 0) {
    
    length <- length(joint_bcs)
    print(paste("There are", length, "cells for", sample, "that can be saved!"))
    print(paste("Starting saviour pipeline"))
    
    #Load gen_freemux data 
    filename3 <- paste0(sample, "_freemux_metadata_all.csv")
    full_gen_freemux_obj <- read.csv(filename3, header = T, row.names = 1)
    
    # Ensure matching cell barcodes between RNA and adt assays
    keep_bcs <- intersect(colnames(seurat_object), rownames(full_gen_freemux_obj))
    seurat_object <- subset(seurat_object, cells = keep_bcs)
    new_metadata <- full_gen_freemux_obj[keep_bcs, ]
    new_metadata <- cbind(new_metadata, hash.ID = seurat_object@meta.data$hash.ID, HTO_maxID = seurat_object@meta.data$HTO_maxID)
    
    # Function to replace Patient1, Patient2, and Patient3 with Pt1, Pt2, and Pt3 respectively
    replace_patient_terms <- function(column) {
      column <- gsub("Pt1", "Patient1", column)
      column <- gsub("Pt2", "Patient2", column)
      column <- gsub("Pt3", "Patient3", column)
      column <- gsub("DBL", "Doublet", column)
      column <- gsub("AMB", "Ambiguous", column)
      return(column)
    }
    
    # Apply the function to every column in new_metadata
    new_metadata[] <- lapply(new_metadata, replace_patient_terms)
    
    # First, prepare the new column based on conditions
    new_metadata$cat_call <- ifelse(
      (new_metadata$hash.ID %in% c('Negative', 'unmapped') & new_metadata$DROPLET.TYPE == 'Ambiguous' | new_metadata$hash.ID == 'Doublet' | new_metadata$DROPLET.TYPE == 'Doublet'),
      'trash',
      paste0(new_metadata$hash.ID, "_", new_metadata$DROPLET.TYPE)
    )
    
    # Create a new column based on conditions
    new_metadata$final_call <- ifelse(
      new_metadata$hash.ID %in% c('Negative', 'unmapped'),
      new_metadata$SNG.BEST.GUESS,
      new_metadata$hash.ID
    )
    
    new_metadata$cat_pt <- paste0(new_metadata$HTO_maxID, "_", new_metadata$SNG.BEST.GUESS)
    new_metadata$final_cat <- paste0(new_metadata$final_call, "_", new_metadata$cat_call)
    ## Add new_metadata to seurat_object
    seurat_object <- AddMetaData(seurat_object, metadata = new_metadata)
    assign(paste0(sample, ".hashtag"), seurat_object, envir = .GlobalEnv)
    # Optionally save the singlet object
    saveRDS(seurat_object, file = paste0(sample, ".combined.hashtag.rds"))
    
    print(paste("Cells have been saved for", sample))
    
  } else {
    warning(paste("No common barcodes found for saving in", sample))
  }
  
  # Add the processed singlet object to the list
  combined_objects[[sample]] <- seurat_object
  
}

# Doublet Removal
# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  
  sample <- metadata$sample[i]
  umis <- metadata$umis[i]
  
  # Construct the filename for the Seurat object
  filename <- paste0(sample, ".combined.hashtag.rds")
  
  # Load the Seurat object from the RDS file
  seurat_object <- readRDS(filename)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  seurat_object <- RunUMAP(seurat_object, dims = 1:10)
  
  # Define the conversion mapping
  conversion_mapping <- c(
    "SNG" = "Singlet",
    "Doublet" = "Doublet",
    "Ambiguous" = "Doublet")
  
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
############# CREATE MERGED OBJ

# Initialize a list to store the processed Seurat objects
seurat_objects_list <- list()
# Initialize a vector to store the id_factors for each object
id_factors <- c()

# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  
  # Construct the file name, load the Seurat object, process it
  file_name <- paste0(sample_name, ".combined.hashtag.rds")
  seurat_object <- readRDS(file = file_name)
  
  DefaultAssay(seurat_object) <- "RNA"
  Idents(seurat_object) <- "orig.ident"
  seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  # optional to subset at this point 
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
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

# Save the merged Seurat object
saveRDS(merged_seurat_object, file = "merged_seurat_object_all.rds")

hto.all <- readRDS(file = paste0("merged_seurat_object_all.rds"))

#### Comparison bar graphs
# Counts of each stim condition for all singlets by patient
stim.cond.counts.merged.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.merged.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.merged.array.final <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=final_call)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  labs(fill = "Final Assignment") +
  theme_classic2(base_size = 22)
print(stim.cond.counts.merged.array.final)

stim.cond.counts.merged.array.hto <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  labs(fill = "Antibody Assignment") +
  theme_classic2(base_size = 22)
print(stim.cond.counts.merged.array.hto)

stim.cond.counts.merged.array.genetic <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=SNG.BEST.GUESS)) + 
  scale_y_continuous(breaks=seq(0,9000,2000)) +
  geom_bar() +
  labs(fill = "Genetic Assignment") +
  theme_classic2(base_size = 22)
print(stim.cond.counts.merged.array.genetic)

pdf(paste0("bar_stim_cond_all_merged_compare.pdf"), height = 8, width = 25)
ggarrange(stim.cond.counts.merged.array.hto, stim.cond.counts.merged.array.genetic, 
          stim.cond.counts.merged.array.final, common.legend = FALSE, nrow = 1)
dev.off()

pdf(paste0("bar_stim_cond_all_merged_compare.pdf"), height = 8, width = 25)
stim.cond.counts.merged.array.hto | stim.cond.counts.merged.array.genetic | stim.cond.counts.merged.array.final
dev.off()

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.merged.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.merged.array)

pdf(paste0("bar_stim_cond_all_merged_sample.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.merged.array, stim.cond.proportion.merged.array, common.legend = TRUE)
dev.off()

pdf(paste0("bar_stim_cond_count_merged_sample.pdf"), height = 8, width = 8)
stim.cond.counts.merged.array
dev.off()
########## Subset 
# Loop over each sample
for (i in 1:nrow(metadata)) {
  
  # Construct the filename for the Seurat object
  sample_name <- metadata$sample[i]
  filename <- paste0(sample_name, ".combined.hashtag.rds")
  
  # Load the Seurat object from the RDS file
  sample_hashtag <- readRDS(filename)
  
  Idents(sample_hashtag) <- "cat_call"
  # Keep only singlet cells
  sample_hashtag.subset <- subset(sample_hashtag, idents = "trash", invert = T)
  
  # Save the singlet object in the local environment with a dynamic name
  assign(paste0(sample_name, ".singlet"), sample_hashtag.subset)
  
  # Optionally save the singlet object
  saveRDS(sample_hashtag.subset, file = paste0(sample_name, ".hashtag.singlet.rds"))
}

# Initialize a list to store singlet objects
HTO_singlet_objects <- list()

# Loop over each sample
for (i in 1:nrow(metadata)) {
  
  # Construct the filename for the Seurat object
  sample_name <- metadata$sample[i]
  filename <- paste0(sample_name, ".hashtag.singlet.rds")
  
  #Read in Seurat object
  sample_hashtag <- readRDS(filename)
  
  # Perform the analysis
  Idents(sample_hashtag) <- "HTO_maxID"
  print(RidgePlot(sample_hashtag, assay = "HTO", features = rownames(sample_hashtag[["HTO"]]), ncol = 1))
  
  # FeatureScatter plots
  DefaultAssay(sample_hashtag) <- "HTO"
  
  Idents(sample_hashtag) <- "HTO_classification.global"
  print(VlnPlot(sample_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE))
  
  # First, we will remove negative cells from the object
  sample_hashtag.subset <- subset(sample_hashtag, idents = "Negative", invert = TRUE)
  
  # Calculate a tSNE embedding of the HTO data
  DefaultAssay(sample_hashtag.subset) <- "HTO"
  sample_hashtag.subset <- ScaleData(sample_hashtag.subset, features = rownames(sample_hashtag.subset),
                                     verbose = FALSE)
  
  sample_hashtag.subset <- RunPCA(sample_hashtag.subset, features = rownames(sample_hashtag.subset), approx = FALSE)
  
  
  sample_hashtag.subset <- RunTSNE(sample_hashtag.subset, dims = 1:4, perplexity = 10, check_duplicates = FALSE)
  print(DimPlot(sample_hashtag.subset))
  
  # To increase the efficiency of plotting, you can subsample cells using the num.cells argument
  print(HTOHeatmap(sample_hashtag, assay = "HTO", ncells = 5000))
  
  # Calculate a tSNE embedding of the HTO data
  DefaultAssay(sample_hashtag) <- "HTO"
  sample_hashtag.singlet <- subset(sample_hashtag, idents = "Singlet")
  
  # Save the singlet object in the local environment with a dynamic name
  # assign(paste0(sample_name, ".singlet"), sample_hashtag.singlet)
  
  # Add the processed singlet object to the list
  HTO_singlet_objects[[sample_name]] <- sample_hashtag.singlet
  
  # Optionally save the singlet object
  saveRDS(sample_hashtag.singlet, file = paste0(sample_name, ".HTO.hashtag.singlet.rds"))
  # Save the modified Seurat object
  #saveRDS(sample_hashtag.subset, file = paste0(sample_name, "_processed.hashtag.rds"))
}

# Initialize a list to store the processed Seurat objects
seurat_objects_list <- list()
# Initialize a vector to store the id_factors for each object
id_factors <- c()

# Loop over each sample to process
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$sample[i]
  
  # Construct the file name, load the Seurat object, process it
  file_name <- paste0(sample_name, ".hashtag.singlet.rds")
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
saveRDS(merged_seurat_object, file = "merged_seurat_object.rds")

hto.all <- readRDS(file = paste0("merged_seurat_object.rds"))

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_all_sample.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

#### Same as above, but looking when patient is called based on hash.ID, then SNG.BEST.GUESS
# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=final_call)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=final_call)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_final_call.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.cat.call <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=cat_call)) + 
  geom_bar() 
print(stim.cond.counts.sing.cat.call)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.final.cat <- ggplot(hto.all@meta.data, aes(x=sample_id, fill=final_cat)) + 
  geom_bar() 
print(stim.cond.counts.sing.final.cat)


pdf(paste0("bar_stim_cond_final_call.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.cat.call, stim.cond.counts.sing.final.cat, common.legend = FALSE)
dev.off()
# Normalization and identification of variable features
hto.all <- NormalizeData(hto.all, verbose = FALSE)
hto.all <- FindVariableFeatures(hto.all, verbose = FALSE)

# Remove low quality cells and join layers 
Idents(hto.all) <- "final_call"
# Visualize QC metrics as a violin plot
VlnPlot(hto.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

hto.id <- subset(hto.all, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hto.id <- JoinLayers(hto.id)
hto.id

# Visualize QC metrics as a violin plot
VlnPlot(hto.id, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalization and identification of variable features
hto.id <- NormalizeData(hto.id, verbose = FALSE)
hto.id <- FindVariableFeatures(hto.id, verbose = FALSE)

# Run the standard workflow for visualization and clustering
hto.id <- ScaleData(hto.id, verbose = FALSE)
hto.id <- RunPCA(hto.id, npcs = 30, verbose = FALSE, reduction.name = "pca", seed.use = 1997)
hto.id <- RunUMAP(hto.id, reduction = "pca", dims = 1:15, reduction.name = "umap", seed.use = 1997)
hto.id <- RunTSNE(hto.id, reduction = "pca", dims = 1:15, check_duplicates = FALSE, reduction.name = "tsne", seed.use = 1997)
hto.id <- FindNeighbors(hto.id, reduction = "pca", dims = 1:15)
hto.id <- FindClusters(hto.id, resolution = 0.6)

saveRDS(hto.id, file = paste0("hto.sing.rds"))
hto.id <- readRDS("hto.sing.rds")
# Visualization
p1 <- DimPlot(hto.id, reduction = "umap", group.by = "hash.ID", 
              cols = c('#F1BB7B','#FD6467','#7294D4','#E6A0C4','#5B1A18'))
p2 <- DimPlot(hto.id, reduction = "umap", group.by = "SNG.BEST.GUESS",
              cols = c('#808080','#FD6467','#7294D4','#E6A0C4'))
p3 <- DimPlot(hto.id, reduction = "umap", group.by = "final_call",
              cols = c('#FD6467','#7294D4','#E6A0C4'))
p4 <- DimPlot(hto.id, reduction = "umap", label = TRUE, repel = TRUE)

p1 

p2 

p3

p4

pdf("umap_bymapping.pdf", width = 15, height = 6)
ggarrange(p1, p2, p3, nrow =1, ncol=3, common.legend = T, align = "v", widths = c(1,1,1),
          label.x = NULL, label.y = NULL)
dev.off()

p5 <- DimPlot(hto.id, reduction = "umap", group.by = "hash.ID", 
              cols = c('#F1BB7B','#FD6467','#7294D4','#E6A0C4','#5B1A18')) & NoLegend()
p6 <- DimPlot(hto.id, reduction = "umap", group.by = "SNG.BEST.GUESS",
              cols = c('#808080','#FD6467','#7294D4','#E6A0C4')) & NoLegend()
p7 <- DimPlot(hto.id, reduction = "umap", group.by = "final_call",
              cols = c('#FD6467','#7294D4','#E6A0C4')) & NoLegend()
p8 <- DimPlot(hto.id, reduction = "umap", label = TRUE, repel = TRUE)

p5 

p6 

p7

p8
pdf("umap_bymapping_nolabel.pdf", width = 15, height = 6)
ggarrange(p5, p6, p7, nrow =1, ncol=3, align = "v", widths = c(1,1,1),
          label.x = NULL, label.y = NULL)
dev.off()


# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_all_sample_mapped.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

# Visualization
p1 <- DimPlot(hto.id, reduction = "umap", group.by = "sample_id")
p2 <- DimPlot(hto.id, reduction = "umap", group.by = "stim")
p3 <- DimPlot(hto.id, reduction = "umap", group.by = "exp")
p4 <- DimPlot(hto.id, reduction = "umap", label = TRUE, repel = TRUE)

p1 

p2 

p3

p4

pdf("umap_unintegrated_mapped.pdf", width = 10)
ggarrange(p1, p2, p3, p4, nrow =2, ncol=2)
dev.off()
############### DOWNSTREAM ANALYSIS
# Perform the analysis
library(ggplot2)
library(cowplot) # for plot_grid to combine plots

# Ensure that the identity variable is a factor and its levels are ordered alphabetically
hto.id@meta.data$HTO_maxID <- factor(hto.id@meta.data$HTO_maxID, levels = rev(sort(unique(hto.id@meta.data$HTO_maxID))))

Idents(hto.id) <- "HTO_maxID"
ridge_list <- RidgePlot(hto.id, assay = "HTO", features = rownames(hto.id[["HTO"]]), 
                        ncol = 1, same.y.lims = TRUE, combine = FALSE)


# Apply the x-axis limits to each plot in the list
ridge_list <- lapply(ridge_list, function(p) {
  p + coord_cartesian(xlim = c(0, 750), expand = FALSE)
})

# Combine the plots into one plot
ridge_combined <- plot_grid(plotlist = ridge_list, align = 'v', ncol = 1)

# Print the combined plot
ridge_combined

pdf("hto_ridge_plot.pdf", height = 10, width = 8)
ridge_combined
dev.off()

# Assuming singlet_objects is a list of Seurat objects
# Merge all singlet objects into one Seurat object
# The names for add.cell.ids should be the same length as the singlet_objects list

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=stim)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Counts of each stim condition for all singlets by patient
stim.cond.counts.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar() 
print(stim.cond.counts.sing.array)

# Proportion of each stim condition for all singlets by patient
stim.cond.proportion.sing.array <- ggplot(hto.id@meta.data, aes(x=sample_id, fill=hash.ID)) + 
  geom_bar(position = "fill") 
print(stim.cond.proportion.sing.array)

pdf(paste0("bar_stim_cond_all_sample.pdf"), height = 15, width = 15)
ggarrange(stim.cond.counts.sing.array, stim.cond.proportion.sing.array, common.legend = TRUE)
dev.off()

# Normalization and identification of variable features
Idents(hto.id) <- "final_call"
#hto.id <- subset(hto.id, idents = 'unmapped', invert = TRUE)
hto.id <- subset(hto.id, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
hto.id <- JoinLayers(hto.id)

DefaultAssay(hto.id) <- "RNA"
hto.id <- NormalizeData(hto.id, verbose = FALSE)
hto.id <- FindVariableFeatures(hto.id, verbose = FALSE)

# Run the standard workflow for visualization and clustering
hto.id <- ScaleData(hto.id, verbose = FALSE)
hto.id <- RunPCA(hto.id, npcs = 30, verbose = FALSE, seed.use = 1997)
hto.id <- RunUMAP(hto.id, reduction = "pca", dims = 1:15, seed.use = 1997)
hto.id <- RunTSNE(hto.id, reduction = "pca", dims = 1:15, check_duplicates = FALSE, seed.use = 1997)
hto.id <- FindNeighbors(hto.id, reduction = "pca", dims = 1:15)
hto.id <- FindClusters(hto.id, resolution = 1.2)

# Visualization
p1 <- DimPlot(hto.id, reduction = "umap", group.by = "sample_id")
p2 <- DimPlot(hto.id, reduction = "umap", group.by = "final_call")
p3 <- DimPlot(hto.id, reduction = "umap", group.by = "exp")
p4 <- DimPlot(hto.id, reduction = "umap", label = TRUE, repel = TRUE)

umap_hash <- DimPlot(hto.id, reduction = "umap", group.by = "hash.ID", repel = TRUE)
umap_genetic <- DimPlot(hto.id, reduction = "umap", group.by = "SNG.BEST.GUESS", repel = TRUE)

pdf('umap_compare.pdf', height = 6, width = 12)
umap_hash | umap_genetic
dev.off()

p1 

p2 

p3

p4

pdf("umap_bysample.pdf", width = 10)
ggarrange(p1, p2, p3, p4, nrow =2, ncol=2)
dev.off()

pdf("umap_1.2_RNA_all.pdf", width = 10)
p2 | p4
dev.off()

ggarrange(p2, p4, nrow = 1)

heatmap.RNA <- DoHeatmap(hto.id, features = c("CD86", "FCGR3A", "CD14", "CD1C", "CD19"), 
                         size = 3) + scale_fill_gradientn(colors = c("lightgrey", "blue")) + guides(colour=FALSE)
heatmap.RNA

################ FIND ALL MARKERS #####################
hto.id.markers <- FindAllMarkers(hto.id, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

hto.id.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(hto.id.markers,file="combined_markers_snn1.2.csv", quote=FALSE)

library(Azimuth)
hto.id <- RunAzimuth(hto.id, reference = "pbmcref")

p1 <- DimPlot(hto.id, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3)
p1 | p4

pdf('umap_predicted.celltype.l2.pdf', height = 6, width = 16)
p2 | p1
dev.off()

p5 <- DimPlot(hto.id, reduction = "umap", group.by = "final_call")
p6 <- DimPlot(hto.id, reduction = "umap", group.by = "predicted.celltype.l2")

pdf('umap_anno_hashID.pdf', height = 8, width = 15)
p5 | p6
dev.off()

#### MULTIMODAL ANALYSIS
DefaultAssay(hto.id) <- "RNA"
hto.id <- SCTransform(hto.id, vars.to.regress = "id_factor") %>%  RunPCA(reduction.name = 'pca', seed.use = 1997)

DefaultAssay(hto.id) <- "ADT"
VariableFeatures(hto.id) <- rownames(hto.id[["ADT"]])
hto.id <- SCTransform(hto.id, assay = "ADT", vars.to.regress = "id_factor", new.assay.name = "SCT_ADT") %>% RunPCA(reduction.name = 'apca', seed.use = 1997)

DefaultAssay(hto.id) <- "RNA"
hto.id <- FindMultiModalNeighbors(hto.id, reduction.list = list("pca", "apca"), dims.list = list(1:30, 1:20), modality.weight.name = "SCT.weight")
hto.id <- FindClusters(hto.id, graph.name = "wsnn", algorithm = 3, resolution = 1.6, verbose = FALSE)
hto.id <- RunUMAP(hto.id, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 1997)
hto.id.markers <- FindAllMarkers(hto.id, only.pos = TRUE)
hto.id.markers <- hto.id.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)

write_csv(hto.id.markers, "wnn_markers_seed1997_RNA.csv")
p7 <- DimPlot(hto.id, reduction = "wnn.umap", group.by = "final_call")
p8 <- DimPlot(hto.id, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5, shuffle = T)
p8
p7 | p8

p9 <-  DimPlot(hto.id, reduction = "wnn.umap", group.by = "predicted.celltype.l2") 
p9
p10 <-  DimPlot(hto.id, reduction = "wnn.umap", group.by = "id_factor")
p10

pdf('cluster_wnn_toanno_1.2_1997.pdf', height = 6, width = 10)
ggarrange(p8, p9, nrow =2, ncol=1)
dev.off()

p9 <- DimPlot(hto.id, reduction = "wnn.umap", group.by = "final_call") + NoLegend()
p10 <- DimPlot(hto.id, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5, shuffle = T) + NoLegend()

p9 | p10

pdf('multi_cluster_wnn_umap_NoLeg_1.pdf', height = 5, width = 5)
p9
dev.off()

pdf('multi_cluster_wnn_umap_NoLeg_2.pdf', height = 5, width = 5)
p10
dev.off()

pdf('multi_cluster_wnn_umap_NoLeg.pdf', height = 6, width = 16)
p9 | p10
dev.off()

new.cluster.ids <- c("CD4 T-Cell", "CD14 Monocyte", "CD4 T-Cell", "NK Cell", "CD14 Monocyte", "Treg", "CD14 Monocyte", 
                     "Treg", "CD14 Monocyte", "CD8 T-Cell", "CD4 T-Cell", "B-Cell", "NK Cell", "CD4 T-Cell", "MAIT Cell", "CD4 T-Cell", 
                     "gdT-Cell", "CD8 T-Cell", "N/A", "CD8 T-Cell", "Dendritic Cell", "CD16 Monocyte", "CD8 T-Cell", "Dendritic Cell", 
                     "Neutrophil", "NK Cell", "N/A", "NK Cell", "Plasmacytoid Dendritic Cell", "NK Cell", "N/A", "N/A", 
                     "Granulocytes")
names(new.cluster.ids) <- levels(hto.id)
hto.id <- RenameIdents(hto.id, new.cluster.ids)
hto.id <- subset(hto.id, idents = "N/A", invert = T)

p19 <- DimPlot(hto.id, reduction = "wnn.umap", group.by = "final_call") 
p20 <- DimPlot(hto.id, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5, shuffle = T) 
p21 <- DimPlot(hto.id, reduction = 'wnn.umap', label = FALSE, shuffle = T) 
p19 | p20
p21
pdf('umap_anno_seed1997.pdf', height = 4, width = 12)
ggarrange(p19, p21, nrow =1, ncol=2)
dev.off()

pdf('umap_anno2_seed1997.pdf', height = 4, width = 12)
p19 | p21
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a1 <- FeaturePlot(hto.id, "CD4", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD4 protein")
DefaultAssay(hto.id) <- "RNA"
a2 <- FeaturePlot(hto.id, "CD4", order = T, reduction = 'wnn.umap') + ggtitle("CD4 RNA")

# place plots side-by-side
a1 | a2

pdf('cd4_featureplots_seed1997.pdf', height = 6, width = 10)
a1 | a2
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a21 <- FeaturePlot(hto.id, "CD8A", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD8 protein")
a21
DefaultAssay(hto.id) <- "RNA"
a22 <- FeaturePlot(hto.id, "CD8A", order = T, reduction = 'wnn.umap') + ggtitle("CD8 RNA")

# place plots side-by-side
a21 | a22

pdf('cd4_cd8_featureplots_seed1997_300.pdf', height = 6, width = 10)
ggarrange(a1, a2, a21, a22, nrow =2, ncol=2, common.legend = T)
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a23 <- FeaturePlot(hto.id, c("CD3D"), cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap', max.cutoff = 300) + ggtitle("CD3 protein")
DefaultAssay(hto.id) <- "RNA"
a24 <- FeaturePlot(hto.id, "CD3D", order = T, reduction = 'wnn.umap', max.cutoff = 4) + ggtitle("CD3 RNA")

# place plots side-by-side
a23 | a24

pdf('tcell_featureplots_seed1997.pdf', height = 6, width = 10)
ggarrange(a1, a2, a21, a22, a23, a24, nrow =3, ncol=2, common.legend = T)
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a25 <- FeaturePlot(hto.id, "CD69", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD69 protein")
DefaultAssay(hto.id) <- "RNA"
a26 <- FeaturePlot(hto.id, "CD69", order = T, reduction = 'wnn.umap') + ggtitle("CD69 RNA")

# place plots side-by-side
a25 | a26

pdf('cd69_featureplots_seed1997.pdf', height = 6, width = 10)
a25 | a26
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a27 <- FeaturePlot(hto.id, "FCGR3A", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD16 protein")
DefaultAssay(hto.id) <- "RNA"
a28 <- FeaturePlot(hto.id, "FCGR3A", order = T, reduction = 'wnn.umap') + ggtitle("CD16 RNA")

# place plots side-by-side
a27 | a28

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a29 <- FeaturePlot(hto.id, "CD14", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD14 protein")
DefaultAssay(hto.id) <- "RNA"
a30 <- FeaturePlot(hto.id, "CD14", order = T, reduction = 'wnn.umap') + ggtitle("CD14 RNA")

# place plots side-by-side
a29 | a30

pdf('cd16_cd14_featureplots_seed1997.pdf', height = 6, width = 10)
ggarrange(a27, a28, a29, a30, nrow =2, ncol=2)
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a31 <- FeaturePlot(hto.id, "CD1C", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD1C protein")
DefaultAssay(hto.id) <- "RNA"
a32 <- FeaturePlot(hto.id, "CD1C", order = T, reduction = 'wnn.umap') + ggtitle("CD1C RNA")

# place plots side-by-side
a31 | a32

pdf('cd14_cd1c_dendritic_featureplots_seed1997.pdf', height = 6, width = 10)
ggarrange(a29, a30, a31, a32, nrow =2, ncol=2, common.legend = T)
dev.off()

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a33 <- FeaturePlot(hto.id, "CD19", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD19 protein")
DefaultAssay(hto.id) <- "RNA"
a34 <- FeaturePlot(hto.id, "CD19", order = T, reduction = 'wnn.umap') + ggtitle("CD19 RNA")

# place plots side-by-side
a33 | a34

# visualize one or the other
DefaultAssay(hto.id) <- "ADT"
a35 <- FeaturePlot(hto.id, "CD86", cols = c("lightgrey", "darkgreen"), order = T, reduction = 'wnn.umap') + ggtitle("CD86 protein")
DefaultAssay(hto.id) <- "RNA"
a36 <- FeaturePlot(hto.id, "CD86", order = T, reduction = 'wnn.umap') + ggtitle("CD86 RNA")

# place plots side-by-side
a35 | a36

pdf('cd19_cd86_featureplots_seed1997.pdf', height = 6, width = 10)
ggarrange(a33, a34, a35, a36, nrow =2, ncol=2)
dev.off()

pdf('supp_featureplots_seed1997.pdf', height = 16, width = 8)
ggarrange(a35, a36, a27, a28, a29, a30, a31, a32, a33, a34, nrow =5, ncol=2)
dev.off()

DefaultAssay(hto.id) <- "RNA"
heatmap.RNA <- DoHeatmap(hto.id, features = c("CD86", "FCGR3A", "CD14", "CD1C", "CD19"), slot = "data", 
                         size = 3) + scale_fill_gradientn(colors = c("lightgrey", "blue")) + guides(colour=FALSE)
heatmap.RNA

heatmap.ADT <- DoHeatmap(hto.id, features = c("CD86", "FCGR3A", "CD14", "CD1C", "CD19"), size = 3, slot = "data",
                     assay = "ADT") + scale_fill_gradientn(colors = c("lightgrey", "darkgreen")) + guides(colour=FALSE)
heatmap.ADT

heatmap.all.rows <- ggarrange(heatmap.ADT, heatmap.RNA, nrow = 1)
heatmap.all.rows

heatmap.all.cols <- ggarrange(heatmap.ADT, heatmap.RNA, ncol = 1)
heatmap.all.cols

pdf('heatmap_markergenes_rows.pdf', height = 4, width = 16)
heatmap.all.rows
dev.off()

pdf('heatmap_markergenes_cols.pdf', height = 8, width = 8)
heatmap.all.cols
dev.off()
