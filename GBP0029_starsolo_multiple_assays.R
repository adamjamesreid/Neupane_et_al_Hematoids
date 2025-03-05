# GBP0029_starsolo_multiple_assays.R

# AUTHOR: Adam Reid
# Copyright (C) 2023 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# Re-process all the data, incorporate spliced, unspliced 'Velocyto' and 'Gene' assays from STARsolo
# Do not integrate the data, but filter each dataset in the same way.

library(Seurat)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)

setwd("/mnt/home3/reid/ajr236/projects/surani/GBP0029/Neupane_peEBs_code/")
set.seed(123)

## FUNCTIONS
# Read in STARsolo data
# N.b. this version uses the 'Gene' assay
read_solo <- function (s, dir = "../starsolo/data/") {
  # Data from STARSolo Gene assay - for main analysis
  gene_matrix <- ReadMtx(
    mtx = paste(dir, s, "/Gene/filtered/matrix.mtx", sep=""), 
    features = paste("features.edit.underscores.tsv", sep=""),
    cells = paste(dir, s, "/Velocyto/filtered/barcodes.tsv", sep="")
  )
  # Data from STARSolo Velocyto assay - only for velocity analysis
  spliced_matrix <- ReadMtx(
    mtx = paste(dir, s, "/Velocyto/filtered/spliced.mtx", sep=""), 
    features = paste("features.edit.underscores.tsv", sep=""),
    cells = paste(dir, s, "/Velocyto/filtered/barcodes.tsv", sep="")
  )
  unspliced_matrix <- ReadMtx(
    mtx = paste(dir, s, "/Velocyto/filtered/unspliced.mtx", sep=""), 
    features = paste("features.edit.underscores.tsv", sep=""),
    cells = paste(dir, s, "/Velocyto/filtered/barcodes.tsv", sep="")
  )
  seurat_object <- CreateSeuratObject(counts = gene_matrix, assay = "RNA", project=s)
  seurat_object[['spliced']] <- CreateAssayObject(counts = spliced_matrix)
  seurat_object[['unspliced']] <- CreateAssayObject(counts = unspliced_matrix)
  
  DefaultAssay(seurat_object) <- "RNA"
  
  # Rename cells to avoid merging problems due to non-unique barcodes
  name = regmatches(x,regexpr("[C-H][0-9]+",x))
  seurat_object <- RenameCells(object = seurat_object, add.cell.id = name)
  
  return(seurat_object)
}

# Convert gene id to gene name
n2e <- function(obj, ids) {
  name_list = c()
  for (x in ids){
    name_list = append(name_list, row.names(subset(obj@assays$RNA[[]], 
                                                   obj@assays$RNA[[]]$gene_name == x)))
  } 
  return(name_list)
}

## Colour settings
custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)

################
# Read in data and filter out cells (previously defined based on doublet analysis, low counts and high MT)
so_list = list()

# Get list of the cells we want to keep (identified in the cellQC.R script)
filtered_meta_data <- read.table("heo_filt_meta.tsv", header=TRUE, row.names=1)
cells_to_keep <- rownames(filtered_meta_data)

# Loop through each sample and create a Seurat object with feature metdata
# added. N.b. the '_' characters which get replaced with '-' are for genes in
# the Y chromosome PAR region (this now is an error and so has been fixed
# in the adjusted features file - features.edit.underscores.tsv)
for (x in c("run_SITTE1_star_Solo.out_copy", "run_SITTF1_star_Solo.out_copy", 
            "run_SITTG1_star_Solo.out_copy", "run_SITTH10_star_Solo.out_copy")) {
    s <- read_solo(x)
  
  # Filter out cells we don't want
  s <- s[,cells_to_keep]
  
  # Append each object to a list
  so_list <- append(so_list, s)
}

# Define list of normalised Seurat objects
sct_list <- c()

# Normalise Seurat objects and run PCA
for (x in so_list) {
  sct_list <- append(sct_list, SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
                       RunPCA(npcs = 50, verbose = TRUE))
}

#############
# Add Meta Data
##############
# Get gene names and sample meta data
feature_meta <- read.table("features.edit.underscores.tsv", header=FALSE, row.names=2, sep="\t")
feature_meta$V3 <- c()
sample_meta <- read.csv("sample.meta.csv", header=TRUE, row.names=1)

# Loop over Seurat objects 
for (i in 1:length(sct_list)) {
  
  DefaultAssay(sct_list[[i]]) <- 'RNA'
  
  # Add gene names
  sct_list[[i]][["RNA"]] <- AddMetaData(sct_list[[i]][["RNA"]], feature_meta, col.name = "gene_name")
  
  # Add cell meta data
  sample_name <- c()
  timepoint <- c()
  treatment <- c()
  
  for (j in 1:length(sct_list[[i]]@meta.data$orig.ident)) {
    
    sample_name[j] <- sample_meta[as.character(sct_list[[i]]@meta.data$orig.ident[j]), 'sample']
    timepoint[j] <- sample_meta[as.character(sct_list[[i]]@meta.data$orig.ident[j]), 'timepoint']
    treatment[j] <- sample_meta[as.character(sct_list[[i]]@meta.data$orig.ident[j]), 'treatment']
  }
  
  sct_list[[i]][['sample_name']] <- sample_name
  sct_list[[i]][['timepoint']] <- timepoint
  sct_list[[i]][['treatment']] <- treatment
  sct_list[[i]][['orig.ident']] <- sct_list[[i]][['sample_name']]
  
  # mitochondrial transcripts
  mt_features = rownames(subset(sct_list[[i]]@assays$RNA[[]], subset=grepl('^MT-', sct_list[[i]]@assays$RNA[[]]$gene_name)))
  sct_list[[i]][["percent.mt"]] <- PercentageFeatureSet(sct_list[[i]], features = mt_features)
                                                          
  # ribosomal gene transcripts
  sct_list[[i]][["percent.ribo"]] <- PercentageFeatureSet(sct_list[[i]], features = 
                                                  row.names(sct_list[[i]]@assays$RNA[[]])[grep('^RP[SL]', 
                                                                                               sct_list[[i]]@assays$RNA[[]]$gene_name)])
  # hemoglobin transcripts
  sct_list[[i]][["percent.hem"]] <- PercentageFeatureSet(sct_list[[i]], features = 
                                                 row.names(sct_list[[i]]@assays$RNA[[]])[grep("^HB[^(P)]", 
                                                                                              sct_list[[i]]@assays$RNA[[]]$gene_name)])
}

################
# Run dimension reduction and clustering
################

for (i in 1:length(sct_list)) { 
  DefaultAssay(sct_list[[i]]) <- 'SCT'
  
  # UMAPs - different n_neighbor parameters give similar results
  sct_list[[i]] <- RunUMAP(sct_list[[i]], dims = 1:30, n.neighbors = 30, reduction.name = "umap30dim30")
  #DimPlot(sct_list[[i]], label = TRUE, reduction = "umap30dim30", cols=custom_colors$discrete)
  
  # Clustering - different resolutions give similar results 
  sct_list[[i]] <- FindNeighbors(sct_list[[i]], reduction = 'pca', dims = 1:30, k.param = 20, 
                       verbose = FALSE, graph.name=c("NN20", "SNN20"))
  
  # get clusters
  sct_list[[i]] <- FindClusters(sct_list[[i]], verbose = TRUE, graph.name="NN20", algorithm=2,
                      resolution = 0.8)
}

# Save objects as individual rds files
for (i in 1:length(sct_list)) { 
  sample_name = sct_list[[i]]@meta.data$sample_name[1]
  file_name = paste('heo_', sample_name, ".rds", sep="")
  print(file_name)
  saveRDS(sct_list[[i]], file=file_name)
}

