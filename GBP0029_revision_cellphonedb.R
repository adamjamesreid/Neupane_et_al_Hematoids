# GBP0029_revision_cellphonedb.R

# AUTHOR: Adam Reid
# Copyright (C) 2025 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# This is the refined CellPhoneDB script for producing figures for the manuscript revision

###########
# LIBRARIES
############

library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

library(SeuratDisk)

################
# SETUP
################

set.seed(123)

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

##############
# FUNCTIONS FOR CONVERTING GENE NAMES
################

# Convert gene id to gene name
n2e <- function(sobj, ids) {
  name_list = c()
  for (x in ids){
    name_list = append(name_list, row.names(subset(sobj@assays$RNA[[]], 
                                                   sobj@assays$RNA[[]]$gene_name == x)))
  } 
  return(name_list)
}

# Convert gene name to gene id
e2n <- function(sobj, ids) {
  name_list = c()
  for (x in ids){
    name_list = append(name_list, sobj@assays$RNA[[]][x,])
  }
  return(name_list)
}

####################
# PREPARE DATA FOR INPUT TO CELLPHONEDB
#####################
heo_f1 <- readRDS("~/projects/surani/GBP0029/starsolo/130723/heo_F1_labelled.rds")
heo_g1 <- readRDS("~/projects/surani/GBP0029/starsolo/130723/heo_G1_labelled.rds")

f1u <- DimPlot(heo_f1, label=TRUE, cols=custom_colors$discrete)

f1u

# Plotting for presentation
DimPlot(heo_f1, label=TRUE, cols=custom_colors$discrete, repel=TRUE) +NoLegend()
ggsave('F1_umap_labelled.pdf', height = 5, width = 5)
FeaturePlot(subset(heo_f1, idents=c("Blood progenitors 1", "Blood progenitors 2", "Mesenchyme 1",
                                    "Mesenchyme 2", "Mesenchyme 3", "Intermediate mesoderm 1",
                                    "Intermediate mesoderm 2", "Intermediate mesoderm 3")), 
            features=n2e(heo_f1, c("APOE", "TREM2", "TYROBP") ))
ggsave('APOE_TREM2_umaps.pdf', height = 5, width = 5)

VlnPlot(subset(heo_f1, idents=c("Blood progenitors 1", "Blood progenitors 2", "Mesenchyme 1",
                                    "Mesenchyme 2", "Mesenchyme 3", "Intermediate mesoderm 1",
                                    "Intermediate mesoderm 2", "Intermediate mesoderm 3")), 
            features=n2e(heo_f1, c("APOE", "TREM2", "TYROBP") ))

FeaturePlot(heo_f1, features=n2e(heo_f1, c("APOE", "TREM2", "TYROBP") ))

VlnPlot(heo_f1, features=n2e(heo_f1, c("APOE", "TREM2", "TYROBP") ))

FeaturePlot(subset(heo_f1, idents=c("Blood progenitors 1", "Blood progenitors 2", "Mesenchyme 1",
                                    "Mesenchyme 2", "Mesenchyme 3", "Intermediate mesoderm 1",
                                    "Intermediate mesoderm 2", "Intermediate mesoderm 3")), 
            features=n2e(heo_f1, c("BMP2", "BMPR1A", "BMPR2") ))

FeaturePlot(subset(heo_f1, idents=c("Blood progenitors 1", "Blood progenitors 2", "Mesenchyme 1",
                                    "Mesenchyme 2", "Mesenchyme 3", "Intermediate mesoderm 1",
                                    "Intermediate mesoderm 2", "Intermediate mesoderm 3")), 
            features=n2e(heo_f1, c("BMP2", "BMR1A", "AVR2B") ))

## F1
# Write out heo_f1 h5ad for reading into Python
# Get rid of assays we don't want (we only want RNA)
DefaultAssay(heo_f1) <- 'RNA' # N.b. there are no normalised counts here - 'data' and 'counts' layers are the same
heo_f1@assays$SCT <- c()
heo_f1@assays$spliced <- c()
heo_f1@assays$unspliced <- c()
# Re-do normalisation without SCTransform
heo_f1 <- NormalizeData(heo_f1)
heo_f1 <- FindVariableFeatures(heo_f1, selection.method = "vst", nfeatures = 2000)
# If we scale, this data gets used as the normalised counts after the conversion to h5ad!!!

# Rename idents for simplified pairwise analysis (remove numbers from cell types)
old_names <- levels(heo_f1@active.ident)
old_names[13] <- 'HCs' #Blood progenitors 1
old_names[16] <- 'ECs' #Blood progenitors 2
# Merge clusters with same cell type (e.g. 'Erythroids 1' and 'Erythroids 2' -> 'Erythroids')
new_names <- gsub(' [0-9]+$', '', old_names, perl=TRUE)
new_names_list <- setNames( new_names, levels(heo_f1@active.ident))
heo_f1_rename <- RenameIdents(object = heo_f1, new_names_list)

# Check naming looks like it should
DimPlot(heo_f1_rename, label=TRUE, cols=custom_colors$discrete)

# Convert to H5AD so we can use the data in Scanpy
SaveH5Seurat(heo_f1_rename, filename = "heo_f1.h5Seurat", overwrite = TRUE)
Convert("heo_f1.h5Seurat", dest = "h5ad", overwrite = TRUE)

## G1
DefaultAssay(heo_g1) <- 'RNA' # N.b. there are no normalised counts here - 'data' and 'counts' layers are the same
heo_g1@assays$SCT <- c()
heo_g1@assays$spliced <- c()
heo_g1@assays$unspliced <- c()
# Re-do normalisation without SCTransform
heo_g1 <- NormalizeData(heo_g1)
heo_g1 <- FindVariableFeatures(heo_g1, selection.method = "vst", nfeatures = 2000)

# Rename idents for simplified pairwise analysis (remove numbers from cell types)
old_names <- levels(heo_g1@active.ident)
old_names[19] <- 'ECs'
new_names <- gsub(' [0-9]+$', '', old_names, perl=TRUE)
new_names_list <- setNames( new_names, levels(heo_g1@active.ident))
heo_g1_rename <- RenameIdents(object = heo_g1, new_names_list)

# Check naming makes sense
DimPlot(heo_g1_rename, label=TRUE, cols=custom_colors$discrete)

# Convert to H5AD for use with Scanpy
SaveH5Seurat(heo_g1_rename, filename = "heo_g1.h5Seurat", overwrite = TRUE)
Convert("heo_g1.h5Seurat", dest = "h5ad", overwrite = TRUE)


### Generate additional files for CPDB (cell type/cluster names for each cell barcode)
# Make metadata files describing cell-> cell type
write.table(data.frame(cell_type = heo_f1_rename@active.ident), file="heo_f1_meta.tsv", sep='\t', quote=FALSE)
write.table(data.frame(cell_type = heo_g1_rename@active.ident), file="heo_g1_meta.tsv", sep='\t', quote=FALSE)

# Make microenvironments file for F1
microenv_f1 <- data.frame(cell_type=c("Mesenchyme", "Intermediate mesoderm", "HCs", "ECs"))
microenv_f1$microenvironment <- rep("env", nrow(microenv_f1))
write.table(microenv_f1, file = "heo_f1_microenviroment.tsv", quote=FALSE, row.names=FALSE, sep='\t')

# Make microenvironments file for G1
microenv_g1 <- data.frame(cell_type=c("Mesenchyme","Mixed/cardiac mesoderm", 
                      "Amnion", "Mixed mesoderm", "Presomitic mesoderm","ECs"))
microenv_g1$microenvironment <- rep("env", nrow(microenv_g1))
write.table(microenv_g1, file = "heo_g1_microenviroment.tsv", quote=FALSE, row.names=FALSE, sep='\t')


###############
# The CPDB heatmap function
## Additional list is a vector of interaction pairs which should be highlighted
cpdb_heatmap <- function(goi_list, cpdb_sig_means_file, additional_list=c(), row_font_size = 10) {
  
  library(pheatmap)
  library(RColorBrewer)
  
  cpdb_sig_means <- read.table(cpdb_sig_means_file, header=TRUE, sep='\t')
  
  sig_means_niche_genes <- subset(cpdb_sig_means, cpdb_sig_means$interacting_pair %in% goi_list)
  sig_means_niche_genes_table <- sig_means_niche_genes[,c(2,15:length(sig_means_niche_genes))]
  rownames(sig_means_niche_genes_table) <- sig_means_niche_genes_table$interacting_pair
  sig_means_niche_genes_table$interacting_pair <- c()
  sig_means_niche_genes_table[is.na(sig_means_niche_genes_table)] <- 0
  
  # Remove rows with no significant values
  sig_means_niche_genes_table <- subset(sig_means_niche_genes_table, 
                                        !(rowSums(sig_means_niche_genes_table == 0) == dim(sig_means_niche_genes_table)[2]))
  
  rank <- sig_means_niche_genes %>% select(interacting_pair, rank)
  rownames(rank) <- rank$interacting_pair
  rank$interacting_pair <- c()
  newCols <- colorRampPalette(grDevices::rainbow(length(unique(rank$rank))))
  mycolors <-colorRampPalette(brewer.pal(7,"YlGn"))(length(unique(rank$rank)))
  mycolors = rev(mycolors)
  names(mycolors) <- unique(rank$rank)
  mycolors <- list(rank = mycolors)
  
  if (length(additional_list) > 0) {
    rank$present_in_g1 <- rownames(rank) %in% additional_list
    rank$present_in_g1 <- as.factor(rank$present_in_g1)
  }
  
  pheatmap(sig_means_niche_genes_table, annotation_row=rank,
           color = colorRampPalette(c("white", "red"))(100),
           annotation_colors = mycolors, fontsize_row = row_font_size)
}

####################
# ANALYSIS
####################

# Figure parameters
height = 6
width = 8

# Read in file of interesting interactions (same as above) to subset cellphonedb results
goi_list <- read.table("niche_interactions.txt", header=FALSE)

## Examine WNT/FZD/LRP interactions
# Read in file of interesting interactions (same as above) to subset cellphonedb results
wnt_ints <- read.table("wnt_interactions.txt", header=FALSE)

# Heatmap of GOIs in F1 with only cell types of interest
f1_goi_hmap <- cpdb_heatmap(goi_list$V1, "cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt")
ggsave(f1_goi_hmap, file="f1_goi_hmap.pdf", width=width, height=height)
ggsave(f1_goi_hmap, file="f1_goi_hmap.png", width=width, height=height)

# Heatmap of WNT ints in F1 with only cell types of interest
f1_wnt_hmap <- cpdb_heatmap(wnt_ints$V1, "cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt")
ggsave(f1_wnt_hmap, file="f1_wnt_hmap.pdf", width=width, height=12)
ggsave(f1_wnt_hmap, file="f1_wnt_hmap.png", width=width, height=12)

# Heatmap of GOIs in G1 with all cell types
g1_goi_hmap <- cpdb_heatmap(goi_list$V1, "cpdb_heo_g1/statistical_analysis_significant_means_cpdb_heo_g1.txt")
ggsave(g1_goi_hmap, file="g1_goi_hmap.pdf", width=width+1, height=height+2)
ggsave(g1_goi_hmap, file="g1_goi_hmap.png", width=width+1, height=height+2)

# Heatmap of WNT ints in G1 with all cell types
g1_wnt_hmap <- cpdb_heatmap(wnt_ints$V1, "cpdb_heo_g1/statistical_analysis_significant_means_cpdb_heo_g1.txt", 
                            row_font_size = 6)
ggsave(g1_wnt_hmap, file="g1_wnt_hmap.pdf", width=width, height=14)
ggsave(g1_wnt_hmap, file="g1_wnt_hmap.png", width=width, height=14)

# Identify interactions which are significant for Intermediate.mesoderm.HCs or HCs.Intermediate.mesoderm in F1
# Make a heatmap of these
res_means_f1 <- read.table("cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt", header=TRUE, sep='\t')
# Select only those which are significant between cell types of interest and have good rank
res_means_f1_subset <- res_means_f1 %>% filter((!is.na(Intermediate.mesoderm.HCs) | !is.na(HCs.Intermediate.mesoderm)) &
                                           rank < 0.2)
hc_im_ints <- res_means_f1_subset$interacting_pair
cpdb_heatmap(hc_im_ints, "cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt")

# Identify which ones are present in G1 interactions
res_means_g1 <- read.table("cpdb_heo_g1/statistical_analysis_significant_means_cpdb_heo_g1.txt", header=TRUE, sep='\t')
res_means_g1_sig <- res_means_g1[rowSums(is.na(res_means_g1[,15:50])) < 36, ]
f1_sig_g1_present <- subset(res_means_g1_sig, interacting_pair %in% hc_im_ints)$interacting_pair

cpdb_heatmap(hc_im_ints, "cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt",
             f1_sig_g1_present)

# Show F1 WNT interactions annotated with which are present in G1
f1_wnt_g1_present <- subset(res_means_g1_sig, interacting_pair %in% wnt_ints$V1)$interacting_pair
f1_wnt_hmap_g1res <- cpdb_heatmap(wnt_ints$V1, "cpdb_heo_f1/statistical_analysis_significant_means_cpdb_heo_f1.txt",
             f1_wnt_g1_present)
ggsave(f1_wnt_hmap_g1res, file="f1_wnt_hmap_g1pres.pdf", width=width, height=12)
ggsave(f1_wnt_hmap_g1res, file="f1_wnt_hmap_g1pres.png", width=width, height=12)

#################
# Write out the session info to provide hope for future package nightmares
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
