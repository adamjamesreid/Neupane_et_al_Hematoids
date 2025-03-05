# GBP0029_blood_analysis_f1_cluster12.R

# AUTHOR: Adam Reid
# Copyright (C) 2025 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# Subcluster Blood Progenitors 1 (aka cluster 12) from F1 sample and look at pseudotime
# to examine development from RUNX1- endothelial cells to RUNX1/CD34/CD45+ 
# hematopoetic cells

library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

set.seed(123)

#################
# FUNCTIONS
###############

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


#################
# Read in data and subset cluster of interest
################

heo_f1 <- readRDS("../heo_F1_labelled.rds")

f1u <- DimPlot(heo_f1, label=TRUE, cols=custom_colors$discrete)

f1u

heo_f1_bp1 <- subset(heo_f1, subset=seurat_clusters==12)

################
# QC the cluster
################

# there are a couple(?) of outliers that probably need removing - one was some dirt on my monitor!
FeaturePlot(heo_f1_bp1, features=n2e(heo_f1_bp1, c("RUNX1")))

# No super-strong signals here
FeaturePlot(heo_f1_bp1, features=c("percent.mt", "percent.ribo", "percent.hem"))

# Cell-cycle analysis
s.genes <- n2e(heo_f1_bp1, cc.genes$s.genes)
g2m.genes <- n2e(heo_f1_bp1, cc.genes$g2m.genes)
heo_f1_bp1 <- CellCycleScoring(heo_f1_bp1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(heo_f1_bp1, features = n2e(heo_f1_bp1, c("PCNA", "TOP2A", "MCM6", "MKI67")), ncol = 2)
# The phases look well mixed in what we might imagine will be pseudotime
heo_f1_bp1@meta.data$Phase <- factor(heo_f1_bp1@meta.data$Phase)
DimPlot(heo_f1_bp1, group.by=c("Phase"))

DimPlot(heo_f1_bp1, label=TRUE, reduction="pca", cols=custom_colors$discrete)
FeaturePlot(heo_f1_bp1, reduction="pca", features=n2e(heo_f1_bp1, c("RUNX1")))

# Remove outlier
heo_f1_bp1 <- heo_f1_bp1[,heo_f1_bp1@reductions$umap30dim30@cell.embeddings[, "UMAP_1"] <13.5]

#################
# recluster
#################

# PCA
heo_f1_bp1 <- RunPCA(heo_f1_bp1)
ElbowPlot(heo_f1_bp1, ndims=50)
DimPlot(heo_f1_bp1, label = TRUE, reduction="pca") + NoLegend()

# UMAP
heo_f1_bp1 <- RunUMAP(heo_f1_bp1, dims = 1:20, n.neighbors = 30, reduction.name = "umap30dim20")
DimPlot(heo_f1_bp1, label = TRUE, reduction="umap30dim20") + NoLegend()
FeaturePlot(heo_f1_bp1, reduction="umap30dim20", features=n2e(heo_f1_bp1, c("RUNX1")))

heo_f1_bp1 <- RunUMAP(heo_f1_bp1, dims = 1:20, n.neighbors = 20, reduction.name = "umap20dim20")
DimPlot(heo_f1_bp1, label = TRUE, reduction="umap20dim20") + NoLegend()

# this one looks good, so going to give it the generic name "umap", which is required by monocle3
heo_f1_bp1 <- RunUMAP(heo_f1_bp1, dims = 1:20, n.neighbors = 20)

# Examine genes of interest
FeaturePlot(heo_f1_bp1, reduction="umap20dim20", features=n2e(heo_f1_bp1, c("RUNX1", "CD34", "PTPRC", "MYB", "CDH5", "ANGPT1", "CD44")))

# Different assays seem to tell the same story for these genes
DefaultAssay(heo_f1_bp1) <- 'RNA'
FeaturePlot(heo_f1_bp1, reduction="umap20dim20", features=n2e(heo_f1_bp1, c("RUNX1", "CD34", "PTPRC", "MYB", "CDH5", "ANGPT1", "CD44")))
DefaultAssay(heo_f1_bp1) <- 'spliced'
FeaturePlot(heo_f1_bp1, reduction="umap20dim20", features=n2e(heo_f1_bp1, c("RUNX1", "CD34", "PTPRC", "MYB", "CDH5", "ANGPT1", "CD44")))
DefaultAssay(heo_f1_bp1) <- 'unspliced'
FeaturePlot(heo_f1_bp1, reduction="umap20dim20", features=n2e(heo_f1_bp1, c("RUNX1", "CD34", "PTPRC", "MYB", "CDH5", "ANGPT1", "CD44")))
DefaultAssay(heo_f1_bp1) <- 'SCT'

# Clustering
heo_f1_bp1 <- FindNeighbors(heo_f1_bp1, dims = 1:20, verbose = FALSE)
# Louvain with multi-level refinement
heo_f1_bp1 <- FindClusters(heo_f1_bp1, verbose = TRUE, algorithm=2, resolution=0.5)

DimPlot(heo_f1_bp1, label = TRUE, reduction="umap20dim20") + NoLegend()

####################
# Pseudotime
####################

# Convert Seurat object to Monocle CDS object
cds<-as.cell_data_set(heo_f1_bp1, assay = "RNA")

# Do some preprocessing
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)

# Seems to be causing a problem - complaints about duplicate columns names
cds$sample_name <- c()

cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition=F)
plot_cells(cds,
           color_cells_by = "ident",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "ident",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# Choose node 1 (can't seem to do this with arguments)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 2,
           trajectory_graph_color = "black")

# The docs say that NB is better for smaller datasets
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", 
                              expression_family = "negbinomial", cores=4)

# There are >5000 de genes by q-value < 0.01 !!!
de_genes <- subset(cds_pr_test_res, cds_pr_test_res$q_value == 0)
de_genes$gene_name <- e2n(heo_f1_bp1, rownames(de_genes))

cds_pr_test_res[n2e(heo_f1_bp1, "RUNX1"),]

#We are getting a lot of ribosomal genes, lets check their signal again
FeaturePlot(heo_f1_bp1, features=c("percent.mt", "percent.ribo", "percent.hem"),
            reduction="umap")
# This does not seem to be driven by cell cycle
DimPlot(heo_f1_bp1, group.by=c("Phase"),
            reduction="umap")

# Set up gene names, which are missing in the Monocle CDS feature data
rowData(cds)$gene_name <- rownames(rowData(cds))
rowData(cds)$gene_short_name <- e2n(heo_f1_bp1, rowData(cds)$gene_name)

# Plot some genes of interest
plot_cells(cds, genes=c("RUNX1", "C1QC", "RPL11", "CD34", "SPN", "PTPRC"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           cell_size = 2,
           label_leaves=FALSE)

plot_cells(cds, genes=c("RUNX1", "CD44", "MYB", "CXCR4", "CDH5", "SOX7", "ERG", "SPN", "PTPRC"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           cell_size = 1.5,
           label_leaves=FALSE)

FeaturePlot(heo_f1_bp1, features=n2e(heo_f1_bp1, c("RUNX1", "CD44", "MYB", "CXCR4", "CDH5", "SOX7", "ERG", "SPN", "PTPRC", "CD34")),
            reduction="umap")

#Remove ribosomal genes from DE list and see what is left
cds_pr_test_res$gene_name <- e2n(heo_f1_bp1, rownames(cds_pr_test_res))
cds_pr_test_res_noribo <- cds_pr_test_res[grep("^RP[LS]", cds_pr_test_res$gene_name, invert=TRUE),]
# Select genes with morans I > 0.4 (~350 genes)
de_genes <- subset(cds_pr_test_res_noribo, cds_pr_test_res_noribo$morans_I > 0.4)

dim(de_genes)

####################
# plot individual genes in pseudotime
####################

# This doesn't look to be working - results come out looking unlike other plots
gois <- n2e(heo_f1_bp1, c("C1QC", "SOX7", "RUNX1", "CD34"))
goi_cds <- cds[rowData(cds)$gene_name %in% gois]

plot_genes_in_pseudotime(goi_cds,
                         color_cells_by="ident",
                         min_expr=0.5)

#####################
# Try plotting pseudotime DE genes on a heatmap
######################

library(pheatmap)

# Gets matrix of de genes from SCT assay
as.matrix(heo_f1_bp1@assays$SCT[c(rownames(de_genes))])

# Get pseudotime values for each cell
heo_f1_bp1$pseudotime <- pseudotime(cds, reduction_method = "UMAP")

#heo_f1_bp1[order(heo_f1_bp1$pseudotime)]@assays$SCT[,c(rownames(de_genes))]

pheatmap(as.data.frame(heo_f1_bp1@assays$SCT[c(rownames(de_genes))]), show_colnames = FALSE,
         labels_row=de_genes$gene_name, scale="row", annotation_col=as.data.frame(heo_f1_bp1$pseudotime))

library(dplyr)

de_gene_pt_order <- as.data.frame(heo_f1_bp1@assays$SCT[c(rownames(de_genes))]) %>% select(order(heo_f1_bp1$pseudotime))
pheatmap(de_gene_pt_order, show_colnames = FALSE,
         labels_row=de_genes$gene_name, scale="row", annotation_col=as.data.frame(heo_f1_bp1$pseudotime), 
         fontsize_row=3, cluster_cols = FALSE)

####################
# Make monocle3 modules of DE genes (i.e. look for correlated expression patterns)
####################
gene_module_df <- find_gene_modules(cds[rownames(de_genes),], resolution=c(10^seq(-6,-1)),
                                    reduction_method = c("UMAP"))

# Visualise the clustering of the genes into modules
plot(gene_module_df$dim_1, gene_module_df$dim_2, col=gene_module_df$module)

# Heatmap with gene module annotation included
pheatmap(de_gene_pt_order, show_colnames = FALSE,
         labels_row=de_genes$gene_name, scale="row", 
         annotation_col=as.data.frame(heo_f1_bp1$pseudotime), 
         annotation_row=gene_module_df %>% select(!dim_1:dim_2) %>% column_to_rownames("id"),
         fontsize_row=3, cluster_cols = FALSE)

# Write out files listing genes in each module for subsequent GO term analysis
for (i in 1:max(as.numeric(gene_module_df$module))) {
    fn = paste("f1_cl12_mod", i, "_gene_list.txt", sep="")
    write.table(gene_module_df[gene_module_df$module==i,]$id, file=fn, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Write table of DE genes 
write.table(cds_pr_test_res[cds_pr_test_res$status=="OK" & cds_pr_test_res$q_value < 0.01,], 
            file="heo_f1_bp1_pseudotime_de_genes.txt", quote=FALSE, sep="\t")


######################
# Save object

saveRDS(heo_f1_bp1, "heo_f1_bp1.rds")

######################
# Load object to continue analysis

heo_f1_bp1 <- readRDS("heo_f1_bp1.rds")


####################
# Generate markers for Seurat clusters

DefaultAssay(heo_f1_bp1) <- 'SCT'
all_markers <- FindAllMarkers(heo_f1_bp1, min.pct = 0.5, assay='SCT', 
                              only.pos = TRUE, min.cells.group = 40, test.use="wilcox",
                              logfc.threshold = 0.5)
all_markers$gene_name <- e2n(heo_f1_bp1, rownames(all_markers))

write.table(all_markers, "heo_f1_bp1_markers.txt", quote=FALSE, sep="\t")

# Save cluster UMAP image
DimPlot(heo_f1_bp1, label = TRUE, reduction="umap20dim20", cols = custom_colors$discrete)
ggsave("heo_f1_bp1_umap.pdf")

#####################
# Dotplot of gois in clusters for the "foot"

gois <- c("CD44", "SPI1", "SPINK2", "CD34", "RUNX1", "MYB", "ANGPT1", "CDH5", 
          "SOX7", "ERG", "PTPRC", "ITGA2B", "IGF1B") 

dp_goi_heo_f1_bp1 <- DotPlot(heo_f1_bp1, features=n2e(heo_f1_bp1, gois), assay="SCT", cols="Spectral",
                         dot.scale = 5, scale = FALSE) +
  scale_x_discrete(labels=gois) + 
  theme(axis.text.x = element_text(angle = 90, size = 20), axis.text.y = element_text(size=20), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
dp_goi_heo_f1_bp1
ggsave(dp_goi_heo_f1_bp1, file="heo_f1_bp1_dotplot_goi.pdf", width=7, height=5)


######################
# Explore hypothesis of a bifurcating trajectory
######################

# Set up gene names, which are missing in the Monocle CDS feature data
rowData(cds)$gene_name <- rownames(rowData(cds))
rowData(cds)$gene_short_name <- e2n(heo_f1_bp1, rowData(cds)$gene_name)

# Learn graph, making a new CDS object
# The minimal_branch_lengh parameter allows us to start deeper in the 'heel'
# but it does make the final part of the calf more complicated
cds_bif <- learn_graph(cds, use_partition=F, learn_graph_control=list(minimal_branch_len = 2, nn.k = NULL))

plot_cells(cds_bif,
           color_cells_by = "ident",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 1,
  #         cell_stroke = 1,
  label_cell_groups = FALSE,
           trajectory_graph_segment_size = 1) + scale_color_manual(values = custom_colors$discrete) +
         theme(legend.position = "none")
ggsave("f1_c12_bif_traj.pdf")

# Choose node 1 (can't seem to do this with arguments)
cds_bif <- order_cells(cds_bif)

plot_cells(cds_bif,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 2,
           trajectory_graph_color = "black")

# Subset cells by their trajectory and add labels back to the original cds object
# This avoids losing the dimension reduction and pseudotime info which is lost
# using 'choose_cells'
cds_calf <- choose_graph_segments(cds_bif)
cds_bif$calf_traj <- rownames(colData(cds_bif)) %in% rownames(colData(cds_calf)) 

cds_toe <- choose_graph_segments(cds_bif)
cds_bif$toe_traj <- rownames(colData(cds_bif)) %in% rownames(colData(cds_toe)) 

# Display some genes of interest in pseudotime along the two different trajectories
# 0, 2, 3 
# MYCT1, MPL, SLC37A1, MPA, GATA1, ITGA2B and CXCL2
# 0, 1, 3
# GAS6, CD4, CD36, FABP3, MRC1 and SLC15A3 
#gois = c("RUNX1", "CD44", "MYB", "CXCR4", "CDH5", "SOX7", "ERG", "SPN", "PTPRC", "CD34")
gois_calf = c("RUNX1", "SPI1", "SOX7", "ERG", "CD34", "CDH5", "PECAM1", "KIT", "GJAG", "SPN", 
              "MYB", "GFI1B", "PTPRC", "GAS6", "CD4", "CD36", "FABP3", "MRC1", "SLC15A3")
gois_toe = c("RUNX1", "SPI1", "SOX7", "ERG", "CD34", "CDH5", "PECAM1", "KIT", "GJAG", "SPN", 
              "MYB", "GFI1B", "PTPRC", "MYCT1", "MPL", "SLC37A1", "MPA", "GATA1", "ITGA2B", "CXCL2")
cds_bif_calf <- cds_bif[rowData(cds_bif)$gene_short_name %in% gois_calf,
                       colData(cds_bif)$calf_traj == TRUE]
cds_bif_toe <- cds_bif[rowData(cds_bif)$gene_short_name %in% gois_toe,
                        colData(cds_bif)$toe_traj == TRUE]
plot_genes_in_pseudotime(cds_bif_calf,
                         color_cells_by="ident",
                         min_expr=0.5) + scale_color_manual(values = c("#FFC312", "#C4E538", "#FDA7DF"))
ggsave("f1_c12_calf_traj_genes.pdf", width=5, height=20)

plot_genes_in_pseudotime(cds_bif_toe,
                         color_cells_by="ident",
                         min_expr=0.5) + scale_color_manual(values = c("#FFC312", "#12CBC4", "#FDA7DF"))
ggsave("f1_c12_toe_traj_genes.pdf", width=5, height=20)

# Call DE genes in pseudotime
cds_bif_pr_test_res <- graph_test(cds_bif, neighbor_graph="principal_graph", 
                              expression_family = "negbinomial", cores=4)

# There are >5000 de genes by q-value < 0.01 !!!
de_genes_bif <- subset(cds_bif_pr_test_res, cds_bif_pr_test_res$q_value < 0.01)
de_genes_bif$gene_name <- e2n(heo_f1_bp1, rownames(de_genes_bif))

write.table(de_genes_bif, file="heo_f1_blood_cluster_trajectory_genes.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


cds_bif_pr_test_res[n2e(heo_f1_bp1, "CYBB"),]


# Examine cluster 0 more carefully
FeaturePlot(heo_f1_bp1, features=n2e(heo_f1_bp1, c("C1QB", "C1QC", "C1QA")), reduction="umap")

# Add boolean variables for membership of the trajectories to the Seurat object
heo_f1_bp1$calf_traj <- cds_bif$calf_traj
heo_f1_bp1$toe_traj <- cds_bif$toe_traj

FeaturePlot(heo_f1_bp1, features="calf_traj", reduction="umap20dim20")
FeaturePlot(heo_f1_bp1, features="toe_traj", reduction="umap20dim20")
