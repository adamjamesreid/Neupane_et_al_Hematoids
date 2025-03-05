# GBP0029_blood_analysis_cardiomyocyte_trajectory.R

# AUTHOR: Adam Reid
# Copyright (C) 2025 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# Here I will do a trajectory analysis of cardiomyocytes

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
# Read in data and subset clusters of interest
################

heo_f1 <- readRDS("../heo_F1_labelled.rds")

f1u <- DimPlot(heo_f1, label=TRUE, cols=custom_colors$discrete)

f1u

heo_g1 <- readRDS("../heo_G1_labelled.rds")

g1u <- DimPlot(heo_g1, label=TRUE, cols=custom_colors$discrete)

g1u

# Examine genes of interest
FeaturePlot(heo_f1, features=n2e(heo_f1, c("MESP1", "NKX2-5", "MYH2",
                                         "ACTC1", "TNNT2", "MYL7", "IRX3", "MYH6", "TNNI1")))

cm_f1 <- subset(heo_f1, idents=c("Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", "Cardiomyocytes 4",
                                 "Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3"))

DimPlot(cm_f1, label=TRUE, cols=custom_colors$discrete)

cm_f1@meta.data$old.labels <- cm_f1@active.ident

#################
# Re-process cardiomycoyte and mesenchyme subset
cm_f1 <- SCTransform(cm_f1, verbose = FALSE)
cm_f1 <- RunPCA(cm_f1)
ElbowPlot(cm_f1, ndims=50)

# UMAP
cm_f1 <- RunUMAP(cm_f1, dims = 1:50, n.neighbors = 30, reduction.name = "umap")

# Clustering
cm_f1 <- FindNeighbors(cm_f1, dims = 1:30, verbose = FALSE)
# Louvain with multi-level refinement
cm_f1 <- FindClusters(cm_f1, verbose = TRUE, algorithm=2)

DimPlot(cm_f1, reduction="umap", label = TRUE) + NoLegend()

# Cell-cycle analysis
s.genes <- n2e(cm_f1, cc.genes$s.genes)
g2m.genes <- n2e(cm_f1, cc.genes$g2m.genes)
cm_f1 <- CellCycleScoring(cm_f1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(cm_f1, features = n2e(cm_f1, c("PCNA", "TOP2A", "MCM6", "MKI67")), ncol = 2)
# The phases look well mised in what we might imagine will be pseudotime
cm_f1@meta.data$Phase <- factor(cm_f1@meta.data$Phase)

## PLOT FOR PAPER - at day 14 most cardiomycoytes are cycling
DimPlot(cm_f1, group.by=c("Phase"), reduction="umap")
ggsave(file="cm_f1_umap_phase_plot.png")
ggsave(file="cm_f1_umap_phase_plot.pdf")

## PLOT FOR PAPER - at day 14 there is no continuity between mesenchyme and cardiomyocytes
DimPlot(cm_f1, reduction="umap", label = TRUE, group.by="old.ident")
ggsave(file="cm_f1_umap_plot.png")
ggsave(file="cm_f1_umap_plot.pdf")

## PLOT FOR PAPER - at day 14 basically all the caridomyocytes are expressing contractility genes
# Examine genes of interest
FeaturePlot(cm_f1, reduction = "umap", features=n2e(cm_f1, c("ACTN2", "MYH6", "TTN")))
ggsave(file="cm_f1_umap_contractility_plot.png")
ggsave(file="cm_f1_umap_contractility_plot.pdf")

# Examine genes of interest
FeaturePlot(cm_f1, reduction = "umap", features=n2e(cm_f1, c("MESP1", "NKX2-5", "MYH2",
                                         "ACTC1", "TNNT2", "MYL7", "IRX3", "MYH6", "TNNI1")))

####################
# Pseudotime
####################

# Convert Seurat object to Monocle CDS object
cds<-as.cell_data_set(cm_f1, assay = "RNA")
rowData(cds)$gene_short_name <- cm_f1@assays$RNA[[]]$gene_name

# Do some preprocessing
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 30)

# Seems to be causing a problem - complaints about duplicate columns names
cds$sample_name <- c()

cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition=F)
plot_cells(cds,
           color_cells_by = "old.labels",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "ident",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# Choose node 1
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 2,
           trajectory_graph_color = "black")

# The docs say that expression_family NB is better for smaller datasets,
# trying quasipoisson here
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", 
                              expression_family = "quasipoisson", cores=4)

# There are >5000 de genes by q-value < 0.01 !!!
de_genes <- subset(cds_pr_test_res, cds_pr_test_res$q_value == 0)
de_genes$gene_name <- e2n(cm_f1, rownames(de_genes))

head(de_genes[order(de_genes$morans_I, decreasing=TRUE),], n =5)

cds_pr_test_res[n2e(cm_f1, "TNNC1"),]
rownames(head(de_genes[order(de_genes$morans_I, decreasing=TRUE),], n=5))
# Plot top genes
goi_f1 <- rownames(head(de_genes[order(de_genes$morans_I, decreasing=TRUE),], n =5))
plot_cells(cds, genes=goi_f1,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           cell_size = 2,
           label_leaves=FALSE)

######################
# Let's look at cardiomyocytes in E1 (day 8) sample
#######################

heo_e1 <- readRDS("../heo_E1_labelled.rds")

e1u <- DimPlot(heo_e1, label=TRUE, cols=custom_colors$discrete)

e1u

# Examine genes of interest
FeaturePlot(heo_e1, features=n2e(heo_e1, c("MESP1", "NKX2-5", "MYH2",
                                          "ACTC1", "TNNT2", "MYL7", "IRX3", "MYH6", "TNNI1")))

cm_e1 <- subset(heo_e1, idents=c("Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", 
                                 "Cardiomyocytes 4", "Cardiomyocytes 5", "Cardiomyocytes 6",
                                 "Mesenchyme 2", "Mesenchyme 1", "Mesenchyme 3", "Mesenchyme 4"))

DimPlot(cm_e1, label=TRUE, cols=custom_colors$discrete)

cm_e1@meta.data$old.labels <- cm_e1@active.ident

# Re-process this subset
cm_e1 <- SCTransform(cm_e1, verbose = FALSE)
cm_e1 <- RunPCA(cm_e1)
ElbowPlot(cm_e1, ndims=50)

# UMAP
cm_e1 <- RunUMAP(cm_e1, dims = 1:50, n.neighbors = 30, reduction.name = "umap")

# Clustering
cm_e1 <- FindNeighbors(cm_e1, dims = 1:30, verbose = FALSE)
# Louvain with multi-level refinement
cm_e1 <- FindClusters(cm_e1, verbose = TRUE, algorithm=2)

DimPlot(cm_e1, reduction="umap", label = TRUE) + NoLegend()

# Plot for paper - main UMAP
cm_e1_umap_plot <- DimPlot(cm_e1, reduction="umap", label = TRUE, group.by="old.labels",
                           repel=TRUE)
ggsave(cm_e1_umap_plot, file="cm_e1_umap_plot.png")
ggsave(cm_e1_umap_plot, file="cm_e1_umap_plot.pdf")

# Cell-cycle analysis
DefaultAssay(cm_e1) <- 'RNA'
s.genes <- n2e(cm_e1, cc.genes$s.genes)
g2m.genes <- n2e(cm_e1, cc.genes$g2m.genes)
cm_e1 <- CellCycleScoring(cm_e1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(cm_e1, features = n2e(cm_e1, c("PCNA", "TOP2A", "MCM6", "MKI67")), ncol = 2)
# The phases look well mised in what we might imagine will be pseudotime
cm_e1@meta.data$Phase <- factor(cm_e1@meta.data$Phase)

# Plot for paper - cycle phase UMAP
cm_e1_umap_phase_plot <- DimPlot(cm_e1, group.by=c("Phase"), reduction="umap")
ggsave(cm_e1_umap_phase_plot, file="cm_e1_umap_phase_plot.png")
ggsave(cm_e1_umap_phase_plot, file="cm_e1_umap_phase_plot.pdf")

# Examine genes of interest
DefaultAssay(cm_e1) <- 'SCT'
FeaturePlot(cm_e1, reduction = "umap", features=n2e(cm_e1, c("MESP1", "NKX2-5", "MYH2",
                                         "ACTC1", "TNNT2", "MYL7", "IRX3", "MYH6", "TNNI1")))

# Contraction genes

FeaturePlot(cm_e1, reduction="umap", features=n2e(cm_e1, c("MYH6", "ACTN2", "TTN")))

# Tyser et al. use these markers to identify the more/less advanced CMs 
FeaturePlot(cm_e1, reduction="umap", features=n2e(cm_e1, c("MEF2C", "GATA4", "PDGFRA", "TNNT2", "MYH6", "TTN")), 
            pt.size = 0.1) 

all_gois = c("ISL1", "FOXC2", "TBX1", "FGF8", "HOXA1", "HOXB1", # SHF
             "MAB21L2", "HAND1", "TBX5", "HCN4", "SFRP5", # FHF
             "MYH6", "ACTN2", "TTN", # Contractility
             "MEF2C", "GATA4", "PDGFRA", "TNNT2") # Tyser markers

DotPlot(cm_e1, features=n2e(cm_e1, all_gois), group.by="old.labels", cluster.idents=TRUE)

####################
# Pseudotime
####################

# Convert Seurat object to Monocle CDS object
cds_e1<-as.cell_data_set(cm_e1, assay = "RNA")
rowData(cds_e1)$gene_short_name <- cm_e1@assays$RNA[[]]$gene_name

# Do some preprocessing
cds_e1 <- estimate_size_factors(cds_e1)
cds_e1 <- preprocess_cds(cds_e1, num_dim = 30)

# Seems to be causing a problem - complaints about duplicate columns names
cds_e1$sample_name <- c()

cds_e1 <- cluster_cells(cds_e1, reduction_method = "PCA")
#plot_cells(cds, color_cells_by = "partition")

cds_e1 <- learn_graph(cds_e1, use_partition=F)

plot_cells(cds_e1,
           color_cells_by = "old.labels",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds_e1,
           color_cells_by = "ident",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# subset to just the cardiomyocytes
#cds_sub <- choose_graph_segments(cds)

# Choose node 1 (can't seem to do this with arguments)
cds_e1 <- order_cells(cds_e1)

# Paper figure - pseudotime trajectory
cm_e1_umap_pt_plot <- plot_cells(cds_e1,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 2,
           trajectory_graph_color = "black")
ggsave(cm_e1_umap_pt_plot, file="cm_e1_umap_pt_plot.png")
ggsave(cm_e1_umap_pt_plot, file="cm_e1_umap_pt_plot.pdf")

# The docs say that expression_family NB is better for smaller datasets,
# trying quasipoisson here
cds_pr_test_res_e1 <- graph_test(cds_e1, neighbor_graph="principal_graph", 
                              expression_family = "quasipoisson", cores=16)

# There are >5000 de genes by q-value < 0.01 !!!
de_genes_e1 <- subset(cds_pr_test_res_e1, (cds_pr_test_res_e1$q_value == 0 & cds_pr_test_res_e1$morans_I >0.5))
de_genes_e1$gene_name <- e2n(cm_e1, rownames(de_genes_e1))

goi_e1 <- rownames(head(de_genes_e1[order(de_genes_e1$morans_I, decreasing=TRUE),], n =5))

cds_pr_test_res_e1[n2e(cm_e1, "TNNC1"),]
rownames(head(de_genes_e1[order(de_genes_e1$morans_I, decreasing=TRUE),], n=5))
# Plot top genes
plot_cells(cds_e1, genes=c("TNNC1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           cell_size = 2,
           label_leaves=FALSE)


plot_cells(cds_e1, genes=c("ENSG00000160808"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

plot_genes_in_pseudotime(cds_e1[goi_e1,],
                         min_expr=0.5)

# Cluster pseudotime genes into modules
gene_module_e1_df <- find_gene_modules(cds_e1[rownames(de_genes_e1),], resolution=c(10^seq(-6,-1)))

# Make a heatmap of gene module/cell cluster scores (not really sure how these are calculated)
cell_group_e1_df <- tibble::tibble(cell=row.names(colData(cds_e1)), 
                                cell_group=colData(cds_e1)$old.labels)
agg_mat_e1 <- aggregate_gene_expression(cds_e1, gene_module_e1_df, cell_group_e1_df)
row.names(agg_mat_e1) <- stringr::str_c("Module ", row.names(agg_mat_e1))
# Order columns by pseudotime - 4, 3, 1, 2 : 5, 2, 4, 3, 6, 1
agg_mat_e1 <- agg_mat_e1[,c("Mesenchyme 4", "Mesenchyme 3", "Mesenchyme 1", "Mesenchyme 2",
             "Cardiomyocytes 5", "Cardiomyocytes 2", "Cardiomyocytes 4", "Cardiomyocytes 3",
             "Cardiomyocytes 6", "Cardiomyocytes 1")]

# Plot known contractility genes in pseudotime
e1_plot_contractility_genes_in_pt <- plot_genes_in_pseudotime(cds_e1[gois_ens,],
                         min_expr=0.5)
ggsave(e1_plot_contractility_genes_in_pt, file="e1_plot_contractility_genes_in_pt.png")
