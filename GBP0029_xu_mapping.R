# GBP0029_xu_mapping.R

# AUTHOR: Adam Reid
# Copyright (C) 2025 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# Here I am mapping our datasests to the Xu et al reference in order to show what 
# sort of cells we have, compared to in vivo embryos at a similar timepoint

library(symphony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scmap) # Sankey diagrams
library(pheatmap)

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

#####################
# Build Symphony reference for Xu et al data
#####################

# Read in raw count expression matrix
xu_data <- Read10X("../../../external_data/Xu_et_al/", gene.column=1)

# Read in meta data for cells
xu_meta_data <- read.csv("../../../external_data/Xu_et_al/GSE157329_cell_annotate.txt", header=TRUE, sep='\t')

# Build Seurat object for reference data
xu_seurat <- CreateSeuratObject(counts = xu_data)

# Add meta data to reference Seurat object
xu_seurat@meta.data <- cbind(xu_seurat@meta.data, xu_meta_data)

# Calculate cell cycle stage predictions to include in Harmonization
# Define cell cycle genes
#s.genes <- n2e(heo_f1, cc.genes$s.genes)
#g2m.genes <- n2e(heo_f1, cc.genes$g2m.genes)
#xu_seurat <- CellCycleScoring(xu_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Filter out cell types which are not needed (and not shown in the Xu et al main figure) - n.b. invert = TRUE!
xu_seurat_filter = subset(xu_seurat, subset=developmental.system %in% c('fibroblast', 'PGC', 'miscellaneous'), invert=TRUE)

# Build Symphony reference
## Use embyro as batch factor
## normalise, because we have raw counts to start with
## use same normalisation for the query data - e.g. start with raw counts
## n.b. I need a theta for each harmonization variable - I don't really understand what this value means
reference = symphony::buildReference(
  xu_seurat_filter@assays$RNA@counts,
  xu_seurat_filter@meta.data,
  #vars = c('embryo', 'Phase'),         # variables to integrate over
  vars = c('embryo'),
  #theta = c(2, 2),          # Theta parameter(s) (cluster diversity enforcement?) for Harmony step
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = TRUE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = 'embryo', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = './xu_uwot_model'
)

saveRDS(reference, "xu_symphony_reference.rds")

reference <- readRDS("xu_symphony_reference.rds")

#########################
# Read in and map F1 data
########################
heo_f1 <- readRDS("../heo_F1_labelled.rds")

# Map query using raw counts, allowing Symphony to do the normalisation
# so that reference and query are normalised in the same way
query_f1 = mapQuery(heo_f1@assays$RNA[],             # query gene expression (genes x cells)
                 heo_f1@meta.data,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = TRUE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP

query_f1 = knnPredict(query_f1, reference, reference$meta_data$developmental.system, k = 5)

head(query_f1$meta_data)

# Add label predictions to F1 Seurat object
heo_f1@meta.data$cell_type_pred_knn = query_f1$meta_data$cell_type_pred_knn
heo_f1@meta.data$cell_type_pred_knn_prob = query_f1$meta_data$cell_type_pred_knn_prob

### Here I was trying to set those calls with low prob to 'unknown' or 'NA'
### I think this may have worked but I'm not sure how i did it
#cell_type_pred_knn_filt = heo_f1@meta.data$cell_type_pred_knn
#heo_f1@meta.data[which(heo_f1@meta.data$cell_type_pred_knn_prob < 0.5),] <- 'unknown'
#test <- heo_f1@meta.data$cell_type_pred_knn_prob
#test[which(test< 0.8)] <- 0

# Plot mapped labels onto original UMAP and compare to our calls
p1 <- DimPlot(heo_f1, group.by="cell_type_pred_knn", cols = custom_colors$discrete)
p2 <- DimPlot(heo_f1, cols = custom_colors$discrete)
p1+p2

# Plot probabilities to see if particular clusters have low probability (not really)
FeaturePlot(heo_f1, features=c("cell_type_pred_knn_prob"))

### Plot on reference UMAP
# Add reference-based UMAP coordinates to F1 Seurat object as reduction "ref_umap"
heo_f1[["ref_umap"]] <- CreateDimReducObject(embeddings = query_f1$umap, key = "ref_umap_", 
                        assay = DefaultAssay(heo_f1))

# Plot F1 cells on reference-based UMAP
DimPlot(heo_f1, reduction="ref_umap", group.by="cell_type_pred_knn", cols = custom_colors$discrete)

### Add Symphony UMAP data to the Xu et al Seurat object
ref_umap <- as.data.frame(reference$umap$embedding)
rownames(ref_umap) <- reference$meta_data$cell_id
xu_seurat_filter[['umap_symph']] <- CreateDimReducObject(embeddings = as.matrix(ref_umap),
                                                  key = "umap_symph", assay = DefaultAssay(xu_seurat_filter))

# Plot Xu et al. UMAP and F1 reference-based UMAP side-by-side
q1 <- DimPlot(xu_seurat_filter, reduction="umap_symph", group.by="developmental.system",
              cols = custom_colors$discrete, label=TRUE, raster=FALSE,
              shuffle=TRUE) + 
              xlim(-6.7, 11.7) +ylim(-11.4, 12.9)
q2 <- DimPlot(heo_f1, reduction="ref_umap", group.by="cell_type_pred_knn", 
              cols = custom_colors$discrete, raster=FALSE) + 
              xlim(-6.7, 11.7) +ylim(-11.4, 12.9)
q1 + q2

ggsave(q1, file="xu_ref_symphony_umap.pdf")

# Make combined cell type annotation for F1 (day 14 late)
day14_df <- as.data.frame(heo_f1@active.ident)
day14_df$new <- heo_f1@active.ident

day14_df$new <- gsub('Intermediate mesoderm [1-3]', 'Splanchnic mesoderm', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Mesenchyme [1-3]', 'Mesenchyme', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Cardiomyocytes [1-5]', 'Cardiomyocytes', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Erythroids [1-5]', 'Erythroids', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Hepatocytes [1-5]', 'Hepatocytes', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Hepatocytes [1-5]', 'Hepatocytes', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Blood progenitors 1', 'Blood progenitors (HC)', day14_df$new, perl = TRUE)
day14_df$new <- gsub('Blood progenitors 2', 'Blood progenitors (EC)', day14_df$new, perl = TRUE)

heo_f1@meta.data$cell_type_combined <- day14_df$new
DimPlot(heo_f1, group.by="cell_type_combined")

# Calculaute proportion of cells in each reference class
f1_heatmap <- pheatmap(table(heo_f1@meta.data$cell_type_combined, heo_f1@meta.data$cell_type_pred_knn)/rowSums(table(heo_f1@meta.data$cell_type_combined, heo_f1@meta.data$cell_type_pred_knn)), 
                       scale='none', cluster_rows=TRUE, cluster_cols = TRUE)

f1_heatmap
ggsave(f1_heatmap, file="heo_f1_xu_mapping_heatmap.pdf", width = 5, height = 4)


############
# Map G1
###########
# Read in G1 data 
heo_g1 <- readRDS("../heo_G1_labelled.rds")

# Map query using raw counts, allowing Symphony to do the normalisation
query_g1 = mapQuery(heo_g1@assays$RNA[],             # query gene expression (genes x cells)
                    heo_g1@meta.data,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    do_normalize = TRUE,  # perform log(CP10k+1) normalization on query
                    do_umap = TRUE)        # project query cells into reference UMAP

# Call labels from reference for G1
query_g1 = knnPredict(query_g1, reference, reference$meta_data$developmental.system, k = 5)

head(query_g1$meta_data)

# Add predictions to G1 Seurat object
heo_g1@meta.data$cell_type_pred_knn = query_g1$meta_data$cell_type_pred_knn
heo_g1@meta.data$cell_type_pred_knn_prob = query_g1$meta_data$cell_type_pred_knn_prob
 
# Make combined cell type annotation for G1 (day 14 early)
day14_early_df <- as.data.frame(heo_g1@active.ident)
day14_early_df$new <- heo_g1@active.ident

day14_early_df$new <- gsub('Mesenchyme [1-4]', 'Mesenchyme', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Neural progenitors [1-2]', 'Neural progenitors', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Neural crest [1-2]', 'Neural crest', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Blood progenitors', 'Blood progenitors (EC)', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Mesenchyme [1-4]', 'Mesenchyme', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Hepatocytes [1-2]', 'Hepatocytes', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Musculoskeletal progenitors/cardiomyocytes', 'Musculoskeletal progenitors', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Hepatocytes [1-2]', 'Hepatocytes', day14_early_df$new, perl = TRUE)
day14_early_df$new <- gsub('Mixed/cardiac mesoderm', 'Mixed mesoderm', day14_early_df$new, perl = TRUE)

heo_g1@meta.data$cell_type_combined <- day14_early_df$new
DimPlot(heo_g1, group.by="cell_type_combined")

# Calculaute proportion of cells in each reference class
g1_heatmap <- pheatmap(table(heo_g1@meta.data$cell_type_combined, heo_g1@meta.data$cell_type_pred_knn)/rowSums(table(heo_g1@meta.data$cell_type_combined, heo_g1@meta.data$cell_type_pred_knn)), 
                       scale='none', cluster_rows=TRUE, cluster_cols = TRUE)
g1_heatmap
ggsave(g1_heatmap, file="heo_g1_xu_mapping_heatmap.pdf", height=5, width=6)

# Plot mapped labels onto original UMAP and compare to our calls
p1_g1 <- DimPlot(heo_g1, group.by="cell_type_pred_knn", cols = custom_colors$discrete)
p2_g1 <- DimPlot(heo_g1, cols = custom_colors$discrete, label = TRUE)
p1_g1+p2_g1

# Add reference UMAP coords to G1 seurat object
heo_g1[["ref_umap"]] <- CreateDimReducObject(embeddings = query_g1$umap, 
                                             key = "ref_umap_", assay = DefaultAssay(heo_g1))

# Plot predicted labels on reference-based UMAP for G1
q2_g1 <- DimPlot(heo_g1, reduction="ref_umap", group.by="cell_type_pred_knn", 
                 cols = custom_colors$discrete)  +xlim(-7, 12.3) +ylim(-10.71, 11.64)
q1 + q2_g1

# Call more detailed annotations for G1
query_g1_finalann = knnPredict(query_g1, reference, reference$meta_data$final_annotation, k = 5)

# Add predictions to Seurat object
heo_g1@meta.data$cell_type_pred_knn_finalann = query_g1_finalann$meta_data$cell_type_pred_knn

# Add our original annotations to the metadata table
heo_g1@meta.data$original_labels <- heo_g1@active.ident

# Save Xu-annotated G1
saveRDS(heo_g1, "heo_g1_xu_labels.rds")

#####################
# Combine all three Seurat objects so we can display the mapped cells on top of the reference
###################
# Merge Seurat objects
embryo.combined <- merge(xu_seurat_filter, y = c(heo_f1, heo_g1), project = "embryo")

# Set ident for the reference
embryo.combined@meta.data[which(!(embryo.combined@meta.data$orig.ident %in% c('F1', 'G1'))),]$orig.ident <- 'ref'

# Add merged UMAP reduction
comb_emb <- rbind(xu_seurat_filter[['umap_symph']]@cell.embeddings, heo_f1@reductions$ref_umap@cell.embeddings, 
                  heo_g1@reductions$ref_umap@cell.embeddings)
embryo.combined[["ref_umap"]] <- CreateDimReducObject(embeddings = comb_emb, 
                                                      key = "ref_umap_", assay = DefaultAssay(embryo.combined))

## Plot F1 and G1 cells on integrated UMAP
r1 <- DimPlot(embryo.combined, group.by="orig.ident", cols=c('red', 'green', 'grey'),
              pt.size=0.01)
q1 + r1

# Plot just F1
umap_int_f1 <- DimPlot(subset(embryo.combined, subset=orig.ident %in% c("F1", "ref")), 
                       group.by="orig.ident", cols=c('red', 'grey'), pt.size=0.01,
                       raster=FALSE, shuffle=TRUE)
umap_int_f1
ggsave(umap_int_f1, file="heo_f1_xu_mapping_umap.pdf")

# Plot just G1
umap_int_g1 <- DimPlot(subset(embryo.combined, subset=orig.ident %in% c("G1", "ref")), 
                       group.by="orig.ident", cols=c('blue', 'grey'), pt.size=0.01,
                       raster=FALSE, shuffle=TRUE)
umap_int_g1
ggsave(umap_int_g1, file="heo_g1_xu_mapping_umap.pdf")


####
# Make a Sankey plot
plot(getSankey(heo_f1@active.ident, heo_f1@meta.data$cell_type_pred_knn))

plot(getSankey(heo_g1@active.ident, heo_g1@meta.data$cell_type_pred_knn))


########
# Try subsetted, detailed blood->blood mapping for F1
# Predict cell labels using "final_annotation"
query_f1_finalann = knnPredict(query_f1, reference, reference$meta_data$final_annotation, k = 5)

# Add predictions to Seurat object
heo_f1@meta.data$cell_type_pred_knn_finalann = query_f1_finalann$meta_data$cell_type_pred_knn
heo_f1@meta.data$cell_type_pred_knn_finalann_prob = query_f1_finalann$meta_data$cell_type_pred_knn_prob

# Add our original annotations to the metadata table
heo_f1@meta.data$original_labels <- heo_f1@active.ident

# Count the frequency of detailed blood labels
final_ann_counts <- table(subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn == 'blood')$cell_type_pred_knn_finalann)

# Get only those cells with sensible labels
final_ann_filt <- subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn_finalann %in% 
                           names(subset(final_ann_counts, final_ann_counts > 10)))

# plot Sankey diagram - save the HTML as a PDF for use in Figures
plot(getSankey(final_ann_filt$original_labels, final_ann_filt$cell_type_pred_knn_finalann))

# UMAP plot of detailed reference blood cell labels
DimPlot(heo_f1[,rownames(final_ann_filt)], group.by="cell_type_pred_knn_finalann")

# Plot only those cells from blood progenitors 1 in F1 with the Xu labels
f1_bp_xu_umap <- DimPlot(subset(heo_f1, idents='Blood progenitors 1'), group.by="cell_type_pred_knn_finalann",
        cols=custom_colors$discrete)
ggsave(f1_bp_xu_umap, file="heo_f1_xu_bloodprog_umap.pdf")


saveRDS(heo_f1, "heo_f1_xu_labels.rds")

########################
# Call markers for Xu-based F1 blood cell types
########################
heo_f1_hc_subset <- subset(heo_f1, idents='Blood progenitors 1')

# set active ident to the one we are interested in
heo_f1_hc_subset <- SetIdent(heo_f1_hc_subset, value="cell_type_pred_knn_finalann") 

# Save this as it is important for figures
saveRDS(heo_f1_hc_subset, file="heo_f1_hc_subset.rds")

# Find cell types of interest with at least n cells
n_cells = 20
cois <- subset(levels(heo_f1_hc_subset@meta.data$cell_type_pred_knn_finalann), 
               table(heo_f1_hc_subset@meta.data$cell_type_pred_knn_finalann) > n_cells)

all_top5_markers = list()

# Make a gene names dataframe
gene_names_df = heo_f1@assays$RNA[[]]
gene_names_df$gene <- rownames(heo_f1@assays$RNA[[]])

for (x in cois) {
  markers <- FindMarkers(heo_f1_hc_subset, ident.1 = x, group.by="cell_type_pred_knn_finalann")
  
  # Add gene names
  markers_named <- merge(markers, gene_names_df, by='row.names')
  markers_named$gene <- c()
  
  # Remove dangerous characters from cell type name for output filename
  output_name_string = gsub('[/ -]', '_', x, perl=TRUE)
  output_fn <- paste("heo_f1_hc_subset_", output_name_string, "_markers.txt", sep="")
  
  write.table(markers_named, file=output_fn, quote=FALSE, sep='\t', row.names=FALSE)
  
  all_top5_markers <- append(all_top5_markers, rownames(markers[1:5,]))
}

# Make and write out HC Xu-based dotplot
dp_heo_f1_hc_subset <- DotPlot(heo_f1_hc_subset, features=unique(unlist(all_top5_markers)), 
                        group.by = "cell_type_pred_knn_finalann", assay="SCT", cols="Spectral",
                     idents = cois, dot.scale = 5) +
  scale_x_discrete(labels=e2n(heo_f1, unique(unlist(all_top5_markers)))) + 
  theme(axis.text.x = element_text(angle = 90, size = 20), axis.text.y = element_text(size=20), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
dp_heo_f1_hc_subset
ggsave(dp_heo_f1_hc_subset, file="heo_f1_hc_subset_dotplot.pdf", width=13, height=6)

p1 <- DimPlot(heo_f1_hc_subset)
p2 <- FeaturePlot(heo_f1_hc_subset, features=c(n2e(heo_f1, c("PTPRC"))))
p1 + p2

# Generate a new dotplot with hand-picked markers
hc_markers <- c("KIT","GATA2","GATA3","RUNX1","SOX7","CD34","CDH5","PECAM1",
                "PTPRC","GATA1","SPI1","SOX4","GATA2","IL18","IL1B","TGFB1","STAT5A",
                "CD69","CD4","CD52","CD69","CD3D","MYC","SPINK2","HOPX","IL18","IL7R",
                "NKG7","RUNX1T1","MYH9","CD226","MPL","NFE2","FLI1","GP1BA","GP1BB",
                "ITGA2B","PF4","CD68","CD14","CD86","CD163","MRC1","CXCL8","NOTCH2",
                "IL6R","MERTK","TREM2")

dp_heo_f1_hc_subset_handpicked_markers <- DotPlot(heo_f1_hc_subset, features=unique(n2e(heo_f1_hc_subset, hc_markers)), 
                               group.by = "cell_type_pred_knn_finalann", assay="SCT", cols="Spectral",
                               idents = cois, dot.scale = 5) +
  scale_x_discrete(labels=unique(hc_markers)) + 
  theme(axis.text.x = element_text(angle = 90, size = 20), axis.text.y = element_text(size=20), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
dp_heo_f1_hc_subset_handpicked_markers
ggsave(dp_heo_f1_hc_subset_handpicked_markers, file="heo_f1_hc_subset_dotplot_handpicked_markers.pdf", width=17, height=6)

# Investigate certainty of annotation calls
DimPlot(subset(heo_f1, idents='Blood progenitors 1'), group.by="cell_type_pred_knn_finalann_prob")

coi_metdata <- subset(subset(heo_f1, idents='Blood progenitors 1')@meta.data, 
       subset(heo_f1, idents='Blood progenitors 1')@meta.data$cell_type_pred_knn_finalann %in% cois)

# Convert from factor to character for plotting
coi_metdata$cell_type_pred_knn_finalann <- as.character(coi_metdata$cell_type_pred_knn_finalann)

# Plot distributions of probabilities for each cell type
boxplot(cell_type_pred_knn_finalann_prob ~ cell_type_pred_knn_finalann, data=coi_metdata)

#########################
# Examine cardiomyocytes/heart cells
#########################
######### F1
# Load xu-labelled heo_f1 object
heo_f1 <- readRDS("heo_f1_xu_labels.rds")

## Do all splnchnic LPM cells
# Count the frequency of detailed splanchnic LPM labels
final_ann_counts <- table(subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn == 'splanchnic LPM')$cell_type_pred_knn_finalann)

# Get only those cells with sensible labels
final_ann_filt <- subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn_finalann %in% 
                           names(subset(final_ann_counts, final_ann_counts > 50)))

# plot Sankey diagram - save the HTML as a PDF for use in Figures
plot(getSankey(final_ann_filt$original_labels, final_ann_filt$cell_type_pred_knn_finalann))

### Just splanchnic lpm: heart cells
f1_splanchnic_heart_counts <- table(subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn_finalann %in% c("atria cardiomyocyte−1",
                              "atria cardiomyocyte−2","ventricle cardiomyocyte−1","ventricle cardiomyocyte−2",
                              "atrioventricular canal","sinoatrial node (SAN)","epicardium","epicardial derived cell−1","epicardial derived cell−2",
                               "cardiomyocyte−1","pericyte (myocardium)","Second heart field (SHF)"))$cell_type_pred_knn_finalann)

f1_splanchnic_heart_counts_filt <- subset(heo_f1@meta.data, heo_f1@meta.data$cell_type_pred_knn_finalann %in% 
                                            names(subset(f1_splanchnic_heart_counts, f1_splanchnic_heart_counts > 10)))

# Subset to just cardiomyocytes and mesenchyme from our data
f1_splanchnic_heart_counts_filt_refilt <- subset(f1_splanchnic_heart_counts_filt, f1_splanchnic_heart_counts_filt$original_labels %in%
                                                   c("Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", "Cardiomyocytes 4",
                                                     "Cardiomyocytes 5", "Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3"))

# plot Sankey diagram - save the HTML as a PDF for use in Figures
plot(getSankey(f1_splanchnic_heart_counts_filt_refilt$original_labels, f1_splanchnic_heart_counts_filt_refilt$cell_type_pred_knn_finalann))

# Calculate percentages of cells from each of our cell types still included
table(f1_splanchnic_heart_counts_filt_refilt$original_labels) * 100 / table(heo_f1@meta.data$original_labels)

DimPlot(heo_f1)

cardio_f1 <- subset(heo_f1, idents = c("Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", "Cardiomyocytes 4"))

DimPlot(cardio_f1, group.by = "cell_type_pred_knn_finalann", label = TRUE)

# Subset cells with major new annotations
cardio_f1_major <- subset(cardio_f1 , subset = cell_type_pred_knn_finalann %in% subset(levels(cardio_f1@meta.data$cell_type_pred_knn_finalann), table(cardio_f1@meta.data$cell_type_pred_knn_finalann) > 20))

f1_cardio_xu_transfer <- DimPlot(cardio_f1_major, group.by = "cell_type_pred_knn_finalann", label = TRUE,
        repel=TRUE, cols = custom_colors$discrete)
ggsave(f1_cardio_xu_transfer, file="heo_f1_cardio_xu_mapping.pdf")

#################################
###### Heart cell mapping for G1
# Read in Xu-labelled Seurat object
heo_g1 <- readRDS("heo_g1_xu_labels.rds")
#########################################
### Just splanchnic lpm: heart cells
g1_splanchnic_heart_counts <- table(subset(heo_g1@meta.data, heo_g1@meta.data$cell_type_pred_knn_finalann %in% c("atria cardiomyocyte−1",
                             "atria cardiomyocyte−2","ventricle cardiomyocyte−1","ventricle cardiomyocyte−2",
                             "atrioventricular canal","sinoatrial node (SAN)","epicardium","epicardial derived cell−1","epicardial derived cell−2",
                               "cardiomyocyte−1","pericyte (myocardium)","Second heart field (SHF)"))$cell_type_pred_knn_finalann)

g1_splanchnic_heart_counts_filt <- subset(heo_g1@meta.data, heo_g1@meta.data$cell_type_pred_knn_finalann %in% 
                                            names(subset(f1_splanchnic_heart_counts, f1_splanchnic_heart_counts > 10)))

# Subset to just cardiomyocytes and mesenchyme from our data
g1_splanchnic_heart_counts_filt_refilt <- subset(g1_splanchnic_heart_counts_filt, g1_splanchnic_heart_counts_filt$original_labels %in%
                                                   c("Cardiomyocytes",
                                                     "Cardiomyocytes 5", "Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3",
                                                     "Mesenchyne 4"))

# plot Sankey diagram - save the HTML as a PDF for use in Figures
plot(getSankey(g1_splanchnic_heart_counts_filt_refilt$original_labels, g1_splanchnic_heart_counts_filt_refilt$cell_type_pred_knn_finalann))

# Calculate percentages of cells from each of our cell types still included
table(g1_splanchnic_heart_counts_filt_refilt$original_labels) * 100 / table(heo_g1@meta.data$original_labels)


###### Cardiomyocyte mapping for E1 
#########################
# Read in and map E1 data
########################
heo_e1 <- readRDS("../heo_E1_labelled.rds")

# Map query using raw counts, allowing Symphony to do the normalisation
# so that reference and query are normalised in the same way
query_e1 = mapQuery(heo_e1@assays$RNA[],             # query gene expression (genes x cells)
                    heo_e1@meta.data,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    do_normalize = TRUE,  # perform log(CP10k+1) normalization on query
                    do_umap = TRUE)        # project query cells into reference UMAP

query_e1 = knnPredict(query_e1, reference, reference$meta_data$developmental.system, k = 5)

head(query_e1$meta_data)

# Add label predictions to E1 Seurat object
heo_e1@meta.data$cell_type_pred_knn = query_e1$meta_data$cell_type_pred_knn
heo_e1@meta.data$cell_type_pred_knn_prob = query_e1$meta_data$cell_type_pred_knn_prob

query_e1_finalann = knnPredict(query_e1, reference, reference$meta_data$final_annotation, k = 5)

# Add predictions to Seurat object
heo_e1@meta.data$cell_type_pred_knn_finalann = query_e1_finalann$meta_data$cell_type_pred_knn

# Add our original annotations to the metadata table
heo_e1@meta.data$original_labels <- heo_e1@active.ident

# Count the frequency of detailed splanchnic LPM labels
final_ann_counts_cardio_e1 <- table(subset(heo_e1@meta.data, heo_e1@meta.data$cell_type_pred_knn == 'splanchnic LPM')$cell_type_pred_knn_finalann)

# Get only those cells with sensible labels
final_ann_filt_cardio_e1 <- subset(heo_e1@meta.data, heo_e1@meta.data$cell_type_pred_knn_finalann %in% 
                           names(subset(final_ann_counts_cardio_e1, final_ann_counts_cardio_e1 > 50)))

# plot Sankey diagram - save the HTML as a PDF for use in Figures
plot(getSankey(final_ann_filt_cardio_e1$original_labels, final_ann_filt_cardio_e1$cell_type_pred_knn_finalann))

DimPlot(heo_e1)

cardio_e1 <- subset(heo_e1, idents = c("Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", "Cardiomyocytes 4",
                                       "Cardiomyocytes 4", "Cardiomyocytes 5", "Cardiomyocytes 6"))

DimPlot(cardio_e1, group.by = "cell_type_pred_knn_finalann", label = TRUE)

# Subset cells with major new annotations
cardio_e1_major <- subset(cardio_e1 , subset = cell_type_pred_knn_finalann %in% subset(levels(cardio_e1@meta.data$cell_type_pred_knn_finalann), table(cardio_e1@meta.data$cell_type_pred_knn_finalann) > 20))

e1_cardio_xu_transfer <- DimPlot(cardio_e1_major, group.by = "cell_type_pred_knn_finalann", label = TRUE,
                                 repel=TRUE, cols = custom_colors$discrete)
e1_cardio_xu_transfer
