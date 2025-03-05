# Comparison of endothelial blood progenitor cells bwtween -/+ and +/+ at days 8 and 14
# We want DEGs, GO terms and GO term plots
# We have established the equivlaence of the endothelial blood cell clusters (at least at day 14)
# from Symphony mapping to Xu et al.

library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(DropletUtils)
library(scuttle)
library(EnhancedVolcano)

setwd('/mnt/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/paper_figures')

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

##################
# Examine potential trajectories
##################

heo_f1 <- readRDS("../heo_F1_labelled.rds")
heo_e1 <- readRDS("../heo_E1_labelled.rds")
heo_h10 <- readRDS("../heo_H10_labelled.rds")
heo_g1 <- readRDS("../heo_G1_labelled.rds")

f1u <- DimPlot(heo_f1, label=TRUE, cols=custom_colors$discrete)
e1u <- DimPlot(heo_e1, label=TRUE, cols=custom_colors$discrete)
h10u <- DimPlot(heo_h10, label=TRUE, cols=custom_colors$discrete)
g1u <- DimPlot(heo_g1, label=TRUE, cols=custom_colors$discrete)

e1u + h10u

f1u + g1u

####################
# Merge day 8 datasets
####################
# Add cluster names to metadata, so hopefully they get added integrated object
heo_e1$orig_cluster_label <- heo_e1@active.ident
heo_h10$orig_cluster_label <- heo_h10@active.ident

# Try simply merging and then doing DGE
heo_e1_h10 <- merge(heo_e1, y = heo_h10)

# Prep find markers
heo_e1_h10 <- PrepSCTFindMarkers(heo_e1_h10, assay = "SCT", verbose = TRUE)

###################
# Look for Blood progenitor markers between conditions at day 8
####################
# Add feature meta data (gene names)
feature_meta <- read.table("../features.edit.underscores.tsv", header=FALSE, row.names=2, sep="\t")
feature_meta$V3 <- c()
heo_e1_h10[["RNA"]] <- AddMetaData(heo_e1_h10[["RNA"]], feature_meta, col.name = "gene_name")

# If we try and compare the clusters as they are they have the same name, but different sample names
# this causes an error telling us to PrepSCTFindMarkers(), although we already did that.
# Instead, rename the clusters from the different samples so that we just compare the two different clusters
heo_e1_h10@meta.data$sample_cluster <- paste(heo_e1_h10@meta.data$sample_name, heo_e1_h10@meta.data$orig_cluster_label)

# Get condition-specific markers
blood_markers <- FindMarkers(heo_e1_h10, group.by='sample_cluster', ident.1 = "E1 Blood progenitors", ident.2 = 'H10 Blood progenitors', assay='SCT')#, recorrect_umi=TRUE)

# Add gene names
blood_markers$gene_name <- e2n(heo_e1_h10, rownames(blood_markers))
blood_markers_sig <- subset(blood_markers, blood_markers$p_val_adj < 0.05 & (blood_markers$avg_log2FC >0.5 | blood_markers$avg_log2FC < -0.5))

# Write out DEGs to a file
write.table(blood_markers, file="blood_progenitors_de_genes_day8.txt", sep="\t", quote=FALSE)
# Write out all day 8 genes to a file
write.table(blood_markers_sig, file="blood_progenitors_de_genes_day8_all.txt", sep="\t", quote=FALSE)

# Read in markers, rather than recreate them
blood_markers <- read.table("blood_progenitors_de_genes_day8.txt", sep="\t")

# Very basic Volcano plot
plot(blood_markers_sig$avg_log2FC, -log(blood_markers_sig$p_val_adj))

# Load volcano plotting package

eh_day8 <- EnhancedVolcano(blood_markers,
                lab = blood_markers$gene_name,
                selectLab = c("ITGB1", "CXCR4", "NRARP", "GJA4", "IGFBP4", "MYCT1", "MECOM"),
                labSize = 5,
                labFace = "bold",
                subtitle = "",
                x = 'avg_log2FC',
                y = 'p_val_adj',
                axisLabSize = 20,
                #legendLabSize = 20,
                legendIconSize = 10,
                FCcutoff = 0.5,
                pCutoff = 0.05,
                xlim = c(min(blood_markers$avg_log2FC) - 0.1, max(blood_markers$avg_log2FC) + 0.1),
                title = "Endothelial blood progenitors DGEs day 8 -/+ versus +/+") + theme_classic() +
                theme(axis.text=element_text(size=20), axis.title = element_text(size = 16),
                      legend.text=element_text(size=16), legend.position = 'none')
eh_day8
ggsave(eh_day8, file="blood_progenitors_de_genes_day8.volcano.pdf", height=5, width=5)


VlnPlot(heo_e1_h10, features=n2e(heo_f1, c("CALCRL", "NRP1", "CD40")), split.by="orig.ident")

# N.b. we have 2-3x as many blood progenitor cells in H10 (+/+) as E1 (-/+)
table(subset(heo_e1_h10@meta.data, orig_cluster_label == "Blood progenitors")$sample_name)

########################
# Analyse ambient contamination i.e. hemoglobins at day 14
########################
# We need to take account of ambient RNA (hemoglobin) in order to compare F1 and G1

heo_f1$orig_cluster_label <- heo_f1@active.ident
heo_g1$orig_cluster_label <- heo_g1@active.ident

f1_raw_dir = "~/projects/surani/GBP0029/starsolo/data/run_SITTF1_star_Solo.out_copy/Gene/raw/"

f1_raw <- ReadMtx(
  mtx = paste(f1_raw_dir, "matrix.mtx", sep=""), 
  features = "../features.edit.underscores.tsv",
  cells = paste(f1_raw_dir, "barcodes.tsv", sep="")
)

g1_raw_dir = "~/projects/surani/GBP0029/starsolo/data/run_SITTG1_star_Solo.out_copy/Gene/raw/"

g1_raw <- ReadMtx(
  mtx = paste(g1_raw_dir, "matrix.mtx", sep=""), 
  features = "../features.edit.underscores.tsv",
  cells = paste(g1_raw_dir, "barcodes.tsv", sep="")
)

# Determine the (pseudo) gene expression profile for ambient-only cells
f1_ambient <- ambientProfileEmpty(f1_raw, 
                                  good.turing=FALSE, round=FALSE)

g1_ambient <- ambientProfileEmpty(g1_raw, 
                                  good.turing=FALSE, round=FALSE)

# Convert Seurat object to SCE so that we can get counts summed across clusters
##!!! just for for the bp subset
heo_f1_sce <- as.SingleCellExperiment(heo_f1)
heo_g1_sce <- as.SingleCellExperiment(heo_g1)
# Sum counts across clusters
summed_f1 <- aggregateAcrossCells(heo_f1_sce, 
                                  ids=DataFrame(sample=heo_f1_sce$sample_name,
                                                label=heo_f1_sce$orig_cluster_label)
)
summed_g1 <- aggregateAcrossCells(heo_g1_sce, 
                                  ids=DataFrame(sample=heo_g1_sce$sample_name,
                                                label=heo_g1_sce$orig_cluster_label)
)

#summed_f1 has had low count genes filtered out (and subsetted to bp), so need to subset ambient as well
#ambient_filtered <- subset(ambient, rownames(ambient) %in% rownames(summed_f1))
f1_ambient_filtered <- subset(f1_ambient, names(f1_ambient) %in% rownames(summed_f1))
g1_ambient_filtered <- subset(g1_ambient, names(g1_ambient) %in% rownames(summed_g1))

max.ambient.f1 <- ambientContribMaximum(counts(summed_f1), 
                                        f1_ambient_filtered, mode="proportion")
max.ambient.g1 <- ambientContribMaximum(counts(summed_g1), 
                                        g1_ambient_filtered, mode="proportion")

# We want a function which returns the mean contamination and mean count for 
# all genes within a particular sample and across certain clusters
# Proportions of ambient contamination per cluster per genes from ambientContribMaximum()
get_ambient_stats <- function(ambient_prop, pseudobulk_sce, cois, label) {
  #ambient_prop <- max.ambient.f1
  #pseudobulk_sce <- summed_f1
  #cois <- c(6, 8)
  
  # Calculate mean propN for each gene in clusters of interest
  prop_means <- c()
  if(length(cois) > 1) {
    prop_means <- rowMeans(ambient_prop[,cois])
  }
  else{
    prop_means <- ambient_prop[,cois]
  }
  # Calculate mean count for each gene across clusters of interest
  counts_means <- c()
  if(length(cois) > 1) {
    count_means <- rowMeans(counts(pseudobulk_sce)[,cois])
    
  }
  else {
    count_means <- counts(pseudobulk_sce)[,cois]
  }
  # Check that the names all match up
  if(!all(names(prop_means) == names(count_means))) {
    errorCondition("Names do not match up!")
  }
  
  pro_count_res <- cbind(prop_means, count_means)
  colnames(pro_count_res) <- c(paste(label, "_prop_means", sep=""), paste(label, "_count_means", sep=""))
  
  return (pro_count_res)
}

# Run funciton to get mean proportions and counts across genes of interest
ambient.stats.f1 <- get_ambient_stats(max.ambient.f1, summed_f1, c(1), 'f1')
ambient.stats.g1 <- get_ambient_stats(max.ambient.g1, summed_g1, c(1), 'g1')

########################
# Merge day 14 samples and do DGE
#########################

# Try simply merging and then doing DGE
heo_f1_g1 <- merge(heo_f1, y = heo_g1)

# Prep find markers
heo_f1_g1 <- PrepSCTFindMarkers(heo_f1_g1, assay = "SCT", verbose = TRUE)

# Add feature meta data (gene names)
heo_f1_g1[["RNA"]] <- AddMetaData(heo_f1_g1[["RNA"]], feature_meta, col.name = "gene_name")

# If we try and compare the clusters as they are they have the same name, but different sample names
# this causes an error telling us to PrepSCTFindMarkers(), although we already did that.
# Instead, rename the clusters from the different samples so that we just compare the two different clusters
heo_f1_g1@meta.data$sample_cluster <- paste(heo_f1_g1@meta.data$sample_name, heo_f1_g1@meta.data$orig_cluster_label)

# Get condition-specific markers
blood_markers_day14 <- FindMarkers(heo_f1_g1, group.by='sample_cluster', ident.1 = "F1 Blood progenitors 2", ident.2 = 'G1 Blood progenitors', assay='SCT')#, recorrect_umi=TRUE)

# Add gene names
blood_markers_day14$gene_name <- e2n(heo_f1_g1, rownames(blood_markers_day14))
blood_markers_day14_sig <- subset(blood_markers_day14, blood_markers_day14$p_val_adj < 0.05 & 
                      (blood_markers_day14$avg_log2FC >0.5 | blood_markers_day14$avg_log2FC < -0.5))

# Write out DEGs to a file
write.table(blood_markers_day14_sig, file="blood_progenitors_de_genes_day14.txt", sep="\t", quote=FALSE)

# Volcano plot
plot(blood_markers_day14_sig$avg_log2FC, -log(blood_markers_day14_sig$p_val_adj))

VlnPlot(heo_f1_g1, features=c("ENSG00000118271"), split.by="orig.ident")

# Merge data together
blood_markers_day14_sig_ambient <- merge(blood_markers_day14_sig, ambient.stats.f1, all.x= TRUE, by = "row.names")
blood_markers_day14_sig_ambient <- merge(blood_markers_day14_sig_ambient, ambient.stats.g1, all.x= TRUE, by.x = "Row.names", by.y= "row.names")
rownames(blood_markers_day14_sig_ambient) <- blood_markers_day14_sig_ambient$Row.names
blood_markers_day14_sig_ambient$Row.names <- c()
# Add gene names
#blood_markers_day14_sig_ambient <- merge(blood_markers_day14_sig_ambient, heo_f1[["RNA"]]@meta.features, all.x= TRUE, by = "row.names")
# Order by adjusted pvalue again
blood_markers_day14_sig_ambient <- blood_markers_day14_sig_ambient[order(blood_markers_day14_sig_ambient$p_val_adj),]
# Label likely contamination
blood_markers_day14_sig_ambient$contamination <- (blood_markers_day14_sig_ambient$f1_prop_means > 0.3 & 
                                                    blood_markers_day14_sig_ambient$f1_count_means > 100) | 
  (blood_markers_day14_sig_ambient$g1_prop_means > 0.3 & blood_markers_day14_sig_ambient$g1_count_means > 100)

head(blood_markers_day14_sig_ambient, n=20)

write.table(blood_markers_day14_sig_ambient, file="endo_blood_progenitors_de_genes_day14_ambient.txt", quote=FALSE, sep="\t")

# Merge ambient RNA results with full set of DGEs for volcano plotting
blood_markers_day14_ambient <- merge(blood_markers_day14, ambient.stats.f1, all.x= TRUE, by = "row.names")
blood_markers_day14_ambient <- merge(blood_markers_day14_ambient, ambient.stats.g1, all.x= TRUE, by.x = "Row.names", by.y= "row.names")
rownames(blood_markers_day14_ambient) <- blood_markers_day14_ambient$Row.names
blood_markers_day14_ambient$Row.names <- c()
blood_markers_day14_ambient <- blood_markers_day14_ambient[order(blood_markers_day14_ambient$p_val_adj),]
# Label likely contamination
blood_markers_day14_ambient$contamination <- (blood_markers_day14_ambient$f1_prop_means > 0.3 & 
                                                blood_markers_day14_ambient$f1_count_means > 100) | 
  (blood_markers_day14_ambient$g1_prop_means > 0.3 & blood_markers_day14_ambient$g1_count_means > 100)

blood_markers_day14_ambient_filtered = subset(blood_markers_day14_ambient, blood_markers_day14_ambient$contamination==FALSE)

# Write out a file of all DEGs for day 14 with ambient annotation
write.table(blood_markers_day14_ambient, file="endo_blood_progenitors_de_genes_day14_ambient_all.txt", sep="\t", quote=FALSE)

# Read data to recreate volcano plot (rather than regenerating the results again from scratch)
blood_markers_day14_ambient <- read.table("endo_blood_progenitors_de_genes_day14_ambient_all.txt", sep="\t")
blood_markers_day14_ambient_filtered = subset(blood_markers_day14_ambient, blood_markers_day14_ambient$contamination==FALSE)

eh_day14 <- EnhancedVolcano(blood_markers_day14_ambient_filtered,
                           lab = blood_markers_day14_ambient_filtered$gene_name,
                           selectLab = c("TMEM88", "TPM1", "PRTG", "HOXB7", "MEIS2", "RNF19A", "CAV1", 
                                         "REN", "BNIP3", "STC1", "FGF23"),
                           labSize = 5,
                           labFace = "bold",
                           subtitle = "",
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           axisLabSize = 20,
                           #legendLabSize = 20,
                           legendIconSize = 10,
                           FCcutoff = 0.5,
                           pCutoff = 0.05,
                           xlim = c(min(blood_markers_day14_ambient_filtered$avg_log2FC) - 0.1, max(blood_markers_day14_ambient_filtered$avg_log2FC) + 0.1),
                           title = "Endothelial blood progenitors DGEs day 14 -/+ versus +/+") + theme_classic() +
                            theme(axis.text=element_text(size=20), axis.title = element_text(size = 16),
                            legend.text=element_text(size=16), legend.position = 'none')
eh_day14
ggsave(eh_day14, file="blood_progenitors_de_genes_day14.volcano.pdf", height=5, width=8)


##############
# Plot GO term results
###############
# I ran runTopGO3.pl here: ~/projects/surani/GBP0029/starsolo/130723/paper_figures

# Function to filter GO terms and plot as a bar chart and PDF
go_bar_plot <- function(fn, outfile, h, w, q=0.001) {

  day8_up_terms <- read.table(fn, header=FALSE, sep="\t")
  day8_up_terms_filt <- subset(day8_up_terms, (day8_up_terms$V1=='bp' & day8_up_terms$V7 < q))
  names(day8_up_terms_filt) <- c("Part", "Term", "Name", "Total", "Expected", "Observed", "padj", "genes")
  
  day8_up_terms_filt_fact <- day8_up_terms_filt %>% mutate(Name = fct_reorder(Name, padj, .desc=TRUE))
  
  day8_up_term_plot <- ggplot(day8_up_terms_filt_fact, aes(x=Name, y=-log10(padj), fill=Part)) + 
      geom_bar(stat="identity") +
      coord_flip() +
      ylab("-log10(P value)") +
      xlab("") + theme_classic() +
    theme(axis.text = element_text(size = 14))
  day8_up_term_plot
  ggsave(day8_up_term_plot, file=outfile, height=h,
         width =w)  
  return(day8_up_term_plot)
}

go_bar_plot("blood_progenitors_de_genes_day8.up_final.dat", "blood_progenitors_de_genes_day8.up_final.pdf", 4, 7)
go_bar_plot("blood_progenitors_de_genes_day8.down_final.dat", "blood_progenitors_de_genes_day8.down_final.pdf", 4, 7)

#go_term_day14_negpos_filter = c("hepatocyte differentiation", "branching involved in blood vessel morph...",
#                                "maternal process involved in parturition", "regulation of actin cytoskeleton organiz...",
#                                "regulation of actin cytoskeleton organiz...")

# N.b. I have manually filtered the results files here as Jitesh asked for plotting
go_bar_plot("endo_blood_progenitors_de_genes_day14_ambient.up_final.filter.dat", "endo_blood_progenitors_de_genes_day14_ambient.up_final.pdf", 
            4, 7, q=0.001)
go_bar_plot("endo_blood_progenitors_de_genes_day14_ambient.down_final.filter.dat", "endo_blood_progenitors_de_genes_day14_ambient.down_final.pdf", 
            4, 7, q=0.001)

########
# Some checking
VlnPlot(heo_e1_h10, features=n2e(heo_f1, c("PECAM1")), split.by="sample_name")
VlnPlot(heo_f1_g1, features=n2e(heo_f1, c("REN")), split.by="sample_name")
