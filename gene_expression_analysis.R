# Detailed investigation of gene expression in hEO datasets

# Load libraries
library(Seurat)

# Set the working directory
setwd('/mnt/bioinfo_sharing/sharing/surani/GBP0029/')

set.seed(123)

########################
# This some set-up code you need to run
########################
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


############################
# Load data
###########################
# Load the Seurat objects
heo_h10 <- readRDS("heo_H10_labelled.rds")
heo_e1 <- readRDS("heo_E1_labelled.rds")
heo_f1 <- readRDS("heo_F1_labelled.rds")
heo_c8 <- readRDS("heo_C8.rds")
heo_g1 <- readRDS("heo_G1.rds")
heo_g9 <- readRDS("heo_G9_labelled.rds")

###########################
# Plot clusters
# e.g. DimPlot(heo_h10, labels = TRUE, col=custom_colors$discrete)
# change "heo_h10" to the dataset object you want

DimPlot(heo_f1, label = TRUE, cols=custom_colors$discrete)

###########################
# Plot the expression values on the UMAP plot
# e.g. FeaturePlot(heo_h10, features=n2e(heo_h10, c("SOX2")))
# change the Seurat object name to the one you are interested in (two places)
# change the gene name for the gene of interest

FeaturePlot(heo_f1, features=n2e(heo_f1, c("STC1")))

###########################
# Violin plot of gene expression across clusters
# E.g. VlnPlot(heo_h10, features=n2e(heo_h10,c("WNT6")))
# change the Seurat object name to the one you are interested in (two places)
# change the gene name for the gene of interest

VlnPlot(heo_g1, features=n2e(heo_g1,c("STC1")), cols=custom_colors$discrete)


##########################
# Plots for Geraldine

library(ggpubr)

DimPlot(heo_f1, label = TRUE, cols=custom_colors$discrete)
VlnPlot(heo_g1, features=c(n2e(heo_g1, "ID2")), assay = "RNA")
VlnPlot(heo_g1, features=c(n2e(heo_f1, "FGF23")))

setwd('~/projects/surani/GBP0029/starsolo/130723/paper_figures/')

for (m in c(heo_f1, heo_g1)) {
  for (x in c("FGF23", "STC1", "ID2", "MIF")) {
    g1 <- plot_gene(m, x)
    outname <- paste(m@meta.data$sample_name[1], '_', x, '.pdf', sep="")
    ggsave(g1, file=outname, height = 15, width = 5)
    #print(outname)
  }  
}

plot_gene <- function(early, gene) {
  e1 <- DimPlot(early, label = TRUE, cols=custom_colors$discrete) + theme(legend.position = "none")
  #l1 <- DimPlot(late, label = TRUE, cols=custom_colors$discrete)
  
  e2 <- FeaturePlot(early, features=c(n2e(early, gene)))
  #l2 <- FeaturePlot(late, features=c(n2e(late, gene)))
  
  e3 <- VlnPlot(early, features=c(n2e(early, gene))) + theme(legend.position = "none")
  #l3 <- VlnPlot(late, features=c(n2e(late, gene)))
  
  #g <- ggarrange(e1, l1, e2, l2, e3, l3, ncol=2, nrow=3)
  g <- ggarrange(e1,  e2, e3, ncol=1, nrow=3, heights = c(5, 5, 8))
  
  return(g)
}

########################
# Barplots for cell type comparison between SB-early and SB-late at day 8 and day 14
########################

library(ggplot2)
library(forcats) # recording bar plots

setwd('/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/paper_figures')

# Load the Seurat objects
heo_h10 <- readRDS("../heo_H10_labelled.rds")
heo_e1 <- readRDS("../heo_E1_labelled.rds")
heo_f1 <- readRDS("../heo_F1_labelled.rds")
heo_g1 <- readRDS("../heo_G1_labelled.rds")

# Massage data into correct format day14 SBlate (F1)
day14_df <- as.data.frame(table(heo_f1@active.ident))
day14_df$treatment <- rep("SB-late", dim(day14_df)[1])

sum(day14_df$Freq)

day14_df <- rbind(day14_df, data.frame(Var1 = "Splanchnic mesoderm",
                           Freq = as.numeric(table(heo_f1@active.ident)["Intermediate mesoderm 1"] + 
                                               table(heo_f1@active.ident)["Intermediate mesoderm 2"] +
                                               table(heo_f1@active.ident)["Intermediate mesoderm 3"]),
                           treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Mesenchyme",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Mesenchyme 1"] + table(heo_f1@active.ident)["Mesenchyme 2"] +
                                                           table(heo_f1@active.ident)["Mesenchyme 3"]),
                                       treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Cardiomyocytes",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Cardiomyocytes 1"] + table(heo_f1@active.ident)["Cardiomyocytes 2"] +
                                                           table(heo_f1@active.ident)["Cardiomyocytes 3"] + table(heo_f1@active.ident)["Cardiomyocytes 4"] +
                                                           table(heo_f1@active.ident)["Cardiomyocytes 5"]),
                                       treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Erythroids",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Erythroids 1"] + table(heo_f1@active.ident)["Erythroids 2"]),
                                       treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Hepatocytes",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Hepatocytes 1"] + table(heo_f1@active.ident)["Hepatocytes 2"]),
                                       treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Blood progenitors (EC)",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Blood progenitors 1"]),
                                       treatment = "SB-late"))
day14_df <- rbind(day14_df, data.frame(Var1 = "Blood progenitors (HC)",
                                       Freq = as.numeric(table(heo_f1@active.ident)["Blood progenitors 2"]),
                                       treatment = "SB-late"))

day14_df <- day14_df[!day14_df$Var1 %in% c("Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3", "Intermediate mesoderm 1",
                               "Cardiomyocytes 1", "Erythroids 1", "Cardiomyocytes 2",
                               "Erythroids 2", "Intermediate mesoderm 2", "Intermediate mesoderm 3",
                               "Cardiomyocytes 3", "Cardiomyocytes 4", "Hepatocytes 1", "Hepatocytes 2",
                               "Cardiomyocytes 5", "Blood progenitors 1", "Blood progenitors 2"),]

sum(day14_df$Freq)
# Add percentages
day14_df$Percentage = day14_df$Freq * 100 / sum(day14_df$Freq)

# Massage data into correct format day14 SBearly (G1)
day14_early_df <- as.data.frame(table(heo_g1@active.ident))
day14_early_df$treatment <- rep("SB-early", dim(day14_early_df)[1])

sum(day14_early_df$Freq)

day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Mesenchyme",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Mesenchyme 1"] + table(heo_g1@active.ident)["Mesenchyme 2"] +
                                                                       table(heo_g1@active.ident)["Mesenchyme 3"] +
                                                                       table(heo_g1@active.ident)["Mesenchyme 4"]),
                                                   treatment = "SB-early"))

day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Neural progenitors",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Neural progenitors 1"] + 
                                                                       table(heo_g1@active.ident)["Neural progenitors 2"]),
                                                 treatment = "SB-early"))

day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Neural crest",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Neural crest 1"] + 
                                                                       table(heo_g1@active.ident)["Neural crest 2"]),
                                                   treatment = "SB-early"))

day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Blood progenitors (EC)",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Blood progenitors"]),
                                                   treatment = "SB-early"))
day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Hepatocytes",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Hepatocytes 1"] + 
                                                                       table(heo_g1@active.ident)["Hepatocytes 2"]),
                                                   treatment = "SB-early"))

day14_early_df <- rbind(day14_early_df, data.frame(Var1 = "Musculoskeletal progenitors",
                                                   Freq = as.numeric(table(heo_g1@active.ident)["Musculoskeletal progenitors/cardiomyocytes"]),
                                                   treatment = "SB-early"))

day14_early_df <- day14_early_df[!day14_early_df$Var1 %in% c("Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3", 
                                           "Mesenchyme 4", "Neural progenitors 1", "Neural progenitors 2",
                                           "Hepatocytes 1", "Hepatocytes 2", "Neural crest 1", "Neural crest 2",
                                           "Blood progenitors", "Musculoskeletal progenitors/cardiomyocytes"),]

# More tricky combination with pre-existing name
mesoderm_temp <- data.frame(Var1 = "Mixed mesoderm",
           Freq = as.numeric(table(heo_g1@active.ident)["Mixed mesoderm"] + 
                               table(heo_g1@active.ident)["Mixed/cardiac mesoderm"]),
           treatment = "SB-early")

day14_early_df <- day14_early_df[!day14_early_df$Var1 %in% c("Mixed mesoderm", "Mixed/cardiac mesoderm"),]
day14_early_df <- rbind(day14_early_df, mesoderm_temp)


sum(day14_early_df$Freq)

day14_early_df$Percentage = day14_early_df$Freq * 100 / sum(day14_early_df$Freq)

# set-up dataframe for plotting
day14_df_combined <- rbind(day14_df, day14_early_df)
day14_df_combined$Var1 <- as.character(day14_df_combined$Var1)
day14_df_combined <- day14_df_combined[order(day14_df_combined$Percentage),]
day14_df_combined$Var1 <- factor(day14_df_combined$Var1, levels=unique(day14_df_combined$Var1))

day14_barplot <- ggplot(day14_df_combined, aes(fill=treatment, y=Percentage, x=Var1)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_minimal() + 
  theme(axis.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size=14), axis.title.y = element_blank(),
        axis.title.x = element_text(size=16)) +labs(fill = "Treatment", y = "Percentage of cells")
day14_barplot
ggsave(day14_barplot, file="day14_earlyvslate_celltype_counts_barplot.pdf", height = 5, width = 6)

### Day 8

# Massage data into correct format day8 SBlate (E1)
day8_late_df <- as.data.frame(table(heo_e1@active.ident))
day8_late_df$treatment <- rep("SB-late", dim(day8_late_df)[1])

sum(day8_late_df$Freq)

day8_late_df <- rbind(day8_late_df, data.frame(Var1 = "Cardiomyocytes",
                                       Freq = as.numeric(table(heo_e1@active.ident)["Cardiomyocytes 1"] + table(heo_e1@active.ident)["Cardiomyocytes 2"] +
                                                           table(heo_e1@active.ident)["Cardiomyocytes 3"] + table(heo_e1@active.ident)["Cardiomyocytes 4"] +
                                                           table(heo_e1@active.ident)["Cardiomyocytes 5"] + table(heo_e1@active.ident)["Cardiomyocytes 6"]),
                                       treatment = "SB-late"))

day8_late_df <- rbind(day8_late_df, data.frame(Var1 = "Mesenchyme",
                                       Freq = as.numeric(table(heo_e1@active.ident)["Mesenchyme 1"] + table(heo_e1@active.ident)["Mesenchyme 2"] +
                                                           table(heo_e1@active.ident)["Mesenchyme 3"] + table(heo_e1@active.ident)["Mesenchyme 4"]),
                                       treatment = "SB-late"))

day8_late_df <- day8_late_df[!day8_late_df$Var1 %in% c("Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3", "Mesenchyme 4",
                                           "Cardiomyocytes 1", "Cardiomyocytes 2", "Cardiomyocytes 3", "Cardiomyocytes 4", 
                                           "Cardiomyocytes 5", "Cardiomyocytes 6"),]

sum(day8_late_df$Freq)
day8_late_df$Percentage = day8_late_df$Freq * 100 / sum(day8_late_df$Freq)

# Massage data into correct format day8 SBearly (H10)
day8_early_df <- as.data.frame(table(heo_h10@active.ident))
day8_early_df$treatment <- rep("SB-early", dim(day8_early_df)[1])

sum(day8_early_df$Freq)

day8_early_df <- rbind(day8_early_df, data.frame(Var1 = "Cardiomyocytes",
                                               Freq = as.numeric(table(heo_h10@active.ident)["Cardiomyocytes 1"] + table(heo_h10@active.ident)["Cardiomyocytes 2"]),
                                               treatment = "SB-early"))

day8_early_df <- rbind(day8_early_df, data.frame(Var1 = "Mesenchyme",
                                               Freq = as.numeric(table(heo_h10@active.ident)["Mesenchyme 1"] + table(heo_h10@active.ident)["Mesenchyme 2"] +
                                                                   table(heo_h10@active.ident)["Mesenchyme 3"] + table(heo_h10@active.ident)["Mesenchyme 4"]),
                                               treatment = "SB-early"))

day8_early_df <- rbind(day8_early_df, data.frame(Var1 = "Paraxial mesoderm",
                                                 Freq = as.numeric(table(heo_h10@active.ident)["Paraxial mesoderm 1"] + table(heo_h10@active.ident)["Paraxial mesoderm 2"]),
                                                 treatment = "SB-early"))

day8_early_df <- day8_early_df[!day8_early_df$Var1 %in% c("Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3", "Mesenchyme 4",
                                                       "Cardiomyocytes 1", "Cardiomyocytes 2", "Paraxial mesoderm 1",
                                                       "Paraxial mesoderm 2"),]

sum(day8_early_df$Freq)

day8_early_df$Percentage = day8_early_df$Freq * 100 / sum(day8_early_df$Freq)

# set-up dataframe for plotting
day8_df_combined <- rbind(day8_late_df, day8_early_df)
day8_df_combined$Var1 <- as.character(day8_df_combined$Var1)
day8_df_combined <- day8_df_combined[order(day8_df_combined$Percentage),]
day8_df_combined$Var1 <- factor(day8_df_combined$Var1, levels=unique(day8_df_combined$Var1))

day8_barplot <- ggplot(day8_df_combined, aes(fill=treatment, y=Percentage, x=Var1)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_minimal() + 
  theme(axis.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size=14), axis.title.y = element_blank(),
        axis.title.x = element_text(size=16)) +labs(fill = "Treatment", y = "Percentage of cells")
day8_barplot

ggsave(day8_barplot, file="day8_earlyvslate_celltype_counts_barplot.pdf", height = 5, width = 6)
                       