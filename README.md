# Neupane_et_al_Hematoids
R code invovled in processing, analysing and generating figures from single-cell RNA-seq data associated with Neupane et al. - "Three-dimensional multilineage organogenesis model captures early hematopoietic niche"

The sample names used here are not used in the paper, the relationship is:

* heo_E1 - peEB day 8, SB-late
* heo_F1 - peEB day 14, SB-late
* heo_H10 - peEB day 8, SB-early
* heo_G1 - peEB day 14, SB-early

The initial reads were processed using STARSolo. 

Generated index using human genome GRCh38 gencode v43. 
Fasta: GRCh38.primary_assembly.genome.fa, 
GTF: gencode.v43.basic.annotation.gtf

STAR solo parameters:
STAR version STAR-2.7.10b
--soloFeatures Gene Velocyto 
--soloType CB_UMI_Simple
--soloCBwhitelist 3M-february-2018.txt
--soloUMIlen 12 
--soloCellFilter EmptyDrops_CR 
--soloMultiMappers EM 
--runThreadN 10 

processed data location:
/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/unzip_fastqs/

## Initial identification of cells to filter out

This was done with the inclusion of several datasets which we did not end up publishing and therefore cannot be precisely recreated using our published data.

### Scrublet
Scrublet was run using the “spliced” assay, initially to identify suitable cutoffs for each dataset, then with those cutoffs to generate lists of likely doublets. These were read into the objects in the next step to help generate the list of filtered cells.

### Cell filtering QC
The data were processed using Seurat (cellQC.R) to identify cells to filter out. The data were normalised, integrated and explored to determine filtering cutoffs. The dataset was then filtered and the metadata for retained cells written to a file (heo_filt_meta.tsv).

* This needs the spliced and unpspliced assays, features.edit.tsv and sample.meta.csv
* Uses the scrublet output

## Creation of normalised Seurat objects, dimred, clustering

The data were then reprocessed from scratch, removing those cells identified in the filtering step to produce filtered, normalised Seurat objects with PCA/UMAP (umap30dim30) dimension reductions and clustering (create_seurat.R).

Requires
* features.edit.underscores.tsv

Objects written out:
heo_E1.rds
heo_F1.rds
heo_G1.rds
heo_H10.rds

These objects have several assays: ‘RNA’ - based on the ‘Gene’ counts from STARSolo, ‘spliced’/’unspliced’ based on the Velocyto counts from STARSolo, ‘SCT’ normalised counts using the ‘RNA’ assay.

## Identification of cell types

The manually chosen labels were added to the objects using this script (add_cell_labels.R)

The final set of objects is:

heo_E1_labelled.rds
heo_F1_labelled.rds
heo_G1_labelled.rds
heo_H10_labelled.rds

n.b. currently this does not label the data properly because the clustering has changed (and the UMAPs). Might need to use a cell id -> label meta data file generated from the published Seurat objects

## Mapping to CS12-16 cells from Xu et al. (2023)

This step was actually performed in two different ways in the paper. Here we describe the use of Symphony.

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/paper_figures/GBP0029_xu_mapping.R

Trajectory analysis:
/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/blood_analysis/GBP0029_blood_analysis_f1_cluster12.R

Initial filtering and QC:

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/GBP0029_starsolo_analysis_sctnew.R
* 

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/040723/GBP0029_starsolo_analysis_040723.R, /mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/GBP0029_starsolo_multiple_assays.R

Xu et al. (2023), mapping:




DEG analysis:
/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/paper_figures/GBP0029_endo_blood_prog_comparison.R

Bar plots of cell counts (Figure S6): /mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/gene_expression_analysis.R

