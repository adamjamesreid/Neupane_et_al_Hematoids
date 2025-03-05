# Neupane_et_al_Hematoids
R code invovled in processing, analysing and generating figures from single-cell RNA-seq data associated with Neupane et al. - "Three-dimensional multilineage organogenesis model captures early hematopoietic niche"

The sample names used here are not used in the paper, the relationship is:

* heo_E1 - SB-late day 8
* heo_F1 - SB-late day 14,
* heo_H10 - SB-early day 8
* heo_G1 - SB-early day 14

The primary data and the *.rds files of containing Seurat objects for each sample are available from [Array Express](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13632)

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


# Initial filtering and QC:

There was an initial processing of the data which ended with a scrublet analysis to idenitfy multiplets. This was somewhat convoluted and is excluded here. We then reprocessed the data from the StarSOLO output, excluding multiplets predicted in the previous round of analysis.

## Meta data from intial QC round with Scrublet analysis
heo_filt_meta.tsv

## Code to generate Seurat objects from STARSolo output and exclude multiplets 
GBP0029_starsolo_multiple_assays.R

# Integration with CS12-16 cells from Xu et al. (2023) using Symphony 

GBP0029_xu_mapping.R

# Trajectory analysis using Mococle

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/blood_analysis/GBP0029_blood_analysis_f1_cluster12.R

# Differential gene expression analysis

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/paper_figures/GBP0029_endo_blood_prog_comparison.R

# Bar plots of cell counts

/mnt/beegfs6/home3/reid/ajr236/projects/surani/GBP0029/starsolo/130723/gene_expression_analysis.R

# CellPhoneDB analysis

# Cardiomyocyte trajectory analysis


