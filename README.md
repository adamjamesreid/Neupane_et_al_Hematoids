# Neupane_et_al_Hematoids
R code invovled in processing, analysing and generating figures from single-cell RNA-seq data associated with Neupane et al. - "A post-implantation embryo model of human development with a definitive hematopoietic niche"

The sample names used here are not used in the paper, the relationship is:

* heo_E1 - SB-late day 8
* heo_F1 - SB-late day 14,
* heo_H10 - SB-early day 8
* heo_G1 - SB-early day 14

The primary data and the *.rds files of containing Seurat objects for each sample are available from [Array Express](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13632)

The initial reads were processed using STARSolo. 

Generated index using human genome GRCh38 gencode v43. 
Fasta: `GRCh38.primary_assembly.genome.fa`, 
GTF: `gencode.v43.basic.annotation.gtf`

STAR solo parameters:

```
STAR version STAR-2.7.10b
--soloFeatures Gene Velocyto 
--soloType CB_UMI_Simple
--soloCBwhitelist 3M-february-2018.txt
--soloUMIlen 12 
--soloCellFilter EmptyDrops_CR 
--soloMultiMappers EM 
--runThreadN 10 
```

The R scripts were developed in R studio are designed to be run interactively. They will require some adjustments to work in alternative environments e.g. packages will need to be installed and paths/filenames changed. If there is anything more complicated please contact ajr236@cam.ac.uk for assistance.

## Initial filtering and QC:

There was an initial processing of the data which ended with a scrublet analysis to idenitfy multiplets. This was somewhat convoluted and is excluded here. We then reprocessed the data from the STARSolo output, excluding multiplets predicted in the previous round of analysis.

Meta data from intial QC round with Scrublet analysis - heo_filt_meta.tsv

Edited features file for STARSolo - features.edit.underscores.tsv

Code to generate Seurat objects from STARSolo output and exclude multiplets  - GBP0029_starsolo_multiple_assays.R

## Integration with CS12-16 cells from Xu et al. (2023) using Symphony 

GBP0029_xu_mapping.R

## Trajectory analysis using Mococle

GBP0029_blood_analysis_f1_cluster12.R

## Differential gene expression analysis

GBP0029_endo_blood_prog_comparison.R

run_topGO3.pl

## CellPhoneDB analysis

Prepare data for CellPhoneDB analysis - GBP0029_revision_cellphonedb.R

run cellphoneDB  - i.e. `python run_cellphonedb.py -d ~/projects/surani/GBP0029/starsolo/130723/cellphonedb/v5.0.0/cellphonedb.zip -m heo_g1_meta.tsv -c heo_g1.h5ad -e heo_g1_microenviroment.tsv -o cpdb_he
o_g1`

## Cardiomyocyte trajectory analysis

GBP0029_blood_analysis_cardiomyocyte_trajectory.R


