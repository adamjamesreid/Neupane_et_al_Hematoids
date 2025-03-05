# Run CellPhoneDB 

import pandas as pd
import sys
import argparse
import os
import anndata
from cellphonedb.src.core.methods import cpdb_analysis_method # "method 1"
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method # "method 2"

pd.set_option('display.max_columns', 100)

parser = argparse.ArgumentParser(description="Run CellPhoneDB")
parser.add_argument("-d", help="File path for CellPhoneDB database zip file", default='v5.0.0/cellphonedb.zip')
parser.add_argument("-m", help="Metadata file")
parser.add_argument("-c", help="Count data as H5AD file")
parser.add_argument("-e", help="Microenvironment file")
parser.add_argument("-o", help="Output string")

args = parser.parse_args()

# Read in data
metadata = pd.read_csv(args.m, sep = '\t')

adata = anndata.read_h5ad(args.c)

microenv = pd.read_csv(args.e,
                       sep = '\t')

microenv.groupby('microenvironment', group_keys = False)['cell_type'] \
    .apply(lambda x : list(x.value_counts().index))


# Method 2 - the statistical method
cpdb_results = cpdb_statistical_analysis_method.call(
             cpdb_file_path = args.d,           # mandatory: CellphoneDB database zip file.
    meta_file_path = args.m,           # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = args.c,       # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'ensembl',               # defines the gene annotation in counts matrix.
    microenvs_file_path = args.e, # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                 # optional: whether to score interactions or not.
    output_path = args.o,                    # Path to save results    microenvs_file_path = None,
    separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    threads = 5,                               # number of threads to use in the analysis.
    threshold = 0.1,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
    debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
    debug_seed = 123,               # Set seed to prevent variation in runs
    output_suffix = args.o      # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)
