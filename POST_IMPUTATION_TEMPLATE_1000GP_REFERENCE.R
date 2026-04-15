#############################################
##### Post-Imputation Processing and QC #####  
#############################################

## ============================================================================
## METADATA (EDIT FOR YOUR DATASET)
## ============================================================================

## analyst:
## date:
## Name of dataset:
## Location of imputed data:
## Ancestry you're working with:
## NOTE: Use this template for admixed populations or when using 1000G reference

## ============================================================================
## LOAD LIBRARIES
## ============================================================================

library(data.table)
library(dplyr)
library(tidyr)

## ============================================================================
## STEPS 0-2: DOWNLOAD AND UNZIP
## ============================================================================

## Same as POST_IMPUTATION_TEMPLATE.R steps 0-2

## ============================================================================
## STEP 3: PREPARE IMPUTED FILES (1000G VERSION)
## ============================================================================

## Read all chromosome info files (same as previous template)
# ... (code identical to previous template for reading info files)

## Create filters for Rsq thresholds (same as previous template)
# ... (code identical for Rsq filtering)

## Create update-name file (same as previous template)
# ... (code identical for update_names)

## DIFFERENCE: Use 1000G reference instead of HRC
list_1000G <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/1000G_FILES/1000GP_reference_legend_biallelic_SNPs_38million_bim_format_nonames.txt", header = FALSE, sep = " ")

keepers_rsq3 <- list_1000G[list_1000G$V2 %in% chrmpos_r2.3, ]
fwrite(keepers_rsq3[, 2], "keeper_markers_rsq3_data_ancestry.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

keepers_rsq7 <- list_1000G[list_1000G$V2 %in% chrmpos_r2.7, ]
fwrite(keepers_rsq7[, 2], "keeper_markers_rsq7_data_ancestry.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

## Create summary table (same as previous template)
# ... (code identical)

## ============================================================================
## STEPS 4-6: MERGE, UPDATE IDs, SCREEN MARKERS
## ============================================================================

## These steps are identical to POST_IMPUTATION_TEMPLATE.R
## The only difference is the reference panel used for filtering

## ============================================================================
## FINAL FILE
## ============================================================================

## data_ancestry_imputed_rsq3_nodup.{bed,bim,fam} is ready for analysis
