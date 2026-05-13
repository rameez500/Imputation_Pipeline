############################################
##### Pre-Imputation Processing and QC #####
############################################

## ============================================================================
## METADATA (EDIT FOR YOUR DATASET)
## ============================================================================

## Analyst:
## Name of data:
## Location:
## Genome build:
## Platform:
## Strand:
## Number variants:
## Number of samples:
## Original filename:

## ============================================================================
## LOAD LIBRARIES
## ============================================================================

library(data.table)   # Fast file I/O for large genetic data
library(ggplot2)      # Visualization
library(tidyverse)    # Data manipulation
library(foreach)      # Parallel processing


## ============================================================================
## EXAMINE ORIGINAL DATA
## ============================================================================

## Check total number of variants
system("wc -l filepath-to-original-data/data.bim")

## Look at markers in data -- check build using dbSNP: https://www.ncbi.nlm.nih.gov/snp/
bim <- fread("filepath-to-original-data/data.bim")
bim[3000:3010,]  # Look in middle of file, not just head

## Run BIM file through Chipendium to determine platform:
## http://mccarthy.well.ox.ac.uk/chipendium/ui/

## Look at FAM file structure
fam <- fread("filepath-to-original-data/data.fam")
head(fam)
nrow(fam)  # Number of samples


## ============================================================================
## STEP 0: LIFT DATA TO NCBI37 (if needed)
## ============================================================================

## Download strand files from: http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html
## Usage: update_build.sh <bed-stem> <strand-file> <output-stem>

system("sh /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/update_build.sh original_data PATH-TO-STRAND-FILE data_37")

## If no lift needed, run through PLINK to standardize
system("plink --bfile filepath-to-original-data/data --make-bed --out data_37")

system("wc -l data_37.fam")  # Verify sample count
system("wc -l data_37.bim")  # Verify variant count


## ============================================================================
## STEP 1: QC THE SAMPLE DATA
## ============================================================================

## Remove SNPs with genotyping rate < 95% (missing in >5% of samples)
system("plink --bfile data_37 --geno 0.05 --threads 20 --make-bed --out temp_data_37_geno")

## Keep only SNPs with MAF > 10% (pre-imputation threshold)
system("plink --bfile temp_data_37_geno --maf 0.1 --threads 20 --make-bed --out temp_data_37_maf")

system("wc -l temp_data_37_maf.fam")  # Remaining samples
system("wc -l temp_data_37_maf.bim")  # Remaining variants


## ============================================================================
## STEP 2: MERGE REFERENCE AND SAMPLE DATA TO IDENTIFY OVERLAPPING SNPs
## ============================================================================

## 1000 Genomes LD pruned reference file
kg_bim <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_temp_1kg_Reference_All_Biallelic_SNPS.bim")
dim(kg_bim)  # 2,050,463 variants
kg_bim_snps <- as.matrix(kg_bim$V2)

## Read sample BIM
data_bim <- fread("temp_data_37_maf.bim")
data_bim_snps <- as.matrix(data_bim$V2)

## Determine overlap between 1KG and sample
table(kg_bim_snps %in% data_bim_snps)
overlap <- kg_bim_snps[kg_bim_snps %in% data_bim_snps]
length(overlap)

## Save overlapping SNPs list
fwrite(as.list(overlap), "data_updated_temp_1kg_ref_snps.txt", sep = " ", quote = FALSE, col.names = FALSE)


## ============================================================================
## STEP 3: RESTRICT DATA TO OVERLAPPING SNPs
## ============================================================================

## Create 1KG reference with only overlapping SNPs
system("plink --bfile /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_temp_1kg_Reference_All_Biallelic_SNPS --extract data_updated_temp_1kg_ref_snps.txt --make-bed --out temp_1kg_reference_for_data")

## Trim sample data to match reference
system("plink --bfile temp_data_37_maf --extract data_updated_temp_1kg_ref_snps.txt --make-bed --out temp_data_match_1kgRef")


## ============================================================================
## STEP 4: RAYNER TOOL (STRAND ALIGNMENT)
## ============================================================================

## Generate frequency file (required by Rayner)
system("plink --bfile temp_data_match_1kgRef --freq --out temp_data_match_1kgRef")

## Run HRC-1000G-check-bim.pl
## Options:
##   -b: BIM file
##   -f: Frequency file
##   -r: Reference legend (1000GP_Phase3_combined.legend)
##   -g: Generate PLINK script
##   -p: Population (ALL for 1000G)
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.2.11_Oct2019/HRC-1000G-check-bim.pl -b temp_data_match_1kgRef.bim -f temp_data_match_1kgRef.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/1000GP_Phase3_combined.legend -g -p ALL")

## Run the generated PLINK script
system("mv Run-plink.sh Run_plink_data.sh")
system("chmod u+x Run_plink_data.sh")
system("./Run_plink_data.sh")

## Check results
system("wc -l Exclude-temp_data_match_1kgRef-1000G.txt")  # Excluded variants
system("wc -l temp_data_match_1kgRef-updated.bim")        # Retained variants


## ============================================================================
## STEP 5: PREPARE FILE FOR PCA
## ============================================================================

## Create 1KG reference excluding problematic variants
system("plink --bfile temp_1kg_reference_for_data --exclude Exclude-temp_data_match_1kgRef-1000G.txt --make-bed --out temp_data_match_1kgRef_raynor")

## Merge updated data with 1KG for PCA
system("plink --bfile temp_data_match_1kgRef-updated --bmerge temp_data_match_1kgRef_raynor --threads 20 --make-bed --out temp_data_match_1kgRef-updated_1KG")


## ============================================================================
## STEP 6: FLASHPCA
## ============================================================================

## Run PCA (20 dimensions)
system("flashpca --bfile temp_data_match_1kgRef-updated_1KG --ndim 20 --suffix _data.txt --numthreads 20 > pca_data.log &")

## Read PCA results
pcs <- fread("pcs_data.txt")
dim(pcs)  # n_samples+1kg x 22

## Load 1000 Genomes super-population labels
thousand_index <- fread("/projects/bga_lab/DATA_REPOSITORIES/1000Genome/1000Genome_labels.csv")
sub_thous <- thousand_index[, c(1,3)]
colnames(sub_thous)[1] <- "IID"

## Merge population labels with PCA data
PC_20_merged <- merge(pcs, sub_thous, by = "IID", all.x = TRUE)

## Label your data (samples not in 1KG)
PC_20_merged$GROUP[is.na(PC_20_merged$GROUP)] <- "DATA"

## Plot ancestry
ggplot(PC_20_merged, aes(x = PC1, y = PC2, color = GROUP)) +
  geom_point(size = 0.2) +
  theme_classic()
ggsave("data_pc2vpc1.png", width = 5, height = 5, units = "in")

ggplot(PC_20_merged, aes(x = PC2, y = PC3, color = GROUP)) +
  geom_point(size = 0.2) +
  theme_classic()
ggsave("data_pc2vpc3.png", width = 5, height = 5, units = "in")


## ============================================================================
## STEP 7: SELECT ANCESTRAL GROUPS
## ============================================================================

## Extract only your data samples
data <- PC_20_merged[PC_20_merged$GROUP == "DATA", ]

## Load ancestry selection functions
source("select_ancestral.R")
source("multiscale_ancestral.R")

## Calculate population statistics and thresholds
dat_select_ancestral <- select_ancestral(PC_20_merged)

## View thresholds for each ancestry
dat_select_ancestral$afr_thresholds
dat_select_ancestral$eur_thresholds
dat_select_ancestral$sas_thresholds
dat_select_ancestral$eas_thresholds
dat_select_ancestral$amr_thresholds

## Combine all thresholds
all_thresholds <- cbind(
  dat_select_ancestral$afr_thresholds,
  dat_select_ancestral$eur_thresholds,
  dat_select_ancestral$sas_thresholds,
  dat_select_ancestral$eas_thresholds,
  dat_select_ancestral$amr_thresholds
)
colnames(all_thresholds) <- c("afr", "eur", "sas", "eas", "amr")
rownames(all_thresholds) <- c("upper_pc1", "lower_pc1", "upper_pc2", "lower_pc2", "upper_pc3", "lower_pc3")

## Apply thresholds to select European ancestry samples
data_eur_1up <- data[data$PC1 <= all_thresholds[1, 2], ]
data_eur_1lo <- data_eur_1up[data_eur_1up$PC1 >= all_thresholds[2, 2], ]
data_eur_2up <- data_eur_1lo[data_eur_1lo$PC2 <= all_thresholds[3, 2], ]
data_eur_2lo <- data_eur_2up[data_eur_2up$PC2 >= all_thresholds[4, 2], ]
data_eur_3up <- data_eur_2lo[data_eur_2lo$PC3 <= all_thresholds[5, 2], ]
data_eur_3lo <- data_eur_3up[data_eur_3up$PC3 >= all_thresholds[6, 2], ]

## Save keep lists
data_eur_3lo_iids <- data_eur_3lo[, c(2, 1)]  # FID, IID format
fwrite(data_eur_3lo_iids, "data_eur_3pcs_keeplist.txt", col.names = FALSE, sep = " ")

## Repeat for other ancestries (AFR, SAS, EAS, AMR) as needed


## ============================================================================
## STEP 8: MULTIDIMENSIONAL OUTLIER DETECTION
## ============================================================================

## Get PCA files for each ancestry
eur_pcs <- data[data$IID %in% data_eur_3lo_iids$IID, ]
fwrite(eur_pcs, "data_eur_pcs.txt")

## Run multidimensional scaling outlier detection
multidimensional_scaling <- multidimensional_scaling(afr_pcs, eur_pcs, sas_pcs, eas_pcs, amr_pcs)

## Check outliers
head(multidimensional_scaling$m_eur_pcs)
table(multidimensional_scaling$m_eur_pcs$outlier)

## Keep non-outliers
eur_mds_keep <- multidimensional_scaling$m_eur_pcs[
  multidimensional_scaling$m_eur_pcs$outlier == 0,
]
fwrite(eur_mds_keep[, c(2, 1)], "data_eur_3pcs_keeplist_mds.txt", col.names = FALSE, sep = " ")


## ============================================================================
## STEP 9: EXTRACT SUPER-POPULATIONS FROM SAMPLE DATA
## ============================================================================

## Create PLINK subsets for each ancestry
system("plink --bfile temp_data_37 --keep data_eur_3pcs_keeplist_mds.txt --make-bed --out data_37_eur")


## ============================================================================
## STEP 10: POPULATION-SPECIFIC RAYNER TOOL
## ============================================================================

## For EUROPEAN (EUR) - Use HRC reference
system("plink --bfile data_37_eur --freq --out data_37_eur")
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b data_37_eur.bim -f data_37_eur.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -t 0.1")
system("mv Run-plink.sh Run_plink_data_37_eur.sh")
system("chmod u+x Run_plink_data_37_eur.sh")
system("./Run_plink_data_37_eur.sh")

## For AFRICAN (AFR) - Use CAAPA reference (build 37)
system("plink --bfile data_37_afr --freq --out data_37_afr")
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b data_37_afr.bim -f data_37_afr.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/CAAPA_FILES/all.caapa.sorted.txt -h -t .1")

## For AFRICAN build 38 - Use TOPMED reference
system("sh /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/update_build.sh data_37_afr /scratch/silo1/BGA_LAB/dbGaP/HumanOmni1-Quad_v1-0_H-b38.Source.strand data_38_afr")
system("plink --bfile data_38_afr --freq --out data_38_afr")
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b data_38_afr.bim -f data_38_afr.frq -h -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/TOPMED_REF/PASS.Variantsbravo-dbsnp-all.tab -t .1")

## For SOUTH ASIAN (SAS) and EAST ASIAN (EAS) - Use Asian 100K reference
system("plink --bfile data_37_sas --freq --out data_37_sas")
system("perl /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b data_37_sas.bim -f data_37_sas.frq -r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/asian_100K/ASIA.Genome.Reference.legend -g -p ALL -t 0.1")


## ============================================================================
## STEP 11: QC BY CHROMOSOME
## ============================================================================

## Split by chromosome (Rayner outputs per-chromosome files)
## For each chromosome, apply:
##   --geno 0.05: Remove SNPs missing in >5%
##   --maf 0.01:  Keep MAF > 1% (post-imputation threshold)
##   --mind 0.1:  Remove samples missing >10%

chrm <- 1:22

## Example for EUR
for(i in chrm) {
  cmd_geno <- paste0("plink --bfile data_37_eur-updated-chr", i, 
                     " --geno 0.05 --threads 20 --make-bed --out data_37_eur-updated-chr", i, "_geno")
  system(cmd_geno)
  
  cmd_maf <- paste0("plink --bfile data_37_eur-updated-chr", i, 
                    "_geno --maf 0.01 --threads 20 --make-bed --out data_37_eur-updated-chr", i, "_maf")
  system(cmd_maf)
  
  cmd_mind <- paste0("plink --bfile data_37_eur-updated-chr", i, 
                     "_maf --mind 0.1 --threads 20 --make-bed --out data_37_eur-updated-chr", i, "_mind")
  system(cmd_mind)
}


## ============================================================================
## STEP 12: PREPARE VCF FOR IMPUTATION SERVER
## ============================================================================

## Convert PLINK to VCF
for(i in chrm) {
  cmd <- paste0("plink --bfile data_37_eur-updated-chr", i, 
                "_mind --recode vcf --out data_37_eur-updated-chr", i)
  system(cmd)
}

## Sort and compress VCF
for(i in chrm) {
  cmd <- paste0("vcf-sort data_37_eur-updated-chr", i, 
                ".vcf | bgzip -c > data_37_eur-updated-chr", i, ".vcf.gz")
  system(cmd)
}

## ============================================================================
## READY FOR IMPUTATION SERVER
## ============================================================================

## Upload VCFs to: https://imputationserver.sph.umich.edu
## Select appropriate reference panel:
##   - EUR: HRC.r1-1.2016 (build 37)
##   - AFR: CAAPA (build 37) or TOPMED (build 38)
##   - SAS/EAS: Asia-100K
