################################################################################
#                                                                              #
#                    PRE-IMPUTATION PROCESSING AND QC PIPELINE                 #
#                                                                              #
################################################################################

# ==============================================================================
# METADATA (EDIT FOR YOUR DATASET)
# ==============================================================================

# Analyst:          [Your Name]
# Data name:        [Dataset Name]
# Location:         [Storage Path]
# Genome build:     hg19/GRCh37
# Platform:         Illumina/Affymetrix
# Strand:           Plus/Minus
# Variants:         [Number]
# Samples:          [Number]
# Original file:    [Filename]

# ==============================================================================
# 0. SETUP - Load libraries & create directories
# ==============================================================================

# Load required packages
library(data.table)
library(ggplot2)
library(tidyverse)

# Create directory structure
dirs <- c("logs", "temp", "output", "figures")
for(d in dirs) {
  if(!dir.exists(d)) dir.create(d)
}

# Set paths (MODIFY THESE)
INPUT_PREFIX <- "filepath-to-original-data/data"
OUTPUT_PREFIX <- "data_37"
REF_1KG <- "/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/IMPUTE2_temp_1kg_Reference_All_Biallelic_SNPS"
HRC_SCRIPT <- "/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl"

# QC thresholds
GENO_THRESH <- 0.05   # Remove SNPs missing in >5% of samples
MAF_PRE <- 0.10       # Pre-imputation MAF threshold (keep SNPs with >10% MAF)
MAF_POST <- 0.01      # Post-imputation MAF threshold
MIND_THRESH <- 0.10   # Remove samples missing >10% of genotypes

# ==============================================================================
# 1. EXAMINE ORIGINAL DATA
# ==============================================================================

cat("\n========== STEP 1: Examining original data ==========\n")

# Check file sizes and counts
system("wc -l filepath-to-original-data/data.bim")
system("wc -l filepath-to-original-data/data.fam")

# Look at first few variants
bim <- fread("filepath-to-original-data/data.bim")
cat("\nFirst 5 variants:\n")
print(head(bim, 5))
cat("\nVariants 3000-3010:\n")
print(bim[3000:3010, ])

# Check sample info
fam <- fread("filepath-to-original-data/data.fam")
cat("\nFirst 5 samples:\n")
print(head(fam))
cat("\nTotal samples:", nrow(fam), "\n")

# ==============================================================================
# 2. LIFT DATA TO NCBI37 (if needed)
# ==============================================================================

cat("\n========== STEP 2: Lifting to build 37 ==========\n")

# Option A: If you need to lift from build 38 to 37
# system("sh /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/update_build.sh original_data /path/to/strand.file data_37")

# Option B: Just standardize existing build 37 files
system(paste("plink --bfile", INPUT_PREFIX, "--make-bed --out", OUTPUT_PREFIX))

# Verify output
cat("\nAfter lift/standardization:\n")
system(paste("wc -l", paste0(OUTPUT_PREFIX, ".fam")))
system(paste("wc -l", paste0(OUTPUT_PREFIX, ".bim")))

# ==============================================================================
# 3. BASIC QUALITY CONTROL
# ==============================================================================

cat("\n========== STEP 3: Basic QC ==========\n")

# Step 3.1: Remove SNPs with low genotyping rate (>5% missing)
system(paste("plink --bfile", OUTPUT_PREFIX, 
             "--geno", GENO_THRESH, 
             "--threads 20", 
             "--make-bed --out temp_geno"))

# Step 3.2: Keep only common variants (MAF > 10%)
system(paste("plink --bfile temp_geno", 
             "--maf", MAF_PRE, 
             "--threads 20", 
             "--make-bed --out temp_maf"))

# Check results
cat("\nAfter QC:\n")
system(paste("wc -l temp_maf.fam"))  # Samples remaining
system(paste("wc -l temp_maf.bim"))  # Variants remaining

# ==============================================================================
# 4. FIND OVERLAPPING SNPS WITH 1000 GENOMES REFERENCE
# ==============================================================================

cat("\n========== STEP 4: Finding overlapping SNPs with 1000G ==========\n")

# Load reference SNPs
kg_bim <- fread(paste0(REF_1KG, ".bim"))
kg_snps <- kg_bim$V2
cat("Reference SNPs:", length(kg_snps), "\n")

# Load data SNPs
data_bim <- fread("temp_maf.bim")
data_snps <- data_bim$V2
cat("Data SNPs:", length(data_snps), "\n")

# Find overlap
overlap <- kg_snps[kg_snps %in% data_snps]
cat("Overlapping SNPs:", length(overlap), 
    sprintf("(%.1f%% of reference)", 100 * length(overlap)/length(kg_snps)), "\n")

# Save overlapping SNP list
fwrite(list(overlap), "overlapping_snps.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ==============================================================================
# 5. RESTRICT DATA TO OVERLAPPING SNPS
# ==============================================================================

cat("\n========== STEP 5: Restricting to overlapping SNPs ==========\n")

# Subset reference panel
system(paste("plink --bfile", REF_1KG, 
             "--extract overlapping_snps.txt", 
             "--make-bed --out temp_1kg_subset"))

# Subset data
system(paste("plink --bfile temp_maf", 
             "--extract overlapping_snps.txt", 
             "--make-bed --out temp_data_matched"))

cat("Data and reference now have matching SNPs\n")

# ==============================================================================
# 6. STRAND ALIGNMENT USING RAYNER TOOL
# ==============================================================================

cat("\n========== STEP 6: Strand alignment (Rayner) ==========\n")

# Generate frequency file (required by Rayner)
system("plink --bfile temp_data_matched --freq --out temp_data_matched")

# Run HRC check script
system(paste(
  "perl", HRC_SCRIPT,
  "-b temp_data_matched.bim",
  "-f temp_data_matched.frq", 
  "-r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/1000GP_Phase3_combined.legend",
  "-g -p ALL"
))

# Run the generated PLINK script
system("mv Run-plink.sh Run_plink_data.sh")
system("chmod u+x Run_plink_data.sh")
system("./Run_plink_data.sh")

# Check results
cat("\nRayner results:\n")
system("wc -l Exclude-temp_data_matched-1000G.txt")  # Excluded variants
system("wc -l temp_data_matched-updated.bim")        # Retained variants

# ==============================================================================
# 7. MERGE WITH 1000 GENOMES FOR PCA
# ==============================================================================

cat("\n========== STEP 7: Merging with 1000G for PCA ==========\n")

# Remove problematic variants from reference
system(paste("plink --bfile temp_1kg_subset", 
             "--exclude Exclude-temp_data_matched-1000G.txt", 
             "--make-bed --out temp_1kg_clean"))

# Merge data with reference
system(paste("plink --bfile temp_data_matched-updated", 
             "--bmerge temp_1kg_clean", 
             "--threads 20", 
             "--make-bed --out merged_1kg"))

# If merge fails due to strand issues
if(!file.exists("merged_1kg.bim")) {
  cat("\nMerge failed, flipping mismatched SNPs...\n")
  system("plink --bfile temp_data_matched-updated --flip merged_1kg-merge.missnp --make-bed --out temp_flipped")
  system(paste("plink --bfile temp_flipped", 
               "--bmerge temp_1kg_clean", 
               "--make-bed --out merged_1kg"))
}

# ==============================================================================
# 8. RUN PCA
# ==============================================================================

cat("\n========== STEP 8: Running PCA ==========\n")

# Run FlashPCA
system("flashpca --bfile merged_1kg --ndim 20 --suffix _data.txt --numthreads 20 > pca.log 2>&1 &")

# Wait for completion (check every 10 seconds)
while(!file.exists("pcs_data.txt")) {
  Sys.sleep(10)
  cat(".")
}
cat("\nPCA complete!\n")

# Read PCA results
pcs <- fread("pcs_data.txt")
cat("PCA dimensions:", dim(pcs), "\n")

# ==============================================================================
# 9. VISUALIZE PCA WITH POPULATION LABELS
# ==============================================================================

cat("\n========== STEP 9: Visualizing ancestry ==========\n")

# Load 1000 Genomes population labels
thousand_labels <- fread("/projects/bga_lab/DATA_REPOSITORIES/1000Genome/1000Genome_labels.csv")
colnames(thousand_labels)[1] <- "IID"
colnames(thousand_labels)[3] <- "POPULATION"

# Merge labels with PCA
pcs_with_pop <- merge(pcs, thousand_labels[, .(IID, POPULATION)], by = "IID", all.x = TRUE)

# Label our data (samples not in 1000G)
pcs_with_pop[, GROUP := ifelse(is.na(POPULATION), "OUR_DATA", POPULATION)]

# Create PCA plots
p1 <- ggplot(pcs_with_pop, aes(x = PC1, y = PC2, color = GROUP)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "PCA: PC1 vs PC2", 
       x = paste0("PC1 (", round(100 * var(pcs$PC1)/sum(apply(pcs[,2:21], 2, var)), 1), "%)"),
       y = paste0("PC2 (", round(100 * var(pcs$PC2)/sum(apply(pcs[,2:21], 2, var)), 1), "%)")) +
  theme(legend.position = "right")

p2 <- ggplot(pcs_with_pop, aes(x = PC2, y = PC3, color = GROUP)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "PCA: PC2 vs PC3") +
  theme(legend.position = "right")

# Save plots
ggsave("figures/pca_pc1vpc2.png", p1, width = 8, height = 6, dpi = 300)
ggsave("figures/pca_pc2vpc3.png", p2, width = 8, height = 6, dpi = 300)

cat("PCA plots saved to 'figures/' directory\n")

# ==============================================================================
# 10. SELECT EUROPEAN ANCESTRY SAMPLES
# ==============================================================================

cat("\n========== STEP 10: Selecting European ancestry samples ==========\n")

# Source ancestry selection functions
source("select_ancestral.R")
source("multiscale_ancestral.R")

# Calculate ancestry thresholds
ancestry_thresholds <- select_ancestral(pcs_with_pop)

# View thresholds
cat("\nEuropean ancestry thresholds:\n")
print(ancestry_thresholds$eur_thresholds)

# Apply thresholds to select European samples
data_only <- pcs_with_pop[GROUP == "OUR_DATA", ]

# Filter using PC1, PC2, PC3 thresholds
eur_selected <- data_only %>%
  filter(PC1 <= ancestry_thresholds$eur_thresholds[1,1],   # Upper PC1
         PC1 >= ancestry_thresholds$eur_thresholds[2,1],   # Lower PC1
         PC2 <= ancestry_thresholds$eur_thresholds[3,1],   # Upper PC2
         PC2 >= ancestry_thresholds$eur_thresholds[4,1],   # Lower PC2
         PC3 <= ancestry_thresholds$eur_thresholds[5,1],   # Upper PC3
         PC3 >= ancestry_thresholds$eur_thresholds[6,1])   # Lower PC3

cat("\nSelected", nrow(eur_selected), "European samples out of", nrow(data_only), "total\n")

# Save keep list
eur_keep <- eur_selected[, .(FID = V1, IID)]
fwrite(eur_keep, "output/eur_samples_keep.txt", sep = " ", col.names = FALSE)

# ==============================================================================
# 11. MULTIDIMENSIONAL OUTLIER DETECTION
# ==============================================================================

cat("\n========== STEP 11: Removing outliers ==========\n")

# Run multidimensional outlier detection
mds_result <- multidimensional_scaling(eur_pcs = eur_selected)

# Check outlier status
cat("\nOutlier summary:\n")
print(table(mds_result$m_eur_pcs$outlier))

# Keep only non-outliers
eur_clean <- mds_result$m_eur_pcs[mds_result$m_eur_pcs$outlier == 0, ]

# Save final keep list
eur_final_keep <- eur_clean[, .(FID = V1, IID)]
fwrite(eur_final_keep, "output/eur_final_keep.txt", sep = " ", col.names = FALSE)

cat("Final European samples after outlier removal:", nrow(eur_clean), "\n")

# ==============================================================================
# 12. EXTRACT EUROPEAN SUBSET
# ==============================================================================

cat("\n========== STEP 12: Extracting European subset ==========\n")

system(paste("plink --bfile", OUTPUT_PREFIX, 
             "--keep output/eur_final_keep.txt", 
             "--make-bed --out data_37_eur"))

# ==============================================================================
# 13. POPULATION-SPECIFIC STRAND ALIGNMENT (EUROPEAN)
# ==============================================================================

cat("\n========== STEP 13: European-specific alignment ==========\n")

# Generate frequency file
system("plink --bfile data_37_eur --freq --out data_37_eur")

# Run Rayner with HRC reference
system(paste(
  "perl", HRC_SCRIPT,
  "-b data_37_eur.bim",
  "-f data_37_eur.frq",
  "-r /projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab",
  "-h -t 0.1"
))

# Apply corrections
system("mv Run-plink.sh Run_plink_eur.sh")
system("chmod u+x Run_plink_eur.sh")
system("./Run_plink_eur.sh")

# ==============================================================================
# 14. PER-CHROMOSOME QC
# ==============================================================================

cat("\n========== STEP 14: Per-chromosome QC ==========\n")

# Process chromosomes 1-22
for(chr in 1:22) {
  cat("Processing chromosome", chr, "...\n")
  
  # Apply QC filters per chromosome
  system(paste("plink --bfile data_37_eur-updated-chr", chr, 
               "--geno", GENO_THRESH, 
               "--threads 20", 
               "--make-bed --out data_37_eur-chr", chr, "_geno"))
  
  system(paste("plink --bfile data_37_eur-chr", chr, "_geno",
               "--maf", MAF_POST, 
               "--threads 20", 
               "--make-bed --out data_37_eur-chr", chr, "_maf"))
  
  system(paste("plink --bfile data_37_eur-chr", chr, "_maf",
               "--mind", MIND_THRESH, 
               "--threads 20", 
               "--make-bed --out data_37_eur-chr", chr, "_final"))
}

cat("Per-chromosome QC complete\n")

# ==============================================================================
# 15. CONVERT TO VCF FOR IMPUTATION SERVER
# ==============================================================================

cat("\n========== STEP 15: Creating VCF files ==========\n")

# Convert each chromosome to VCF
for(chr in 1:22) {
  cat("Converting chromosome", chr, "to VCF...\n")
  
  # Convert to VCF
  system(paste("plink --bfile data_37_eur-chr", chr, "_final",
               "--recode vcf --out data_37_eur-chr", chr))
  
  # Sort and compress
  system(paste("vcf-sort data_37_eur-chr", chr, ".vcf | bgzip -c > data_37_eur-chr", chr, ".vcf.gz"))
  
  # Index
  system(paste("tabix -p vcf data_37_eur-chr", chr, ".vcf.gz"))
}

cat("\n========================================\n")
cat("✓ PIPELINE COMPLETE!\n")
cat("========================================\n")
cat("\nOutput files:\n")
cat("  - European subset: data_37_eur\n")
cat("  - VCF files: data_37_eur-chr*.vcf.gz\n")
cat("  - PCA plots: figures/pca_*.png\n")
cat("  - Sample lists: output/eur_*_keep.txt\n")
cat("\nReady for upload to imputation server:\n")
cat("  https://imputationserver.sph.umich.edu\n")
cat("  Recommended reference panel: HRC.r1-1.2016 (build 37)\n")
cat("\n========================================\n")
