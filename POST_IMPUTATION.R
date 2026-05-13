################################################################################
#                                                                              #
#                   POST-IMPUTATION PROCESSING AND QC PIPELINE                 #
#                                                                              #
################################################################################

# ==============================================================================
# METADATA (EDIT FOR YOUR DATASET)
# ==============================================================================

# Analyst:              [Your Name]
# Date:                 [YYYY-MM-DD]
# Dataset name:         [Dataset Name]
# Location of imputed data: [Path to Michigan Server downloads]
# Ancestry:             AFR / EUR / SAS / EAS / AMR
# Imputation password:  [YOUR PASSWORD]

# ==============================================================================
# 0. SETUP - Load libraries & create directories
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)

# Create output directories
dirs <- c("logs", "output", "temp", "summary_tables")
for(d in dirs) {
  if(!dir.exists(d)) dir.create(d)
}

# Set paths (MODIFY THESE)
IMPUTATION_PASSWORD <- "your_password_here"
ANCESTRY <- "eur"  # Change to afr, sas, eas, amr as needed
DATA_NAME <- "your_dataset"

# QC thresholds
RSQ_THRESHOLDS <- c(0.3, 0.5, 0.7, 0.9)
HRC_REF_FILE <- "/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC_reference_biallelic_snps_bim_format_CHRPOS.txt"

# ==============================================================================
# 1. UNZIP ALL FILES FROM MICHIGAN SERVER
# ==============================================================================

cat("\n========== STEP 1: Unzipping chromosome files ==========\n")

# Create unzip script for chromosome zip files
chrm <- 1:22
unzip_commands <- paste0("7z x chr_", chrm, ".zip -p'", IMPUTATION_PASSWORD, "'")
write.table(unzip_commands, file = "unzip_data.sh", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

system("dos2unix unzip_data.sh")
system("chmod +x unzip_data.sh")
system("nohup ./unzip_data.sh > logs/unzip_data.log 2>&1 &")
cat("Unzipping chromosomes in background...\n")

# ==============================================================================
# 2. UNZIP INFO FILES
# ==============================================================================

cat("\n========== STEP 2: Unzipping .info files ==========\n)

# Create unzip script for info files
unzip_info <- paste0("gunzip chr", chrm, ".info.gz")
write.table(unzip_info, file = "unzip_info.sh", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

system("dos2unix unzip_info.sh")
system("chmod +x unzip_info.sh")
system("nohup ./unzip_info.sh > logs/unzip_info.log 2>&1 &")
cat("Unzipping info files in background...\n")

# Wait for files to unzip (check every 30 seconds)
cat("Waiting for files to unzip...\n")
all_unzipped <- FALSE
while(!all_unzipped) {
  Sys.sleep(30)
  info_files_exist <- all(file.exists(paste0("chr", chrm, ".info")))
  if(info_files_exist) {
    all_unzipped <- TRUE
    cat("All files unzipped!\n")
  }
}

# ==============================================================================
# 3. COMBINE INFO FILES AND CALCULATE IMPUTATION STATISTICS
# ==============================================================================

cat("\n========== STEP 3: Processing imputation quality scores ==========\n")

# Read all chromosome info files
cat("Reading info files...\n")
all_chromosomes <- list()

for(i in 1:22) {
  cat("  Reading chr", i, "...\n", sep = "")
  all_chromosomes[[i]] <- fread(paste0("chr", i, ".info"), 
                                header = TRUE, sep = "\t", na.strings = "-")
}

# Combine into single data table
all_imputed <- rbindlist(all_chromosomes)
cat("\nTotal variants imputed:", nrow(all_imputed), "\n")

# Summary of genotyped vs imputed
geno_summary <- table(all_imputed$Genotyped)
cat("\nGenotyped vs Imputed variants:\n")
print(geno_summary)

# Imputation quality statistics
cat("\nImputation quality (Rsq) statistics:\n")
cat("  Mean Rsq:", mean(all_imputed$Rsq, na.rm = TRUE), "\n")
cat("  Median Rsq:", median(all_imputed$Rsq, na.rm = TRUE), "\n")
cat("  Min Rsq:", min(all_imputed$Rsq, na.rm = TRUE), "\n")
cat("  Max Rsq:", max(all_imputed$Rsq, na.rm = TRUE), "\n")
cat("  Five-number summary:\n")
print(fivenum(all_imputed$Rsq))

# Create SNP identifier (chr:pos format)
all_imputed[, chrpos := substr(SNP, 1, nchar(SNP) - 4)]  # Remove :A1:A2 suffix

# Filter by different Rsq thresholds
cat("\nCreating Rsq-filtered datasets...\n")
rsq_filters <- list()

for(r in RSQ_THRESHOLDS) {
  threshold_name <- paste0("rsq", r * 10)
  rsq_filters[[threshold_name]] <- all_imputed[Rsq > r, ]
  cat("  Rsq >", r, ":", nrow(rsq_filters[[threshold_name]]), "variants\n")
  
  # Save keep lists
  fwrite(rsq_filters[[threshold_name]][, .(SNP)], 
         file = paste0("output/keep_", threshold_name, "_", DATA_NAME, ".txt"),
         col.names = FALSE, quote = FALSE)
}

# ==============================================================================
# 4. FILTER TO HRC BIALLELIC VARIANTS
# ==============================================================================

cat("\n========== STEP 4: Filtering to HRC biallelic variants ==========\n")

# Load HRC reference
cat("Loading HRC reference...\n")
hrc_ref <- fread(HRC_REF_FILE, header = FALSE, sep = " ")
hrc_snps <- hrc_ref$V2
cat("HRC reference variants:", length(hrc_snps), "\n")

# Filter each Rsq threshold to HRC variants
biallelic_counts <- c()

for(r in RSQ_THRESHOLDS) {
  threshold_name <- paste0("rsq", r * 10)
  
  # Get chr:pos identifiers for this threshold
  snps_to_keep <- rsq_filters[[threshold_name]]$chrpos
  
  # Find which ones are in HRC
  keep_in_hrc <- hrc_ref[hrc_ref$V2 %in% snps_to_keep, ]
  
  cat("  Rsq >", r, ": Keeping", nrow(keep_in_hrc), "biallelic HRC variants\n")
  
  # Save final keep list
  fwrite(keep_in_hrc[, .(V2)], 
         file = paste0("output/keepers_", threshold_name, "_", DATA_NAME, ".txt"),
         col.names = FALSE, quote = FALSE)
  
  biallelic_counts <- c(biallelic_counts, nrow(keep_in_hrc))
}

# ==============================================================================
# 5. CREATE SUMMARY TABLE
# ==============================================================================

cat("\n========== STEP 5: Creating summary table ==========\n")

# Build comprehensive summary table
summary_table <- data.frame(
  QC_Step = c(
    "Total imputed variants",
    "Genotyped variants",
    "Imputed variants",
    paste0("Rsq > ", RSQ_THRESHOLDS[1]),
    paste0("Rsq > ", RSQ_THRESHOLDS[2]),
    paste0("Rsq > ", RSQ_THRESHOLDS[3]),
    paste0("Rsq > ", RSQ_THRESHOLDS[4]),
    paste0("Biallelic HRC (Rsq > ", RSQ_THRESHOLDS[1], ")"),
    paste0("Biallelic HRC (Rsq > ", RSQ_THRESHOLDS[2], ")"),
    paste0("Biallelic HRC (Rsq > ", RSQ_THRESHOLDS[3], ")"),
    paste0("Biallelic HRC (Rsq > ", RSQ_THRESHOLDS[4], ")")
  ),
  Variants = c(
    nrow(all_imputed),
    geno_summary["1"],
    geno_summary["0"],
    nrow(rsq_filters$rsq3),
    nrow(rsq_filters$rsq5),
    nrow(rsq_filters$rsq7),
    nrow(rsq_filters$rsq9),
    biallelic_counts[1],
    biallelic_counts[2],
    biallelic_counts[3],
    biallelic_counts[4]
  )
)

# Save summary table
write.csv(summary_table, 
          file = paste0("summary_tables/", DATA_NAME, "_", ANCESTRY, "_imputation_summary.csv"),
          row.names = FALSE)

cat("\nSummary table saved!\n")
print(summary_table)

# ==============================================================================
# 6. CONVERT VCF TO PLINK FORMAT
# ==============================================================================

cat("\n========== STEP 6: Converting VCF to PLINK ==========\n")

# Convert each chromosome from VCF to PLINK
for(i in 1:22) {
  cat("Converting chromosome", i, "...\n")
  
  cmd <- paste0("plink --vcf chr", i, ".dose.vcf.gz ",
                "--const-fid ",
                "--make-bed ",
                "--out temp/", DATA_NAME, "_", ANCESTRY, "_imputed_chr", i)
  system(cmd)
}

cat("VCF to PLINK conversion complete\n")

# ==============================================================================
# 7. MERGE ALL CHROMOSOMES
# ==============================================================================

cat("\n========== STEP 7: Merging chromosomes ==========\n")

# Create merge list
merge_list <- paste0("temp/", DATA_NAME, "_", ANCESTRY, "_imputed_chr", 1:22)
writeLines(merge_list, "temp/chrmerge_list.txt")

# Merge all chromosomes
cat("Merging all chromosomes...\n")
system(paste("plink --bfile temp/", DATA_NAME, "_", ANCESTRY, "_imputed_chr1 ",
             "--merge-list temp/chrmerge_list.txt ",
             "--make-bed ",
             "--out output/", DATA_NAME, "_", ANCESTRY, "_imputed_allchr"))

cat("Chromosome merge complete\n")

# ==============================================================================
# 8. UPDATE VARIANT NAMES TO CHR:POS FORMAT
# ==============================================================================

cat("\n========== STEP 8: Updating variant names ==========\n")

# Create update name file (original SNP -> chr:pos)
update_names <- all_imputed[, .(SNP, chrpos)]
fwrite(update_names, "temp/update_names.txt", 
       sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Apply name update
system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_imputed_allchr ",
             "--update-name temp/update_names.txt ",
             "--make-bed ",
             "--out output/", DATA_NAME, "_", ANCESTRY, "_chrpos"))

cat("Variant names updated to chr:pos format\n")

# ==============================================================================
# 9. FIX SAMPLE IDs (IMPUTATION SERVER ADDS EXTRA TEXT)
# ==============================================================================

cat("\n========== STEP 9: Fixing sample IDs ==========\n")

# Read original FAM file (your pre-imputation samples)
original_fam <- fread("filepath-to-original-data/data.fam", header = FALSE)

# Read imputed FAM file
imp_fam <- fread(paste0("output/", DATA_NAME, "_", ANCESTRY, "_chrpos.fam"), 
                 header = FALSE)

# Extract original IDs (imputation server adds _IID suffix)
imp_fam[, IID_clean := gsub(".*_", "", V2)]  # Keep part after underscore
imp_fam[, FID := IID_clean]                  # Use same as IID

# Create new FAM file with correct IDs
new_fam <- imp_fam[, .(FID, IID_clean, V3, V4, V5, V6)]
colnames(new_fam) <- c("FID", "IID", "Pat", "Mat", "Sex", "Pheno")

# Save updated FAM file
fwrite(new_fam, "temp/corrected_ids.fam", 
       sep = " ", col.names = FALSE, quote = FALSE)

# Apply corrected IDs to PLINK files
system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_chrpos ",
             "--fam temp/corrected_ids.fam ",
             "--make-bed ",
             "--out output/", DATA_NAME, "_", ANCESTRY, "_ids_fixed"))

cat("Sample IDs corrected\n")

# ==============================================================================
# 10. APPLY RSQ AND HRC FILTERS
# ==============================================================================

cat("\n========== STEP 10: Applying quality filters ==========\n)

# For Rsq > 0.3 (recommended minimum)
cat("\nCreating Rsq > 0.3 filtered dataset...\n")
system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_ids_fixed ",
             "--extract output/keepers_rsq3_", DATA_NAME, ".txt ",
             "--make-bed ",
             "--out output/", DATA_NAME, "_", ANCESTRY, "_rsq3"))

# For Rsq > 0.7 (higher quality, recommended for most analyses)
cat("Creating Rsq > 0.7 filtered dataset...\n")
system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_ids_fixed ",
             "--extract output/keepers_rsq7_", DATA_NAME, ".txt ",
             "--make-bed ",
             "--out output/", DATA_NAME, "_", ANCESTRY, "_rsq7"))

# ==============================================================================
# 11. REMOVE DUPLICATE VARIANTS (IF ANY)
# ==============================================================================

cat("\n========== STEP 11: Removing duplicate variants ==========\n")

# Check for duplicates in Rsq3 dataset
bim_rsq3 <- fread(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3.bim"))
duplicate_check <- table(duplicated(bim_rsq3$V2))

cat("Duplicates in Rsq3 dataset:\n")
print(duplicate_check)

if(any(duplicated(bim_rsq3$V2))) {
  # Keep only unique variants
  unique_snps <- bim_rsq3[!duplicated(V2), .(V2)]
  fwrite(unique_snps, "temp/unique_snps.txt", 
         col.names = FALSE, quote = FALSE)
  
  # Remove duplicates
  system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_rsq3 ",
               "--extract temp/unique_snps.txt ",
               "--make-bed ",
               "--out output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup"))
  
  cat("Duplicates removed from Rsq3 dataset\n")
} else {
  cat("No duplicates found\n")
  # Copy file if no duplicates
  file.copy(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3.bim"),
            paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.bim"))
  file.copy(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3.bed"),
            paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.bed"))
  file.copy(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3.fam"),
            paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.fam"))
}

# Repeat for Rsq7 dataset
bim_rsq7 <- fread(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq7.bim"))
if(any(duplicated(bim_rsq7$V2))) {
  unique_snps_rsq7 <- bim_rsq7[!duplicated(V2), .(V2)]
  fwrite(unique_snps_rsq7, "temp/unique_snps_rsq7.txt", 
         col.names = FALSE, quote = FALSE)
  
  system(paste("plink --bfile output/", DATA_NAME, "_", ANCESTRY, "_rsq7 ",
               "--extract temp/unique_snps_rsq7.txt ",
               "--make-bed ",
               "--out output/", DATA_NAME, "_", ANCESTRY, "_rsq7_nodup"))
}

# ==============================================================================
# 12. FINAL QUALITY CONTROL CHECKS
# ==============================================================================

cat("\n========== STEP 12: Final QC checks ==========\n")

# Function to get file info
get_file_info <- function(prefix) {
  bim_file <- paste0(prefix, ".bim")
  fam_file <- paste0(prefix, ".fam")
  
  if(file.exists(bim_file) & file.exists(fam_file)) {
    variants <- system(paste("wc -l", bim_file), intern = TRUE)
    samples <- system(paste("wc -l", fam_file), intern = TRUE)
    return(list(variants = as.numeric(strsplit(variants, " ")[[1]][1]),
                samples = as.numeric(strsplit(samples, " ")[[1]][1])))
  } else {
    return(list(variants = NA, samples = NA))
  }
}

# Check final datasets
cat("\nFinal dataset statistics:\n")
cat("  Rsq > 0.3 (with duplicates):\n")
rsq3_info <- get_file_info(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3"))
cat("    Variants:", rsq3_info$variants, "\n")
cat("    Samples:", rsq3_info$samples, "\n")

cat("\n  Rsq > 0.3 (no duplicates):\n")
rsq3_nodup_info <- get_file_info(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup"))
cat("    Variants:", rsq3_nodup_info$variants, "\n")
cat("    Samples:", rsq3_nodup_info$samples, "\n")

cat("\n  Rsq > 0.7 (no duplicates):\n")
rsq7_nodup_info <- get_file_info(paste0("output/", DATA_NAME, "_", ANCESTRY, "_rsq7_nodup"))
cat("    Variants:", rsq7_nodup_info$variants, "\n")
cat("    Samples:", rsq7_nodup_info$samples, "\n")

# ==============================================================================
# 13. CREATE FINAL RENAMED OUTPUT
# ==============================================================================

cat("\n========== STEP 13: Creating final output files ==========\n)

# For Rsq > 0.3 (more variants, lower quality)
final_name_rsq3 <- paste0(DATA_NAME, "_", ANCESTRY, "_imputed_rsq3")
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.bed ", final_name_rsq3, ".bed", sep = ""))
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.bim ", final_name_rsq3, ".bim", sep = ""))
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq3_nodup.fam ", final_name_rsq3, ".fam", sep = ""))

# For Rsq > 0.7 (fewer variants, higher quality)
final_name_rsq7 <- paste0(DATA_NAME, "_", ANCESTRY, "_imputed_rsq7")
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq7_nodup.bed ", final_name_rsq7, ".bed", sep = ""))
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq7_nodup.bim ", final_name_rsq7, ".bim", sep = ""))
system(paste("cp output/", DATA_NAME, "_", ANCESTRY, "_rsq7_nodup.fam ", final_name_rsq7, ".fam", sep = ""))

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

cat("\n")
cat("========================================\n")
cat("✓ POST-IMPUTATION PIPELINE COMPLETE!\n")
cat("========================================\n")
cat("\nOutput files:\n")
cat("  - High coverage (Rsq>0.3):", final_name_rsq3, "{bed,bim,fam}\n")
cat("  - High quality (Rsq>0.7):", final_name_rsq7, "{bed,bim,fam}\n")
cat("  - Summary table: summary_tables/", DATA_NAME, "_", ANCESTRY, "_imputation_summary.csv\n", sep = "")
cat("\nRecommended Rsq threshold for analysis:\n")
cat("  - GWAS: Rsq > 0.7\n")
cat("  - Fine-mapping: Rsq > 0.9\n")
cat("  - Exploratory: Rsq > 0.3\n")
cat("\n========================================\n")
