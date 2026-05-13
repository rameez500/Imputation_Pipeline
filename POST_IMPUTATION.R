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
## change: "data" to your dataset and "ancestry" to AFR or EUR (or admixed)

## ============================================================================
## LOAD LIBRARIES
## ============================================================================

library(data.table)
library(dplyr)
library(tidyr)

## ============================================================================
## STEP 0: DOWNLOAD DATA FROM MICHIGAN SERVER
## ============================================================================

## Imputation password: [YOUR PASSWORD]

## ============================================================================
## STEP 1: UNZIP EACH CHROMOSOME FOLDER
## ============================================================================

chrm <- 1:22
x <- NULL
for(i in 1:22) {
  x <- rbind(x, paste("7z x chr_", chrm[i], ".zip -p'password'", sep = ""))
}
write.table(x, file = "unzip_data.sh", row.names = FALSE, col.names = FALSE, quote = FALSE)
system("dos2unix unzip_data.sh")
system("chmod +x unzip_data.sh")
system("nohup ./unzip_data.sh > unzip_data.log &")

## ============================================================================
## STEP 2: UNZIP EACH .INFO FILE
## ============================================================================

x2 <- NULL
for(i in 1:22) {
  x2 <- rbind(x2, paste("gunzip chr", chrm[i], ".info.gz", sep = ""))
}
write.table(x2, file = "unzip_info_data.sh", row.names = FALSE, col.names = FALSE, quote = FALSE)
system("dos2unix unzip_info_data.sh")
system("chmod +x unzip_info_data.sh")
system("nohup ./unzip_info_data.sh > unzip_info_data.log &")

## ============================================================================
## STEP 3: PREPARE IMPUTED FILES FOR MERGING
## ============================================================================

## Read in all chromosome info files
chr1 <- fread("chr1.info", header = TRUE, sep = "\t", na.strings = "-")
chr2 <- fread("chr2.info", header = TRUE, sep = "\t", na.strings = "-")
# ... repeat for chr3-22

## Bind all together
allchrm <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10,
                 chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19,
                 chr20, chr21, chr22)

nrow(allchrm)  # Total number imputed
table(allchrm$Genotyped)  # Genotyped vs Imputed count
mean(allchrm$Rsq)  # Mean imputation quality
fivenum(allchrm$Rsq)  # Five-number summary of Rsq

## Create filters for different Rsq thresholds
allchrm_r2.3 <- allchrm[allchrm$Rsq > 0.3, ]
allchrm_r2.5 <- allchrm[allchrm$Rsq > 0.5, ]
allchrm_r2.7 <- allchrm[allchrm$Rsq > 0.7, ]
allchrm_r2.9 <- allchrm[allchrm$Rsq > 0.9, ]

## Save keep lists
fwrite(allchrm_r2.3, "allchrm_keep_rsq3_data.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
fwrite(allchrm_r2.5, "allchrm_keep_rsq5_data.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
fwrite(allchrm_r2.7, "allchrm_keep_rsq7_data.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
fwrite(allchrm_r2.9, "allchrm_keep_rsq9_data.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Create chr:pos identifiers for matching with HRC
library(tidyr)
allchrm_r2.3_chrmpos <- separate(allchrm_r2.3, SNP, c("chr", "pos"), sep = ":", extra = 'drop')
chrmpos_r2.3 <- paste(allchrm_r2.3_chrmpos$chr, allchrm_r2.3_chrmpos$pos, sep = ":")

## Create update-name file for BIM file
## Converts from "chr:pos:A1:A2" to "chr:pos" format
allchrm$chrpos <- substr(allchrm$SNP, 1, nchar(allchrm$SNP) - 4)  # Remove last 4 chars (:A1:A2)
update_names <- allchrm[, c("SNP", "chrpos")]
fwrite(update_names, "update_names.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Filter to HRC biallelic variants
list_HRC <- fread("/projects/bga_lab/DATA_REPOSITORIES/IMP_PIPELINE_FILES/HRC_FILES/HRC_reference_biallelic_snps_bim_format_CHRPOS.txt", header = FALSE, sep = " ")

keepers_rsq3 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.3, ]
fwrite(keepers_rsq3[, 2], "keeper_markers_rsq3_data.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

keepers_rsq7 <- list_HRC[list_HRC$V2 %in% chrmpos_r2.7, ]
fwrite(keepers_rsq7[, 2], "keeper_markers_rsq7_data.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

## Create summary table
all_imputed <- as.data.frame(table(allchrm$Genotyped))
names(all_imputed) <- c("QC", "SNPs")

rsq <- c(length(chrmpos_r2.3), length(chrmpos_r2.5), length(chrmpos_r2.7), length(chrmpos_r2.9))
rsq_label <- c("Imputation Quality: R2 > .3", "Imputation Quality: R2 > .5", 
               "Imputation Quality: R2 > .7", "Imputation Quality: R2 > .9")
rsq2_table <- as.data.frame(cbind(rsq_label, rsq))
names(rsq2_table) <- c("QC", "SNPs")

biallelic <- c(nrow(keepers_rsq3), nrow(keepers_rsq5), nrow(keepers_rsq7), nrow(keepers_rsq9))
biallelic_label <- c("Biallelic: R2 > .3", "Biallelic: R2 > .5", 
                     "Biallelic: R2 > .7", "Biallelic: R2 > .9")
biallelic_table <- as.data.frame(cbind(biallelic_label, biallelic))
names(biallelic_table) <- c("QC", "SNPs")

imputed_table <- rbind(all_imputed, rsq2_table, biallelic_table)
write.csv(imputed_table, "data_imputed_ancestry_summary_table.csv", row.names = FALSE, quote = TRUE)

## ============================================================================
## STEP 4: MERGE IMPUTED FILES
## ============================================================================

## Convert VCF to PLINK
chrm <- 1:22
for(i in chrm) {
  cmd <- paste0("plink --vcf chr", i, ".dose.vcf.gz --const-fid --make-bed --out data_ancestry_imputed_chr", i)
  system(cmd)
}

## Create merge list
for(i in chrm) {
  write(paste0("data_ancestry_imputed_chr", i), file = "chrmerge_list.txt", append = TRUE)
}

## Merge all chromosomes
system("nohup plink --bfile data_ancestry_imputed_chr1 --merge-list chrmerge_list.txt --make-bed --out data_ancestry_imputed_allchrm --threads 20 &")

## Update variant names to match HRC format
system("plink --bfile data_ancestry_imputed_allchrm --update-name update_names.txt --make-bed --out allchrm_data_ancestry_chrpos")

## ============================================================================
## STEP 5: UPDATE PARTICIPANT IDs TO MATCH ORIGINAL IDs
## ============================================================================

## Read original FAM file
original_fam <- fread("filepath-to-original-data/data.fam")

## Read imputed FAM file
imp_fam <- fread("allchrm_data_ancestry_chrpos.fam", header = FALSE, sep = " ")

## Split combined IDs (format: FID_IID)
imp_fam$newV2 <- gsub("_", ",", imp_fam$V2)
imp_fam_new <- imp_fam %>% separate(col = newV2, into = c("FID", "IID"), sep = ",")

## Write new FAM file
write.table(imp_fam_new[, c(7:8, 3:6)], "imp_fam_IDadjust.fam", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Update IDs in PLINK file
system("plink --bfile allchrm_data_ancestry_chrpos --fam imp_fam_IDadjust.fam --make-bed --out allchrm_data_ancestry_chrpos_ids")

## ============================================================================
## STEP 6: SCREEN MARKERS
## ============================================================================

## Keep well-imputed biallelic markers
system("plink --bfile allchrm_data_ancestry_chrpos_ids --extract keeper_markers_rsq3_data.txt --make-bed --out data_ancestry_imputed_rsq3")
system("plink --bfile allchrm_data_ancestry_chrpos_ids --extract keeper_markers_rsq7_data.txt --make-bed --out data_ancestry_imputed_rsq7")

## Check for duplicate markers
bim <- fread("data_ancestry_imputed_rsq3.bim")
table(duplicated(bim$V2))

## Keep only unique variants
keepers <- bim[!(duplicated(bim$V2) | duplicated(bim$V2, fromLast = TRUE)), ]
dups <- bim[(duplicated(bim$V2) | duplicated(bim$V2, fromLast = TRUE)), ]

fwrite(keepers[, 2], "unique_markers_rsq3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

## Remove duplicates
system("plink --bfile data_ancestry_imputed_rsq3 --extract unique_markers_rsq3.txt --make-bed --out data_ancestry_imputed_rsq3_nodup")

## ============================================================================
## FINAL FILE
## ============================================================================

## data_ancestry_imputed_rsq3_nodup.{bed,bim,fam} is ready for analysis
system("plink --bfile data_ancestry_imputed_rsq3 --make-bed --out FINAL_NAME_TBD")
