##### Step 7: Select Ancestral Groups ####

## ============================================================================
## FUNCTION: select_ancestral
## ============================================================================
## Purpose: Calculate population statistics and thresholds for ancestry assignment
## Input: PC_20_merged - Data frame with PC1-20 and GROUP column (AFR/EUR/etc.)
## Output: List containing means, SDs, and ±2 SD thresholds for each ancestry
## ============================================================================

select_ancestral <- function(PC_20_merged) {
  
  ## --------------------------------------------------------------------------
  ## AFRICAN (AFR)
  ## --------------------------------------------------------------------------
  just_afr <- PC_20_merged[PC_20_merged$GROUP == "AFR", ]
  print(head(just_afr))
  print(dim(just_afr))  # 661 samples in 1000G AFR
  
  afr_m1 <- mean(just_afr$PC1)
  afr_m2 <- mean(just_afr$PC2)
  afr_m3 <- mean(just_afr$PC3)
  afr_sd1 <- sd(just_afr$PC1)
  afr_sd2 <- sd(just_afr$PC2)
  afr_sd3 <- sd(just_afr$PC3)
  
  afr_means <- rbind(afr_m1, afr_m2, afr_m3)
  afr_sds <- rbind(afr_sd1, afr_sd2, afr_sd3)
  afr_stats <- cbind(afr_means, afr_sds)
  colnames(afr_stats) <- c("mean", "stand_dev")
  
  ## --------------------------------------------------------------------------
  ## EUROPEAN (EUR)
  ## --------------------------------------------------------------------------
  just_eur <- PC_20_merged[PC_20_merged$GROUP == "EUR", ]
  print(dim(just_eur))  # 503 samples in 1000G EUR
  
  eur_m1 <- mean(just_eur$PC1)
  eur_m2 <- mean(just_eur$PC2)
  eur_m3 <- mean(just_eur$PC3)
  eur_sd1 <- sd(just_eur$PC1)
  eur_sd2 <- sd(just_eur$PC2)
  eur_sd3 <- sd(just_eur$PC3)
  
  eur_means <- rbind(eur_m1, eur_m2, eur_m3)
  eur_sds <- rbind(eur_sd1, eur_sd2, eur_sd3)
  eur_stats <- cbind(eur_means, eur_sds)
  colnames(eur_stats) <- c("mean", "stand_dev")
  
  ## --------------------------------------------------------------------------
  ## SOUTH ASIAN (SAS)
  ## --------------------------------------------------------------------------
  just_sas <- PC_20_merged[PC_20_merged$GROUP == "SAS", ]
  print(dim(just_sas))  # 489 samples in 1000G SAS
  
  sas_m1 <- mean(just_sas$PC1)
  sas_m2 <- mean(just_sas$PC2)
  sas_m3 <- mean(just_sas$PC3)
  sas_sd1 <- sd(just_sas$PC1)
  sas_sd2 <- sd(just_sas$PC2)
  sas_sd3 <- sd(just_sas$PC3)
  
  sas_means <- rbind(sas_m1, sas_m2, sas_m3)
  sas_sds <- rbind(sas_sd1, sas_sd2, sas_sd3)
  sas_stats <- cbind(sas_means, sas_sds)
  colnames(sas_stats) <- c("mean", "stand_dev")
  
  ## --------------------------------------------------------------------------
  ## EAST ASIAN (EAS)
  ## --------------------------------------------------------------------------
  just_eas <- PC_20_merged[PC_20_merged$GROUP == "EAS", ]
  print(dim(just_eas))  # 504 samples in 1000G EAS
  
  eas_m1 <- mean(just_eas$PC1)
  eas_m2 <- mean(just_eas$PC2)
  eas_m3 <- mean(just_eas$PC3)
  eas_sd1 <- sd(just_eas$PC1)
  eas_sd2 <- sd(just_eas$PC2)
  eas_sd3 <- sd(just_eas$PC3)
  
  eas_means <- rbind(eas_m1, eas_m2, eas_m3)
  eas_sds <- rbind(eas_sd1, eas_sd2, eas_sd3)
  eas_stats <- cbind(eas_means, eas_sds)
  colnames(eas_stats) <- c("mean", "stand_dev")
  
  ## --------------------------------------------------------------------------
  ## ADMIXED AMERICAN (AMR)
  ## --------------------------------------------------------------------------
  just_amr <- PC_20_merged[PC_20_merged$GROUP == "AMR", ]
  print(dim(just_amr))  # 347 samples in 1000G AMR
  
  amr_m1 <- mean(just_amr$PC1)
  amr_m2 <- mean(just_amr$PC2)
  amr_m3 <- mean(just_amr$PC3)
  amr_sd1 <- sd(just_amr$PC1)
  amr_sd2 <- sd(just_amr$PC2)
  amr_sd3 <- sd(just_amr$PC3)
  
  amr_means <- rbind(amr_m1, amr_m2, amr_m3)
  amr_sds <- rbind(amr_sd1, amr_sd2, amr_sd3)
  amr_stats <- cbind(amr_means, amr_sds)
  colnames(amr_stats) <- c("mean", "stand_dev")
  
  ## --------------------------------------------------------------------------
  ## COMPUTE THRESHOLDS (±2 STANDARD DEVIATIONS)
  ## --------------------------------------------------------------------------
  
  ## AFR thresholds
  afr_UP_thresh_1 <- afr_stats[1, 1] + (2 * afr_stats[1, 2])
  afr_LO_thresh_1 <- afr_stats[1, 1] - (2 * afr_stats[1, 2])
  afr_UP_thresh_2 <- afr_stats[2, 1] + (2 * afr_stats[2, 2])
  afr_LO_thresh_2 <- afr_stats[2, 1] - (2 * afr_stats[2, 2])
  afr_UP_thresh_3 <- afr_stats[3, 1] + (2 * afr_stats[3, 2])
  afr_LO_thresh_3 <- afr_stats[3, 1] - (2 * afr_stats[3, 2])
  
  afr_thresholds <- rbind(afr_UP_thresh_1, afr_LO_thresh_1, 
                          afr_UP_thresh_2, afr_LO_thresh_2,
                          afr_UP_thresh_3, afr_LO_thresh_3)
  
  ## EUR thresholds
  eur_UP_thresh_1 <- eur_stats[1, 1] + (2 * eur_stats[1, 2])
  eur_LO_thresh_1 <- eur_stats[1, 1] - (2 * eur_stats[1, 2])
  eur_UP_thresh_2 <- eur_stats[2, 1] + (2 * eur_stats[2, 2])
  eur_LO_thresh_2 <- eur_stats[2, 1] - (2 * eur_stats[2, 2])
  eur_UP_thresh_3 <- eur_stats[3, 1] + (2 * eur_stats[3, 2])
  eur_LO_thresh_3 <- eur_stats[3, 1] - (2 * eur_stats[3, 2])
  
  eur_thresholds <- rbind(eur_UP_thresh_1, eur_LO_thresh_1,
                          eur_UP_thresh_2, eur_LO_thresh_2,
                          eur_UP_thresh_3, eur_LO_thresh_3)
  
  ## SAS thresholds
  sas_UP_thresh_1 <- sas_stats[1, 1] + (2 * sas_stats[1, 2])
  sas_LO_thresh_1 <- sas_stats[1, 1] - (2 * sas_stats[1, 2])
  sas_UP_thresh_2 <- sas_stats[2, 1] + (2 * sas_stats[2, 2])
  sas_LO_thresh_2 <- sas_stats[2, 1] - (2 * sas_stats[2, 2])
  sas_UP_thresh_3 <- sas_stats[3, 1] + (2 * sas_stats[3, 2])
  sas_LO_thresh_3 <- sas_stats[3, 1] - (2 * sas_stats[3, 2])
  
  sas_thresholds <- rbind(sas_UP_thresh_1, sas_LO_thresh_1,
                          sas_UP_thresh_2, sas_LO_thresh_2,
                          sas_UP_thresh_3, sas_LO_thresh_3)
  
  ## EAS thresholds
  eas_UP_thresh_1 <- eas_stats[1, 1] + (2 * eas_stats[1, 2])
  eas_LO_thresh_1 <- eas_stats[1, 1] - (2 * eas_stats[1, 2])
  eas_UP_thresh_2 <- eas_stats[2, 1] + (2 * eas_stats[2, 2])
  eas_LO_thresh_2 <- eas_stats[2, 1] - (2 * eas_stats[2, 2])
  eas_UP_thresh_3 <- eas_stats[3, 1] + (2 * eas_stats[3, 2])
  eas_LO_thresh_3 <- eas_stats[3, 1] - (2 * eas_stats[3, 2])
  
  eas_thresholds <- rbind(eas_UP_thresh_1, eas_LO_thresh_1,
                          eas_UP_thresh_2, eas_LO_thresh_2,
                          eas_UP_thresh_3, eas_LO_thresh_3)
  
  ## AMR thresholds
  amr_UP_thresh_1 <- amr_stats[1, 1] + (2 * amr_stats[1, 2])
  amr_LO_thresh_1 <- amr_stats[1, 1] - (2 * amr_stats[1, 2])
  amr_UP_thresh_2 <- amr_stats[2, 1] + (2 * amr_stats[2, 2])
  amr_LO_thresh_2 <- amr_stats[2, 1] - (2 * amr_stats[2, 2])
  amr_UP_thresh_3 <- amr_stats[3, 1] + (2 * amr_stats[3, 2])
  amr_LO_thresh_3 <- amr_stats[3, 1] - (2 * amr_stats[3, 2])
  
  amr_thresholds <- rbind(amr_UP_thresh_1, amr_LO_thresh_1,
                          amr_UP_thresh_2, amr_LO_thresh_2,
                          amr_UP_thresh_3, amr_LO_thresh_3)
  
  ## --------------------------------------------------------------------------
  ## RETURN RESULTS
  ## --------------------------------------------------------------------------
  ancestral_list <- list(
    "afr_stats" = afr_stats,
    "eur_stats" = eur_stats,
    "sas_stats" = sas_stats,
    "eas_stats" = eas_stats,
    "amr_stats" = amr_stats,
    "afr_thresholds" = afr_thresholds,
    "eur_thresholds" = eur_thresholds,
    "sas_thresholds" = sas_thresholds,
    "eas_thresholds" = eas_thresholds,
    "amr_thresholds" = amr_thresholds
  )
  
  return(ancestral_list)
}
