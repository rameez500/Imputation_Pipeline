## ============================================================================
## FUNCTION: multidimensional_scaling
## ============================================================================
## Purpose: Identify outliers using multidimensional Euclidean distance
## Method: 
##   1. Calculate quartiles for PC1-3 for each reference population
##   2. Compute distance from population median for each sample
##   3. Threshold = sqrt(Q3² + Q3² + Q3²) + 1.5 × sqrt(IQR² + IQR² + IQR²)
##   4. Classify outliers as distance ≥ threshold
## ============================================================================

multidimensional_scaling <- function(afr_pcs, eur_pcs, sas_pcs, eas_pcs, amr_pcs) {
  
  ## --------------------------------------------------------------------------
  ## AFRICAN (AFR)
  ## --------------------------------------------------------------------------
  
  ## Calculate quartiles for each PC
  afr_quart_pc1 <- as.matrix(quantile(afr_pcs$PC1))
  afr_quart_pc2 <- as.matrix(quantile(afr_pcs$PC2))
  afr_quart_pc3 <- as.matrix(quantile(afr_pcs$PC3))
  
  ## Interquartile ranges
  afr_iqr_pc1 <- afr_quart_pc1[4, 1] - afr_quart_pc1[2, 1]  # Q3 - Q1
  afr_iqr_pc2 <- afr_quart_pc2[4, 1] - afr_quart_pc2[2, 1]
  afr_iqr_pc3 <- afr_quart_pc3[4, 1] - afr_quart_pc3[2, 1]
  
  ## Combine statistics
  afr_pc1_stats <- rbind(afr_iqr_pc1, afr_quart_pc1)
  afr_pc2_stats <- rbind(afr_iqr_pc2, afr_quart_pc2)
  afr_pc3_stats <- rbind(afr_iqr_pc3, afr_quart_pc3)
  
  afr_stats2 <- cbind(afr_pc1_stats, afr_pc2_stats, afr_pc3_stats)
  rownames(afr_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
  colnames(afr_stats2) <- c("PC1", "PC2", "PC3")
  
  ## Calculate Euclidean distance from median
  m_afr_pcs <- afr_pcs[, 1:5]  # Keep FID, IID, PC1, PC2, PC3
  m_afr_pcs$distance <- 999  # Placeholder
  m_afr_pcs$outlier <- 999   # Placeholder
  
  for (i in 1:nrow(m_afr_pcs)) {
    m_afr_pcs[i, 6] <- sqrt(
      ((m_afr_pcs[i, 3] - afr_stats2[4, 1])^2) +
      ((m_afr_pcs[i, 4] - afr_stats2[4, 2])^2) +
      ((m_afr_pcs[i, 5] - afr_stats2[4, 3])^2)
    )
  }
  
  ## Calculate outlier threshold
  threshold_afr <- sqrt(
    (afr_stats2[5, 1]^2) + (afr_stats2[5, 2]^2) + (afr_stats2[5, 3]^2)
  ) + 1.5 * sqrt(
    (afr_stats2[1, 1]^2) + (afr_stats2[1, 2]^2) + (afr_stats2[1, 3]^2)
  )
  
  ## Classify outliers
  for (i in 1:nrow(m_afr_pcs)) {
    m_afr_pcs[i, 7] <- ifelse(m_afr_pcs[i, 6] >= threshold_afr, 1, 0)
  }
  
  
  ## --------------------------------------------------------------------------
  ## EUROPEAN (EUR)
  ## --------------------------------------------------------------------------
  
  eur_quart_pc1 <- as.matrix(quantile(eur_pcs$PC1))
  eur_quart_pc2 <- as.matrix(quantile(eur_pcs$PC2))
  eur_quart_pc3 <- as.matrix(quantile(eur_pcs$PC3))
  
  eur_iqr_pc1 <- eur_quart_pc1[4, 1] - eur_quart_pc1[2, 1]
  eur_iqr_pc2 <- eur_quart_pc2[4, 1] - eur_quart_pc2[2, 1]
  eur_iqr_pc3 <- eur_quart_pc3[4, 1] - eur_quart_pc3[2, 1]
  
  eur_pc1_stats <- rbind(eur_iqr_pc1, eur_quart_pc1)
  eur_pc2_stats <- rbind(eur_iqr_pc2, eur_quart_pc2)
  eur_pc3_stats <- rbind(eur_iqr_pc3, eur_quart_pc3)
  
  eur_stats2 <- cbind(eur_pc1_stats, eur_pc2_stats, eur_pc3_stats)
  rownames(eur_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
  colnames(eur_stats2) <- c("PC1", "PC2", "PC3")
  
  m_eur_pcs <- eur_pcs[, 1:5]
  m_eur_pcs$distance <- 999
  m_eur_pcs$outlier <- 999
  
  for (i in 1:nrow(m_eur_pcs)) {
    m_eur_pcs[i, 6] <- sqrt(
      ((m_eur_pcs[i, 3] - eur_stats2[4, 1])^2) +
      ((m_eur_pcs[i, 4] - eur_stats2[4, 2])^2) +
      ((m_eur_pcs[i, 5] - eur_stats2[4, 3])^2)
    )
  }
  
  threshold_eur <- sqrt(
    (eur_stats2[5, 1]^2) + (eur_stats2[5, 2]^2) + (eur_stats2[5, 3]^2)
  ) + 1.5 * sqrt(
    (eur_stats2[1, 1]^2) + (eur_stats2[1, 2]^2) + (eur_stats2[1, 3]^2)
  )
  
  for (i in 1:nrow(m_eur_pcs)) {
    m_eur_pcs[i, 7] <- ifelse(m_eur_pcs[i, 6] >= threshold_eur, 1, 0)
  }
  
  
  ## --------------------------------------------------------------------------
  ## SOUTH ASIAN (SAS)
  ## --------------------------------------------------------------------------
  
  sas_quart_pc1 <- as.matrix(quantile(sas_pcs$PC1))
  sas_quart_pc2 <- as.matrix(quantile(sas_pcs$PC2))
  sas_quart_pc3 <- as.matrix(quantile(sas_pcs$PC3))
  
  sas_iqr_pc1 <- sas_quart_pc1[4, 1] - sas_quart_pc1[2, 1]
  sas_iqr_pc2 <- sas_quart_pc2[4, 1] - sas_quart_pc2[2, 1]
  sas_iqr_pc3 <- sas_quart_pc3[4, 1] - sas_quart_pc3[2, 1]
  
  sas_pc1_stats <- rbind(sas_iqr_pc1, sas_quart_pc1)
  sas_pc2_stats <- rbind(sas_iqr_pc2, sas_quart_pc2)
  sas_pc3_stats <- rbind(sas_iqr_pc3, sas_quart_pc3)
  
  sas_stats2 <- cbind(sas_pc1_stats, sas_pc2_stats, sas_pc3_stats)
  rownames(sas_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
  colnames(sas_stats2) <- c("PC1", "PC2", "PC3")
  
  m_sas_pcs <- sas_pcs[, 1:5]
  m_sas_pcs$distance <- 999
  m_sas_pcs$outlier <- 999
  
  for (i in 1:nrow(m_sas_pcs)) {
    m_sas_pcs[i, 6] <- sqrt(
      ((m_sas_pcs[i, 3] - sas_stats2[4, 1])^2) +
      ((m_sas_pcs[i, 4] - sas_stats2[4, 2])^2) +
      ((m_sas_pcs[i, 5] - sas_stats2[4, 3])^2)
    )
  }
  
  threshold_sas <- sqrt(
    (sas_stats2[5, 1]^2) + (sas_stats2[5, 2]^2) + (sas_stats2[5, 3]^2)
  ) + 1.5 * sqrt(
    (sas_stats2[1, 1]^2) + (sas_stats2[1, 2]^2) + (sas_stats2[1, 3]^2)
  )
  
  for (i in 1:nrow(m_sas_pcs)) {
    m_sas_pcs[i, 7] <- ifelse(m_sas_pcs[i, 6] >= threshold_sas, 1, 0)
  }
  
  
  ## --------------------------------------------------------------------------
  ## EAST ASIAN (EAS)
  ## --------------------------------------------------------------------------
  
  eas_quart_pc1 <- as.matrix(quantile(eas_pcs$PC1))
  eas_quart_pc2 <- as.matrix(quantile(eas_pcs$PC2))
  eas_quart_pc3 <- as.matrix(quantile(eas_pcs$PC3))
  
  eas_iqr_pc1 <- eas_quart_pc1[4, 1] - eas_quart_pc1[2, 1]
  eas_iqr_pc2 <- eas_quart_pc2[4, 1] - eas_quart_pc2[2, 1]
  eas_iqr_pc3 <- eas_quart_pc3[4, 1] - eas_quart_pc3[2, 1]
  
  eas_pc1_stats <- rbind(eas_iqr_pc1, eas_quart_pc1)
  eas_pc2_stats <- rbind(eas_iqr_pc2, eas_quart_pc2)
  eas_pc3_stats <- rbind(eas_iqr_pc3, eas_quart_pc3)
  
  eas_stats2 <- cbind(eas_pc1_stats, eas_pc2_stats, eas_pc3_stats)
  rownames(eas_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
  colnames(eas_stats2) <- c("PC1", "PC2", "PC3")
  
  m_eas_pcs <- eas_pcs[, 1:5]
  m_eas_pcs$distance <- 999
  m_eas_pcs$outlier <- 999
  
  for (i in 1:nrow(m_eas_pcs)) {
    m_eas_pcs[i, 6] <- sqrt(
      ((m_eas_pcs[i, 3] - eas_stats2[4, 1])^2) +
      ((m_eas_pcs[i, 4] - eas_stats2[4, 2])^2) +
      ((m_eas_pcs[i, 5] - eas_stats2[4, 3])^2)
    )
  }
  
  threshold_eas <- sqrt(
    (eas_stats2[5, 1]^2) + (eas_stats2[5, 2]^2) + (eas_stats2[5, 3]^2)
  ) + 1.5 * sqrt(
    (eas_stats2[1, 1]^2) + (eas_stats2[1, 2]^2) + (eas_stats2[1, 3]^2)
  )
  
  for (i in 1:nrow(m_eas_pcs)) {
    m_eas_pcs[i, 7] <- ifelse(m_eas_pcs[i, 6] >= threshold_eas, 1, 0)
  }
  
  
  ## --------------------------------------------------------------------------
  ## ADMIXED AMERICAN (AMR)
  ## --------------------------------------------------------------------------
  
  amr_quart_pc1 <- as.matrix(quantile(amr_pcs$PC1))
  amr_quart_pc2 <- as.matrix(quantile(amr_pcs$PC2))
  amr_quart_pc3 <- as.matrix(quantile(amr_pcs$PC3))
  
  amr_iqr_pc1 <- amr_quart_pc1[4, 1] - amr_quart_pc1[2, 1]
  amr_iqr_pc2 <- amr_quart_pc2[4, 1] - amr_quart_pc2[2, 1]
  amr_iqr_pc3 <- amr_quart_pc3[4, 1] - amr_quart_pc3[2, 1]
  
  amr_pc1_stats <- rbind(amr_iqr_pc1, amr_quart_pc1)
  amr_pc2_stats <- rbind(amr_iqr_pc2, amr_quart_pc2)
  amr_pc3_stats <- rbind(amr_iqr_pc3, amr_quart_pc3)
  
  amr_stats2 <- cbind(amr_pc1_stats, amr_pc2_stats, amr_pc3_stats)
  rownames(amr_stats2) <- c("iqr", "min", "q1", "median", "q3", "max")
  colnames(amr_stats2) <- c("PC1", "PC2", "PC3")
  
  m_amr_pcs <- amr_pcs[, 1:5]
  m_amr_pcs$distance <- 999
  m_amr_pcs$outlier <- 999
  
  for (i in 1:nrow(m_amr_pcs)) {
    m_amr_pcs[i, 6] <- sqrt(
      ((m_amr_pcs[i, 3] - amr_stats2[4, 1])^2) +
      ((m_amr_pcs[i, 4] - amr_stats2[4, 2])^2) +
      ((m_amr_pcs[i, 5] - amr_stats2[4, 3])^2)
    )
  }
  
  threshold_amr <- sqrt(
    (amr_stats2[5, 1]^2) + (amr_stats2[5, 2]^2) + (amr_stats2[5, 3]^2)
  ) + 1.5 * sqrt(
    (amr_stats2[1, 1]^2) + (amr_stats2[1, 2]^2) + (amr_stats2[1, 3]^2)
  )
  
  for (i in 1:nrow(m_amr_pcs)) {
    m_amr_pcs[i, 7] <- ifelse(m_amr_pcs[i, 6] >= threshold_amr, 1, 0)
  }
  
  
  ## --------------------------------------------------------------------------
  ## RETURN RESULTS
  ## --------------------------------------------------------------------------
  
  multi_scale_list <- list(
    "afr_stats2" = afr_stats2,
    "eur_stats2" = eur_stats2,
    "sas_stats2" = sas_stats2,
    "eas_stats2" = eas_stats2,
    "amr_stats2" = amr_stats2,
    "threshold_afr" = threshold_afr,
    "m_afr_pcs" = m_afr_pcs,
    "threshold_eur" = threshold_eur,
    "m_eur_pcs" = m_eur_pcs,
    "threshold_sas" = threshold_sas,
    "m_sas_pcs" = m_sas_pcs,
    "threshold_eas" = threshold_eas,
    "m_eas_pcs" = m_eas_pcs,
    "threshold_amr" = threshold_amr,
    "m_amr_pcs" = m_amr_pcs
  )
  
  return(multi_scale_list)
}
