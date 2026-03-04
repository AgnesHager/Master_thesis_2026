# Plotting the DEGs to see the expression across the conditions

# Save_plot function


# Plotting the genes (DEGs), individually
dds_ARC <- dds_full[, dds_full$Tissue == "ARC"]
dds_ARC <- estimateSizeFactors(dds_ARC)

dds_LC  <- dds_full[, dds_full$Tissue == "LC"]
dds_LC <- estimateSizeFactors(dds_LC)

dds_VMH <- dds_full[, dds_full$Tissue == "VMH"]
dds_VMH <- estimateSizeFactors(dds_VMH)


# ---------------------------------------
# Plotting the DEGs in ARC, females(change the code to account for only LC, and only VMH) 
# F - KI vs WT on WD, the DEGs are: Chrna6 + Chrnb3

plot_gene_boxplot_ARC_F <- function(gene_id, gene_name) {
  
  # Extract normalized counts from ARC
  plot_df <- as.data.frame(colData(dds_ARC)) %>%
    dplyr::mutate(NormCount = counts(dds_ARC, normalized = TRUE)[gene_id, ]) %>%
    dplyr::filter(Sex == "F") %>%
    dplyr::mutate(
      Diet_Genotype = factor(paste(Genotype, Diet, sep = ", "),
                       levels = c("WT, CD", "KI, CD", "WT, WD", "KI, WD"))
    )
  
  ggplot(plot_df, aes(x = Diet_Genotype, y = NormCount, fill = Diet_Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c("WT, CD" = "grey60", "KI, CD" = "#f4a0c0",
                                 "WT, WD" = "grey30", "KI, WD" = "red")) +
    labs(title = paste(gene_name, "- ARC Females"),
         x = "", y = "Normalized counts") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

#  Usage:
save_plot("Chrna6_ARC_F_KIvsWT_WD.pdf", width = 6, height = 5) # Using my already made code to save my plots
plot_gene_boxplot_ARC_F("ENSMUSG00000031491", "Chrna6")
dev.off()

# Plotting the DEGs in ARC, males (Change the code to account for only LC, and only VMH) 
plot_gene_boxplot_ARC_M <- function(gene_id, gene_name) {
  
  # Extract normalized counts from ARC
  plot_df <- as.data.frame(colData(dds_ARC)) %>%
    dplyr::mutate(NormCount = counts(dds_ARC, normalized = TRUE)[gene_id, ]) %>%
    dplyr::filter(Sex == "M") %>%
    dplyr::mutate(
      Diet_Genotype = factor(paste(Genotype, Diet, sep = ", "),
                       levels = c("WT, CD", "KI, CD", "WT, WD", "KI, WD"))
    )
  
  ggplot(plot_df, aes(x = Diet_Genotype, y = NormCount, fill = Diet_Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c("WT, CD" = "grey60", "KI, CD" = "lightblue",
                                 "WT, WD" = "grey30", "KI, WD" = "darkgreen")) +
    labs(title = paste(gene_name, "- ARC Males"),
         x = "", y = "Normalized counts") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# Usage:
save_plot("Myef2l_ARC_M_KIvsWT_CD.pdf", width = 6, height = 5) # Using my already made code to save my plots
plot_gene_boxplot_ARC_M("ENSMUSG00000049230", "Myef2l")
dev.off()
