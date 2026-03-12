# My DEG path, since I put them in another folder
deg_path <- "C:/Users/agnes/OneDrive/BIOMEDISIN MASTER/GOI_2025_Featurecounts/Master_project_Recode_05.02/DEGs_13.02/"
# Load DEGs
degs_ARC_M_WD <- read.csv(paste0(deg_path, "DEGs_ARC_M_WD_KIvsWT_sig_clean.csv"))

# Extract WD male samples
meta_ARC     <- as.data.frame(colData(dds_ARC_filt))
WD_M_samples <- rownames(meta_ARC[meta_ARC$Diet == "WD" & meta_ARC$Sex == "M", ])

# Subset dds to WD males only and apply VST
dds_ARC_M_WD <- dds_ARC_filt[, WD_M_samples]
vsd_ARC_M_WD <- vst(dds_ARC_M_WD, blind = FALSE)

# Extract expression for DEGs only
deg_genes          <- degs_ARC_M_WD$ensembl_gene_id
expr_mat           <- assay(vsd_ARC_M_WD)[deg_genes, ]
rownames(expr_mat) <- degs_ARC_M_WD$mgi_symbol

# Annotation
annotation_col <- data.frame(Genotype = meta_ARC[WD_M_samples, "Genotype"])
rownames(annotation_col) <- WD_M_samples
ann_colors <- list(Genotype = c(WT = "darkorange", KI = "darkgreen"))

# Plotting the heatmap
pheatmap(
  expr_mat,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 10,
  color             = colorRampPalette(c("blue", "white", "red"))(100),
  main              = "Differentially expressed genes: ARC Males WD (KI vs WT)"
)
