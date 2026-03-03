# QC_analysis.R
# Quality control analysis for Slc7a10 KI RNA-seq data

# LIBRARIES
library(DESeq2)
library(RNAseqQC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(biomaRt)
library(pheatmap)

# Helper function to save plots as PDF with Arial font
save_plot <- function(filename, width = 8, height = 6) {
  cairo_pdf(filename, width = width, height = height, family = "Arial")
}

# Colour schemes
female_genotype_colors <- c("KI" = "#E41A1C", "WT" = "#E69F00")
male_genotype_colors   <- c("KI" = "#4DAF4A", "WT" = "#A6CEE3")
sex_colors             <- c("F" = "#984EA3", "M" = "#F0E442")
diet_colors            <- c("CD" = "#3C5488", "WD" = "#E64B35")
tissue_colors          <- c("ARC" = "#F8766D", "LC" = "#00BA38", "VMH" = "#619CFF")
pt_size <- 4.5

pca_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    axis.text  = element_text(size = 6.5),
    axis.title = element_text(size = 7)
  )


# 1) TECHNICAL ANNOTATION PCA (all 150 samples)
# Defining sample categories
control_samples <- c(
  "21_arc_S242",
  "21_arc_test100ng_S286",
  "21_arc_test35ng_S289",
  "21_arc_test65ng_S290",
  "24_arc_S244",
  "24_arc_test100ng_S288",
  "24_arc_test35ng_S292",
  "24_arc_test65ng_S291",
  "24_arc_testFrag_S287"
)

true_low_input_samples <- c(
  "13_arc_S239",
  "29_arc_S247",
  "40_arc_S256",
  "55_arc_S263",
  "67_arc_S268",
  "13_vmh_S143"
)

large_library_samples <- c(
  "24_arc_testFrag_S287",
  "11_vmh_S183",
  "45_vmh_S163"
)

low_RIN_samples <- c(
  "70_arc_S270",
  "39_vmh_S159",
  "31_arc_S249",
  "66_vmh_S171"
)

# Build QC annotation table
qc_annotation <- data.frame(
  Sample = colnames(dds_full),
  Sample_Category = "Standard",
  stringsAsFactors = FALSE
)
rownames(qc_annotation) <- qc_annotation$Sample

qc_annotation$Sample_Category[qc_annotation$Sample %in% control_samples]        <- "Control"
qc_annotation$Sample_Category[qc_annotation$Sample %in% true_low_input_samples] <- "True_Low_Input"
qc_annotation$Sample_Category[qc_annotation$Sample %in% large_library_samples]  <- "Large_Library"
qc_annotation$Sample_Category[qc_annotation$Sample %in% low_RIN_samples]        <- "Low_RIN"

colData(dds_bio)$Sample_Category <- qc_annotation[colnames(dds_bio), "Sample_Category"]

# VST and PCA with technical annotations
vsd_bio <- vst(dds_bio, blind = TRUE)

pca_data <- plotPCA(vsd_bio, intgroup = "Tissue", returnData = TRUE)
pca_data$Sample_Category <- qc_annotation[rownames(pca_data), "Sample_Category"]
percentVar <- round(100 * attr(pca_data, "percentVar"))

# PCA plot - colored by Tissue, shaped by Sample_Category
# No outlier removal: flagged samples clustered correctly within
# tissue groups, indicating biological rather than technical variation
save_plot("PCA_QC_technical_annotation.pdf")
ggplot(pca_data, aes(x = PC1, y = PC2, color = Tissue, shape = Sample_Category)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_shape_manual(
    values = c(
      "Standard"       = 16,
      "Control"        = 17,
      "True_Low_Input" = 15,
      "Large_Library"  = 18,
      "Low_RIN"        = 3
    ),
    name = "Sample category"
  ) +
  theme_bw() +
  ggtitle("PCA of biological samples with technical annotations") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
dev.off()

# 2) PRE-DESeq2 QC ON CLEANED 140 SAMPLES
# Library sizes
save_plot("library_sizes.pdf")
libsize <- colSums(counts(dds_full))
boxplot(libsize, las = 2, main = "Library sizes", ylab = "Total counts")
dev.off()

# Raw count distribution - global
log_counts <- log10(as.matrix(raw_counts) + 1)

save_plot("raw_counts_global_hist.pdf")
ggplot(data.frame(x = as.vector(log_counts)), aes(x = x)) +
  geom_histogram(bins = 200) +
  xlab("log10(raw counts + 1)") +
  ylab("Frequency") +
  ggtitle("Global raw count distribution")
dev.off()

# Raw count distribution - per sample
log_counts_long <- tidyr::pivot_longer(
  as.data.frame(log_counts),
  cols = everything(),
  names_to = "Sample",
  values_to = "log_counts"
)

save_plot("per_sample_count_distribution.pdf", width = 20, height = 6)
ggplot(log_counts_long, aes(x = Sample, y = log_counts)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample") +
  ylab("log10(counts + 1)") +
  ggtitle("Per-sample count distribution")
dev.off()

# Total counts, library complexity, gene detection
save_plot("total_counts.pdf")
plot_total_counts(dds_full)
dev.off()

save_plot("library_complexity.pdf")
plot_library_complexity(dds_full)
dev.off()

save_plot("gene_detection.pdf")
plot_gene_detection(dds_full)
dev.off()

# 3) GENE BIOTYPE QC (per tissue)
dds_biotype <- dds_full

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes <- rownames(dds_biotype)

biotypes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = genes,
  mart       = ensembl
)

rowData(dds_biotype)$gene_biotype <- biotypes$gene_biotype[
  match(rownames(dds_biotype), biotypes$ensembl_gene_id)
]

dds_arc_biotype <- dds_biotype[, dds_biotype$Tissue == "ARC"]
dds_lc_biotype  <- dds_biotype[, dds_biotype$Tissue == "LC"]
dds_vmh_biotype <- dds_biotype[, dds_biotype$Tissue == "VMH"]

save_plot("bio_types_ARC.pdf", width = 12, height = 6)
p <- plot_biotypes(dds_arc_biotype)
p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

save_plot("bio_types_LC.pdf", width = 12, height = 6)
p <- plot_biotypes(dds_lc_biotype)
p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

save_plot("bio_types_VMH.pdf", width = 12, height = 6)
p <- plot_biotypes(dds_vmh_biotype)
p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

# 4) MEAN-SD PLOTS (before and after VST)
save_plot("mean_SD_plot_before_vst.pdf")
mean_sd_plot(dds_full)
dev.off()

vsd_full <- vst(dds_full, blind = TRUE)

save_plot("mean_SD_plot_after_vst.pdf")
mean_sd_plot(vsd_full)
dev.off()

# 5) SUBSET VSD OBJECTS BY TISSUE AND SEX
# By tissue
vsd_arc <- vst(dds_arc_biotype, blind = TRUE)
vsd_lc  <- vst(dds_lc_biotype,  blind = TRUE)
vsd_vmh <- vst(dds_vmh_biotype, blind = TRUE)

# By tissue and sex
vsd_arc_f <- vsd_arc[, vsd_arc$Sex == "F"]
vsd_arc_m <- vsd_arc[, vsd_arc$Sex == "M"]
vsd_lc_f  <- vsd_lc[,  vsd_lc$Sex  == "F"]
vsd_lc_m  <- vsd_lc[,  vsd_lc$Sex  == "M"]
vsd_vmh_f <- vsd_vmh[, vsd_vmh$Sex == "F"]
vsd_vmh_m <- vsd_vmh[, vsd_vmh$Sex == "M"]

# By diet
vsd_CD <- vsd_full[, vsd_full$Diet == "CD"]
vsd_WD <- vsd_full[, vsd_full$Diet == "WD"]

# By tissue, sex and diet - females
vsd_arc_f_CD <- vsd_CD[, vsd_CD$Tissue == "ARC" & vsd_CD$Sex == "F"]
vsd_arc_f_WD <- vsd_WD[, vsd_WD$Tissue == "ARC" & vsd_WD$Sex == "F"]
vsd_lc_f_CD  <- vsd_CD[, vsd_CD$Tissue == "LC"  & vsd_CD$Sex == "F"]
vsd_lc_f_WD  <- vsd_WD[, vsd_WD$Tissue == "LC"  & vsd_WD$Sex == "F"]
vsd_vmh_f_CD <- vsd_CD[, vsd_CD$Tissue == "VMH" & vsd_CD$Sex == "F"]
vsd_vmh_f_WD <- vsd_WD[, vsd_WD$Tissue == "VMH" & vsd_WD$Sex == "F"]

# By tissue, sex and diet - males
vsd_arc_m_CD <- vsd_CD[, vsd_CD$Tissue == "ARC" & vsd_CD$Sex == "M"]
vsd_arc_m_WD <- vsd_WD[, vsd_WD$Tissue == "ARC" & vsd_WD$Sex == "M"]
vsd_lc_m_CD  <- vsd_CD[, vsd_CD$Tissue == "LC"  & vsd_CD$Sex == "M"]
vsd_lc_m_WD  <- vsd_WD[, vsd_WD$Tissue == "LC"  & vsd_WD$Sex == "M"]
vsd_vmh_m_CD <- vsd_CD[, vsd_CD$Tissue == "VMH" & vsd_CD$Sex == "M"]
vsd_vmh_m_WD <- vsd_WD[, vsd_WD$Tissue == "VMH" & vsd_WD$Sex == "M"]

# -------------------------------------------------------------------
# 6) PCA DATA OBJECTS

# Global
pca_global <- plotPCA(vsd_full, intgroup = "Tissue", returnData = TRUE)
percentVar_global <- round(100 * attr(pca_global, "percentVar"))

# By tissue and sex
pca_arc_f <- plotPCA(vsd_arc_f, intgroup = c("Genotype", "Diet"), returnData = TRUE)
pca_arc_m <- plotPCA(vsd_arc_m, intgroup = c("Genotype", "Diet"), returnData = TRUE)
pca_lc_f  <- plotPCA(vsd_lc_f,  intgroup = c("Genotype", "Diet"), returnData = TRUE)
pca_lc_m  <- plotPCA(vsd_lc_m,  intgroup = c("Genotype", "Diet"), returnData = TRUE)
pca_vmh_f <- plotPCA(vsd_vmh_f, intgroup = c("Genotype", "Diet"), returnData = TRUE)
pca_vmh_m <- plotPCA(vsd_vmh_m, intgroup = c("Genotype", "Diet"), returnData = TRUE)

percentVar_arc_f <- round(100 * attr(pca_arc_f, "percentVar"))
percentVar_arc_m <- round(100 * attr(pca_arc_m, "percentVar"))
percentVar_lc_f  <- round(100 * attr(pca_lc_f,  "percentVar"))
percentVar_lc_m  <- round(100 * attr(pca_lc_m,  "percentVar"))
percentVar_vmh_f <- round(100 * attr(pca_vmh_f, "percentVar"))
percentVar_vmh_m <- round(100 * attr(pca_vmh_m, "percentVar"))

# By sex (for sex comparison plots)
pca_arc_sex <- plotPCA(vsd_arc, intgroup = "Sex", returnData = TRUE)
pca_lc_sex  <- plotPCA(vsd_lc,  intgroup = "Sex", returnData = TRUE)
pca_vmh_sex <- plotPCA(vsd_vmh, intgroup = "Sex", returnData = TRUE)

percentVar_arc_sex <- round(100 * attr(pca_arc_sex, "percentVar"))
percentVar_lc_sex  <- round(100 * attr(pca_lc_sex,  "percentVar"))
percentVar_vmh_sex <- round(100 * attr(pca_vmh_sex, "percentVar"))

# Diet-split by tissue, sex and genotype
pca_arc_f_CD <- plotPCA(vsd_arc_f_CD, intgroup = "Genotype", returnData = TRUE)
pca_arc_f_WD <- plotPCA(vsd_arc_f_WD, intgroup = "Genotype", returnData = TRUE)
pca_lc_f_CD  <- plotPCA(vsd_lc_f_CD,  intgroup = "Genotype", returnData = TRUE)
pca_lc_f_WD  <- plotPCA(vsd_lc_f_WD,  intgroup = "Genotype", returnData = TRUE)
pca_vmh_f_CD <- plotPCA(vsd_vmh_f_CD, intgroup = "Genotype", returnData = TRUE)
pca_vmh_f_WD <- plotPCA(vsd_vmh_f_WD, intgroup = "Genotype", returnData = TRUE)
pca_arc_m_CD <- plotPCA(vsd_arc_m_CD, intgroup = "Genotype", returnData = TRUE)
pca_arc_m_WD <- plotPCA(vsd_arc_m_WD, intgroup = "Genotype", returnData = TRUE)
pca_lc_m_CD  <- plotPCA(vsd_lc_m_CD,  intgroup = "Genotype", returnData = TRUE)
pca_lc_m_WD  <- plotPCA(vsd_lc_m_WD,  intgroup = "Genotype", returnData = TRUE)
pca_vmh_m_CD <- plotPCA(vsd_vmh_m_CD, intgroup = "Genotype", returnData = TRUE)
pca_vmh_m_WD <- plotPCA(vsd_vmh_m_WD, intgroup = "Genotype", returnData = TRUE)

# 7) SHARED AXIS LIMITS
# Females (genotype + diet)
all_pc1_f <- c(pca_arc_f$PC1, pca_lc_f$PC1, pca_vmh_f$PC1)
all_pc2_f <- c(pca_arc_f$PC2, pca_lc_f$PC2, pca_vmh_f$PC2)
pc1_lim_f <- range(all_pc1_f) * 1.1
pc2_lim_f <- range(all_pc2_f) * 1.1

# Males (genotype + diet)
all_pc1_m <- c(pca_arc_m$PC1, pca_lc_m$PC1, pca_vmh_m$PC1)
all_pc2_m <- c(pca_arc_m$PC2, pca_lc_m$PC2, pca_vmh_m$PC2)
pc1_lim_m <- range(all_pc1_m) * 1.1
pc2_lim_m <- range(all_pc2_m) * 1.1

# Sex
all_pc1_sex <- c(pca_arc_sex$PC1, pca_lc_sex$PC1, pca_vmh_sex$PC1)
all_pc2_sex <- c(pca_arc_sex$PC2, pca_lc_sex$PC2, pca_vmh_sex$PC2)
pc1_lim_sex <- range(all_pc1_sex) * 1.1
pc2_lim_sex <- range(all_pc2_sex) * 1.1

# Diet-split females
all_pc1_f_diet <- c(pca_arc_f_CD$PC1, pca_arc_f_WD$PC1,
                    pca_lc_f_CD$PC1,  pca_lc_f_WD$PC1,
                    pca_vmh_f_CD$PC1, pca_vmh_f_WD$PC1)
all_pc2_f_diet <- c(pca_arc_f_CD$PC2, pca_arc_f_WD$PC2,
                    pca_lc_f_CD$PC2,  pca_lc_f_WD$PC2,
                    pca_vmh_f_CD$PC2, pca_vmh_f_WD$PC2)
pc1_lim_f_diet <- range(all_pc1_f_diet) * 1.1
pc2_lim_f_diet <- range(all_pc2_f_diet) * 1.1

# Diet-split males
all_pc1_m_diet <- c(pca_arc_m_CD$PC1, pca_arc_m_WD$PC1,
                    pca_lc_m_CD$PC1,  pca_lc_m_WD$PC1,
                    pca_vmh_m_CD$PC1, pca_vmh_m_WD$PC1)
all_pc2_m_diet <- c(pca_arc_m_CD$PC2, pca_arc_m_WD$PC2,
                    pca_lc_m_CD$PC2,  pca_lc_m_WD$PC2,
                    pca_vmh_m_CD$PC2, pca_vmh_m_WD$PC2)
pc1_lim_m_diet <- range(all_pc1_m_diet) * 1.1
pc2_lim_m_diet <- range(all_pc2_m_diet) * 1.1

# 8) PCA PLOTS
# Global PCA, showing all tissues (ARC, LC, and VMH)
p_global_ellipse <- ggplot(pca_global, aes(PC1, PC2, color = Tissue)) +
  geom_point(size = pt_size, alpha = 0.9) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  scale_color_manual(values = tissue_colors) +
  xlab(paste0("PC1: ", percentVar_global[1], "%")) +
  ylab(paste0("PC2: ", percentVar_global[2], "%")) +
  ggtitle("Global PCA – All tissues") +
  pca_theme +
  coord_fixed()
save_plot("PCA_Global_Tissue_Ellipse.pdf", width = 6, height = 5)
print(p_global_ellipse)
dev.off()

# Females: Genotype + Diet
p_arc_f_gd <- ggplot(pca_arc_f, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = female_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_arc_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_arc_f[2], "%")) +
  ggtitle("ARC – Females") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_ARC_Female_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_arc_f_gd)
dev.off()

p_lc_f_gd <- ggplot(pca_lc_f, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = female_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_lc_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_lc_f[2], "%")) +
  ggtitle("LC – Females") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_LC_Female_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_lc_f_gd)
dev.off()

p_vmh_f_gd <- ggplot(pca_vmh_f, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = female_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_vmh_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_vmh_f[2], "%")) +
  ggtitle("VMH – Females") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_VMH_Female_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_vmh_f_gd)
dev.off()

# Males: Genotype + Diet
p_arc_m_gd <- ggplot(pca_arc_m, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = male_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_arc_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_arc_m[2], "%")) +
  ggtitle("ARC – Males") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_ARC_Male_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_arc_m_gd)
dev.off()

p_lc_m_gd <- ggplot(pca_lc_m, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = male_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_lc_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_lc_m[2], "%")) +
  ggtitle("LC – Males") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_LC_Male_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_lc_m_gd)
dev.off()

p_vmh_m_gd <- ggplot(pca_vmh_m, aes(PC1, PC2, color = Genotype, shape = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = male_genotype_colors) +
  xlab(paste0("PC1: ", percentVar_vmh_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_vmh_m[2], "%")) +
  ggtitle("VMH – Males") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_VMH_Male_GenoColor_DietShape.pdf", width = 6, height = 5)
print(p_vmh_m_gd)
dev.off()

# Only sex plots 
p_arc_sex <- ggplot(pca_arc_sex, aes(PC1, PC2, color = Sex)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = sex_colors) +
  xlab(paste0("PC1: ", percentVar_arc_sex[1], "%")) +
  ylab(paste0("PC2: ", percentVar_arc_sex[2], "%")) +
  ggtitle("ARC – Sex") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_sex, ylim = pc2_lim_sex)
save_plot("PCA_ARC_Sex.pdf", width = 6, height = 5)
print(p_arc_sex)
dev.off()

p_lc_sex <- ggplot(pca_lc_sex, aes(PC1, PC2, color = Sex)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = sex_colors) +
  xlab(paste0("PC1: ", percentVar_lc_sex[1], "%")) +
  ylab(paste0("PC2: ", percentVar_lc_sex[2], "%")) +
  ggtitle("LC – Sex") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_sex, ylim = pc2_lim_sex)
save_plot("PCA_LC_Sex.pdf", width = 6, height = 5)
print(p_lc_sex)
dev.off()

p_vmh_sex <- ggplot(pca_vmh_sex, aes(PC1, PC2, color = Sex)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = sex_colors) +
  xlab(paste0("PC1: ", percentVar_vmh_sex[1], "%")) +
  ylab(paste0("PC2: ", percentVar_vmh_sex[2], "%")) +
  ggtitle("VMH – Sex") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_sex, ylim = pc2_lim_sex)
save_plot("PCA_VMH_Sex.pdf", width = 6, height = 5)
print(p_vmh_sex)
dev.off()

# Females: Diet only
p_arc_f_diet <- ggplot(pca_arc_f, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_arc_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_arc_f[2], "%")) +
  ggtitle("ARC – Females (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_ARC_Female_Diet.pdf", width = 6, height = 5)
print(p_arc_f_diet)
dev.off()

p_lc_f_diet <- ggplot(pca_lc_f, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_lc_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_lc_f[2], "%")) +
  ggtitle("LC – Females (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_LC_Female_Diet.pdf", width = 6, height = 5)
print(p_lc_f_diet)
dev.off()

p_vmh_f_diet <- ggplot(pca_vmh_f, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_vmh_f[1], "%")) +
  ylab(paste0("PC2: ", percentVar_vmh_f[2], "%")) +
  ggtitle("VMH – Females (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_f, ylim = pc2_lim_f)
save_plot("PCA_VMH_Female_Diet.pdf", width = 6, height = 5)
print(p_vmh_f_diet)
dev.off()

# Males: Diet only
p_arc_m_diet <- ggplot(pca_arc_m, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_arc_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_arc_m[2], "%")) +
  ggtitle("ARC – Males (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_ARC_Male_Diet.pdf", width = 6, height = 5)
print(p_arc_m_diet)
dev.off()

p_lc_m_diet <- ggplot(pca_lc_m, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_lc_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_lc_m[2], "%")) +
  ggtitle("LC – Males (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_LC_Male_Diet.pdf", width = 6, height = 5)
print(p_lc_m_diet)
dev.off()

p_vmh_m_diet <- ggplot(pca_vmh_m, aes(PC1, PC2, color = Diet)) +
  geom_point(size = pt_size) +
  scale_color_manual(values = diet_colors) +
  xlab(paste0("PC1: ", percentVar_vmh_m[1], "%")) +
  ylab(paste0("PC2: ", percentVar_vmh_m[2], "%")) +
  ggtitle("VMH – Males (Diet)") +
  pca_theme +
  coord_fixed(xlim = pc1_lim_m, ylim = pc2_lim_m)
save_plot("PCA_VMH_Male_Diet.pdf", width = 6, height = 5)
print(p_vmh_m_diet)
dev.off()

# Helper function for diet-split genotype plots ( not sure if this is relevant, but I still made them, but did not include them in the presentation)
plot_diet_geno <- function(pca_data, title, lim_pc1, lim_pc2, colors) {
  pv <- round(100 * attr(pca_data, "percentVar"))
  ggplot(pca_data, aes(PC1, PC2, color = Genotype)) +
    geom_point(size = pt_size) +
    scale_color_manual(values = colors) +
    xlab(paste0("PC1: ", pv[1], "%")) +
    ylab(paste0("PC2: ", pv[2], "%")) +
    ggtitle(title) +
    pca_theme +
    coord_fixed(xlim = lim_pc1, ylim = lim_pc2)
}

# Females: diet-split genotype
save_plot("PCA_ARC_Female_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_arc_f_CD, "ARC – Females (CD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

save_plot("PCA_LC_Female_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_lc_f_CD, "LC – Females (CD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

save_plot("PCA_VMH_Female_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_vmh_f_CD, "VMH – Females (CD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

save_plot("PCA_ARC_Female_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_arc_f_WD, "ARC – Females (WD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

save_plot("PCA_LC_Female_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_lc_f_WD, "LC – Females (WD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

save_plot("PCA_VMH_Female_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_vmh_f_WD, "VMH – Females (WD)", pc1_lim_f_diet, pc2_lim_f_diet, female_genotype_colors))
dev.off()

# Males: diet-split genotype
save_plot("PCA_ARC_Male_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_arc_m_CD, "ARC – Males (CD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()

save_plot("PCA_LC_Male_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_lc_m_CD, "LC – Males (CD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()

save_plot("PCA_VMH_Male_CD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_vmh_m_CD, "VMH – Males (CD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()

save_plot("PCA_ARC_Male_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_arc_m_WD, "ARC – Males (WD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()

save_plot("PCA_LC_Male_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_lc_m_WD, "LC – Males (WD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()

save_plot("PCA_VMH_Male_WD_GenoColor.pdf", width = 6, height = 5)
print(plot_diet_geno(pca_vmh_m_WD, "VMH – Males (WD)", pc1_lim_m_diet, pc2_lim_m_diet, male_genotype_colors))
dev.off()
