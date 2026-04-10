# Loading packages
library(WGCNA)
library(matrixStats)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(survival)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(gridExtra)
library(pheatmap)
library(tidyverse)
library(matrixStats)

save_plot <- function(filename, width = 8, height = 6) {
  cairo_pdf(filename, width = width, height = height, family = "Arial")
}

# Phenotype tables, using to correlate afterwards
pheno_F <- Reduce(function(a, b) merge(a, b, by = "Mouse_ID"),
                  list(
                    df_BWgain_female[, c("Mouse_ID","BW_gain")],
                    df_iWAT_size_female[, c("Mouse_ID","iWAT_size")],
                    df_oWAT_size_female[, c("Mouse_ID","oWAT_size")],
                    df_auc_insulin_female[, c("Mouse_ID","auc_insulin")],
                    df_matsuda_index_female[, c("Mouse_ID","matsuda_index")],
                    df_food_intake_daily_female[, c("Mouse_ID","food_intake_daily")],
                    df_RER_female[, c("Mouse_ID","Tot_AUC_RER")],
                    df_WAT_female[, c("Mouse_ID", "WAT")]
                  ))
# Male
pheno_M <- Reduce(function(a, b) merge(a, b, by = "Mouse_ID"),
                  list(
                    df_BWgain_male[, c("Mouse_ID","BW_gain")],
                    df_iWAT_size_male[, c("Mouse_ID","iWAT_size")],
                    df_eWAT_size_male[, c("Mouse_ID","eWAT_size")],
                    df_auc_insulin_male[, c("Mouse_ID","auc_insulin")],
                    df_matsuda_index_male[, c("Mouse_ID","matsuda_index")],
                    df_food_intake_daily_male[, c("Mouse_ID","food_intake_daily")],
                    df_RER_male[, c("Mouse_ID","Tot_AUC_RER")],
                    df_WAT_male[, c("Mouse_ID", "WAT")]
                  ))

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

TISSUE <- "Tissue_X"  
SEX    <- "F"        # Change to what you're looking at 
CONDITION <- "CD"    # Change to what you're looking at

RUN_ID <- paste(TISSUE, SEX, CONDITION, sep = "_")              
OUT_DIR <- file.path("Out", RUN_ID)                    # making a folder where the results will be saved, example "Out/ARC_F_WD"
# sanity check BEFORE loading
if (!file.exists("dds_input.rds")) {
  stop("dds_input.rds not found in working directory")
}
DDS_OBJ <- readRDS("dds_input.rds") # Must be a DESeq2 object (DESeqDataSet), required colData columns: Sex, Diet, Genotype, MouseID

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# STEP 1: Data preparation: subset samples and normalize
metaData <- as.data.frame(colData(DDS_OBJ))
keep      <- metaData$Diet == CONDITION & metaData$Sex == SEX

dds_sub  <- DDS_OBJ[, keep]
metaData <- metaData[keep, ]

cat("Samples retained (DIET +", SEX, "):", ncol(dds_sub), "\n")

# VST on the subset only
vsd         <- vst(dds_sub, blind = FALSE)
norm.counts <- assay(vsd) %>% t()   # samples (rows) x genes (columns)

cat("Matrix dimensions (samples x genes):", dim(norm.counts), "\n")

# STEP 2: Build trait table for module-trait associations, Genotype is coded as binary: WT=0, KI=1
# Match phenotype to expression samples via MouseID in colData
# rownames of metaData are sample IDs 
pheno <- if (SEX == "F") pheno_F else pheno_M

traitData <- merge(
  data.frame(MouseID  = metaData$MouseID,
             SampleID = rownames(metaData)),
  pheno,
  by.x = "MouseID", by.y = "Mouse_ID",
  all.x = TRUE
)
rownames(traitData) <- traitData$SampleID
traitData$SampleID  <- NULL
traitData$MouseID   <- NULL

# Add genotype as numeric trait (WT=0, KI=1)
traitData$genotype_num <- ifelse(metaData[rownames(traitData), "Genotype"] == "WT", 0, 1)

cat("Traits included:", paste(colnames(traitData), collapse = ", "), "\n")
cat("Trait table dimensions:", dim(traitData), "\n")


# STEP 3: Quality check — remove bad samples / genes
gsg <- goodSamplesGenes(norm.counts, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)   > 0)
    cat("Removing genes:",   paste(colnames(norm.counts)[!gsg$goodGenes],   collapse = ", "), "\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Removing samples:", paste(rownames(norm.counts)[!gsg$goodSamples], collapse = ", "), "\n")
  
  norm.counts <- norm.counts[gsg$goodSamples, gsg$goodGenes]
  traitData   <- traitData[gsg$goodSamples, , drop = FALSE]
}

# STEP 4: Sample clustering (coloured by genotype)
sampleTree     <- hclust(dist(norm.counts), method = "average")
genotypeColors <- ifelse(traitData$genotype_num == 0, "blue", "red")

save_plot(file.path(OUT_DIR, "sampleClustering.pdf"), width = 12, height = 9)
par(cex = 0.6, mar = c(4, 4, 2, 0))
plot(sampleTree,
     main = paste("Sample clustering: ", RUN_ID),
     sub  = "Blue = WT  |  Red = KI",
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
n <- length(sampleTree$order)
rect(seq_len(n) - 0.5,
     par("usr")[3] - diff(par("usr")[3:4]) * 0.06,
     seq_len(n) + 0.5,
     par("usr")[3],
     col    = genotypeColors[sampleTree$order],
     border = NA, xpd = TRUE)
dev.off()

# STEP 5: Filter to top 50% most variable genes, Reduces the pairwise correlation matrix to a smaller sizekeeping only the most variable genes
# For n < 30 samples, retaining top 50% is recommended
# For large datasets (> 50 samples), keeping the top 10-25% by variance is recommended
# Source: Nguyen P, Zeng E. A Protocol for Weighted Gene Co-expression Network Analysis With Module Preservation and Functional Enrichment Analysis for Tumor and Normal Transcriptomic Data. Bio-Protoc. 2025 Sep 20;15(18):e5447. doi:10.21769/BioProtoc.5447 PubMed PMID: 41000162; PubMed Central PMCID: PMC12457846
gene_vars      <- matrixStats::colVars(norm.counts)
n_keep         <- floor(length(gene_vars) * 0.50)   # <-- adjust if needed
high_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_keep]
norm.counts    <- norm.counts[, high_var_genes]

cat("Genes retained after variance filtering:", ncol(norm.counts), "\n")


# STEP 6: Calculating pairwise correlations + soft thresholding
power <- c(1:20)
sft   <- pickSoftThreshold(norm.counts,
                           powerVector = power,
                           networkType = "signed hybrid",
                           verbose     = 5)

sft.data <- sft$fitIndices

p1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05) +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(x = "Soft threshold (power)",
       y = "Scale free topology model fit (R^2)",
       title = paste("Scale independence: ", RUN_ID)) +
  theme_classic()

p2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = "Soft threshold (power)",
       y = "Mean connectivity",
       title = paste("Mean connectivity: ", RUN_ID)) +
  theme_classic()

save_plot(file.path(OUT_DIR, "softThreshold.pdf"), width = 10, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# Choose the LOWEST power where signed R² >= 0.8
soft_power <- x # Insert number after checking softThreshold.pdf


# STEP 7: Topological overlap matrix (TOM) + network construction, 
# corType = "bicor": biweight midcorrelation, robust to outliers
bwnet <- blockwiseModules(
  norm.counts,
  maxBlockSize   = 14000,
  networkType = "signed hybrid", # recommended by Langfelder over pure "signed"
  TOMType        = "signed",  
  power          = soft_power,
  corType        = "bicor",
  maxPOutliers   = 0.05,         # max 5% of samples treated as outliers per gene
  # stabilizes bicor estimates in small datasets
  deepSplit      = 2,
  mergeCutHeight = 0.25,
  minModuleSize  = 30,
  numericLabels  = FALSE,
  randomSeed     = 1234,
  verbose        = 3
)

# Save immediately 
save(bwnet, file = file.path(OUT_DIR, paste0(RUN_ID, "_bwnet.RData")))


# STEP 8: Module eigengene calculation, inspect modules
module_eigengenes <- bwnet$MEs

cat("\nNumber of genes per module:\n")
print(table(bwnet$colors))

save_plot(file.path(OUT_DIR, "moduleTree.pdf"), width = 12, height = 9)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  cbind(bwnet$unmergedColors, bwnet$colors),
  c("Unmerged", "Merged"),
  dendroLabels = FALSE,
  addGuide     = TRUE,
  hang         = 0.03,
  guideHang    = 0.05,
  main         = paste("Gene dendrogram: ", RUN_ID)
)
dev.off()


# STEP 9: Module-trait associations
# Correlate module eigengenes to external phenotypic traits
nSamples <- nrow(norm.counts)
nGenes   <- ncol(norm.counts)

# Explicitly align rows before cbind
traitData    <- traitData[rownames(module_eigengenes), , drop = FALSE]
heatmap.data <- cbind(module_eigengenes, traitData)

# Calculate bicor correlations and p-values
module.trait.corr       <- bicor(module_eigengenes, traitData,
                                 maxPOutliers = 0.05,
                                 use = "pairwise.complete.obs")
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Add significance stars (* p<0.05, ** p<0.01, *** p<0.001)
sig_stars <- matrix("", nrow = nrow(module.trait.corr), ncol = ncol(module.trait.corr))
sig_stars[module.trait.corr.pvals < 0.001] <- "***"
sig_stars[module.trait.corr.pvals >= 0.001 & module.trait.corr.pvals < 0.01]  <- "**"
sig_stars[module.trait.corr.pvals >= 0.01  & module.trait.corr.pvals < 0.05]  <- "*"

display_numbers <- matrix(
  paste0(round(module.trait.corr, 2), sig_stars),
  nrow = nrow(module.trait.corr),
  dimnames = dimnames(module.trait.corr)
)

# Plot heatmap using pre-calculated bicor values
save_plot(file.path(OUT_DIR, "moduleTraitHeatmap.pdf"), width = 12, height = 14)
pheatmap(
  module.trait.corr,
  display_numbers = display_numbers,
  color           = colorRampPalette(c("blue1", "skyblue", "white", "pink", "red"))(100),
  breaks          = seq(-1, 1, length.out = 101),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  fontsize_number = 7,
  fontsize_row    = 8,
  fontsize_col    = 9,
  angle_col       = 45,
  main            = paste("Module-trait relationships:", RUN_ID)
)
dev.off()

# STEP 10: Rank modules by genotype association
genotype_cor <- module.trait.corr[, "genotype_num"]
genotype_p   <- module.trait.corr.pvals[, "genotype_num"]

genotype_hits <- data.frame(
  module  = names(genotype_cor),
  r       = round(genotype_cor, 3),
  p_value = signif(genotype_p, 3)
) %>% arrange(p_value)

write.table(genotype_hits,
            file = file.path(OUT_DIR, "genotype_module_correlations.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nTop modules correlated with genotype (WT vs KI):\n")
print(head(genotype_hits, 10))


# STEP 11: Module membership (kME) and gene significance (GS)
# kME = correlation of each gene to its module eigengene, (how central the gene is to its module)
geneModuleMembership <- bicor(norm.counts, module_eigengenes,
                              maxPOutliers = 0.05,
                              use = "pairwise.complete.obs")
MMPvalue             <- corPvalueStudent(as.matrix(geneModuleMembership), nSamples)

geneTraitSignificance <- bicor(norm.counts, traitData$genotype_num,
                               maxPOutliers = 0.05,
                               use = "pairwise.complete.obs")
GSPvalue              <- corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)

colnames(geneTraitSignificance) <- "GS.genotype"
colnames(GSPvalue)              <- "p.GS.genotype"

# Use match() for safe indexing
hub_stats <- data.frame(
  gene   = rownames(geneModuleMembership),
  module = bwnet$colors[match(rownames(geneModuleMembership), names(bwnet$colors))],
  GS     = round(geneTraitSignificance[, 1], 3),
  GS.p   = signif(GSPvalue[, 1], 3)
)

hub_stats$kME <- mapply(function(gene, mod) {
  me_col <- paste0("ME", mod)
  if (me_col %in% colnames(geneModuleMembership))
    round(geneModuleMembership[gene, me_col], 3)
  else NA
}, hub_stats$gene, hub_stats$module)

hub_stats <- hub_stats %>% arrange(module, desc(abs(kME)))

# Standard hub gene thresholds: |kME| > 0.8 and |GS| > 0.2
hub_genes <- hub_stats %>% filter(abs(kME) > 0.8, abs(GS) > 0.2)
cat("\nHub genes identified (|kME| > 0.8 & |GS| > 0.2):", nrow(hub_genes), "\n")

write.table(hub_stats,
            file = file.path(OUT_DIR, "hub_gene_stats.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(hub_genes,
            file = file.path(OUT_DIR, "hub_genes_filtered.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)



# STEP 12: Save gene -> module assignment
module.gene.mapping <- as.data.frame(bwnet$colors)
colnames(module.gene.mapping) <- "module"

write.table(module.gene.mapping,
            file = file.path(OUT_DIR, "gene2module.txt"),
            sep = "\t", col.names = NA, quote = FALSE)

save(module_eigengenes, module.gene.mapping, traitData,
     geneModuleMembership, geneTraitSignificance,
     file = file.path(OUT_DIR, paste0(RUN_ID, "_networkConstruction.RData")))


# STEP 13: Cytoscape export for top genotype module
# Skip grey module since it contains unassigned genes that are not co-expressed 
# Grey is excluded from network export as it has no real biological module structure: 
# Source: Langfelder P, Luo R, Oldham MC, Horvath S. Is My Network Module Preserved and Reproducible? PLoS Comput Biol. 2011 Jan 20;7(1):e1001057. doi:10.1371/journal.pcbi.1001057 PubMed PMID: 21283776; PubMed Central PMCID: PMC3024255.
top_color <- sub("^ME", "", genotype_hits$module[genotype_hits$module != "MEgrey"][1])
cat("\nExporting top genotype module to Cytoscape:", top_color, "\n")

inModule  <- bwnet$colors == top_color
modProbes <- names(bwnet$colors)[inModule]

mod.counts <- norm.counts[, modProbes]
mod.TOM    <- TOMsimilarityFromExpr(
  mod.counts,
  power       = soft_power,
  networkType = "signed hybrid"
)
dimnames(mod.TOM) <- list(modProbes, modProbes)

exportNetworkToCytoscape(
  mod.TOM,
  edgeFile  = file.path(OUT_DIR, paste0("Cytoscape-edges-", top_color, ".txt")),
  nodeFile  = file.path(OUT_DIR, paste0("Cytoscape-nodes-", top_color, ".txt")),
  weighted  = TRUE,
  threshold = 0.1,
  nodeNames = modProbes,
  nodeAttr  = bwnet$colors[inModule]
)
