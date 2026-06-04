# Performing Differential expressed genes (DEG) analysis between tissues, to confirm/check if the tissue is "correct"
# Have a list of "biology" genes, taken from literature, which will be used to check after the DEG analysis, using VMH as the reference tissue

# Region enriched genes
arc_biology_genes <- c( 
  "Pomc", "Agrp", "Npy", "Cartpt", "Ghrh", "Th", "Kiss1", 
  
  "Lepr", "Ghsr", "Htr2c", "Calb1", "Acvr1c", "Mc4r", 
  
  "Glp1r", "Cnr1", "Calcr", "Crhr2", "Otp", "Crhr1", 
  
  "Sst", "Oxtr", "Prlr", "Npy2r", "Mc3r" 
) 

vmh_biology_genes <- c(
  "Myl9", "Ly6c1", "Mrc1", "Ctss", "Aqp4", "Opalin", 

"Slc17a6", "Adcyap1", "Fezf1", "Inpp4b", 

"Chst9", "Gldn", "Gpc3", "Tnfaip8l3",  

"Ctxn3", "Nup62cl", "Satb2", 

"Rspo3", "Slc35f4", "Gpr101", "Arhgap36", 

"Egflam", "Nr5a1", "Foxp2", "Scgn", 

"C1ql2", "Six3", "Klhl1", "Sh3bgrl2"
)

lc_biology_genes <- c(
  "Agtr1a", "Asb4", "Calr3", "Chodl", "Chrna6", "Cilp", 

"Cyb561", "Dbh", "Ddc", "Adgre1", "Eya2", "Shisal2b", 

"Cfap144", "Fibcd1", "Gal", "Gch1", "Glra2", "Gng4", 

"Gpx3", "Gtf2a1l", "Hcrtr1", "Igsf5", "Maoa", "Mrap2", 

"Myom2", "Neurog2", "Slc9b2", "Nxph4", "Ovgp1", "Pcbd1", 

"Phox2a", "Phox2b", "Pla2g4d", "Ptger2", "Slc18a2", 

"Slc31a1", "Slc6a2", "Stbd1", "Syt17", "Th", 

"Tm4sf1", "Tm4sf5", "Traf3ip2" 
)


# Splitting the dataset by sex
dds_tissue_F <- dds_full_filt[, dds_full_filt$Sex == "F"]
dds_tissue_M <- dds_full_filt[, dds_full_filt$Sex == "M"]

# Setting VMH as the reference tissue
dds_tissue_F$Tissue <- relevel(dds_tissue_F$Tissue, ref = "VMH")
dds_tissue_M$Tissue <- relevel(dds_tissue_M$Tissue, ref = "VMH")

# Setting tissue to design
design(dds_tissue_F) <- ~ Tissue
design(dds_tissue_M) <- ~ Tissue

# Running the DESeq2
dds_tissue_F <- DESeq(dds_tissue_F)
dds_tissue_M <- DESeq(dds_tissue_M)


# Getting the results
res_ARC_vs_VMH_F <- results(dds_tissue_F, contrast = c("Tissue", "ARC", "VMH"), lfcThreshold = 0.263, altHypothesis = "greaterAbs", alpha = 0.05)
res_LC_vs_VMH_F  <- results(dds_tissue_F, contrast = c("Tissue", "LC",  "VMH"), lfcThreshold = 0.263, altHypothesis = "greaterAbs", alpha = 0.05)
res_ARC_vs_VMH_M <- results(dds_tissue_M, contrast = c("Tissue", "ARC", "VMH"), lfcThreshold = 0.263, altHypothesis = "greaterAbs", alpha = 0.05)
res_LC_vs_VMH_M  <- results(dds_tissue_M, contrast = c("Tissue", "LC",  "VMH"), lfcThreshold = 0.263, altHypothesis = "greaterAbs", alpha = 0.05)


# Check overlap with ARC/LC biology genes females
overlapp_ARC_F <- intersect(arc_bio_ens,
                             rownames(res_ARC_vs_VMH_F[
                               which(res_ARC_vs_VMH_F$padj < 0.05 &
                                     res_ARC_vs_VMH_F$log2FoldChange > 0), ]))

overlapp_LC_F <- intersect(lc_bio_ens,
                            rownames(res_LC_vs_VMH_F[
                              which(res_LC_vs_VMH_F$padj < 0.05 &
                                    res_LC_vs_VMH_F$log2FoldChange > 0), ]))

# VMH: genes downregulated in both ARC and LC vs VMH
vmh_down_both_F <- intersect(
  rownames(res_ARC_vs_VMH_F[which(res_ARC_vs_VMH_F$padj < 0.05 & res_ARC_vs_VMH_F$log2FoldChange < 0), ]),
  rownames(res_LC_vs_VMH_F[which(res_LC_vs_VMH_F$padj < 0.05  & res_LC_vs_VMH_F$log2FoldChange < 0), ])
)
overlapp_VMH_F <- intersect(vmh_bio_ens, vmh_down_both_F)


# Check overlap with ARC/LC biology genes males
overlapp_ARC_M <- intersect(arc_bio_ens,
                             rownames(res_ARC_vs_VMH_M[
                               which(res_ARC_vs_VMH_M$padj < 0.05 &
                                     res_ARC_vs_VMH_M$log2FoldChange > 0), ]))

overlapp_LC_M <- intersect(lc_bio_ens,
                            rownames(res_LC_vs_VMH_M[
                              which(res_LC_vs_VMH_M$padj < 0.05 &
                                    res_LC_vs_VMH_M$log2FoldChange > 0), ]))

# VMH: genes downregulated in both ARC and LC vs VMH
vmh_down_both_M <- intersect(
  rownames(res_ARC_vs_VMH_M[which(res_ARC_vs_VMH_M$padj < 0.05 & res_ARC_vs_VMH_M$log2FoldChange < 0), ]),
  rownames(res_LC_vs_VMH_M[which(res_LC_vs_VMH_M$padj < 0.05  & res_LC_vs_VMH_M$log2FoldChange < 0), ])
)
overlapp_VMH_M <- intersect(vmh_bio_ens, vmh_down_both_F)
