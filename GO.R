# Gene ontology analysis on relevant modules after performing WGCNA
# Load WGCNA results
load(file.path(OUT_DIR, paste0(RUN_ID, "_bwnet.RData")))

# Select module of interest
MODULE <- "pink"   # change as needed

genes_module <- names(bwnet$colors[bwnet$colors == MODULE])
cat("Number of genes in", MODULE, "module:", length(genes_module), "\n")

# GO enrichment
go_res <- enrichGO(
  gene          = genes_module,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Save plot
save_plot(file.path(OUT_DIR, paste0("GO_BP_", MODULE, "_dotplot.pdf")),
          width = 6, height = 10)

dotplot(go_res, showCategory = 20,
        title = paste("GO Biological Process:", RUN_ID, MODULE, "module")) +
  theme(axis.text.y = element_text(size = 8))

dev.off()
