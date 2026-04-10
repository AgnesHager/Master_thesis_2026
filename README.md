# RNA-Seq Analysis: Adipose Gene Of Interest (GOI) KI Effects on Brain Transcriptome 
----------------------------
Project Overview:
------------------------
Analysis of RNA-seq data from hypothalamic and brainstem regions (ARC, LC, VMH) in mice with adipose-specific GOI overexpression under Carbohydrate diet (CD) and Western diet (WD) conditions.

Experimental Design:
* 3 tissues (ARC, LC, VMH)
* 2 sexes (Male/female)
* 2 diets (CD, WD)
* 2 genotypes (WT, KI)
* n = 5-6 per group

## Pipeline Overview:

### 1. Raw Data Processing (Performed by Genomics Core Facility)
**Alignment:**
* HISAT2 v2.2.2.1
* Reference: GRCm38.p5
* Parameters: -p 140 --new-summary

**Quantification:**
* featureCounts
* Annotation: Gencode vM13
* Parameters: -t exon -g gene_id -p

--------------------
### 2. Quality Control ('QC_analysis.R')
**Purpose:** Assess global data quality and sample relationships
**Input:** DESeqDataSet with all 150 samples ('dds_full'), including test fragments and different concentrations 

Scripts
*'QC_analysis.R'

**Steps:**
1. Define sample categories (Standard, Control, True_Low_Input, Large_Library, Low_RIN)
2. VST transformation ('blind = TRUE') for visualization
3. PCA with technical annotations to visually assess whether samples with low RIN, large libraries, or low input clustered abnormally
          -PCA plot colored by Tissue, shaped by Sample_category
4. Remove mouse 600 (sick) and test fragment samples
       -Reduces from 150 to 140 samples
5. Global low-count filtering:'rowSums(counts >= 10) >=5'
       -Filters genes with <10 counts in <5 samples across all samples
       -Reduces genes from ~50,000 to 20,614 genes
6. Pre-DESeq2 QC on cleaned 140 samples:
       -Library size distribution
       -Raw count distribution
       -Gene biotype composition (per tissue: ARC, LC, VMH)
       -Mean-SD plots before and after VST
7. VST transformation ('blind = TRUE') for visualization
8. Global PCA analysis to assess:
       -Tissue separation 
       -Sex effects 
       -Diet effects
       -Genotype effects

**Output:**
-'dds_full': Filtered count matrix (20,614 genes x 140 samples)
-'vsd_full': VST-transformed data for visualization
-QC plots (PDFs): library sizes, count distributions, biotypes, PCA's

**Note:** No samples were removed based on their QC metrics as all samples clustered correctly within their tissue groups. Only mouse 600 (sick) and test fragment samples were removed.  'dds_full' is the base dataset. For each differential expression contrast, samples are subset and DESeq2 is re-run on the subset to ensure proper dispersion estimates for that specific comparison

--------------------
### 3. Differential Expression Analysis ('DE_analysis_DESeq2.R')
**Purpose:** Identify differentially expressed genes across 24 contrasts
**Input:** 
- Raw count matrix (after removing test samples and mouse 600, and all-zero genes)
- Metadata table

Software & packages
* DESeq2 (v. 1.46.0)
* apeglm (v. 1.28.0)


Scripts
* DE_analysis_DESeq2.R

**Analysis Strategy:**
For each of 24 contrasts:
1. **Subset samples** by tissue, sex and condition
2. **Contrast-specific filtering**:
     -Calculate minimum group size (n_min)
     -Keep genes with >= 10 counts in >=n_min samples
     -Ensures only expressed genes in that specific comparison are tested
3. **Drop unused factor levels** with 'droplevels()'
4. **Set reference level** (WT for genotype, CD for diet)
5. **Run DESeq2** normalization and dispersion estimation
6. **Extract results** and apply apeglm LFC shrinkage
7. **Filter significant DEGs**: padj < 0.05 AND |log2FC| >= 0.263

**Total Contrasts:** 24
-Genotype effects: KI vs WT (per tissue x sex x diet)
-Diet effects: WD vs CD (per tissue x sex x genotype)

**Gene Exclusions:**
- Gt(ROSA)26Sor: Technical artifact (KI construct insertion site)
- Adipoq: Technical artifact (KI promoter construct)

**Output:**
- `*_stats.csv`: All genes with statistics
- `*_sig.csv`: Significant DEGs only


-------------------
## 4. Weighted Gene Co-expression Network Analysis (WGCNA)

### Script
`WGCNA_pipeline.R`

### User parameters
- TISSUE (ARC / LC / VMH)
- SEX (M / F)
- CONDITION (CD / WD)

### Workflow

1. Variance filtering (top 50% genes)
2. Sample clustering QC
3. Soft-threshold selection
4. Network construction (blockwiseModules)
5. Module eigengene calculation
6. Module–trait correlation (bicor)
7. Hub gene detection (kME + GS)
8. Cytoscape export (top module)

------------------
### 6. Functional Enrichment Analysis (GO)

Gene Ontology enrichment is performed on selected WGCNA modules using `clusterProfiler`.

**Input:**
- Module gene lists from WGCNA output

**Analysis:**
- enrichGO (Biological Process)
- Multiple testing correction (BH)

**Output:**
- Dotplot visualizations

----------------------------
### 5. Plotting 

**Plot types**
* PCA
* Boxplots (individual DEG expression)
* Heatmap (DEGs)
* Venn diagrams
* MA and volcano
* Module–trait heatmaps
* Enrichment dotplots

**Scripts**
* QC_control.R (PCA plots)
* Plotting_DEGs.R
* Heatmap_DEGs.R
* WGCNA.R


  
