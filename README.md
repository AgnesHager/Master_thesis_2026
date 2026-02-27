# RNA-Seq Analysis: Adipose Slc7a10 KI Effects on Brain Transcriptome 
----------------------------
Project Overview:
------------------------
Analysis of RNA-seq data from hypothalamic and brainstem regions (ARC, LC, VMH) in mice with adipose-specific Slc7a10 overexpression under Carbohydrate diet (CD) and Western diet (WD) conditions.

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
**Input:** FeatureCounts matrix (all 140 samples, 'counts.txt')

**Steps:**
1. Remove problematic samples ( mouse 600 (sick), test fragments)
2. Global low-count filtering:'rowSums(counts >= 10) >=5'
       -Filters genes with <10 counts in <5 samples across all samples
       -Reduces genes from ~50,000 to 20,614 genes
3. Create global DESeqDataSet ('dds_full')
4. VST transformation ('blind = TRUE') for visualization
5. Global PCA analysis to assess:
       -Tissue separation (primary driver)
       -Sex effects (secondary driver)
       -Outlier detection

**Output:**
-'dds_full': Filtered count matrix (20,614 genes x 140 samples)
-'vsd_full': VST-transformed data for visualization
-Global QC plots (PCAs)

**Note:** 'dds_full' is the base dataset. For each differential expression contrast, samples are subset and DESeq2 is re-run on the subset to ensure proper dispersion estimates for that specific comparison

--------------------
### 3. Differential Expression Analysis ('DE_analysis_DESeq2.R')
**Purpose:** Identify differentially expressed genes across 24 contrasts
**Input:** 
- Raw count matrix (after removing problematic samples and all-zero genes)
- Metadata table

Software & packages
* DESeq2


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
### 4. Plotting 

**Software & packages**
* PCA
* Venn diagrams
* MA and volcano
* Enrichment dotplots

**Scripts**
* pca_DESeq2.R

  
