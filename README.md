Project Title:
Exploratory Transcriptomic Analysis of Renal Tubular Cells in Rats Exposed to Ochratoxin A

Project Description:
This project involves an exploratory analysis of RNA-seq data derived from rat renal proximal tubular cells following 13-week exposure to Ochratoxin A (OTA) a known nephrotoxic and karyomegaly-inducing compound. The dataset was obtained from the Gene Expression Omnibus (GEO Accession: GSE231379), and includes expression profiles from control and chemically treated rats.

The analysis focuses on:

Loading and preprocessing raw gene expression data

Filtering out lowly expressed genes using edgeR

Normalizing expression data using log<sub>2</sub> Counts Per Million (CPM)

Performing summary statistics and data visualization to assess sample quality and expression distribution

Identifying the top 100 most variable genes

Visualizing expression patterns using heatmaps and Principal Component Analysis (PCA)

This project serves as an exploratory step toward understanding OTA-induced transcriptional changes in renal tissue and sets the stage for downstream differential expression or pathway enrichment analyses.

 Technologies Used:
R (v4.5.0)

edgeR

ggplot2

pheatmap

tibble

Exploratory Analysis Summary
- Raw count data from rat kidney samples were loaded and inspected.
- Boxplots and histograms were generated to visualize sample distribution.
- Lowly expressed genes were filtered out using edgeR’s `filterByExpr()`.
- Normalized expression values were computed using log2 CPM.
- PCA and heatmap visualizations revealed sample clustering and high-variance genes.
