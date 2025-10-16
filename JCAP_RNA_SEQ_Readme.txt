🧬 JCAP RNA-SEQ APP
Summary of the Main Functions

This Shiny application accepts a counts matrix and a phenotype file and performs a complete exploratory RNA-seq workflow through point-and-click modules.

The app now supports multiple species and can automatically map Ensembl gene IDs to the correct annotation before enrichment.

Supported Species

Human (hsapiens)

Mouse (mmusculus)

Zebrafish (drerio)

Fruit Fly (dmelanogaster)

You select the species from a drop-down menu before analysis. All downstream annotation and pathway enrichment automatically use the appropriate organism databases.

Core Functions

Differential Expression + Feature Selection

Conducted with limma-voom after QC (low-count filtering + outlier removal).

Top 200 DE genes → feature selection via random forest.

The 20 most influential genes (based on permutation importance) are listed with:

Ensembl ID 

HGNC symbol / orthologous symbol

log₂ fold change

adjusted p-value (FDR)

protein product name (if available)

Table exportable as CSV.

Principal Component Analysis (PCA)

Colored by phenotype; hover to identify samples.

Uniform Manifold Approximation and Projection (UMAP)

Visualizes nonlinear sample relationships; hover for sample IDs.

Volcano Plot

Shows log₂ fold change vs adjusted p-value.

Hover reveals gene symbol; save via camera icon.

Pathway Enrichment Analysis

Runs both g:Profiler2 and EnrichR behind the scenes.

Choose among Reactome, KEGG, GO Biological Process, WikiPathways, etc.

Separate buttons for:

All top genes

Upregulated genes

Downregulated genes

Species-specific databases selected automatically.

Parallel computing accelerates differential expression and dimensionality-reduction steps.

Prerequisites

A counts matrix derived from read quantification (.csv or .txt).

A phenotype table describing sample groups (also .csv or .txt).

Only one phenotype variable is analyzed per run.

You can re-analyze for different phenotypes by re-uploading.

How to Use
Inputs

Input Counts Matrix → upload file → “Upload complete” ribbon.

Input Phenotype Data → upload file → “Upload complete.”

Select Species from the drop-down (Human, Mouse, Zebrafish, Fly).

Click “Perform Differential Expression Analysis.”

You’ll receive start/finish notifications.

Results appear under “Differential Expression Results.”

Plots

Display PCA Plot → colored by phenotype.

Display UMAP Plot → alternative nonlinear projection.

Display Volcano Plot → after DE table generation.

Pathway Enrichment

Enrich Pathways (All Genes)

Enrich Pathways (Upregulated)

Enrich Pathways (Downregulated)
Each returns an interactive enrichment table and bar plot.

Sliders

PCA components displayed

UMAP neighbors

Volcano p-value threshold

Volcano log FC threshold

Downloads

Export DE Results as CSV

Export All-Genes Enrichment as CSV

Export Upregulated Enrichment as CSV

Export Downregulated Enrichment as CSV

Notes

Analysis buttons remain disabled until required inputs are present.

Download buttons activate only after their respective results are displayed.

All visualizations support zoom, pan, and image export via Plotly controls.

The README tab reproduces this documentation inside the app.

Example Test Dataset

A helper script make_mouse_cancer_like_data.R can generate a biologically realistic mouse dataset enriched for cell-cycle and DNA replication pathways.
Run:

source("make_mouse_cancer_like_data.R")


Then load the resulting counts_mmusculus_cancer.csv and phenotype_mmusculus_cancer.csv into the app with Species = Mouse to test the full pipeline.

Acknowledgments

Developed by J. Caperella, combining standard Bioconductor methods with modern machine-learning-based feature selection and cross-species pathway annotation.
