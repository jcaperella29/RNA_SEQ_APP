🧪 RNA-SEQ APP – Quick Start Guide
📂 1. Upload Your Data
Counts Matrix (.csv or .txt)
→ Rows = genes, Columns = samples

Phenotype File (.csv or .txt)
→ One column = condition/group labels

✅ Once uploaded, select the correct phenotype column from the dropdown.

🧬 2. Differential Expression
Click "Run Differential Expression"

App runs voom + limma pipeline

Output = Top 50 DE genes (adjusted p-value sorted)

DEGs used for downstream Random Forest modeling

🧾 View: Differential Expression Results tab

🌲 3. Feature Selection + Random Forest
From top 50 DEGs, top 20 genes are selected by Random Forest importance

RF model is trained on those top 20 predictors

Performance metrics calculated (Train/Test split)

📊 Check Classification > Metrics tab for:

Metric	Description
Accuracy	Overall classification accuracy
Sensitivity	True positive rate
Specificity	True negative rate
AUC (ROC)	Area under ROC curve (2-class)

📊 4. Visualization
Tool	What It Shows
PCA Plot	Sample clustering (Principal Components)
UMAP Plot	Non-linear low-dim embedding
Volcano Plot	DE results: log2FC vs -log10 adj.P
Heatmap	Clustered Z-score heatmap of DEGs
👉 All plots are interactive via plotly.

🧠 5. Pathway Enrichment
Use EnrichR to analyze gene set enrichment:

Click one of:

Enrich All DE Genes

Enrich Upregulated

Enrich Downregulated

Choose from:

KEGG

GO BP

Reactome

WikiPathways

📈 Results shown as:

Table (adjusted p, combined score)

Barplot (-log10 p-adjusted)

📈 6. Power Analysis
Click "Run Power Analysis"

Estimates power based on your phenotype data

Supports:

t-test (2 groups)

ANOVA (>2 groups)

Adjustable effect size & sample size range

Plots power curves using pwr package

📤 7. Downloadable Outputs
Each major result has a download button:

Output	Button Label
DEGs	Download DE Results
Enrichment tables	Download Enrichment
Classifier predictions	Download RF Predictions
RF metrics	Download RF Metrics
Power summary	Download Power Summary
⚠️ Error Logging
All errors are saved to error_log.txt

Developers can enable scheduled error emails via:

Rscript email_log.R
Built with 💻 + 🧬 by JCaperella
MIT Licensed — Reproducible science > black boxes
