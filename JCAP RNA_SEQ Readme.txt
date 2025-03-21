🧪 RNA-SEQ APP – Quick Start Guide
📂 1. Upload Your Data
Upload Counts Matrix (CSV or TXT)
Upload Phenotype Data (CSV or TXT)
✅ Counts = rows as genes, columns as samples
✅ Phenotype = 1 column with condition labels

🧬 2. Differential Expression
Click "Perform Differential Expression Analysis"
View the top 50 DEGs in the Differential Expression Results tab
These DEGs are used for feature selection + modeling.

🌲 3. Feature Selection + Random Forest
App selects the top 20 genes using Random Forest importance
Then builds a new RF model with those 20 as predictors
Check Random Forest Metrics tab for:
Accuracy
Sensitivity
Specificity
AUC
📊 4. Visualization
PCA Plot and UMAP Plot for sample clustering
Tune UMAP neighbors & PCA components
Volcano Plot: logFC vs adjusted p-value
🧠 5. Enrichment Analysis
Click any of:

Enrich Pathways (All Genes)
Enrich (Upregulated Genes)
Enrich (Downregulated Genes)
→ Uses Enrichr to return pathway/gene set enrichment results
📈 6. Power Calculation
Click "Calculate Statistical Power"
Uses your phenotype data + effect size to estimate power using pwr
📤 7. Download Results
Every major table supports CSV export
Look for buttons like:
"Download Results", "Export Enrichment", "Download Performance Metrics"

⚠️ Errors?
All errors are auto-logged to error_log.txt
Developers can enable daily email reporting via email_log.R

Built with 💻 + 🧬 by JCaperella
MIT Licensed. Reproducible science > black boxes.
