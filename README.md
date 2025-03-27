# 🧬 JCAP RNA-SEQ ANALYZER

An interactive **R Shiny** web app for RNA-seq analysis — integrating differential expression, power estimation, feature selection (Random Forest), dimensionality reduction, and pathway enrichment in one clean UI.


![R](https://img.shields.io/badge/R-%3E%3D4.1-blue?logo=r&logoColor=white)
![Shiny](https://img.shields.io/badge/Shiny-Interactive%20App-orange?logo=rstudio)
![Dockerized](https://img.shields.io/badge/Docker-Ready-2496ED?logo=docker&logoColor=white)
![Singularity](https://img.shields.io/badge/Singularity-HPC%20Compatible-2A2A2A)
![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20Windows-lightgrey)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
![Sample Data](https://img.shields.io/badge/Sample%20Data-Included-blueviolet)
![Status](https://img.shields.io/badge/Status-Production%20Ready-brightgreen)

> A complete R Shiny app for RNA-seq differential expression, enrichment, and ML analysis.

> **Built for bioinformatics teams, core facilities, and high-performance computing clusters.**  
> **No coding required. Fully Docker/Singularity compatible.**

---

## 🚀 Features

- 📂 Upload gene expression counts & phenotype data
- 🧬 Differential expression (limma + voom)
- 🌋 Interactive volcano plot (logFC vs adjusted P)
- 🌌 PCA and UMAP with dynamic tooltips
- 🌲 Random Forest feature selection + classification
- 📊 Enrichment analysis with **EnrichR**
- 📈 Statistical power calculation and power curves
- 🔥 Downloadable results, ROC curves, and metrics
- 🛠️ Built-in error logging system (`error_log.txt`)
- 📦 Dockerized for local + HPC deployment
- 📄 In-app readme tab for offline guidance

---

## 📁 Folder Structure

```text
RNA_SEQ_APP/
├── app.R                    # Main Shiny app
├── run.sh                  # Docker launcher script
├── Dockerfile              # Docker image definition
├── Singularity.def         # HPC container definition
├── JCAP RNA_SEQ Readme.txt # In-app readme content
├── error_log.txt           # Auto-generated server logs
├── email_log.R             # Cron-compatible email notifier
└── www/
    └── Arcane_Alchemy_Theme.css  # Optional UI theme

⚙️ Deployment Options
Option 1: Run Locally with Docker
git clone https://github.com/your-user/RNA_SEQ_APP.git
cd RNA_SEQ_APP
bash run.sh
Open your browser at: http://localhost:8787 

Option 2: Run on HPC with Singularity / Apptainer
Build the container:

singularity build rna-seq.sif Singularity.def
# or with Apptainer
apptainer build rna-seq.sif Singularity.def

Run the app:

singularity run --bind $(pwd):/mnt rna-seq.sif
Remote access:

ssh -L 8080:localhost:8080 youruser@cluster

Open: http://localhost:8080

📊 UI Tabs & Outputs

Tab	Description
Differential Expression Results	Top 50 DEGs with annotations + RF importance
PCA Plot / UMAP Plot	Interactive sample projections
Volcano Plot	logFC vs -log10(adj.P) w/ tooltips
Heatmap	Z-scored clustered heatmap
Pathway Enrichment	Enrichr: All / Up / Down-regulated genes
Power Summary / Curve	pwr-based sample size power curves
Classification	Random Forest: ROC, metrics, predictions
Read Me	In-app embedded guidance text

Counts + Phenotype
        ↓
Differential Expression (limma-voom)
        ↓
Feature Selection (RF importance top 20)
        ↓
Classification (train/test split + performance metrics)
        ↓
Optional: Enrichment (All / Up / Down)
        ↓
Visualization (PCA / UMAP / Volcano / Heatmap)
        ↓
Statistical Power Estimation (pwr)

🧠 Tech Stack
R, Shiny, limma, edgeR, biomaRt, enrichR

randomForest, caret, pROC, pwr

plotly, ggplot2, umap

Docker + Singularity for portability

Error handling + custom CSS theme support



📥 Downloadable Outputs


Output	Description
DE_results_YYYY-MM-DD.csv	Final differentially expressed genes
enrichment_all.csv	All DEGs enriched pathways
enrichment_upregulated.csv	Upregulated gene enrichment
enrichment_downregulated.csv	Downregulated gene enrichment
rf_predictions.csv	RF classification output
rf_metrics.csv	Sensitivity, specificity, AUC
power_summary.csv	Power analysis summary

📧 Automated Error Monitoring (Optional)
All errors logged to error_log.txt

Schedule email_log.R via cron (Linux) or taskscheduleR (Windows)

Example Cron Job:

0 0 * * * Rscript /path/to/email_log.R

📢 Feedback / Support
🐛 GitHub Issues

💡 Contributions welcome! Fork, PR, or ideas anytime.


RNA-SEQ processing shouldn’t feel like black magic.
🧙‍♂️ JCAP RNA-SEQ ANALYZER turns your data into insights — fast, reproducibly, and interactively.



