# 🧬 RNA_SEQ_APP: Differential Expression + Enrichment + ML

An interactive **R Shiny app** for RNA-seq data analysis — combining statistical testing, power estimation, feature selection (Random Forest), dimensionality reduction (PCA/UMAP), and functional enrichment (Enrichr).

No coding required. Built for **bio labs, genomics teams, and HPC environments**.

---

## 🚀 Key Features

- 📂 Upload count + phenotype data (CSV or TXT)
- 🧬 Differential expression (voom + limma)
- 🧠 Pathway enrichment (via `enrichR`)
- 🧪 Power analysis (via `pwr`)
- 🌌 PCA + UMAP plots (interactive via `plotly`)
- 🌋 Volcano plots (logFC vs adj.P.Val)
- 🌲 Random Forest classification + performance metrics
- 🛠️ Error logging system + scheduled email alerts
- 📦 Dockerized for local deployment
- 🖥️ HPC-ready via **Singularity / Apptainer**
- 📄 README tab embedded in-app

---

## 📁 Folder Structure

RNA_SEQ_APP/ ├── app.R # Main Shiny app ├── run.sh # 🔁 Docker launcher ├── Dockerfile # 🐳 Full container definition ├── Singularity.def # 🧬 HPC build file ├── JCAP RNA_SEQ Readme.txt # Embedded in-app ├── error_log.txt # Auto-generated error logs ├── email_log.R # 📨 Scheduled log emailer ├── www/ │ └── custom.css # Optional UI styles


---

## 🧪 Option 1: Run Locally with Docker

### 🐳 Quick Start

```bash
git clone https://github.com/your-user/RNA_SEQ_APP.git
cd RNA_SEQ_APP
bash run.sh

Then open: http://localhost:8787

run.sh builds + runs a Docker container mapped to your system port.


🧠 Option 2: Run on HPC via Singularity

📦 Build

singularity build rna-seq.sif Singularity.def

with Apptainer:

apptainer build rna-seq.sif Singularity.def

🚀 Run


singularity run --bind $(pwd):/mnt rna-seq.sif

Then port-forward if remote:

ssh -L 8080:localhost:8080 youruser@cluster

Open: http://localhost:8080

📧 Error Logging + Email Monitoring
All server-side exceptions are logged to error_log.txt.

🛠 Schedule email_log.R via cron or taskscheduleR
Emails logs daily at midnight
Archives old logs into logs/
Internal-only; users are never exposed

Example (Linux cron):

0 0 * * * Rscript /path/to/RNA_SEQ_APP/email_log.R

📊 Outputs & Tabs


📊 Outputs & Tabs
Tab	Description
Differential Expression Results	Shows the top 50 differentially expressed genes using voom + limma. These genes are annotated with gene symbols and descriptions (biomaRt) and are used for downstream feature selection and modeling.
PCA Plot	Interactive PCA plot based on the expression matrix; samples colored by phenotype
UMAP Plot	Interactive UMAP using the expression data; number of neighbors is tunable
Volcano Plot	Log2 Fold Change vs Adjusted P-value plot with threshold sliders and tooltips
Pathway Enrichment	Enrichment results from Enrichr using:
All top 20 genes selected via RF
Upregulated genes
Downregulated genes | | Random Forest Performance Metrics | A second Random Forest model is trained using only the top 20 genes (selected via importance ranking from the top 50 DEGs). This tab shows the performance of that classifier: accuracy, sensitivity, specificity, AUC, and prevalence. | | Power Calculation | Uses pwr to calculate statistical power based on group sizes and desired effect size | | Read Me | Live-rendered content from  JCAP RNA_SEQ Readme.txt, embedded directly into the app |

#General processing  workflow
Counts + Phenotype
      ↓
 Differential Expression (top 50 genes)
      ↓
  Feature Selection (RF importance → top 20)
      ↓
   New RF model trained on top 20 genes
      ↓
   Performance Evaluation (Metrics Tab)
      ↓
    Optional Pathway Enrichment 



