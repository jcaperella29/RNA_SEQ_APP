FROM rocker/r-ver:4.3.1

# ⛔ Prevents interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV MAKEFLAGS="-j4"

# 🛠️ Install system libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libxt-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libncurses-dev \
    libicu-dev \
    pandoc \
    curl \
    wget && apt-get clean

# 📦 Install CRAN packages
RUN R -e "install.packages(c( \
  'shiny', 'shinyjs', 'shinythemes', 'DT', 'plotly', 'pROC', 'pwr', \
  'ggplot2', 'ggfortify', 'caret', 'randomForest', 'varImp', \
  'dplyr', 'tidyverse', 'foreach', 'doParallel', 'umap' \
), repos='https://cloud.r-project.org')"

# 📦 Install Bioconductor + enrichR
RUN R -e "install.packages('enrichR', repos='https://cloud.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(version='3.18', ask=FALSE)"
RUN R -e "BiocManager::install(c('limma', 'edgeR', 'biomaRt'), ask=FALSE, update=FALSE, dependencies=TRUE)"

# 📁 Copy everything into /app
COPY . /app
WORKDIR /app

# 🌐 Expose port for Cloud Run/local Docker
EXPOSE 8080

# 🚀 Start Shiny app manually
CMD ["Rscript", "-e", "options(shiny.port=as.integer(Sys.getenv('PORT')), shiny.host='0.0.0.0'); shiny::runApp('/app')"]
