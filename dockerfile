FROM rocker/r-ver:4.3.1

ENV DEBIAN_FRONTEND=noninteractive
ENV MAKEFLAGS="-j4"
ENV PORT=8080

RUN apt-get update && apt-get install -y \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    libgit2-dev libxt-dev libpng-dev libjpeg-dev libfreetype6-dev \
    libbz2-dev liblzma-dev libz-dev libncurses-dev libicu-dev \
    pandoc curl wget && apt-get clean

RUN R -e "install.packages(c( \
  'shiny', 'shinyjs', 'shinythemes', 'DT', 'plotly', 'pROC', 'pwr', \
  'ggplot2', 'caret', 'randomForest', 'dplyr', \
  'foreach', 'doParallel', 'umap', 'enrichR' \
), repos='https://cloud.r-project.org')"

RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
          BiocManager::install(version='3.18', ask=FALSE); \
          BiocManager::install(c('limma', 'edgeR', 'biomaRt'), ask=FALSE, update=FALSE, dependencies=TRUE)"

COPY . /app
WORKDIR /app

EXPOSE 8080

CMD ["Rscript", "-e", "options(shiny.port=as.integer(Sys.getenv('PORT')), shiny.host='0.0.0.0'); shiny::runApp('/app')"]
