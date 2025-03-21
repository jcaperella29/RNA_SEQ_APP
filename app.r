# === RNA-SEQ APP: Differential Expression & Enrichment ===
library(shiny)
library(shinythemes)
library(shinyjs)
library(DT)
library(plotly)
library(limma)
library(edgeR)
library(biomaRt)
library(enrichR)
library(randomForest)
library(caret)
library(pROC)
library(pwr)
library(doParallel)
library(foreach)
library(umap)
library(ggplot2)
library(dplyr)

# === Error Logging ===
log_file <- "error_log.txt"
options(shiny.error = function() {
  err <- geterrmessage()
  timestamp <- Sys.time()
  msg <- paste0("[", timestamp, "] ", err, "\n\n")
  cat(msg, file = log_file, append = TRUE)
})

# === UI ===
ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("cyborg"),
  titlePanel("JCAP RNA-SEQ Analyzer"),

  sidebarLayout(
    sidebarPanel(
      fileInput("counts_input", "Input Counts", accept = c(".csv", ".txt")),
      fileInput("phenotype_input", "Input Phenotype Data", accept = c(".csv", ".txt")),
      actionButton("running_Dif", "Run Differential Expression"),
      actionButton("plot_pca", "PCA Plot"),
      actionButton("plot_umap", "UMAP Plot"),
      actionButton("plot_volcano", "Volcano Plot"),
      actionButton("run_rf", "Run Random Forest"),
      actionButton("pathway_all", "Enrich All Genes"),
      actionButton("pathway_up", "Enrich Upregulated"),
      actionButton("pathway_down", "Enrich Downregulated"),
      actionButton("calculate_power", "Calculate Power"),
      sliderInput("n_components_pca", "PCA Components", min = 2, max = 10, value = 2),
      sliderInput("n_neighbors_umap", "UMAP Neighbors", min = 5, max = 50, value = 15),
      sliderInput("p_value_threshold", "P-Value Threshold", min = 0, max = 0.1, value = 0.05, step = 0.001),
      sliderInput("log_fc_threshold", "Log Fold Change Threshold", min = 0, max = 5, value = 1, step = 0.1),
      downloadButton("output", "Download DE Results"),
      downloadButton("All_genes", "Download All Enrichment"),
      downloadButton("Up_genes", "Download Up Enrichment"),
      downloadButton("Down_genes", "Download Down Enrichment"),
      downloadButton("download_metrics", "Download RF Metrics")
    ),

    mainPanel(
      tabsetPanel(id = "main_tabset",
        tabPanel("Differential Expression Results", DTOutput("Dif_expr_results")),
        tabPanel("PCA Plot", plotlyOutput("pcaplot")),
        tabPanel("UMAP Plot", plotlyOutput("umapplot")),
        tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),
        tabPanel("Pathway Enrichment",
          tabsetPanel(
            tabPanel("All Genes", DTOutput("enrich_table_all")),
            tabPanel("Upregulated Genes", DTOutput("enrich_table_up")),
            tabPanel("Downregulated Genes", DTOutput("enrich_table_down"))
          )
        ),
        tabPanel("Random Forest Performance Metrics", DTOutput("rf_metrics")),
        tabPanel("Power Calculation", verbatimTextOutput("power_results")),
        tabPanel("Read Me", verbatimTextOutput("readme"))
      )
    )
  )
)

# === Secure File Reading Helper ===
read_uploaded_file <- function(file_input) {
  req(file_input)
  ext <- tools::file_ext(file_input$name)
  if (!(ext %in% c("csv", "txt"))) {
    stop("Unsupported file format. Only .csv or .txt allowed.")
  }
  data <- if (ext == "csv") {
    read.csv(file_input$datapath, row.names = 1, check.names = FALSE)
  } else {
    read.table(file_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
  }
  validate(
    need(nrow(data) > 1, "File must contain more than one row."),
    need(ncol(data) > 1, "File must contain more than one column."),
    need(all(sapply(data, is.numeric)), "All values must be numeric.")
  )
  return(data)
}

# === Run the App ===
shinyApp(ui = ui, server = server)
