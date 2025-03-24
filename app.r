
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
      uiOutput("phenotype_column_selector"),
      actionButton("running_Dif", "Run Differential Expression"),
      actionButton("plot_pca", "PCA Plot"),
      actionButton("plot_umap", "UMAP Plot"),
      actionButton("plot_volcano", "Volcano Plot"),
      sliderInput("n_neighbors_umap", "UMAP Neighbors", min = 5, max = 50, value = 15),
      downloadButton("output", "Download DE Results")
    ),
    
    mainPanel(
      tabsetPanel(id = "main_tabset",
                  tabPanel("Differential Expression Results", DTOutput("Dif_expr_results")),
                  tabPanel("PCA Plot", plotlyOutput("pcaplot")),
                  tabPanel("UMAP Plot", plotlyOutput("umapplot")),
                  tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),
                  tabPanel("Read Me", verbatimTextOutput("readme"))
      )
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  
  final_results <- reactiveVal()
  
  log_error <- function(e, context = "unknown") {
    msg <- paste0("[", Sys.time(), "] [", context, "] ", conditionMessage(e), "\n")
    write(msg, file = "error_log.txt", append = TRUE)
  }
  
  # === UI: Dynamic phenotype column dropdown ===
  output$phenotype_column_selector <- renderUI({
    req(input$phenotype_input)
    ext <- tools::file_ext(input$phenotype_input$name)
    pheno <- if (ext == "csv") {
      read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
    } else {
      read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
    }
    selectInput("phenotype_column", "Choose Phenotype Column", choices = colnames(pheno))
  })
  
  # === Download Handler ===
  output$output <- downloadHandler(
    filename = function() {
      paste0("DE_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(final_results())
      write.csv(final_results(), file, row.names = FALSE)
    }
  )
  
  # === DE + RF Pipeline ===
  observeEvent(input$running_Dif, {
    tryCatch({
      withProgress(message = 'Running DE + Feature Selection...', value = 0, {
        req(input$counts_input, input$phenotype_input, input$phenotype_column)
        
        incProgress(0.05, detail = "Reading data...")
        ext_counts <- tools::file_ext(input$counts_input$name)
        counts <- if (ext_counts == "csv") {
          read.csv(input$counts_input$datapath, row.names = 1, check.names = FALSE)
        } else {
          read.table(input$counts_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
        }
        
        ext_pheno <- tools::file_ext(input$phenotype_input$name)
        pheno <- if (ext_pheno == "csv") {
          read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
        } else {
          read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
        }
        
        incProgress(0.10, detail = "Filtering low-expression genes...")
        counts_filtered <- counts[rowSums(counts) > 10, ]
        
        if (ncol(counts_filtered) >= 10) {
          sample_sums <- colSums(counts_filtered)
          q1 <- quantile(sample_sums, 0.25)
          q3 <- quantile(sample_sums, 0.75)
          iqr <- q3 - q1
          inlier_samples <- names(sample_sums[sample_sums > (q1 - 1.5 * iqr) & sample_sums < (q3 + 1.5 * iqr)])
          counts_final <- counts_filtered[, inlier_samples]
          pheno_final <- pheno[inlier_samples, , drop = FALSE]
        } else {
          counts_final <- counts_filtered
          pheno_final <- pheno[colnames(counts_final), , drop = FALSE]
        }
        
        incProgress(0.35, detail = "Running limma/voom DE analysis...")
        df <- data.frame(t(counts_final))
        df$phenotype <- as.factor(pheno_final[[input$phenotype_column]])
        
        dge <- DGEList(counts = t(df[, -ncol(df)]))
        dge <- dge[filterByExpr(dge, group = df$phenotype), , keep.lib.sizes = FALSE]
        dge <- calcNormFactors(dge)
        
        design <- model.matrix(~0 + df$phenotype)
        colnames(design) <- levels(df$phenotype)
        
        v <- voom(dge, design, normalize.method = "quantile")
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        
        res <- topTable(fit, coef = 1, number = Inf)
        res$Ensembl_IDs <- rownames(res)
        
        top50 <- head(res[order(res$adj.P.Val), ], 50)
        
        incProgress(0.55, detail = "Running Random Forest...")
        rf_data <- data.frame(t(counts_final[rownames(counts_final) %in% top50$Ensembl_IDs, ]))
        rf_data$phenotype <- as.factor(pheno_final[[input$phenotype_column]])
        
        rf_model <- randomForest(phenotype ~ ., data = rf_data, importance = TRUE, ntree = 500)
        rf_importance <- importance(rf_model, type = 1)
        rf_top <- sort(rf_importance[, 1], decreasing = TRUE)[1:20]
        
        final_res <- res[res$Ensembl_IDs %in% names(rf_top), c("adj.P.Val", "logFC", "Ensembl_IDs")]
        final_res$MeanDecreaseAccuracy <- rf_top[final_res$Ensembl_IDs]
        
        incProgress(0.75, detail = "Annotating genes via subprocess...")
        system2("Rscript", args = c("get_annotation.R", shQuote(paste(final_res$Ensembl_IDs, collapse = ","))))
        gene_info <- read.csv("annotation_output.csv", stringsAsFactors = FALSE)
        
        merged <- merge(final_res, gene_info, by.x = "Ensembl_IDs", by.y = "ensembl_gene_id", all.x = TRUE)
        merged <- merged[, c("Ensembl_IDs", "hgnc_symbol", "description", "adj.P.Val", "logFC", "MeanDecreaseAccuracy")]
        
        final_results(merged)
        output$Dif_expr_results <- renderDT(datatable(merged))
        incProgress(1, detail = "✅ Done.")
      })
    }, error = function(e) {
      log_error(e, "DE + RF Pipeline")
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # === PCA ===
  observeEvent(input$plot_pca, {
    tryCatch({
      req(input$counts_input, input$phenotype_input, input$phenotype_column)
      
      ext_counts <- tools::file_ext(input$counts_input$name)
      counts <- if (ext_counts == "csv") {
        read.csv(input$counts_input$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.table(input$counts_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
      }
      
      ext_pheno <- tools::file_ext(input$phenotype_input$name)
      pheno <- if (ext_pheno == "csv") {
        read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
      }
      
      df <- data.frame(t(counts))
      df$Phenotype <- as.factor(pheno[[input$phenotype_column]])
      
      df_vars <- df[, -ncol(df)]
      zero_var_cols <- apply(df_vars, 2, function(x) var(x) == 0)
      df_clean <- df_vars[, !zero_var_cols]
      
      pca <- prcomp(df_clean, scale. = TRUE)
      pca_df <- data.frame(pca$x, Sample = rownames(df), Phenotype = df$Phenotype)
      
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, text = Sample)) +
        geom_point() + theme_minimal()
      
      output$pcaplot <- renderPlotly({ ggplotly(p, tooltip = "text") })
    }, error = function(e) { log_error(e, "PCA Plot") })
  })
  
  # === UMAP ===
  observeEvent(input$plot_umap, {
    tryCatch({
      req(input$counts_input, input$phenotype_input, input$phenotype_column)
      
      ext_counts <- tools::file_ext(input$counts_input$name)
      counts <- if (ext_counts == "csv") {
        read.csv(input$counts_input$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.table(input$counts_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
      }
      
      ext_pheno <- tools::file_ext(input$phenotype_input$name)
      pheno <- if (ext_pheno == "csv") {
        read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
      }
      
      df <- data.frame(t(counts))
      df$Phenotype <- as.factor(pheno[[input$phenotype_column]])
      
      n_samples <- nrow(df)
      n_neighbors <- input$n_neighbors_umap
      if (n_neighbors >= n_samples) {
        n_neighbors <- max(2, n_samples - 1)
      }
      
      umap_config <- umap.defaults
      umap_config$n_neighbors <- n_neighbors
      
      umap_res <- umap(df[, -ncol(df)], config = umap_config)
      umap_df <- data.frame(umap_res$layout, Sample = rownames(df), Phenotype = df$Phenotype)
      
      p <- ggplot(umap_df, aes(x = X1, y = X2, color = Phenotype, text = Sample)) +
        geom_point() + theme_minimal()
      
      output$umapplot <- renderPlotly({ ggplotly(p, tooltip = "text") })
    }, error = function(e) { log_error(e, "UMAP Plot") })
  })
  
  # === Volcano Plot ===
  observeEvent(input$plot_volcano, {
    tryCatch({
      req(final_results())
      df <- final_results()
      
      df$log10p <- -log10(df$adj.P.Val)
      df$significant <- ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) > 1, "Significant", "Not Significant")
      
      p <- ggplot(df, aes(x = logFC, y = log10p, color = significant, text = hgnc_symbol)) +
        geom_point(alpha = 0.8) +
        scale_color_manual(values = c("gray70", "firebrick")) +
        labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value") +
        theme_minimal()
      
      output$volcano_plot <- renderPlotly({
        ggplotly(p, tooltip = "text")
      })
    }, error = function(e) {
      log_error(e, "Volcano Plot")
      showNotification("Error generating volcano plot.", type = "error")
    })
  })
  
  # === README Tab ===
  output$readme <- renderPrint({
    tryCatch({
      readLines("JCAP RNA_SEQ Readme.txt")
    }, error = function(e) { "README file not found." })
  })
}

# === Run the App ===
shinyApp(ui = ui, server = server)
