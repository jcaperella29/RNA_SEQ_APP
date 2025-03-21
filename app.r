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
# === Server Logic ===
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)

  final_results <- reactiveVal()
  enrich_res_up <- reactiveVal()
  enrich_res_down <- reactiveVal()
  rf_metrics <- reactiveVal()
  combined_res_all <- NULL

  enrichr_dbs <- c("KEGG_2021_Human", "GO Molecular Function 2023", "GO Biological Process 2023", "GWAS Catalog 2023", "UK Biobank GWAS v1")

  # ==== Differential Expression Analysis ====
  observeEvent(input$running_Dif, {
    tryCatch({
      counts <- if (grepl("\\.csv$", input$counts_input$name)) read.csv(input$counts_input$datapath, row.names = 1) else read.table(input$counts_input$datapath, header = TRUE, row.names = 1)
      pheno <- if (grepl("\\.csv$", input$phenotype_input$name)) read.csv(input$phenotype_input$datapath, row.names = 1) else read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1)

      pheno_fact <- as.factor(pheno[,1])
      df <- data.frame(t(counts))
      df$phenotype <- pheno_fact

      df <- df[rowSums(df[, -ncol(df)]) > 10, ]
      df$phenotype <- as.factor(df$phenotype)

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
      
      # Map to gene symbols
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
      gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                         filters = "ensembl_gene_id",
                         values = top50$Ensembl_IDs,
                         mart = mart)
      merged <- merge(top50, gene_info, by.x = "Ensembl_IDs", by.y = "ensembl_gene_id", all.x = TRUE)
      final_results(merged)
      output$Dif_expr_results <- renderDT(datatable(merged))
    }, error = function(e) {
      cat(e$message, file = "error_log.txt", append = TRUE)
    })
  })

  # ==== PCA ====
  observeEvent(input$plot_pca, {
    req(input$counts_input, input$phenotype_input)
    counts <- read.csv(input$counts_input$datapath, row.names = 1)
    pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)
    df <- data.frame(t(counts))
    df$Phenotype <- as.factor(pheno[,1])
    df$Sample <- rownames(df)

    df_numeric <- df[, sapply(df, is.numeric)]
    df_numeric$Phenotype <- df$Phenotype
    df_numeric$Sample <- df$Sample

    pca <- prcomp(df_numeric[, !colnames(df_numeric) %in% c("Phenotype", "Sample")], scale. = TRUE)
    pca_df <- data.frame(pca$x)
    pca_df$Sample <- rownames(pca_df)
    pca_df$Phenotype <- df_numeric$Phenotype[match(rownames(pca_df), rownames(df_numeric))]

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, text = Sample)) +
      geom_point() + theme_minimal()
    output$pcaplot <- renderPlotly({ ggplotly(p, tooltip = "text") })
  })

  # ==== Volcano Plot ====
  observeEvent(input$plot_volcano, {
    req(final_results())
    res <- final_results()
    p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), text = res$hgnc_symbol)) +
      geom_point(aes(color = adj.P.Val < input$p_value_threshold & abs(logFC) > input$log_fc_threshold)) +
      scale_color_manual(values = c("black", "red")) + theme_minimal()
    output$volcano_plot <- renderPlotly({ ggplotly(p, tooltip = "text") })
  })


    # ==== UMAP ====
  observeEvent(input$plot_umap, {
    req(input$counts_input, input$phenotype_input)
    counts <- read.csv(input$counts_input$datapath, row.names = 1)
    pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)

    df <- data.frame(t(counts))
    df$phenotype <- as.factor(pheno[,1])
    df$Sample <- rownames(df)

    n_neighbors <- min(input$n_neighbors_umap, nrow(df) - 1)
    umap_res <- umap(df[, !(colnames(df) %in% c("phenotype", "Sample"))], n_neighbors = n_neighbors)

    umap_df <- data.frame(umap_res$layout)
    umap_df$Sample <- rownames(umap_df)
    umap_df$Phenotype <- df$phenotype

    p <- ggplot(umap_df, aes(x = X1, y = X2, color = Phenotype, text = Sample)) +
      geom_point() + theme_minimal()

    output$umapplot <- renderPlotly({ ggplotly(p, tooltip = "text") })
  })

  # ==== Enrichment Helpers ====
  enrich_genes <- function(genes) {
    enrichr(genes, enrichr_dbs)
  }

  combine_enrichr_results <- function(enrichr_results) {
    combined_results <- do.call(rbind, lapply(names(enrichr_results), function(db) {
      df <- enrichr_results[[db]]
      if (nrow(df) > 0) {
        df$Database <- db
        return(df)
      }
      NULL
    }))
    as.data.frame(combined_results)
  }

  enrich_all_genes <- function(res) {
    genes <- res$hgnc_symbol
    combined_res_all <<- combine_enrichr_results(enrich_genes(genes))
    combined_res_all
  }

  enrich_up_genes <- function(res) {
    up_genes <- res$hgnc_symbol[res$logFC > 0]
    if (length(up_genes) == 0) return(data.frame(Message = "No upregulated genes found"))
    combine_enrichr_results(enrich_genes(up_genes))
  }

  enrich_down_genes <- function(res) {
    down_genes <- res$hgnc_symbol[res$logFC < 0]
    if (length(down_genes) == 0) return(data.frame(Message = "No downregulated genes found"))
    combine_enrichr_results(enrich_genes(down_genes))
  }

  # ==== Enrichment Observers ====
  observeEvent(input$pathway_all, {
    req(final_results())
    output$enrich_table_all <- renderDT({
      datatable(enrich_all_genes(final_results()))
    })
  })

  observeEvent(input$pathway_up, {
    req(final_results())
    up <- enrich_up_genes(final_results())
    enrich_res_up(up)
    output$enrich_table_up <- renderDT({
      if ("Message" %in% colnames(up)) datatable(data.frame(Message = "No upregulated genes found."), options = list(dom = 't'))
      else datatable(up)
    })
  })

  observeEvent(input$pathway_down, {
    req(final_results())
    down <- enrich_down_genes(final_results())
    enrich_res_down(down)
    output$enrich_table_down <- renderDT({
      if ("Message" %in% colnames(down)) datatable(data.frame(Message = "No downregulated genes found."), options = list(dom = 't'))
      else datatable(down)
    })
  })

  # ==== Random Forest ====
  observeEvent(input$run_rf, {
    req(final_results(), input$counts_input, input$phenotype_input)
    selected_genes <- final_results()[, "Ensembl_IDs"]
    counts <- read.csv(input$counts_input$datapath, row.names = 1)
    pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)

    df <- data.frame(t(counts))
    df$phenotype <- as.factor(pheno[, 1])
    df <- df[, c(intersect(colnames(df), selected_genes), "phenotype")]

    train_index <- createDataPartition(df$phenotype, p = 0.7, list = FALSE)
    train_data <- df[train_index, ]
    test_data <- df[-train_index, ]

    model <- randomForest(
      x = train_data[, -which(names(train_data) == "phenotype")],
      y = train_data$phenotype,
      importance = TRUE,
      ntree = 500,
      mtry = 3
    )

    predictions <- predict(model, newdata = test_data[, -which(names(test_data) == "phenotype")])
    cm <- confusionMatrix(predictions, test_data$phenotype)
    auc_val <- roc(as.numeric(test_data$phenotype), as.numeric(predictions))$auc

    metrics <- data.frame(
      Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC", "Prevalence"),
      Value = c(
        cm$overall["Accuracy"],
        cm$byClass["Sensitivity"],
        cm$byClass["Specificity"],
        auc_val,
        mean(test_data$phenotype == levels(test_data$phenotype)[1])
      )
    )
    rf_metrics(metrics)
    output$rf_metrics <- renderDT(datatable(metrics))
  })
  # ==== Power Calculation ====
  observeEvent(input$calculate_power, {
    req(input$phenotype_input)
    pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)
    class_counts <- table(pheno[, 1])
    k <- length(class_counts)

    power_result <- if (k == 2) {
      n1 <- class_counts[1]
      n2 <- class_counts[2]
      pwr.t.test(n = min(n1, n2), d = 0.5, sig.level = 0.05, power = NULL, type = "two.sample")
    } else {
      n_harm <- 1 / mean(1 / class_counts)
      pwr.anova.test(k = k, n = n_harm, f = 0.5, sig.level = 0.05, power = NULL)
    }

    output$power_results <- renderPrint({ print(power_result) })
  })

  # ==== Download Handlers ====
  output$output <- downloadHandler(
    filename = function() "differential_expression_results.csv",
    content = function(file) {
      write.csv(final_results(), file, row.names = FALSE)
    }
  )

  output$All_genes <- downloadHandler(
    filename = function() "enrichment_all_genes.csv",
    content = function(file) {
      write.csv(combined_res_all, file, row.names = FALSE)
    }
  )

  output$Up_genes <- downloadHandler(
    filename = function() "enrichment_upregulated_genes.csv",
    content = function(file) {
      up_data <- enrich_res_up()
      if (!is.null(up_data)) write.csv(up_data, file, row.names = FALSE)
    }
  )

  output$Down_genes <- downloadHandler(
    filename = function() "enrichment_downregulated_genes.csv",
    content = function(file) {
      down_data <- enrich_res_down()
      if (!is.null(down_data)) write.csv(down_data, file, row.names = FALSE)
    }
  )

  output$download_metrics <- downloadHandler(
    filename = function() "rf_performance_metrics.csv",
    content = function(file) {
      write.csv(rf_metrics(), file, row.names = FALSE)
    }
  )

  # ==== README Tab ====
  output$readme <- renderPrint({
    tryCatch({
      readme_path <- normalizePath("JCAP RNA_SEQ Readme.txt", mustWork = TRUE)
      readLines(readme_path)
    }, error = function(e) {
      "README file not found."
    })
  })
  
}# Inside server <- function(...)
session$onSessionEnded(function() {
  tryCatch({
    message("📬 Triggering email_log.R...")
    source("email_log.R")
  }, error = function(e) {
    cat("[EMAIL_LOG ERROR] ", e$message, "\n", file = "error_log.txt", append = TRUE)
  })
})


# ==== Run the App ====
shinyApp(ui = ui, server = server)

  
