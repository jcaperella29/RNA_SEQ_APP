
    

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
library(pheatmap)

# === Error Logging ===
log_file <- "error_log.txt"
options(shiny.error = function() {
  err <- geterrmessage()
  timestamp <- Sys.time()
  msg <- paste0("[", timestamp, "] ", err, "\n\n")
  cat(msg, file = log_file, append = TRUE)
})
#ui


    

      
      
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
      actionButton("plot_heatmap", "Plot Heatmap"),
      downloadButton("output", "Download DE Results"),
      downloadButton("download_heatmap", "Download Heatmap"),
      
      hr(),
      h4("Pathway Enrichment"),
      selectInput("enrich_db", "Select Pathway Database", choices = c(
        "GO_Biological_Process_2021",
        "KEGG_2021_Human",
        "WikiPathway_2021_Human",
        "Reactome_2022"
      ), selected = "KEGG_2021_Human"),
      actionButton("enrich_all_btn", "Enrich All DE Genes"),
      actionButton("enrich_up_btn", "Enrich Upregulated"),
      actionButton("enrich_down_btn", "Enrich Downregulated"),
      
      hr(),
      h4("Power Analysis"),
      actionButton("run_power", "Run Power Analysis"),
      sliderInput("effect_size", "Effect Size (Cohen's d / f)", min = 0.1, max = 1.5, value = 0.8, step = 0.1),
      downloadButton("download_power", "Download Power Summary"),
      selectInput("power_test_type", "Power Curve Type",
                  choices = c("t-test (2 groups)" = "ttest", "ANOVA (>2 groups)" = "anova")),
      sliderInput("curve_n_range", "Sample Size Range", min = 2, max = 100, value = c(2, 30)),
      actionButton("plot_power_curve", "Plot Power Curve"),
      
      hr(),
      h4("Classification"),
      actionButton("run_rf_classifier", "Run Random Forest Classifier")
    ),
    
    mainPanel(
      tabsetPanel(id = "main_tabset",
                  tabPanel("Differential Expression Results", DTOutput("Dif_expr_results")),
                  tabPanel("PCA Plot", plotlyOutput("pcaplot")),
                  tabPanel("UMAP Plot", plotlyOutput("umapplot")),
                  tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),
                  tabPanel("Heatmap", plotOutput("heatmap_plot", height = "800px")),
                  tabPanel("Read Me", verbatimTextOutput("readme")),
                  
                  tabPanel("Pathway Analysis",
                           tabsetPanel(
                             tabPanel("All DE Genes",
                                      tabsetPanel(
                                        tabPanel("Results",
                                                 downloadButton("download_enrich_all", "Download All DE Enrichment"),
                                                 DTOutput("enrich_all_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_all_plot"))
                                      )
                             ),
                             tabPanel("Upregulated",
                                      tabsetPanel(
                                        tabPanel("Results",
                                                 downloadButton("download_enrich_up", "Download Upregulated Enrichment"),
                                                 DTOutput("enrich_up_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_up_plot"))
                                      )
                             ),
                             tabPanel("Downregulated",
                                      tabsetPanel(
                                        tabPanel("Results",
                                                 downloadButton("download_enrich_down", "Download Downregulated Enrichment"),
                                                 DTOutput("enrich_down_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_down_plot"))
                                      )
                             )
                           )
                  ),
                  
                  tabPanel("Power Analysis Table", tableOutput("power_summary")),
                  tabPanel("Power Curve", plotlyOutput("power_curve_plot")),
                  
                  tabPanel("Classification",
                           tabsetPanel(
                             tabPanel("Predictions", DTOutput("rf_predictions_dt")),
                             tabPanel("Metrics", tableOutput("rf_metrics_table")),
                             tabPanel("ROC Curve", plotlyOutput("rf_roc_plot"))
                           )
                  )
      )
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  power_result_table <- reactiveVal()
  
  
  final_results <- reactiveVal()
  heatmap_matrix <- reactiveVal()
  
  log_error <- function(e, context = "unknown") {
    msg <- paste0("[", Sys.time(), "] [", context, "] ", conditionMessage(e), "\n")
    write(msg, file = "error_log.txt", append = TRUE)
  }
  
  # === UI Dynamic ===
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
  
  # === Download DE Results ===
  output$output <- downloadHandler(
    filename = function() {
      paste0("DE_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(final_results())
      write.csv(final_results(), file, row.names = FALSE)
    }
  )
  
  # === Differential Expression + Random Forest Pipeline ===
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
  
  # === PCA Plot ===
  observeEvent(input$plot_pca, {
    showNotification("Generating PCA plot...", type = "message")
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
        geom_point() +
        theme_minimal()
      
      output$pcaplot <- renderPlotly({
        ggplotly(p, tooltip = "text")
      })
      
      showNotification("PCA plot ready ✅", type = "default")
    }, error = function(e) {
      log_error(e, "PCA Plot")
      showNotification("Error generating PCA plot.", type = "error")
    })
  })
  
  # === UMAP Plot ===
  observeEvent(input$plot_umap, {
    tryCatch({
      req(input$counts_input, input$phenotype_input, input$phenotype_column)
      showNotification("Generating UMAP plot...", type = "message")
      
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
      n_neighbors <- 15
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
      showNotification("UMAP plot ready ✅", type = "default")
      
    }, error = function(e) { log_error(e, "UMAP Plot") })
  })
  
  # === ✅ FIXED Volcano Plot ===
  observeEvent(input$plot_volcano, {
    tryCatch({
      req(final_results())
      showNotification("Generating Volcano plot...", type = "message")
      
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
      
      showNotification("Volcano plot ready ✅", type = "default")
    }, error = function(e) {
      log_error(e, "Volcano Plot")
      showNotification("Error generating volcano plot.", type = "error")
    })
  })
 #heatmap
  
  observeEvent(input$plot_heatmap, {
    tryCatch({
      req(input$counts_input, final_results())
      showNotification("Generating Heatmap...", type = "message")
      
      # Load counts
      ext_counts <- tools::file_ext(input$counts_input$name)
      counts <- if (ext_counts == "csv") {
        read.csv(input$counts_input$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.table(input$counts_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
      }
      
      # Clean gene IDs
      rownames(counts) <- sub("\\..*", "", rownames(counts))
      final_df <- final_results()
      final_df$Ensembl_IDs <- sub("\\..*", "", final_df$Ensembl_IDs)
      
      # Match genes
      matching_genes <- intersect(rownames(counts), final_df$Ensembl_IDs)
      if (length(matching_genes) < 2) {
        showNotification("Not enough DE genes matched in count matrix.", type = "error")
        return()
      }
      
      counts <- counts[matching_genes, ]
      
      # Normalize
      norm_counts <- t(scale(t(as.matrix(counts)), center = TRUE, scale = TRUE))
      norm_counts[is.na(norm_counts)] <- 0
      
      # Match gene labels
      matched <- match(rownames(norm_counts), final_df$Ensembl_IDs)
      gene_labels <- final_df$hgnc_symbol[matched]
      
      # Fallbacks
      fallbacks <- is.na(gene_labels) | gene_labels == ""
      gene_labels[fallbacks] <- rownames(norm_counts)[fallbacks]
      
      # Assign rownames safely
      if (length(gene_labels) == nrow(norm_counts)) {
        rownames(norm_counts) <- gene_labels
      } else {
        warning("Gene label count mismatch — fallback to Ensembl IDs")
        rownames(norm_counts) <- rownames(norm_counts)
      }
      
      heatmap_matrix(norm_counts)
      
      output$heatmap_plot <- renderPlot({
        mat <- heatmap_matrix()
        if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
          plot.new()
          text(0.5, 0.5, "Insufficient data for heatmap", cex = 1.5)
        } else {
          p <- pheatmap::pheatmap(
            mat,
            show_rownames = TRUE,
            show_colnames = TRUE,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "complete",
            main = "Heatmap - Final DE Genes (HGNC Labels)",
            silent = TRUE
          )
          grid::grid.newpage()
          grid::grid.draw(p$gtable)
          
          # 🎉 NEW: Notify once rendered
          showNotification("Heatmap successfully generated ✅", type = "message")
        }
      })
      
      
    }, error = function(e) {
      log_error(e, "Heatmap")
      showNotification("Error generating heatmap.", type = "error")
    })
  })
  
    
    # === Heatmap Download ===
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0("heatmap_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(heatmap_matrix())
        png(file, width = 1200, height = 1000, res = 150)
        print(
          pheatmap::pheatmap(
            heatmap_matrix(),
            show_rownames = FALSE,
            show_colnames = TRUE,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "complete",
            main = "Heatmap - Final DE Genes"
          )
        )
        dev.off()
      }
    )
    # === ENRICHMENT UTILS ===
    perform_enrichment <- function(gene_list, db) {
      if (length(gene_list) < 2) return(NULL)
      enrichr(gene_list, databases = db)[[1]]
    }
    
    plot_enrichment_bar <- function(df, title) {
      top <- head(df[order(df$Adjusted.P.value), ], 10)
      top$Term <- factor(top$Term, levels = rev(top$Term))
      ggplot(top, aes(x = Term, y = -log10(Adjusted.P.value))) +
        geom_bar(stat = "identity", fill = "#2c7bb6") +
        coord_flip() +
        theme_minimal() +
        labs(title = title, y = "-log10 Adjusted P-value", x = "Pathway")
    }
    
    # === ENRICHMENT UTILS ===
    perform_enrichment <- function(gene_list, db) {
      if (length(gene_list) < 2) return(NULL)
      enrichr(gene_list, databases = db)[[1]]
    }
    
   
    
    plot_enrichment_bar <- function(df, title) {
      top <- head(df[order(df$Adjusted.P.value), ], 10)
      top$Term <- factor(top$Term, levels = rev(top$Term))
      
      plotly::plot_ly(
        data = top,
        x = ~-log10(Adjusted.P.value),
        y = ~Term,
        type = "bar",
        orientation = "h",
        hoverinfo = "text",
        text = ~paste0("P.adj: ", signif(Adjusted.P.value, 3),
                       "<br>Score: ", round(Combined.Score, 2))
      ) %>%
        layout(
          title = list(text = title),
          xaxis = list(title = "-log10 Adjusted P-value"),
          yaxis = list(title = ""),
          margin = list(l = 200)
        )
    }
    
    
    observeEvent(input$enrich_all_btn, {
      req(final_results())
      showNotification("Enriching all DE genes",type = "message")
      db <- input$enrich_db
      genes <- na.omit(final_results()$hgnc_symbol)
      res <- perform_enrichment(genes, db)
      enrich_all_res(res)
      showNotification("All DE Genes Enriched ✅", type = "message")
    })
    
    observeEvent(input$enrich_up_btn, {
      req(final_results())
      showNotification("Enriching upregulated DE genes",type = "message")
      
      db <- input$enrich_db
      df <- final_results()
      genes <- df$hgnc_symbol[df$logFC > 1 & df$adj.P.Val < 0.05]
      res <- perform_enrichment(genes, db)
      enrich_up_res(res)
      showNotification("Upregulated Genes Enriched ✅", type = "message")
    })
    
    enrich_all_res <- reactiveVal()
    enrich_up_res <- reactiveVal()
    enrich_down_res <- reactiveVal()
    
    
    observeEvent(input$enrich_down_btn, {
      req(final_results())
      
      showNotification("Enriching downregulated DE genes",type = "message")
      
      db <- input$enrich_db
      df <- final_results()
      genes <- df$hgnc_symbol[df$logFC < -1 & df$adj.P.Val < 0.05]
      res <- perform_enrichment(genes, db)
      enrich_down_res(res)
      showNotification("Downregulated Genes Enriched ✅", type = "message")
    })
    
    
    
    
    # === DT Tables ===
    output$enrich_all_dt <- renderDT({
      req(enrich_all_res())
      datatable(enrich_all_res()[, c("Term", "Adjusted.P.value", "Combined.Score")])
    })
    
    output$enrich_up_dt <- renderDT({
      req(enrich_up_res())
      datatable(enrich_up_res()[, c("Term", "Adjusted.P.value", "Combined.Score")])
    })
    
    output$enrich_down_dt <- renderDT({
      req(enrich_down_res())
      datatable(enrich_down_res()[, c("Term", "Adjusted.P.value", "Combined.Score")])
    })
    
    output$enrich_all_plot <- renderPlotly({
      req(enrich_all_res())
      showNotification("Enrichment plot (All) ready ✅", type = "default")
      plot_enrichment_bar(enrich_all_res(), "All DE Genes")
    })
    
    output$enrich_up_plot <- renderPlotly({
      req(enrich_up_res())
      showNotification("Enrichment plot (Upregulated) ready ✅", type = "default")
      plot_enrichment_bar(enrich_up_res(), "Upregulated Genes")
    })
    
    output$enrich_down_plot <- renderPlotly({
      req(enrich_down_res())
      showNotification("Enrichment plot (Downregulated) ready ✅", type = "default")
      plot_enrichment_bar(enrich_down_res(), "Downregulated Genes")
    })
    

    # === Download Handlers for Enrichment Results ===
    output$download_enrich_all <- downloadHandler(
      filename = function() paste0("enrichment_all_", Sys.Date(), ".csv"),
      content = function(file) {
        req(enrich_all_res())
        write.csv(enrich_all_res(), file, row.names = FALSE)
      }
    )
    
    output$download_enrich_up <- downloadHandler(
      filename = function() paste0("enrichment_upregulated_", Sys.Date(), ".csv"),
      content = function(file) {
        req(enrich_up_res())
        write.csv(enrich_up_res(), file, row.names = FALSE)
      }
    )
    
    output$download_enrich_down <- downloadHandler(
      filename = function() paste0("enrichment_downregulated_", Sys.Date(), ".csv"),
      content = function(file) {
        req(enrich_down_res())
        write.csv(enrich_down_res(), file, row.names = FALSE)
      }
    )
##power analysis
    observeEvent(input$run_power, {
      req(input$phenotype_input, input$phenotype_column)
      
      showNotification("Running power analysis...", type = "message")
      
      tryCatch({
        ext <- tools::file_ext(input$phenotype_input$name)
        pheno <- if (ext == "csv") {
          read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
        } else {
          read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
        }
        
        group_col <- input$phenotype_column
        groups <- as.factor(pheno[[group_col]])
        group_sizes <- table(groups)
        k <- length(group_sizes)
        effect_size <- input$effect_size
        sig <- 0.05
        
        result <- NULL
        
        if (k == 2) {
          n1 <- as.numeric(group_sizes[1])
          n2 <- as.numeric(group_sizes[2])
          power_res <- pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = effect_size, sig.level = sig)
          
          result <- data.frame(
            `Test` = "t-test",
            `Effect Size (d)` = effect_size,
            `Groups` = paste(names(group_sizes), collapse = ", "),
            `n1` = n1,
            `n2` = n2,
            `Significance Level` = sig,
            `Estimated Power` = round(power_res$power, 3),
            check.names = FALSE
          )
          
        } else if (k > 2) {
          total_n <- sum(group_sizes)
          power_res <- pwr::pwr.anova.test(k = k, n = total_n / k, f = effect_size, sig.level = sig)
          
          result <- data.frame(
            `Test` = "ANOVA",
            `Effect Size (f)` = effect_size,
            `Groups` = paste(names(group_sizes), collapse = ", "),
            `Samples per Group` = round(total_n / k),
            `Significance Level` = sig,
            `Estimated Power` = round(power_res$power, 3),
            check.names = FALSE
          )
        } else {
          result <- data.frame(Message = "⚠️ Not enough groups for power analysis")
        }
        
        power_result_table(result)  # save it for rendering and download
        showNotification("Power analysis complete ✅", type = "message")
        
      }, error = function(e) {
        power_result_table(data.frame(Error = e$message))
        showNotification("Power analysis failed ❌", type = "error")
      })
    })
    output$power_summary <- renderTable({
      req(power_result_table())
      power_result_table()
    })
    output$download_power <- downloadHandler(
      filename = function() {
        paste0("power_summary_", Sys.Date(), ".csv")
      },
      content = function(file) {
        df <- power_result_table()
        if (is.null(df)) {
          writeLines("No power analysis result available.", file)
        } else {
          write.csv(df, file, row.names = FALSE)
        }
      }
    )
    
    output$power_curve_plot <- renderPlotly({
      req(input$plot_power_curve)  # makes it reactive to button
      showNotification("generating power curve",type ="message")
      isolate({
        effect_size <- input$effect_size
        test_type <- input$power_test_type
        n_seq <- seq(input$curve_n_range[1], input$curve_n_range[2])
        sig <- 0.05
        
        power_vals <- sapply(n_seq, function(n) {
          if (test_type == "ttest") {
            # Equal n per group assumed
            pwr::pwr.t.test(n = n, d = effect_size, sig.level = sig, type = "two.sample")$power
          } else {
            # ANOVA assumes n per group = n
            pwr::pwr.anova.test(k = 3, n = n, f = effect_size, sig.level = sig)$power
          }
        })
        
        df <- data.frame(SampleSize = n_seq, Power = power_vals)
        showNotification("Power curve ready ✅", type = "message")  # 💬 Notification after plot
        
        
        plot_ly(df, x = ~SampleSize, y = ~Power, type = 'scatter', mode = 'lines+markers',
                line = list(color = "#00cc99", width = 3)) %>%
          layout(title = "Power Curve",
                 xaxis = list(title = "Sample Size (per group)"),
                 yaxis = list(title = "Power", range = c(0, 1)),
                 shapes = list(
                   list(type = "line", x0 = min(n_seq), x1 = max(n_seq),
                        y0 = 0.8, y1 = 0.8,
                        line = list(dash = 'dash', color = "red"))
                  
                 )) 
      })
    })
    observeEvent(input$run_rf_classifier, {
      tryCatch({
        req(final_results(), input$counts_input, input$phenotype_input, input$phenotype_column)
        showNotification("Running Random Forest classifier...", type = "message")
        
        # Load counts
        ext_counts <- tools::file_ext(input$counts_input$name)
        counts <- if (ext_counts == "csv") {
          read.csv(input$counts_input$datapath, row.names = 1, check.names = FALSE)
        } else {
          read.table(input$counts_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
        }
        
        # Load phenotype
        ext_pheno <- tools::file_ext(input$phenotype_input$name)
        pheno <- if (ext_pheno == "csv") {
          read.csv(input$phenotype_input$datapath, row.names = 1, check.names = FALSE)
        } else {
          read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1, check.names = FALSE)
        }
        
        pheno_vec <- as.factor(pheno[[input$phenotype_column]])
        
        # Filter DE genes
        de_genes <- final_results()$Ensembl_IDs
        counts <- counts[rownames(counts) %in% de_genes, ]
        counts <- t(counts)
        counts <- as.data.frame(counts)
        counts$Phenotype <- pheno_vec
        
        # Train/Test Split (70/30)
        set.seed(42)
        train_idx <- caret::createDataPartition(counts$Phenotype, p = 0.7, list = FALSE)
        train_data <- counts[train_idx, ]
        test_data <- counts[-train_idx, ]
        
        # Train RF model
        rf_model <- randomForest(Phenotype ~ ., data = train_data, ntree = 500, importance = TRUE)
        probs <- predict(rf_model, newdata = test_data, type = "prob")
        preds <- predict(rf_model, newdata = test_data)
        
        # Predictions table
        pred_table <- data.frame(
          Sample = rownames(test_data),
          Actual = test_data$Phenotype,
          Predicted = preds,
          Prob = apply(probs, 1, max),
          check.names = FALSE
        )
        
        # Metrics
        cm <- caret::confusionMatrix(preds, test_data$Phenotype)
        sens <- cm$byClass['Sensitivity']
        spec <- cm$byClass['Specificity']
        
        # ROC
        pheno_levels <- levels(test_data$Phenotype)
        if (length(pheno_levels) != 2) {
          roc_plot <- plotly::plot_ly() %>%
            layout(title = "ROC only supported for 2-class problems ❌")
          auc_val <- NA
        } else {
          roc_obj <- pROC::roc(test_data$Phenotype, probs[, pheno_levels[2]])
          auc_val <- round(roc_obj$auc, 3)
          roc_df <- data.frame(
            FPR = 1 - roc_obj$specificities,
            TPR = roc_obj$sensitivities
          )
          
          roc_plot <- plot_ly(
            data = roc_df,
            x = ~FPR,
            y = ~TPR,
            type = 'scatter',
            mode = 'lines',
            line = list(color = '#1f77b4', width = 2)
          ) %>%
            layout(
              title = paste("ROC Curve (AUC =", auc_val, ")"),
              xaxis = list(title = "False Positive Rate"),
              yaxis = list(title = "True Positive Rate"),
              showlegend = FALSE
            )
        }
        
        # Outputs
        output$rf_predictions_dt <- renderDT({
          datatable(pred_table)
        })
        
        output$rf_metrics_table <- renderTable({
          data.frame(
            Sensitivity = round(sens, 3),
            Specificity = round(spec, 3),
            `AUC (ROC)` = ifelse(is.na(auc_val), "N/A", auc_val)
          )
        }, rownames = FALSE)
        
        output$rf_roc_plot <- renderPlotly({
          roc_plot
        })
        
        showNotification("Random Forest classification (w/ split) done ✅", type = "message")
      }, error = function(e) {
        showNotification(paste("RF Classification Error:", e$message), type = "error")
      })
    })
    
