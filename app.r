
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
      downloadButton("download_heatmap", "Download Heatmap"),  # ✅ This line was missing a comma before the next group
      
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
      actionButton("enrich_down_btn", "Enrich Downregulated")
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
                                        tabPanel("Results", DTOutput("enrich_all_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_all_plot"))  # ✅ Changed to plotlyOutput
                                      )
                             ),
                             tabPanel("Upregulated",
                                      tabsetPanel(
                                        tabPanel("Results", DTOutput("enrich_up_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_up_plot"))  # ✅ Changed to plotlyOutput
                                      )
                             ),
                             tabPanel("Downregulated",
                                      tabsetPanel(
                                        tabPanel("Results", DTOutput("enrich_down_dt")),
                                        tabPanel("Barplot", plotlyOutput("enrich_down_plot"))  # ✅ Changed to plotlyOutput
                                      )
                             )
                           )
                  )
      )
    )
  )
)


# === Server ===
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  
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
    
    
    
} # <- closes server function

# === Launch App ===
shinyApp(server = server, ui = ui)
