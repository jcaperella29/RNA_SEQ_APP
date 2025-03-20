# Load necessary libraries
library(tidyverse)
library(limma)
library(randomForest)
library(varImp)
library(caret)
library(FactoMineR)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(shiny)
library(shinythemes)
library(DT)
library(umap)
library(doParallel)
library(foreach)
library(enrichR)
library(plotly)
library(edgeR)
library(shinyjs)
library(pROC)  # For AUC calculation
library(pwr)
# Power Calculation Function for Two or More Classes
calculate_power <- function(phenotype_data, effect_size = 0.5, alpha = 0.05) {
  # Count the number of samples in each class
  class_counts <- table(phenotype_data)
  
  # Get the number of groups
  k <- length(class_counts)
  
  if (k == 2) {
    # For two classes, perform a two-sample t-test power calculation
    n_group1 <- class_counts[1]
    n_group2 <- class_counts[2]
    
    # Perform power calculation for a two-sample t-test
    result <- pwr.t.test(n = min(n_group1, n_group2), d = effect_size, sig.level = alpha, power = NULL, type = "two.sample")
    
  } else {
    # For more than two classes, perform a one-way ANOVA power calculation
    n_harmonic <- 1 / mean(1 / class_counts)
    
    # Perform power calculation for a one-way ANOVA
    result <- pwr.anova.test(k = k, n = n_harmonic, f = effect_size, sig.level = alpha, power = NULL)
  }
  
  return(result)
}


# Required Libraries
library(randomForest)
library(caret)
library(pROC)

# Regularized Random Forest Classification Function
regularized_random_forest <- function(counts_data, phenotype_data, selected_genes) {
  set.seed(1)
  
  # Load expression data
  if (grepl("\\.csv$", counts_data)) {
    Exprmat <- read.csv(counts_data, row.names = 1)
  } else {
    Exprmat <- read.table(counts_data, header = TRUE, row.names = 1)
  }
  
  # Load phenotype data
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  # Prepare data frame
  df <- data.frame(t(Exprmat))
  df$phenotype <- as.factor(disease[, 1])
  
  # Ensure selected_genes exist in the data
  selected_genes <- intersect(selected_genes, colnames(df))
  if (length(selected_genes) == 0) {
    stop("None of the selected genes are found in the expression data.")
  }
  
  # Subset the data to include only selected genes and phenotype
  df_subset <- df[, c(selected_genes, "phenotype")]
  
  # Split the data into training and testing sets
  set.seed(123)
  train_index <- createDataPartition(df_subset$phenotype, p = 0.7, list = FALSE)
  train_data <- df_subset[train_index, ]
  test_data <- df_subset[-train_index, ]
  
  # Train Random Forest with regularization parameters
  model <- randomForest(
    x = train_data[, -which(colnames(train_data) == "phenotype")], 
    y = train_data$phenotype,
    importance = TRUE,
    ntree = 500,          # Number of trees
    mtry = 3,             # Number of features considered at each split
    nodesize = 5          # Minimum size of terminal nodes
  )
  
  # Make predictions
  predictions <- predict(model, newdata = test_data[, -which(colnames(test_data) == "phenotype")])
  
  # Evaluate model performance
  confusion_matrix <- confusionMatrix(predictions, test_data$phenotype)
  
  accuracy <- confusion_matrix$overall['Accuracy']
  sensitivity <- confusion_matrix$byClass['Sensitivity']
  specificity <- confusion_matrix$byClass['Specificity']
  
  # Calculate AUC
  auc <- roc(response = as.numeric(test_data$phenotype), predictor = as.numeric(predictions))$auc
  
  prevalence <- mean(test_data$phenotype == levels(test_data$phenotype)[1])
  
  metrics <- data.frame(
    Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC", "Prevalence"),
    Value = c(accuracy, sensitivity, specificity, auc, prevalence)
  )
  
  return(metrics)
}



# Define EnrichR databases to use
enrichr_dbs <- c("KEGG_2021_Human", "GO Molecular Function 2023", "GO Biological Process 2023", "GWAS Catalog 2023", "UK Biobank GWAS v1")

# Remove Outlier Samples Function
remove_outliers <- function(df) {
  df$phenotype <- as.numeric(df$phenotype)  # Convert phenotype to numeric
  sample_sums <- rowSums(df[, -ncol(df)])  # Exclude the phenotype column
  
  Q1 <- quantile(sample_sums, 0.25)
  Q3 <- quantile(sample_sums, 0.75)
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  non_outliers <- rownames(df)[sample_sums >= lower_bound & sample_sums <= upper_bound]
  
  df_clean <- df[rownames(df) %in% non_outliers, ]
  df_clean <- as.data.frame(df_clean)  # Ensure df_clean is a data frame
  
  df_clean$phenotype <- as.factor(df_clean$phenotype)  # Convert phenotype back to factor
  
  return(df_clean)
}

# Differential Expression Analysis Function with Debugging
dif_expr <- function(count_data, phenotype_data) {
  set.seed(1)
  if (grepl("\\.csv$", count_data)) {
    Countsmat <- read.csv(count_data, row.names = 1)
  } else {
    Countsmat <- read.table(count_data, header = TRUE, row.names = 1)
  }
  
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  disease_fact <- as.factor(disease[, 1])
  df <- data.frame(t(Countsmat))
  df$phenotype <- disease_fact
  
  if (nrow(df) != length(disease_fact)) {
    stop("The number of samples in count_data and phenotype_data do not match.")
  }
  
  df <- remove_outliers(df)
  
  good_stuff <- Countsmat[rowSums(Countsmat) > 10, ]
  good_stuff$goodgenes <- row.names(good_stuff)
  
  df <- df[, colnames(df) %in% c(good_stuff$goodgenes, "phenotype")]
  
  if (!"phenotype" %in% colnames(df)) {
    stop("Phenotype column is missing after subsetting.")
  }
  
  counts <- t(df[, colnames(df) != "phenotype"])
  
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = df$phenotype)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~0 + df$phenotype)
  colnames(design) <- levels(df$phenotype)
  
  v <- voom(dge, design, plot = FALSE, normalize.method = "quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  initial_results <- topTable(fit, coef = 1, number = Inf, genelist = fit$genes, adjust.method = "BH",
                              sort.by = "B", resort.by = NULL, p.value = 0.05, lfc = 0, confint = FALSE)
  initial_results$Ensembl_IDs <- row.names(initial_results)
  
  initial_results <- initial_results[, c("Ensembl_IDs", "adj.P.Val",  "logFC"  )]
  initial_results <- initial_results[order(initial_results$adj.P.Val), ]
  top50_results <- head(initial_results, 50)
  
  df <- df[, colnames(df) %in% c(top50_results$Ensembl_IDs, "phenotype")]
  model <- randomForest(phenotype ~ ., data = df, importance = TRUE)
  
  importance_df <- data.frame(importance(model))
  importance_df$genenames <- row.names(importance_df)
  importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  Top20_RFgenes <- importance_df$genenames[1:20]
  top20_results <- top50_results[top50_results$Ensembl_IDs %in% Top20_RFgenes, ]
  
  retry_getBM <- function(genes, mart, retries = 5) {
    for (i in 1:retries) {
      tryCatch({
        return(getBM(filters = "ensembl_gene_id", 
                     attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), 
                     values = genes, mart = mart, bmHeader = TRUE, useCache = FALSE))
      }, error = function(e) {
        if (i == retries) stop(e)
        Sys.sleep(10)
      })
    }
  }
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- top20_results$Ensembl_IDs
  G_list <- retry_getBM(genes, mart)
  
  final_results <- merge(top20_results, G_list, by.x = "Ensembl_IDs", by.y = "Gene stable ID", all.x = TRUE)
  
  colnames(final_results)[5] <- "Protein Name" 
  
  return(final_results)
}

  



# Pathway Analysis Functions with Combined Results
enrich_genes <- function(genes) {
  if (length(genes) == 0) {
    stop("No genes provided for enrichment analysis.")
  }
  enrichr(genes, enrichr_dbs)
}

enrich_all_genes <- function(final_results) {
  genes <- final_results[,4]
  print("Enriching all genes:")
  print(genes)
  enrich_res <- enrich_genes(genes)
  combined_res_all <<- combine_enrichr_results(enrich_res)
  return(combined_res_all)
}
enrich_up_genes <- function(final_results) {
  up_genes <- final_results[,4][final_results$logFC > 0]
  if (length(up_genes) == 0) {
    return(data.frame(Message = "No upregulated genes found"))
  } else {
    enrich_res <- enrich_genes(up_genes)
    combined_res_up <<- combine_enrichr_results(enrich_res)
    return(combined_res_up)
  }
}

enrich_down_genes <- function(final_results) {
  down_genes <- final_results[,4][final_results$logFC < 0]
  if (length(down_genes) == 0) {
    return(data.frame(Message = "No downregulated genes found"))
  } else {
    enrich_res <- enrich_genes(down_genes)
    combined_res_down <<- combine_enrichr_results(enrich_res)
    return(combined_res_down)
  }
}


combine_enrichr_results <- function(enrichr_results) {
  combined_results <- do.call(rbind, lapply(names(enrichr_results), function(db) {
    df <- enrichr_results[[db]]
    if (nrow(df) > 0) {
      df$Database <- db
      return(df)
    } else {
      return(NULL)
    }
  }))
  return(as.data.frame(combined_results))
}

Make_PCA <- function(count_data, phenotype_data, n_components) {
  # Read data
  if (grepl("\\.csv$", count_data)) {
    Countsmat <- read.csv(count_data, row.names = 1)
  } else {
    Countsmat <- read.table(count_data, header = TRUE, row.names = 1)
  }
  
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  disease_fact <- as.factor(disease[, 1])
  
  # Transpose and prepare the data frame
  df <- data.frame(t(Countsmat))
  df$phenotype <- disease_fact
  df$Sample <- rownames(df)
  
  # Debugging prints
  print("Column names of the initial data frame:")
  print(colnames(df))
  
  # Keep only numeric columns for PCA
  numeric_cols <- sapply(df, is.numeric)
  df_numeric <- df[, numeric_cols]
  
  # Debugging prints
  print("Column names of the numeric data frame:")
  print(colnames(df_numeric))
  
  # Remove columns with zero variance
  zero_var_columns <- sapply(df_numeric, function(x) var(x, na.rm = TRUE) == 0)
  print("Columns with zero variance:")
  print(names(zero_var_columns)[zero_var_columns])
  
  df_numeric <- df_numeric[, !zero_var_columns]
  
  # Ensure 'phenotype' and 'Sample' columns are not removed
  df_numeric$phenotype <- disease_fact
  df_numeric$Sample <- rownames(df)
  
  # Debugging prints
  print("Column names after removing zero variance columns:")
  print(colnames(df_numeric))
  
  # Perform PCA
  pca <- prcomp(df_numeric[, !colnames(df_numeric) %in% c("phenotype", "Sample")], scale. = TRUE)
  pca_df <- data.frame(pca$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df$Phenotype <- df_numeric$phenotype[match(rownames(pca_df), rownames(df_numeric))]
  
  # Debugging prints
  print("PCA Data Frame:")
  print(head(pca_df))
  
  # Plot PCA results with hover functionality
  pca_plot <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = "Phenotype", text = "Sample")) +
    geom_point() +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Phenotype") +  # Change legend title
    theme_minimal()
  
  return(ggplotly(pca_plot, tooltip = "text"))
}

# UMAP Function with Parallel Processing, Sliders, and Hover
Make_UMAP <- function(count_data, phenotype_data, n_neighbors) {
  if (grepl("\\.csv$", count_data)) {
    Countsmat <- read.csv(count_data, row.names = 1)
  } else {
    Countsmat <- read.table(count_data, header = TRUE, row.names = 1)
  }
  
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  disease_fact <- as.factor(disease[, 1])
  
  df <- data.frame(t(Countsmat))
  df$phenotype <- disease_fact
  df$Sample <- rownames(df)
  
  if (nrow(df) != length(disease_fact)) {
    stop("The number of samples in count_data and phenotype_data do not match.")
  }
  
  n_neighbors <- min(n_neighbors, nrow(df) - 1)
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  umap_res <- foreach(i = 1, .combine = rbind, .packages = 'umap') %dopar% {
    umap(df[, !(colnames(df) %in% c("phenotype", "Sample"))], n_neighbors = n_neighbors)
  }
  
  stopCluster(cl)
  
  umap_df <- data.frame(umap_res$layout)
  umap_df$Sample <- rownames(umap_df)
  umap_df$Phenotype <- df$phenotype
  
  if (nrow(umap_df) != nrow(df)) {
    stop("The number of rows in UMAP result and input data do not match.")
  }
  
  umap_plot <- ggplot(umap_df, aes(x = X1, y = X2, color = Phenotype, text = Sample)) +
    geom_point() +
    labs(title = "UMAP Plot", x = "UMAP1", y = "UMAP2") +
    theme_minimal()
  
  return(ggplotly(umap_plot, tooltip = "text"))
}

# Volcano Plot Function with Sliders and Hover
Make_Volcano <- function(results, p_value_threshold, log_fc_threshold) {
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), text = results[,4])) +
    geom_point(aes(color = adj.P.Val < p_value_threshold & abs(logFC) > log_fc_threshold)) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal()
  
  return(ggplotly(volcano_plot, tooltip = "text"))
}
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  theme = shinytheme("cyborg"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('notify', function(message) {
        alert(message);
      });
    "))
  ),
  titlePanel("JCAP Differential Expression and Feature Selection on RNA SEQ data"),
  sidebarLayout(
    sidebarPanel(
      fileInput("counts_input", "Input Counts", accept = c(".csv", ".txt")),
      fileInput("phenotype_input", "Input Phenotype Data", accept = c(".csv", ".txt")),
      actionButton("running_Dif", "Perform Differential Expression Analysis"),
      actionButton("plot_pca", "Display PCA Plot"),
      actionButton("plot_umap", "Display UMAP Plot"),
      actionButton("pathway_all", "Enrich Pathways (All Genes)"),
      actionButton("pathway_up", "Enrich Pathways (Upregulated Genes)"),
      actionButton("pathway_down", "Enrich Pathways (Downregulated Genes)"),
      actionButton("plot_volcano", "Display Volcano Plot"),
      actionButton("run_rf", "Run Random Forest Analysis"),  # Button for RF analysis
      actionButton("calculate_power", "Calculate Statistical Power"),
      sliderInput("n_components_pca", "Number of PCA Components", min = 2, max = 10, value = 2),
      sliderInput("n_neighbors_umap", "Number of UMAP Neighbors", min = 5, max = 50, value = 15),
      sliderInput("p_value_threshold", "P-Value Threshold for Volcano Plot", min = 0, max = 0.1, value = 0.05, step = 0.001),
      sliderInput("log_fc_threshold", "Log Fold Change Threshold for Volcano Plot", min = 0, max = 5, value = 1, step = 0.1),
      downloadButton("output", "Export Results as CSV"),
      downloadButton("All_genes","Export Enrichr results for all 20 genes as a CSV"),
      downloadButton("Up_genes","Export Enrichr results for upregulated genes as a CSV"),
      downloadButton("Down_genes","Export Enrichr results for downregulated genes as a CSV")
    ),
    mainPanel(
      tabsetPanel(id = "main_tabset",
                  tabPanel("Differential Expression Results", DTOutput("Dif_expr_results")),
                  tabPanel("PCA Plot", plotlyOutput("pcaplot")),  # Change to plotlyOutput
                  tabPanel("UMAP Plot", plotlyOutput("umapplot")),  # Change to plotlyOutput
                  tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),  # Change to plotlyOutput
                  tabPanel("Pathway Enrichment", 
                           tabsetPanel(
                             tabPanel("All Genes", DTOutput("enrich_table_all")),
                             tabPanel("Upregulated Genes", DTOutput("enrich_table_up")),
                             tabPanel("Downregulated Genes", DTOutput("enrich_table_down"))
                           )
                  ),
                  tabPanel("Random Forest Performance Metrics", 
                           DTOutput("rf_metrics"),
                           downloadButton("download_metrics", "Download Performance Metrics as CSV")  # New download button
                  ),
                  tabPanel("Power Calculation", verbatimTextOutput("power_results")),
                  tabPanel("Read Me", verbatimTextOutput("readme"))  # Add a new tab for README
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  # Reactive values to store enrichment results
  enrich_res_up <- reactiveVal()
  enrich_res_down <- reactiveVal()
  final_results <- reactiveVal()
  rf_metrics <- reactiveVal()  # Store the RF metrics
  
  observeEvent(input$running_Dif, {
    shinyjs::info("Starting Differential Expression Analysis...")
    res <- dif_expr(count_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath)
    final_results(res)
    output$Dif_expr_results <- renderDT({
      datatable(res)
    })
    session$sendCustomMessage("notify", "Differential Expression Analysis Completed!")
  })
  
  output$output <- downloadHandler(
    filename = "differential_expression_results.csv",
    content = function(file) {
      write.csv(final_results(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$plot_pca, {
    req(input$counts_input, input$phenotype_input)
    shinyjs::info("Generating PCA Plot...")
    output$pcaplot <- renderPlotly({
      Make_PCA(count_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath, n_components = input$n_components_pca)
    })
    session$sendCustomMessage("notify", "PCA Plot Generated!")
  })
  
  observeEvent(input$plot_umap, {
    req(input$counts_input, input$phenotype_input)
    shinyjs::info("Generating UMAP Plot...")
    output$umapplot <- renderPlotly({
      Make_UMAP(count_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath, n_neighbors = input$n_neighbors_umap)
    })
    session$sendCustomMessage("notify", "UMAP Plot Generated!")
  })
  
  observeEvent(input$plot_volcano, {
    req(final_results())
    shinyjs::info("Generating Volcano Plot...")
    output$volcano_plot <- renderPlotly({
      Make_Volcano(final_results(), input$p_value_threshold, input$log_fc_threshold)
    })
    session$sendCustomMessage("notify", "Volcano Plot Generated!")
  })
  
  observeEvent(input$pathway_all, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (All Genes)...")
    res <- final_results()
    enrich_res <- enrich_all_genes(res)
    output$enrich_table_all <- renderDT({
      datatable(enrich_res)
    })
    session$sendCustomMessage("notify", "Pathway Enrichment (All Genes) Completed!")
  })
  
  output$All_genes <- downloadHandler(
    filename = "Enrichment_on_all_genes.csv",
    content = function(file) {
      write.csv(combined_res_all, file, row.names = FALSE)
    }
  )

  
  observeEvent(input$pathway_up, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (Upregulated Genes)...")
    res <- final_results()
    up_result <- enrich_up_genes(res)
    enrich_res_up(up_result)
    output$enrich_table_up <- renderDT({
      if ("Message" %in% colnames(up_result)) {
        datatable(data.frame(Message = "No upregulated genes found."), options = list(dom = 't'))
      } else {
        datatable(up_result)
      }
    })  
    session$sendCustomMessage("notify", "Pathway Enrichment (Upregulated Genes) Completed!")
  })
  
  observeEvent(input$pathway_down, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (Downregulated Genes)...")
    res <- final_results()
    down_result <- enrich_down_genes(res)
    enrich_res_down(down_result)
    output$enrich_table_down <- renderDT({
      if ("Message" %in% colnames(down_result)) {
        datatable(data.frame(Message = "No downregulated genes found."), options = list(dom = 't'))
      } else {
        datatable(down_result)
      }
    }) 
    session$sendCustomMessage("notify", "Pathway Enrichment (Downregulated Genes) Completed!")
  })
  
  output$Up_genes <- downloadHandler(
    filename = "Enrichment_on_upregulated_genes.csv",
    content = function(file) {
      up_data <- enrich_res_up()  # Get the reactive value
      if (!is.null(up_data) && !("Message" %in% colnames(up_data))) {
        write.csv(up_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No upregulated genes found."), file, row.names = FALSE)
      }
    }
  )
  
  output$Down_genes <- downloadHandler(
    filename = "Enrichment_on_downregulated_genes.csv",
    content = function(file) {
      down_data <- enrich_res_down()  # Get the reactive value
      if (!is.null(down_data) && !("Message" %in% colnames(down_data))) {
        write.csv(down_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No downregulated genes found."), file, row.names = FALSE)
      }
    }
  )
  
  # Random Forest Button Action
  observeEvent(input$run_rf, {
    req(final_results(), input$counts_input, input$phenotype_input)
    shinyjs::info("Running Random Forest Analysis...")
    selected_genes <- final_results()[,1]
    # Get gene list from final results
    metrics <- regularized_random_forest(counts_data = input$counts_input$datapath, 
                                            phenotype_data = input$phenotype_input$datapath, 
                                            selected_genes = selected_genes)
    rf_metrics(metrics)  # Store the metrics
    output$rf_metrics <- renderDT({
      datatable(metrics)
    })
    session$sendCustomMessage("notify", "Random Forest Analysis Completed!")
  })
  
  # Download handler for the performance metrics
  output$download_metrics <- downloadHandler(
    filename = "random_forest_performance_metrics.csv",
    content = function(file) {
      write.csv(rf_metrics(), file, row.names = FALSE)
    }
  )
  observeEvent(input$calculate_power, {
    # Ensure data is loaded
    print("Power analysis button clicked.")
    
    shinyjs::info("Starting Statistical Power Calculation...")
    
    req(input$counts_input, input$phenotype_input)
    
    # Load pheno data
    if (grepl("\\.csv$", input$phenotype_input$datapath)) {
      pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)
    } else {
      pheno <- read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1)
    }
    
    # Calculate power based on the number of samples in the expression data
    power_result <- calculate_power(phenotype_data = pheno, effect_size = 0.5)
    
    # Display the power result
    output$power_results <- renderPrint({
      print(power_result)
    })
    
    session$sendCustomMessage("notify", "Statistical Power Calculation Completed!")
  })
  
  
  
  # Render the README file content
  output$readme <- renderPrint({
    readme_path <- normalizePath("JCAP RNA_SEQ Readme.txt", mustWork = TRUE)
    readLines(readme_path)  # Make sure to place your README.txt in the same directory as the app
  })
}

# ==== GLOBAL ERROR HANDLER FOR LOGGING ==== #
log_file <- "error_log.txt"

options(shiny.error = function() {
  err <- geterrmessage()
  timestamp <- Sys.time()
  msg <- paste0("[", timestamp, "] ", err, "\n\n")
  cat(msg, file = log_file, append = TRUE)
})

                             
# Run the application
shinyApp(ui = ui, server = server)
