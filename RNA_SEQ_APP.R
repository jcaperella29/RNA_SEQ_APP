
library(limma)
library("edgeR")
library('biomaRt')

library(shiny)
library(shinythemes)
library(DT)
library("FactoMineR")
library(ggfortify)
library(ggplot2)

library(randomForest)
library(varImp)
library(dplyr )
dif_expr <- function(count_data,phenotype_data){
  
  
  Countsmat<-read.table(count_data,header = TRUE)
  
  
  disease<-read.table(phenotype_data,header=TRUE)
  disease_fact<-as.factor(disease[,1])
  df<-data.frame(Countsmat)
  
  df<-Countsmat
  
  df<-data.frame(t(Countsmat))
  
  df$phenotype<-disease_fact
  good_stuff<- Countsmat[rowSums(Countsmat[])>10,]
  good_stuff$goodgenes<-row.names(good_stuff)
  df<-df[,colnames(df) %in% c(good_stuff$goodgenes,"phenotype")]
  
  Countsmat<-subset(df, select = -phenotype )
  
  Countsmat<-t(Countsmat)
  
  disease_fact<-(subset(df,select = "phenotype"))
  
  disease_fact<-as.factor(disease[,1])

  disease_data<-data.frame(sapply(levels(disease_fact) ,function(x) as.integer(x == disease[,1])))
  
  PCA_data<-data.frame(t(Countsmat))
  PCA_data$pheno<-disease_data[,2]
  pca <<- prcomp(PCA_data)
  #dark dots = has the phenotype of interest
  PCA_PLOT<<-autoplot(pca,colour="pheno")
  
  d0 <- DGEList(Countsmat)
  
  d0 <- calcNormFactors(d0)
  
  design <- model.matrix(~0+disease_fact)
  
  v <- voom(Countsmat, design, plot=FALSE, normalize="quantile")
  fit <- lmFit(v, design)
  
  
  fit <- eBayes(fit)
  results<<-topTable(fit, coef=NULL, number=200, genelist=fit$genes, adjust.method="BH",
                    sort.by="B", resort.by=NULL, p.value=.05, lfc=0, confint=FALSE)
  results$Ensembl_IDs<<-row.names(results)
  results$log_fold_change<<-log2(results$disease_factCOVID.19/results$disease_factnormal)  
  
  results<<-results[,c("Ensembl_IDs","adj.P.Val","log_fold_change" )]
  row.names(results)<<-NULL
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- results$Ensembl_IDs
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                 "hgnc_symbol", "entrezgene_description"),values=genes,mart= mart)
  results$hgnc_symbol<<-G_list[1:200,2]
  results$protein_name<<-G_list[1:200,3]
  df<-df[,colnames(df) %in% c(results$Ensembl_IDs,"phenotype")]
  model<-randomForest(phenotype~.,data = df,importance=TRUE)
  importance_df<-data.frame(importance(model))
  importance_df$genenames<-row.names(importance_df)
  importance_df<-importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE),]
  Top20_RFgenes<-importance_df$genenames[1:20]
  final_results<<- results[results$Ensembl_IDs %in% Top20_RFgenes,]
  return(final_results )
}


Make_PCA <-function (count_data,phenotype_data){
  Countsmat<-read.table(count_data,header = TRUE)
  
  
  disease<-read.table(phenotype_data,header=TRUE)
disease_data<-data.frame(sapply(levels(disease_fact) ,function(x) as.integer(x == disease[,1])))

PCA_data<-data.frame(t(Countsmat))
PCA_data$pheno<-disease_data[,2]
pca <<- prcomp(PCA_data)
#dark dots = has the phenotype of  interest
PCA_PLot<<-autoplot(pca,colour="pheno")
return(PCA_PLot)
}


ui <- fluidPage(theme = shinytheme("cyborg") ,titlePanel(" Differential Expression and Feature Selection on RNA SEQ data" ),
                inputPanel(
                  fileInput("counts_input", "input counts"),
                  br(),
                  fileInput("phenotype_input","input phenotype data")),
                br(),
                actionButton("running_Dif", "Perform differential expression analysis and Feature Selection on your RNA_SEQ  data"),
                br(),
                br(),
                actionButton("plot_pca" ,"prepare and display a PCA plot of the data"),
                mainPanel(DTOutput("Dif_expr_results"), plotOutput("pcaplot"),
   br(),
   br(),
   downloadButton(outputId = "output",label = "export as  csv")))
   
  
            

   
    



  






server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)

    observeEvent(input$running_Dif , { 
    output$Dif_expr_results<-DT::renderDT(dif_expr(count_data=input$counts_input$datapath,phenotype=input$phenotype_input$datapath)
                                          ) })
    
    output$output <-downloadHandler(filename = "my_output.csv",
                                     content= function(file){
                                       write.csv(final_results,file)}) 

    observeEvent(input$plot_pca , {                                 
    output$pcaplot<- renderPlot({PCA_PLOT})
                                     } ) }
    
    
    
      

    
     
     



shinyApp(ui=ui ,server = server)

