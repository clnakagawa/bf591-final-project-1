library(shiny)
library(DESeq2)
library(ggplot2)
library(DT)
library(pheatmap)
library(fgsea)
library(tidyr)
library(dplyr)
library(ggbeeswarm)


# Define UI
ui <- fluidPage(
  titlePanel("Huntington's Disease mRNA-seq Data Analysis"),
  
  # Define main tabs
  tabsetPanel(
    tabPanel("Sample Data Exploration",
             sidebarLayout(
               sidebarPanel(
                 fileInput("sampleFile", "Upload Sample Information Matrix", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", DTOutput("dataSummary")),
                   tabPanel("Full Table", DT::dataTableOutput("sampleDataTable")),
                   tabPanel("Plot", 
                            selectInput("column", "Select a column for the violin plot:", choices = NULL),
                            fluidRow(
                              column(4, plotOutput("sampleHist")),
                              column(4, plotOutput("sampleViolin")),
                              column(4, plotOutput("sampleDist"))
                            ))
                 )
               )
             )
    ),
    
    tabPanel("Counts Matrix Exploration",
             sidebarLayout(
               sidebarPanel(
                 fileInput("countsFile", "Upload Normalized Counts Matrix (CSV)", accept = ".csv"),
                 sliderInput("varianceSlider", "Minimum Variance", min = 0, max = 1, value = 0.5),
                 sliderInput("nonZeroSlider", "Minimum Non-Zero Samples", min = 0, max = 69, value = 50),
                 actionButton("applyFilter", "Apply Filter"),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filtering Effects", DTOutput("filterSummaryTable")),
                   tabPanel("Diagnostic Plots",
                            plotOutput("medianVsVariancePlot"),
                            plotOutput("medianVsZerosPlot")),
                   tabPanel("Clustered Heatmap", 
                            plotOutput("heatmapPlot", height = "700px")),
                   tabPanel("PCA",
                            numericInput("numPCs", "Number of Principal Components to Plot:", value = 2, min = 1),
                            plotOutput("pcaBeeswarmPlot")
                   )
                 )
               )
             )
    ),
    
    tabPanel("Differential Expression Analysis",
             sidebarLayout(
               sidebarPanel(
                 fileInput("deFile", "Upload DESeq2 Results DataFrame (CSV)", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("DE Table", DT::dataTableOutput("deResults")),
                   tabPanel("Volcano Plot", 
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("pValueCutoff", "P-value Cutoff (Log Scale):", min = -30, max = -1, value = -5, step = 1),
                                sliderInput("logFCutoff", "Log2 Fold Change Cutoff:", min = 0, max = 5, value = 1, step = 0.1)
                              ),
                              mainPanel(
                                plotOutput("volcanoPlot")
                              )
                            ))
                 )
               )
             )
    ),
    
    tabPanel("Gene Set Enrichment Analysis",
             sidebarLayout(
               sidebarPanel(
                 fileInput("gseaFile", "Upload GSEA Results (CSV)", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   # Tab 1: Barplot of top pathways by NES
                   tabPanel("Top Pathways Barplot",
                            sidebarPanel(
                              sliderInput("topPathways", "Number of top pathways to plot", min = 1, max = 20, value = 5)
                            ),
                            mainPanel(plotOutput("barPlot"))
                   ),
                   # Tab 2: Filtered Table of GSEA Results
                   tabPanel("Filtered Results Table",
                            sidebarPanel(
                              sliderInput("adjustedPValue", "Adjusted P-value Cutoff (Log Scale):", min = -50, max = -1, value = -5, step = 1),
                              radioButtons("nesFilter", "Select NES Direction", choices = c("All" = "all", "Positive" = "positive", "Negative" = "negative"))
                            ),
                            mainPanel(
                              DT::dataTableOutput("filteredGSEAResults"),
                              downloadButton("downloadTable", "Download Filtered Results")
                            )
                   ),
                   # Tab 3: Scatter Plot of NES vs -log10 Adjusted P-value
                   tabPanel("NES vs Adjusted P-value Scatter",
                            sidebarPanel(
                              sliderInput("scatterPValue", "Adjusted P-value Cutoff (Log Scale):", min = -50, max = -1, value = -5, step = 1)
                            ),
                            mainPanel(
                              plotOutput("scatterPlot")
                            )
                   )
                 )
               )
             )
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  # change maximum upload size for later counts data
  options(shiny.maxRequestSize=30*1024^2)
  
  # Sample Data Exploration
  observe({
    # make sure file is loaded
    req(input$sampleFile)
    
    # get sample_data 
    sample_data <- read.csv(input$sampleFile$datapath, row.names = 1)
    
    # update dropdown options for plot tab
    numeric_columns <- names(sample_data)[sapply(sample_data, is.numeric)]
    updateSelectInput(session, "column", choices = numeric_columns)
    
    # Make summary table
    output$dataSummary <- renderDT({
      # Create a summary dataframe
      summary_df <- lapply(names(sample_data), function(col_name) {
        col_data <- sample_data[[col_name]]
        
        # Determine data type
        data_type <- class(col_data)[1]
        
        # Create summary based on data type
        if (is.numeric(col_data)) {
          # Calculate mean and SD
          col_mean <- mean(col_data, na.rm = TRUE)
          col_sd <- sd(col_data, na.rm = TRUE)
          # Format mean/sd for summary table
          summary_value <- sprintf("%.2f (+/-%.2f)", col_mean, col_sd)
        } else {
          # Concatenate distinct values for categorical data
          distinct_values <- unique(col_data)
          summary_value <- paste(distinct_values, collapse = ", ")
        }
        
        # Return a row of the summary
        list(
          Column = col_name,
          DataType = data_type,
          Summary = summary_value
        )
      })
      
      # Convert list to dataframe
      summary_df <- do.call(rbind, lapply(summary_df, as.data.frame))
      
      # Render the summary using DT
      datatable(summary_df)
    })
    
    # Show full data table 
    output$sampleDataTable <- DT::renderDataTable({
      DT::datatable(sample_data, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$sampleHist <- renderPlot({
      req(input$column)  # Ensure a column is selected
      selected_column <- input$column
      
      # use ggplot to make histogram
      ggplot(sample_data, aes(x = .data[[selected_column]])) +
        geom_histogram(fill = "lightblue", color = "black") + 
        labs(
          x = selected_column,
          y = "Frequency",
          title = paste("Histogram for", selected_column)
        ) +
        theme(
          plot.title = element_text(size = 16, face = "bold")  # Larger, bold title
        )
    })
    
    output$sampleViolin <- renderPlot({
      req(input$column)  # Ensure a column is selected
      selected_column <- input$column
      
      # Create the violin plot
      ggplot(sample_data, aes(x = "", y = .data[[selected_column]])) +
        geom_violin(fill = "lightblue", color = "black") +
        labs(
          x = "",
          y = selected_column,
          title = paste("Violin Plot of", selected_column)
        ) +
        theme(
          plot.title = element_text(size = 16, face = "bold")  # Larger, bold title
        )
    })
    
    output$sampleDist <- renderPlot({
      req(input$column)  # Ensure a column is selected
      selected_column <- input$column
      
      # Create the density plot
      ggplot(sample_data, aes(x = .data[[selected_column]])) +
        geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
        labs(
          x = selected_column,
          y = "Density",
          title = paste("Distribution of", selected_column)
        ) +
        theme(
          plot.title = element_text(size = 16, face = "bold")  # Larger, bold title
        )
    })
  })
  
  # Counts Matrix Exploration
  observe({
    req(input$countsFile)
    counts_data <- read.csv(input$countsFile$datapath, row.names = 1)
    
    # Filtered counts matrix based on sliders
    filtered_counts <- eventReactive(input$applyFilter, {
      variance_threshold <- input$varianceSlider
      non_zero_threshold <- input$nonZeroSlider
      variance <- apply(counts_data, 1, var)
      non_zero_samples <- apply(counts_data != 0, 1, sum) 
      counts_data[variance > quantile(variance, variance_threshold) & non_zero_samples > non_zero_threshold, ]
    })
    
    output$filterSummaryTable <- renderDT({
      # filtered data
      filtered <- filtered_counts()
      
      # Calculate summary
      num_samples <- ncol(counts_data)  # Number of columns (samples)
      total_genes <- nrow(counts_data)  # Total number of genes
      passing_genes <- nrow(filtered)  # Genes passing the filter
      not_passing_genes <- total_genes - passing_genes  # Genes not passing the filter
      percent_passing <- (passing_genes / total_genes) * 100
      percent_not_passing <- (not_passing_genes / total_genes) * 100
      
      # Create summary table
      summary_df <- data.frame(
        Metric = c("Number of Samples", 
                   "Total Number of Genes", 
                   "Number of Genes Passing Filter", 
                   "Percentage of Genes Passing Filter", 
                   "Number of Genes Not Passing Filter", 
                   "Percentage of Genes Not Passing Filter"),
        Value = c(num_samples, 
                  total_genes, 
                  passing_genes, 
                  sprintf("%.2f%%", percent_passing), 
                  not_passing_genes, 
                  sprintf("%.2f%%", percent_not_passing))
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$medianVsVariancePlot <- renderPlot({
      # Calculate median and variance for all genes
      filtered <- filtered_counts()
      
      gene_stats <- data.frame(
        Median = apply(counts_data, 1, median),
        Variance = apply(counts_data, 1, var),
        PassedFilter = rownames(counts_data) %in% rownames(filtered)
      )
      
      # Create scatter plot
      ggplot(gene_stats, aes(x = Median, y = Variance, color = PassedFilter)) +
        geom_point(alpha = 0.7) +
        scale_x_log10() + 
        scale_y_log10() +
        labs(
          title = "Median Count vs Variance",
          x = "Median Count (log scale)",
          y = "Variance (log scale)"
        ) +
        scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightgray")) +
        theme_minimal() +
        theme(plot.title = element_text(size = 16, face = "bold"))
    })
    
    # Render median count vs number of zeros scatter plot
    output$medianVsZerosPlot <- renderPlot({
      # Calculate median and number of zeros for all genes
      filtered <- filtered_counts()
      
      gene_stats <- data.frame(
        Median = apply(counts_data, 1, median),
        NumZeros = rowSums(counts_data == 0),
        PassedFilter = rownames(counts_data) %in% rownames(filtered)
      )
      
      # Create scatter plot
      ggplot(gene_stats, aes(x = Median, y = NumZeros, color = PassedFilter)) +
        geom_point(alpha = 0.7) +
        scale_x_log10() +
        labs(
          title = "Median Count vs Number of Zeros",
          x = "Median Count (log scale)",
          y = "Number of Zeros"
        ) +
        scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightgray")) +
        theme_minimal() +
        theme(plot.title = element_text(size = 16, face = "bold"))
    })
    
    output$heatmapPlot <- renderPlot({
      req(filtered_counts())
      # Get the filtered counts matrix
      filtered <- log2(filtered_counts() + 1)
      
      #print(filtered)
      # Create heatmap
      pheatmap(
        mat = as.matrix(filtered),
        cluster_cols = T,
        cluster_rows = T,
        legend = TRUE,
        main = "Clustered Heatmap of Filtered Counts (log2 scaled)",
        fontsize = 12,
        show_rownames = F,
        show_colnames = F
      )
    })
    
    output$pcaBeeswarmPlot <- renderPlot({
      req(filtered_counts())  # Ensure the filtered counts matrix exists
      
      counts <- filtered_counts()
      # Perform PCA
      pca_result <- prcomp(t(counts), center = TRUE, scale. = TRUE)  # PCA on samples (transpose counts)
      
      # Extract PCA variance explained
      var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
      
      # Extract PCA scores
      pca_scores <- as.data.frame(pca_result$x)
      pca_scores$Sample <- rownames(pca_scores)  # Add sample names
      
      # Limit to top N principal components
      num_pcs <- min(input$numPCs, ncol(pca_scores) - 1)  # Ensure valid number of PCs
      pcs_to_plot <- paste0("PC", 1:num_pcs)
      pca_scores_long <- pca_scores %>%
        select(Sample, all_of(pcs_to_plot)) %>%
        tidyr::pivot_longer(
          cols = starts_with("PC"),
          names_to = "PrincipalComponent",
          values_to = "Score"
        )
      
      # Add variance explained as a label
      pca_scores_long$VarianceExplained <- sapply(
        pca_scores_long$PrincipalComponent,
        function(pc) round(var_explained[as.numeric(sub("PC", "", pc))], 2)
      )
      
      # Generate beeswarm plot
      ggplot(pca_scores_long, aes(x = PrincipalComponent, y = Score)) +
        geom_beeswarm(size = 3, alpha = 0.8, color = "blue") +
        labs(
          title = paste("Beeswarm Plot of Top", num_pcs, "Principal Components"),
          x = "Principal Component",
          y = "PCA Score"
        ) +
        scale_x_discrete(labels = function(pc) {
          sapply(pc, function(pc_name) {
            pc_index <- as.numeric(sub("PC", "", pc_name))
            paste0(pc_name, " (", round(var_explained[pc_index], 2), "%)")
          })
        }) +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    })
  })
  
  # Differential Expression Analysis
  observe({
    req(input$deFile)
    de_data <- read.csv(input$deFile$datapath, row.names = 1)
    output$deResults <- DT::renderDataTable({
      DT::datatable(de_data, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$volcanoPlot <- renderPlot({
      # get p-val and logFC cutoff
      pvalue_cutoff <- 10^input$pValueCutoff
      logfc_cutoff <- input$logFCutoff
      
      ggplot(de_data, aes(x = log2FoldChange, y = -log10(pvalue), color = pvalue <= pvalue_cutoff & abs(log2FoldChange) >= logfc_cutoff)) +
        geom_point(alpha = 0.6) +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12)
        ) +
        labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value")
    })
  })
  
  # Gene Set Enrichment Analysis (GSEA)
  gsea_data <- reactive({
    req(input$gseaFile)
    read.csv(input$gseaFile$datapath)
  })
  
  # Tab 1: Barplot of top pathways by NES
  output$barPlot <- renderPlot({
    gsea_results <- gsea_data()
    top_pathways <- gsea_results[order(gsea_results$padj),][1:input$topPathways, ]
    ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Top Pathways by NES", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
      theme_minimal()
  })
  
  # Tab 1: Display table entry of selected pathway
  output$selectedPathwayTable <- renderDataTable({
    gsea_results <- gsea_data()
    top_pathways <- gsea_results[order(gsea_results$padj),][1:input$topPathways, ]
    DT::datatable(top_pathways, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Tab 2: Filtered Table of GSEA Results
  output$filteredGSEAResults <- renderDataTable({
    gsea_results <- gsea_data()
    
    # Filter by adjusted p-value
    filtered_data <- gsea_results[gsea_results$padj <= 10^input$adjustedPValue, ]
    
    # Filter by NES direction (positive or negative)
    if (input$nesFilter == "positive") {
      filtered_data <- filtered_data[filtered_data$NES > 0, ]
    } else if (input$nesFilter == "negative") {
      filtered_data <- filtered_data[filtered_data$NES < 0, ]
    }
    
    DT::datatable(filtered_data, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Tab 2: Download button for filtered GSEA results
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("filtered_gsea_results.csv")
    },
    content = function(file) {
      write.csv(gsea_data(), file)
    }
  )
  
  # Tab 3: Scatter Plot of NES vs -log10 Adjusted P-value
  output$scatterPlot <- renderPlot({
    gsea_results <- gsea_data()
    
    ggplot(gsea_results, aes(x = NES, y = -log10(padj))) +
      geom_point(aes(color = ifelse(padj <= 10^input$scatterPValue, "Significant", "Not Significant")), size = 2) +
      scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
      labs(title = "NES vs -log10 Adjusted P-value", x = "Normalized Enrichment Score (NES)", y = "-log10 Adjusted P-value") +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
