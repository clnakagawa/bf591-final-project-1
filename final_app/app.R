library(shiny)
library(DESeq2)
library(ggplot2)
library(DT)
library(pheatmap)
library(fgsea)

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
                 sliderInput("varianceSlider", "Percentile of Variance", min = 0, max = 100, value = 50),
                 sliderInput("nonZeroSlider", "Minimum Non-Zero Samples", min = 0, max = 100, value = 50)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filtering Effects", DT::dataTableOutput("filteredCounts")),
                   tabPanel("Diagnostic Plots", plotOutput("diagnosticPlot")),
                   tabPanel("Clustered Heatmap", plotOutput("heatmapPlot")),
                   tabPanel("PCA Plot", plotOutput("pcaPlot"))
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
                   tabPanel("Volcano Plot", plotOutput("volcanoPlot"))
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
                            mainPanel(
                              plotOutput("barPlot"),
                              DT::dataTableOutput("selectedPathwayTable")
                            )
                   ),
                   # Tab 2: Filtered Table of GSEA Results
                   tabPanel("Filtered Results Table",
                            sidebarPanel(
                              sliderInput("adjustedPValue", "Filter by Adjusted P-value", min = 0, max = 1, value = 0.05),
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
                              sliderInput("scatterPValue", "Filter by Adjusted P-value", min = 0, max = 1, value = 0.05)
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
    filtered_counts <- reactive({
      variance_threshold <- input$varianceSlider / 100
      non_zero_threshold <- input$nonZeroSlider / 100
      variance <- apply(counts_data, 1, var)
      non_zero_samples <- apply(counts_data != 0, 1, sum) / ncol(counts_data)
      counts_data[variance > quantile(variance, variance_threshold) & non_zero_samples > non_zero_threshold, ]
    })
    
    output$filteredCounts <- DT::renderDataTable({
      DT::datatable(filtered_counts(), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$diagnosticPlot <- renderPlot({
      hist(filtered_counts())
    })
    
    output$heatmapPlot <- renderPlot({
      pheatmap(filtered_counts())
    })
    
    output$pcaPlot <- renderPlot({
      pca <- prcomp(t(filtered_counts()))
      ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        labs(title = "PCA Plot")
    })
  })
  
  # Differential Expression Analysis
  observe({
    req(input$deFile)
    de_data <- read.csv(input$deFile$datapath)
    output$deResults <- DT::renderDataTable({
      DT::datatable(de_data, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$volcanoPlot <- renderPlot({
      ggplot(de_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(alpha = 0.6) +
        theme_minimal() +
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
    filtered_data <- gsea_results[gsea_results$padj <= input$adjustedPValue, ]
    
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
    filtered_results <- gsea_results[gsea_results$padj <= input$scatterPValue, ]
    
    ggplot(filtered_results, aes(x = NES, y = -log10(padj))) +
      geom_point(aes(color = ifelse(padj <= input$scatterPValue, "Significant", "Not Significant")), size = 2) +
      scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
      labs(title = "NES vs -log10 Adjusted P-value", x = "Normalized Enrichment Score (NES)", y = "-log10 Adjusted P-value") +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
