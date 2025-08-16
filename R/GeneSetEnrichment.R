# R Shiny app for Gene Expression Enrichment Analysis

# Load necessary libraries. If you don't have them installed, run:
# install.packages("shiny")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# install.packages("ggplot2")
# install.packages("dplyr")

library(shiny)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# UI for the application
ui <- fluidPage(
  
  # App title
  titlePanel("Gene Expression Enrichment Analysis"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      h4("Gene List Input"),
      
      # Text area for gene list input
      textAreaInput(
        "gene_list_input",
        "Enter your gene list (one gene per line):",
        # Example gene list for demonstration
        value = "RPL3\nRPB1\nRPS6\nRPL10A",
        rows = 10
      ),
      
      # Button to trigger analysis
      actionButton("run_analysis", "Run Analysis", class = "btn-primary")
      
    ),
    
    # Main panel for output
    mainPanel(
      h4("Enrichment Results (Gene Ontology)"),
      
      # Display the enrichment plot
      plotOutput("enrichment_plot"),
      
      # Display a message while processing
      textOutput("processing_message")
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive expression to store enrichment results
  enrichment_result <- eventReactive(input$run_analysis, {
    
    # Show processing message
    output$processing_message <- renderText("Running enrichment analysis...")
    
    # 1. Process the gene list input
    # Split the input string by newlines, remove empty entries, and trim whitespace
    gene_list_text <- input$gene_list_input
    if (is.null(gene_list_text) || gene_list_text == "") {
      showModal(modalDialog(
        title = "No Genes Entered",
        "Please enter a list of genes to analyze."
      ))
      return(NULL)
    }
    
    # Clean the gene list
    genes <- trimws(unlist(strsplit(gene_list_text, "\n")))
    genes <- genes[genes != ""]
    
    # Check if any genes remain
    if (length(genes) == 0) {
      showModal(modalDialog(
        title = "No Genes Found",
        "The entered list is empty or contains only spaces."
      ))
      return(NULL)
    }
    
    # 2. Perform Gene Ontology (GO) enrichment analysis
    # NOTE: `org.Hs.eg.db` uses Entrez IDs. For simplicity, we'll assume the
    # user's list are gene symbols that match the database.
    
    # --- CORRECTED: Use all human genes from the database as the universe ---
    # This is a much more appropriate background for a general enrichment analysis.
    all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    
    tryCatch({
      # Perform the enrichment analysis using enrichGO
      go_enrich_results <- enrichGO(
        gene = genes,
        universe = all_genes, # Use the comprehensive list of all human genes
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP", # Use Biological Process ontology
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
      )
      
      # Hide processing message
      output$processing_message <- renderText("")
      
      return(go_enrich_results)
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Analysis Error",
        paste("An error occurred during analysis:", e$message),
        footer = modalButton("OK")
      ))
      output$processing_message <- renderText("")
      return(NULL)
    })
  })
  
  # Render the enrichment plot
  output$enrichment_plot <- renderPlot({
    
    # Get the results from the reactive expression
    results <- enrichment_result()
    
    # Check if the enrichment result is a valid object with significant hits
    if (is.null(results) || !inherits(results, "enrichResult") || nrow(results) == 0) {
      # Return a blank plot with a message if no significant results are found
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No significant enrichment results found.",
                 size = 5, color = "grey50") +
        theme_void()
    } else {
      # Use dotplot from clusterProfiler to visualize the results
      dotplot(results, showCategory = 15, title = "Gene Ontology Enrichment Analysis")
    }
  })
}

# Run the application
shinyApp(ui, server)
