suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(shinydashboard)
  library(plotly)
  library(shinycssloaders)
  library(aws.s3)          # For s3 functionality
  library(clusterProfiler) # For gene set enrichment
  library(org.Hs.eg.db)    # For gene set enrichment
  library(vroom)           # for reading in orthologs for gene set enrichment
  library(jsonlite)        # for tool tips
})

# Define UI for application that draws a histogram
ui <- function(request) {   # Note that I might need to remove "function(request)" for Google Analytics to work.  Revisit later if this breaks anything.
  dashboardPage(
    
    dashboardHeader(title = "Cell HierARchy Gene Explorer (CHARGE)", titleWidth = 500),
    
    dashboardSidebar(title = tags$img(src='CHARGE_logo.png', width = '225', style= 'display: block; margin-left: auto; margin-right: auto;'), #), #disable = TRUE
    
      width = 250,
                  
      h3("What is CHARGE?"),
      p("Cell HierARchy Gene Explorer (CHARGE) is an interactive tool for identifying genes of potential biological interest in a cell type taxonomy. CHARGE uses fast, cluster-centric approaches to find global or local marker genes, differentially expressed genes, and genes following pre-defined or user selected gradients. Outputs of CHARGE include basic visualization of interesting genes and download of gene names and analysis statistics. CHARGE emphasizes speed and usability over statistical rigor, and we encourage users to explore cell level data in the Allen Brain Cell Atlas, CELLxGENE, or related tools."),

      h3("Get started"),
      
      actionButton(inputId = "usecase", 
                   icon = icon("circle-play", lib = "font-awesome"), 
                   a("COMING SOON!",
                     style="color: #000000; border-color: #2e6da4",
                     target = "_blank", 
                     href="https://portal.brain-map.org/")
      ),
      
      h3("Contribute"),
      
      actionButton(inputId = "email1", 
                   icon = icon("envelope", lib = "font-awesome"), 
                   a("PROVIDE FEEDBACK", 
                     style="color: #000000; border-color: #2e6da4",
                     href="mailto:jeremym@alleninstitute.org?
                                  body=''
                                  &subject='CHARGE' app comments")
      ),
      actionButton(inputId = "GitHub", 
                   icon = icon("code", lib = "font-awesome"), 
                   a("ACCESS SOURCE CODE",
                     style="color: #000000; border-color: #2e6da4",
                     target = "_blank", 
                     href="https://github.com/AllenInstitute/CHARGE/")
      ),

      br(),
      h4("Click the three lines next to the title above to minimize this sidebar."),
      br(),
      p("----------------"),
      br(),
      h3("Acknowledgements"),
      p("App developed by Jeremy Miller using some original code developed by Lucas Graybuck and Cindy van Velthoven, and connects to the AIT format developed in collaboration with Nelson Johansen and Inkar Kapen.  Included tables developed through BICAN."),
      br(),
      p("If you would like to contribute to this app, please reach out via email or GitHub using the links above."),
      br()
    ),
    
    dashboardBody(
    
    # Replace with correct google-analytics code block and html file below if I decided to use it.
    #   tags$head(includeHTML("google-analytics.html"),  # Tag for general Google Analytics!
    #             tags$script('var dimension = [0, 0];
    #                       $(document).on("shiny:connected", function(e) {
    #                           dimension[0] = window.innerWidth;
    #                           dimension[1] = window.innerHeight;
    #                           Shiny.onInputChange("dimension", dimension);
    #                       });
    #                       $(window).resize(function(e) {
    #                           dimension[0] = window.innerWidth;
    #                           dimension[1] = window.innerHeight;
    #                           Shiny.onInputChange("dimension", dimension);
    #                       });
    #                       ')),
      
      # You can add custom CSS here
      tags$head(tags$style(HTML("
      .main-sidebar {
        background-color: #1C2532 !important; /* A shade of blue */
      }
      #enrichment_table td {
        white-space: normal !important;
        height: auto; /* Forces the cell to adopt its content height */
        line-height: normal; /* Ensures line spacing is standard */
      }
      #enrichment_table table.dataTable thead th {
        vertical-align: top;
      }
      "))),
      
      #useShinyjs(),  # shinyjs not currently used
      
      fluidRow(width = 12,
               
               box(title = "Select data set",
                   solidHeader = TRUE, status = "primary", width = 12,
                   collapsible = TRUE, collapsed = FALSE,
                   fluidRow(
                     column(6,
                            uiOutput("select_category")
                     ),
                     column(6,
                            uiOutput("select_textbox")
                     ),
                   ),
                   fluidRow(
                     column(8,
                            uiOutput("database_textbox")
                     ),
                     column(2,
                            fileInput("database_upload", "or UPLOAD")
                     ),
                     column(2,
                            uiOutput("checkInput")
                     )
                   ),
                   fluidRow(
                     column(11,
                            uiOutput("dataset_description")
                     )
                   ),
               ),
               
               box(title = "Select cell types for analysis and analysis type",
                   solidHeader = TRUE, status = "primary", width = 12,
                   collapsible = TRUE, collapsed = FALSE, height= 800,
                   fluidRow(
                     column(5,
                            selectInput("hierarchy_level","Choose level of hierarchy:",choices=NULL),
                            uiOutput("plot_type_selection"),
                            selectInput(
                              inputId = "background_type",
                              label = "Choose analysis type:",
                              choices = c("Foreground vs. local types",
                                          "Foreground vs. custom types",
                                          "Foreground vs. all other types",
                                          "Trajectory analysis"),
                              selected = "Foreground vs. local types"
                            ),
                            h4("Foreground cell types:"),
                            verbatimTextOutput("currentFilterIDs"),  
                            fluidRow(
                              column(6,
                                     actionButton("clearFilter", "Clear Foreground Filter")
                              ),
                              column(6,
                                     actionButton(
                                       "find_degenes",
                                       "Find Relevant Genes",
                                       style = "color: #fff; background-color: #39B54A; border-color: #006838; font-weight: bold;"),
                              )
                            ),
                            
                            p(),
                            uiOutput("conditional_list_selection"),
                            uiOutput("conditional_background_title"),
                            uiOutput("conditional_background_filter"),
                            uiOutput("conditional_background_clear")
                            
                     ),
                     column(7,
                            plotlyOutput("sunburst", height = "740"),
                     ),
                     
                   ),
             ),
      ),
      
      
      conditionalPanel(
        condition = "input.find_degenes > 0",
        hr(style = "border-top: 3px solid #000000;"),
        h3("Gene analysis results"),
        uiOutput("gene_analysis_results_text"), 
        hr(style = "border-top: 1px dotted #777777;")
      ),
      
      fluidRow(column(
        12,
        withSpinner(DT::dataTableOutput("de_table")),
        verbatimTextOutput("click_info")
      )),
      
      fluidRow(column(
        12,
        withSpinner(plotOutput("dotplot", height = "800px"))
      )),   
      
      
      fluidRow(
        div(
          style = "display: flex; gap: 10px;", # Use a flexible container
          div(
            uiOutput("download_table_button")
          ),
          div(
            uiOutput("gene_set_enrichment_button")
          ),
          div(
            conditionalPanel(
              condition = "input.gene_set_enrichment > 0",
              p("Enrichment is a go!")
            ),
          ),
          div(
            textOutput("processing_message")  # This isn't working, but it's also not hurting anything
          ),
        ) 
      ),
      
      conditionalPanel(
        condition = "input.gene_set_enrichment > 0",
        plotOutput("enrichment_plot",height = "400px")
      ),
      
      fluidRow(width = 12, br(), br())
      
    )
  )
  
}