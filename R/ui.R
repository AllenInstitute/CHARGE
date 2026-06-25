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
      
      dashboardHeader(
        title = tags$a(
          href = "https://alleninstitute.org",
          target = "_blank",
          tags$img(src = "allen_institute_logo.svg", alt = "Allen Institute", height = "30")
        ),
        titleWidth = 250,
        tags$li(
          class = "dropdown charge-app-title",
          tags$span("Cell HierARchy Gene Explorer (CHARGE)")
        ),
        tags$li(
          class = "dropdown charge-nav-link",
          tags$a(
            href = "https://alleninstitute.org",
            target = "_blank",
            icon("building"), " Allen Institute"
          )
        ),
        tags$li(
          class = "dropdown charge-nav-link",
          tags$a(
            href = "https://brain-map.org",
            target = "_blank",
            icon("globe"), " Brain Map"
          )
        )
      ),
      
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
        
        tags$head(
          tags$style(HTML("
            /* ==========================================================
               CHARGE CUSTOM HEADER
               ========================================================== */
        
            .skin-blue .main-header .logo,
            .skin-blue .main-header .navbar {
              background-color: #000000 !important;
            }
        
            .skin-blue .main-header .logo:hover,
            .skin-blue .main-header .navbar .sidebar-toggle:hover {
              background-color: rgba(255,255,255,0.15) !important;
            }
        
            .main-header .logo {
              display: flex !important;
              align-items: center;
              justify-content: center;
              height: 60px !important;
              padding: 0 18px !important;
            }
        
            .main-header .navbar,
            .main-header .sidebar-toggle {
              min-height: 60px !important;
              height: 60px !important;
            }
        
            .main-header .sidebar-toggle {
              display: flex !important;
              align-items: center;
              padding: 0 18px !important;
              color: #ffffff !important;
            }
        
            .main-header .navbar .nav > li > a {
              color: #ffffff !important;
              padding-top: 20px !important;
              padding-bottom: 20px !important;
              line-height: 20px !important;
            }
           
            .main-header .navbar-custom-menu .navbar-nav > li.charge-app-title {
              position: fixed !important;
              top: 0;
              left: 50vw;
              transform: translateX(-50%);
              height: 60px;
              line-height: 60px;
              z-index: 2000;
              float: none !important;
              padding: 0;
              margin: 0;
              color: #ffffff;
              font-size: 24px;      /* adjust larger/smaller here */
              font-weight: 700;
              white-space: nowrap;
              pointer-events: none; /* keeps links/buttons clickable underneath */
              max-width: calc(100vw - 520px);
              overflow: hidden;
              text-overflow: ellipsis;
            }
            
            .charge-app-title > span {
              display: block;
            }
        
            .charge-nav-link > a {
              font-size: 14px;
              border-radius: 3px;
            }
        
            .charge-nav-link > a:hover {
              background-color: rgba(255,255,255,0.15) !important;
            }
        
            @media (max-width: 900px) {
              .main-header .navbar-custom-menu .navbar-nav > li.charge-app-title {
                font-size: 16px;
                max-width: calc(100vw - 120px);
              }

              .charge-nav-link {
                display: none !important;
              }
            }
        
            /* ==========================================================
               LAYOUT SPACING
               ========================================================== */
        
            .content-wrapper,
            .right-side {
              padding-top: 0px;
            }
        
            .main-sidebar,
            .left-side {
              padding-top: 70px !important;
            }
        
            /* ==========================================================
               SIDEBAR
               ========================================================== */
        
            .main-sidebar {
              background-color: #000000 !important;
            }
        
            /* ==========================================================
               BOX HEADERS
               ========================================================== */
        
            .box.box-solid.box-primary > .box-header,
            .box.box-solid.box-primary > .box-header .box-title,
            .box.box-solid.box-primary > .box-header a {
              background-color: #000000 !important;
              color: #ffffff !important;
            }
        
            .box.box-solid.box-primary {
              border-color: #000000 !important;
            }
        
            /* ==========================================================
               ENRICHMENT TABLE
               ========================================================== */
        
            #enrichment_table td {
              white-space: normal !important;
              height: auto;
              line-height: normal;
            }
        
            #enrichment_table table.dataTable thead th {
              vertical-align: top;
            }
        
            /* ==========================================================
               SHINY NOTIFICATIONS
               ========================================================== */
        
            .shiny-notification {
              position: fixed;
              top: 30%;
              left: 35%;
              width: 30%;
              font-size: 24px;
              font-weight: bold;
              background-color: #800080 !important;
              color: white !important;
              border: 5px solid #301934;
            }
            
          "))
        ),
        
        #useShinyjs(),  # shinyjs not currently used
        
        fluidRow(width = 12,
                 
                 box(title = "Select data set",
                     solidHeader = TRUE, status = "primary", width = 12,
                     collapsible = TRUE, collapsed = FALSE,
                     fluidRow(
                       column(11,
                              uiOutput("dataset_description")
                       )
                     ),
                     fluidRow(
                       #column(6,
                       #        uiOutput("select_category")
                       # ),
                       column(8,
                              selectizeInput(
                                inputId = "select_textbox",
                                label   = "Choose an existing data set (or select 'Enter your own location')",
                                choices = NULL
                              )
                       ),
                       column(4,
                              uiOutput("checkInput")
                       )
                     ),
                     fluidRow(
                       column(3,
                              fileInput("database_upload", "UPLOAD")
                       ),
                       column(6,
                              uiOutput("database_textbox")
                       ),
                       bookmarkButton(label = "Save Current View") # The button
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
                              p("After adjusting your filters, press the green or purple button and then scroll down."),
                              fluidRow(
                                column(4,
                                       actionButton("clearFilter", "Clear Foreground Filter")
                                ),
                                column(4,
                                       actionButton(
                                         "find_degenes",
                                         "Find Relevant Genes",
                                         style = "color: #fff; background-color: #39B54A; border-color: #006838; font-weight: bold;"),
                                ),
                                column(4,
                                       actionButton(
                                         "known_genes",
                                         "Provide a gene list",
                                         style = "color: #fff; background-color: #8E4585; border-color: #006838; font-weight: bold;"),
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
          hr(style = "border-top: 1px dotted #777777;"),
          
          fluidRow(column(
            12,
            withSpinner(DT::dataTableOutput("de_table")),
            verbatimTextOutput("click_info")
          )),
          
          fluidRow(column(
            12,
            withSpinner(plotOutput("dotplot", height = "800px"))
          )),
          
          # For downloading the dot plot
          fluidRow(
            column(8,
                   strong("Options and download button for the DOT OR TRAJECTORY PLOT:"),
                   fluidRow(
                     column(3,textInput("dlw","Width (in)",12)),
                     column(3,textInput("dlh","Height (in)",8)),
                     column(3,textInput("dlf","Font (pt)",10)),
                     column(3,downloadButton('CHARGE_gene_plot',"Download PDF"))
                   )
            )
          )
        ),
        
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
                p("Enrichment is a go! This process can take up to a minute.")
              ),
            ),
            div(
              textOutput("processing_message")  # This isn't working, but it's also not hurting anything
            ),
          ) 
        ),
        
        conditionalPanel(
          condition = "input.gene_set_enrichment > 0",
          plotOutput("enrichment_plot",height = "400px"),
          
          # For downloading the enrichment plot
          fluidRow(
            column(8,
                   strong("Options and download button for the GO ENRICHMENT PLOT:"),
                   fluidRow(
                     column(3,textInput("enrichment_dlw","Width (in)",12)),
                     column(3,textInput("enrichment_dlh","Height (in)",8)),
                     column(3,textInput("enrichment_dlf","Font (pt)",10)),
                     column(3,downloadButton('CHARGE_enrichment_plot',"Download PDF"))
                   )
            )
          ),
          
        ),
        
        
        fluidRow(width = 12, br(), br())
        
      )
    )
    
  }
  