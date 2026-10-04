# Libraries are loaded in ui.R
options(stringsAsFactors = F)
options(shiny.maxRequestSize = 1.9 * 1024^3)  # 1.9GB max For uploading files, which I believe is just below the browser limit

source("initialization.R")
source("sunburst.R")
source("de_genes_functions.R")
source("group_dot_plots.R")

guess_type <- function(x) {
  if(try(sum(is.na(as.numeric(x))) > 0,silent = T)) {
    "cat"
  } else {
    "num"
  }
}

## DEFAULT VALUES
default_vals <- list(db = "Enter a file path or URL here, or choose from dropdown above.",
                     sf = "Enter a file path or URL here, or choose from dropdown above.",
                     list_selection = "Foreground",
                     plot_selection = "Sunburst"
)

## ORTHOLOGS FOR GENE SET ENRICHMENT
orthologs <- vroom("https://github.com/AllenInstitute/GeneOrthology/raw/05300f8dd1fdbe38db56d318337c682437227974/csv/mammalian_orthologs_20231113.csv")
orthologs <- as.data.frame(orthologs)

## TOOL TIPS FOR THE de_table
tooltip_data <- read.csv("table_column_definitions.csv", stringsAsFactors = FALSE)

######################################################
## Default table information is in initialization.R ##
######################################################

server <- function(input, output, session) {

  ###########################
  ##  State Initialization ##
  ###########################
  
  init <- reactiveValues(vals = list())
  
  # Build initial values list
  # These are used to set the state of the input values for UI elements
  
  # First from default_vals,
  # then dropdown_vals,
  # then from URL parsing

  observe({
    vals <- default_vals
    
    if (nzchar(session$clientData$url_search)) {
      query <- parseQueryString(session$clientData$url_search)
      
      query <- lapply(query, function(x) {
        if (length(x) == 0) return(NULL)
        gsub('^"|"$', '', x)
      })
      
      for (nm in names(query)) {
        vals[[nm]] <- query[[nm]]
      }
    }
    
    init$vals <- vals
  })
  
  
  observe({  
   # Get all current input IDs
   all_inputs <- names(input)

   # Keep only app states related to current selection
   keep_list <- c(
     # "Select data set" box
     "select_textbox",
     "db",
     
     # "Select cell types for analysis and analysis type" box
     "hierarchy_level",
     "local_context_level",
     "plot_selection" ,
     "background_type"
     
   )
  
   exclude_list <- all_inputs[!(all_inputs %in% keep_list)]
   
   setBookmarkExclude(exclude_list)
  })
  
  onBookmark(function(state) {
    state$values$foreground <- paste(rv_sunburst$selected_nodes$foreground, collapse = "|")
    state$values$background <- paste(rv_sunburst$selected_nodes$background, collapse = "|")
  })
  
  onRestore(function(state) {
    fg <- state$values$foreground
    bg <- state$values$background
    
    if (!is.null(fg) && nzchar(fg)) {
      rv_sunburst$selected_nodes$foreground <- strsplit(fg, "\\|")[[1]]
    }
    
    if (!is.null(bg) && nzchar(bg)) {
      rv_sunburst$selected_nodes$background <- strsplit(bg, "\\|")[[1]]
    }
  })
  
  ## NEW BOOKMARKING FUNCTIONALITY
  observeEvent(input$save_view, {
    session$doBookmark()
  })
  
  onBookmarked(function(url) {
    showModal(modalDialog(
      title = "Share this link",
      easyClose = TRUE,
      tags$p("Press Ctrl+A then Ctrl+C to copy the link below:"),
      tags$textarea(
        style = "width:100%; height:80px;",
        url
      )
    ))
  })
  
  
  observeEvent(init$vals, {
    req(table_info$table_name)
    
    choices <- c(
      "Select data set...",
      table_info$table_name,
      "Enter your own location"
    )
    
    selected <- init$vals[["select_textbox"]]
    if (length(selected) == 0) selected <- NULL
    
    if (!is.null(selected) && !(selected %in% choices)) {
      warning("select_textbox value not found in choices: ", selected)
      selected <- NULL
    }
    
    updateSelectizeInput(
      session,
      inputId = "select_textbox",
      choices = choices,
      selected = selected,
      server = TRUE
    )
  }, once = TRUE)
  
  
  output$database_textbox <- renderUI({
    req(init$vals)
    
    id      <- "db"
    label   <- "Location of 'CHARGE'd data set"
    initial <- input$Not_on_list
    upload2 <- input$database_upload
    
    if (length(input$select_textbox)>0){
      if (!is.element(input$select_textbox,c("Select data set...",'Enter your own location'))) {
        initial = table_info[table_info$table_name==input$select_textbox,"table_loc"]
      }
      if((input$select_textbox=='Enter your own location')&(!is.null(upload2))){
        initial <- normalizePath(upload2$datapath)
      }
    }
    
    textInput(inputId = id, 
              label = strong(label), 
              value = initial, 
              width = "100%")
    
  })
  
  
  # This function adds the data set description AND hides irrelevant visualization panels for preset data input
   output$dataset_description <- renderUI({
    req(init$vals)
    
    header_text = "READ ME"
    text_desc = "Select a category and a data set from the boxes above -OR- to compare your own annotation data, choose 'Enter your own location' from the 'Select annotation category' and enter the location of your data file or upload a file yourself. Once a data set is selected, wait for the panels below to refresh."
    
    print("input$select_textbox")
    print(input$select_textbox)
    
    if (length(input$select_textbox)>0){
    
        if (input$select_textbox == 'Enter your own location') {
          header_text = "Upload user-provided data"
          text_desc = "User-provided data set file, created using the 'chargeTaxonomy' R function (see GitHub page for details)."
        } else if (input$select_textbox == 'Select data set...') {
          # Do nothing... text_desc should remain as initialized above
          header_text = "READ ME"
        } else {
          header_text = "Dataset description"
          text_desc = table_info[table_info$table_name==input$select_textbox,"description"]
        }
      #}
    }
    
    # Return the description text
    div(style = "font-size:14px;", strong(header_text),br(),text_desc,p(" "))
    
  })
  

  
  ##################################
  ## Loading tables from input$db ##
  ##################################
  
  # Check the path provided by input$db
  # returns a corrected path
  #
  # note: rv_ prefix stands for reactive value
  #
  # rv_path() - length 1 character vector
  #
  rv_path <- reactive({
    req(input$db)
    write("Checking and setting input$db.", stderr())
    
    input$db
  })
  
  
  
  # Read the CELL annotations table from the dataset
  #
  # rv_anno() - a data.frame
  #
   rv_anno <- reactive({
     req(rv_path())
     
     file <- rv_path()
     write("Reading file.", stderr())
     
     withProgress(
       message = "Loading data set...",
       detail = "Please wait.",
       value = NULL,
       {
         
         # THIS IS WHERE THE DATA GETS READ IN.
         if (substr(file, 1, 2) == "s3") {
           
           ## READ FROM s3 bucket
           file2 <- substr(file, 6, 10000)
           file2 <- strsplit(file2, "/")[[1]]
           bucket <- file2[1]
           filename <- paste(file2[2:length(file2)], collapse = "/")
           
           write("filename:", stderr())
           write(filename, stderr())
           write("bucket:", stderr())
           write(bucket, stderr())
           
           objIn <- objects()
           a <- try({
             s3load(object = filename, bucket = bucket)
           })
           
           if (inherits(a, "try-error")) {
             write(
               paste("s3", file, "does not exist or cannot be accessed."),
               stderr()
             )
             return(NULL)
           }
           
           objOut <- objects()
           objs <- setdiff(objOut, objIn)
           
         } else {
           
           ## READ LOCALLY
           if (file.exists(file)) {
             objs <- load(file)
           } else {
             write(
               paste("Local", file, "does not exist."),
               stderr()
             )
             return(NULL)
           }
         }
         
         eval(parse(
           text = paste0(
             "data=list(",
             paste(objs, collapse = ","),
             ")"
           )
         ))
         
         names(data) <- objs
         data
       }
     )
   })
  
  
  # Build the annotation descriptions table
  # 
  # rv_desc() - a data.frame
  #
  rv_desc <- reactive({
    req(rv_anno())
    write("Building desc.", stderr())
    
    data <- rv_anno()
    anno <- data$cluster_info
    names <- colnames(anno)[grepl("_label$",colnames(anno))]
    names <- substr(names,1,nchar(names)-6)
    desc_table <- data.frame(base=names,name=names)
    
    suppressWarnings({
      desc_table <- desc_table %>%
        rowwise() %>%
        mutate(type = guess_type(anno[[paste0(base,"_label")]]))
    })
    
    return(desc_table)
    
  }) # end of rv_desc()
  
  
  
  #############################
  ##    Define hierarchy     ##
  #############################
  
  
  rv_hierarchy_options <- reactive({
    data <- rv_anno()
    data$hierarchy
  })
  
  observeEvent(rv_hierarchy_options(), {
    hierarchy_options = rv_hierarchy_options()
    updateSelectInput(session, 
                      inputId = "hierarchy_level", 
                      label = "Choose level of hierarchy:", 
                      choices = hierarchy_options,
                      selected = hierarchy_options[1]
    )
    
  })
  
  output$local_context_level_ui <- renderUI({
    req(rv_hierarchy_options())
    req(input$hierarchy_level)
    
    if (input$background_type != "Foreground vs. local types") {
      return(NULL)
    }
    
    hierarchy_options <- rv_hierarchy_options()
    current_level_index <- match(
      input$hierarchy_level,
      hierarchy_options
    )
    
    if (
      is.na(current_level_index) ||
      current_level_index >= length(hierarchy_options)
    ) {
      return(
        helpText(
          "No higher hierarchy level is available for local context."
        )
      )
    }
    
    higher_levels <- hierarchy_options[
      (current_level_index + 1):length(hierarchy_options)
    ]
    
    current_selection <- isolate(input$local_context_level)
    
    if (
      is.null(current_selection) ||
      !(current_selection %in% higher_levels)
    ) {
      current_selection <- higher_levels[1]
    }
    
    selectInput(
      inputId = "local_context_level",
      label = "Choose level of local context:",
      choices = higher_levels,
      selected = current_selection
    )
  })
  
  observeEvent(
    list(rv_anno(), input$hierarchy_level),
    {
      req(rv_anno())
      req(input$hierarchy_level)
      
      data <- rv_anno()
      label_column <- paste0(input$hierarchy_level, "_label")
      
      req(label_column %in% colnames(data$cluster_info))
      
      cell_type_choices <- unique(
        as.character(data$cluster_info[[label_column]])
      )
      
      cell_type_choices <- cell_type_choices[
        !is.na(cell_type_choices) &
          nzchar(cell_type_choices)
      ]
      
      current_foreground <- isolate(input$manual_foreground_types)
      current_comparison <- isolate(input$manual_comparison_types)
      
      updateSelectizeInput(
        session,
        inputId = "manual_foreground_types",
        choices = cell_type_choices,
        selected = intersect(current_foreground, cell_type_choices),
        server = TRUE
      )
      
      updateSelectizeInput(
        session,
        inputId = "manual_comparison_types",
        choices = cell_type_choices,
        selected = intersect(current_comparison, cell_type_choices),
        server = TRUE
      )
    },
    ignoreInit = FALSE
  )
  
  
  
  ######################################################################
  ##      Constellation plots, Sunburst plots, and plot selection     ##
  ######################################################################
  
  
  
  output$plot_type_selection <- renderUI({

    radioButtons(
      inputId = "plot_selection",
      label = "Choose plot selection type:",
      choices = list(
        "Sunburst" = "Sunburst",
        "Constellation" = "Constellation",
        "Manual entry" = "Manual entry"
      ),
      selected = "Sunburst", 
      inline = TRUE # Display buttons side-by-side
    )
    
  })

  rv_sunburst <- reactiveValues(
    selected_nodes = list(foreground = character(0), background = character(0))
  )
  
  # NOTE, THIS FUNCTION CALLS BOTH SUNBURST AND CONSTELLATION PLOTS!
  output$sunburst <- renderPlotly({
    req(rv_anno())
    data <- rv_anno()
    constellation <- data$constellation
    
    if(input$plot_selection=="Sunburst"){
      ## SUNBURST PLOT GENERATION
      write("Building sunburst plot", stderr())
      
      # Define the hierarchy based on input selection
      sunburst_hierarchy = data$hierarchy[length(data$hierarchy):1]    
      level = which(sunburst_hierarchy==input$hierarchy_level)
      if(length(level)==1)
        sunburst_hierarchy = sunburst_hierarchy[1:level]
      
      sunburstDF <- as.sunburstDF(data$cluster_info, sunburst_hierarchy,rootname="all")

      p <- plot_ly() %>%
        add_trace(ids = sunburstDF$ids,
                  labels = sunburstDF$labels,
                  parents =sunburstDF$parent,
                  values = sunburstDF$values,
                  type = 'sunburst',
                  sort=FALSE,
                  marker = list(colors = sunburstDF$color),
                  domain = list(column = 1),
                  branchvalues = 'total'
        )%>%
        layout(grid = list(columns =1, rows = 1),
               margin = list(l = 0, r = 0, b = 0, t = 0)
        )
      
      p <- htmlwidgets::onRender(
        p,
        "
          function(el, x) {
            el.on('plotly_sunburstclick', function(d) {
              if (d.points && d.points.length > 0) {
                Shiny.setInputValue(
                  'sunburst_node_click',
                  {
                    label: d.points[0].label,
                    nonce: Date.now()
                  },
                  {priority: 'event'}
                );
              }
              return false;
            });
          }
        "
      )
      
    } else {
      ## CONSTELLATION PLOT GENERATION
      write("Building constellation plot", stderr())
      p <- constellation[[input$hierarchy_level]]
      
      if(is.null(p)){
        p <- plot_ly() %>%
          add_annotations(
            text = paste("No constellation diagram for",input$hierarchy_level),
            x = 0.5, y = 0.5,          # Center coordinates
            xref = "paper", yref = "paper", # Relative to plot area
            showarrow = FALSE,
            font = list(size = 18, color = "black") # Basic font styling
          ) %>%
          layout(
            xaxis = list(visible = FALSE), # Hide X-axis
            yaxis = list(visible = FALSE), # Hide Y-axis
            # Optional: make background transparent if embedding or don't want default gray
            plot_bgcolor = 'rgba(0,0,0,0)',
            paper_bgcolor = 'white'
          )
      }
    }
    
    # RETURN PLOT
    #p
    event_register(p, "plotly_click")
    
  })
  

  # This function sets the selected nodes
  observeEvent(
    list(
      input$sunburst_node_click,
      event_data("plotly_click")
    ),
    {
    
    req(rv_anno())
    data <- rv_anno()
    constellation <- data$constellation
    
    # Register the event
    d <- NULL
    
    if (input$plot_selection != "Sunburst") {
      d <- event_data("plotly_click")
      req(d)
    }
    
    if (input$plot_selection == "Sunburst") {
      
      req(input$sunburst_node_click$label)
      clicked_node_id <- input$sunburst_node_click$label
      
    } else {
      
      dat <- constellation[[input$hierarchy_level]]$x$layoutAttrs[[1]]$annotations
      xval <- as.numeric(lapply(dat, function(x) x$x))
      yval <- as.numeric(lapply(dat, function(x) x$y))
      kp <- which((xval == d$x) & (yval == d$y))
      
      clicked_node_id <- as.character(
        lapply(dat, function(x) x$text)
      )[kp[1]]
    }
    
    write(clicked_node_id, stderr())
    
    # Determine whether the click modifies foreground or comparison types
    which_list <- "foreground"
    
    if (
      input$background_type == "Foreground vs. custom types" &&
      input$list_selection == "Comparison"
    ) {
      which_list <- "background"
    }
    
    current_selected <- rv_sunburst$selected_nodes[[which_list]]
    
    if (input$plot_selection == "Sunburst") {
      
      target_level <- input$hierarchy_level
      target_column <- paste0(target_level, "_label")
      cluster_info <- data$cluster_info
      
      req(target_column %in% colnames(cluster_info))
      
      hierarchy_columns <- paste0(data$hierarchy, "_label")
      hierarchy_columns <- intersect(
        hierarchy_columns,
        colnames(cluster_info)
      )
      
      if (identical(clicked_node_id, "all")) {
        
        clicked_descendants <- unique(
          as.character(cluster_info[[target_column]])
        )
        
      } else {
        
        clicked_columns <- hierarchy_columns[
          vapply(
            hierarchy_columns,
            function(column_name) {
              clicked_node_id %in%
                as.character(cluster_info[[column_name]])
            },
            logical(1)
          )
        ]
        
        if (length(clicked_columns) == 0) {
          
          clicked_descendants <- character(0)
          
        } else {
          
          descendant_rows <- Reduce(
            `|`,
            lapply(
              clicked_columns,
              function(column_name) {
                as.character(cluster_info[[column_name]]) ==
                  clicked_node_id
              }
            )
          )
          
          clicked_descendants <- unique(
            as.character(
              cluster_info[descendant_rows, target_column]
            )
          )
        }
      }
      
      clicked_descendants <- clicked_descendants[
        !is.na(clicked_descendants) &
          nzchar(clicked_descendants)
      ]
      
      if (
        length(clicked_descendants) > 0 &&
        all(clicked_descendants %in% current_selected)
      ) {
        
        # Clicking an already selected branch removes all its descendants
        selected_nodes <- current_selected[
          !(current_selected %in% clicked_descendants)
        ]
        
      } else {
        
        # Add missing descendants in their cluster_info order
        selected_nodes <- c(
          current_selected,
          clicked_descendants[
            !(clicked_descendants %in% current_selected)
          ]
        )
      }
      
    } else {
      
      # Preserve the existing Constellation behavior
      if (clicked_node_id %in% current_selected) {
        selected_nodes <- setdiff(
          current_selected,
          clicked_node_id
        )
      } else {
        selected_nodes <- unique(
          c(current_selected, clicked_node_id)
        )
      }
    }
    
    # Retain only valid types from the currently selected hierarchy level
    level <- input$hierarchy_level
    
    if (length(level) == 0) {
      level <- data$hierarchy[1]
    }
    
    all_types <- unique(
      as.character(
        data$cluster_info[[paste0(level, "_label")]]
      )
    )
    
    selected_nodes <- selected_nodes[
      selected_nodes %in% all_types
    ]
    
    rv_sunburst$selected_nodes[[which_list]] <- selected_nodes
    
    },
    ignoreInit = TRUE
  )
  
  ## FOREGROUND FILTERS
  
  output$currentFilterIDs <- renderPrint({
    if (length(rv_sunburst$selected_nodes$foreground) == 0) {
      "None selected."
    } else {
      data.frame(cell_type=rv_sunburst$selected_nodes$foreground)
    }
  })
  
  observeEvent(input$clearFilter, {
    rv_sunburst$selected_nodes$foreground <- character(0) # Reset the filter
  })
  
  observeEvent(input$replace_foreground_types, {
    rv_sunburst$selected_nodes$foreground <-
      as.character(input$manual_foreground_types)
  })
  
  ## BACKGROUND FILTERS
  
  
  output$conditional_background_title <- renderUI({
    
    if(input$background_type=="Foreground vs. custom types"){
      h4("Comparison cell types (e.g., the ones used as background):")
    } else if(input$background_type=="Trajectory analysis"){
      return("Trajectory analysis can take a up to about a minute to run. Please be patient!")
    } else {
      return("(Comparison types automatically selected.)")
    }
    
  })
  
  output$conditional_background_filter <- renderUI({
    
    if(input$background_type!="Foreground vs. custom types")
      return(NULL)
    
    verbatimTextOutput("currentBackgroundFilterIDs")
  })
  
  
  output$currentBackgroundFilterIDs <- renderPrint({
    
    if (length(rv_sunburst$selected_nodes$background) == 0) {
      "None selected."
    } else {
      data.frame(cell_type=rv_sunburst$selected_nodes$background)
    }
    
  })
  
  
  output$conditional_background_clear <- renderUI({
    
    if(input$background_type!="Foreground vs. custom types")
      return(NULL)
    
    actionButton("conditional_background_clear", "Clear Comparison Filter")
    
  })
  
  observeEvent(input$conditional_background_clear, {
    rv_sunburst$selected_nodes$background <- character(0) # Reset the filter
  })
  
  observeEvent(input$replace_comparison_types, {
    req(input$background_type == "Foreground vs. custom types")
    
    rv_sunburst$selected_nodes$background <-
      as.character(input$manual_comparison_types)
  })
  
  output$conditional_list_selection <- renderUI({
    
    if(input$background_type!="Foreground vs. custom types")
      return(NULL)
    
    radioButtons(
      inputId = "list_selection",
      label = "Choose cell type for:",
      choices = list(
        "Foreground" = "Foreground",
        "Comparison" = "Comparison"
      ),
      selected = "Foreground", # Blue will be pre-selected (by its value "B")
      inline = TRUE # Display buttons side-by-side
    )
    
  })
 
  

  ##################################################
  #######   DIFFERENTIAL GENE CALULATIONS    #######
  #######              - OR -                #######
  #######   PROVIDING LIST OF KNOWN GENES    #######
  ##################################################
  
  
  # 1. Create a "tracker" to store which button was clicked last
  last_clicked <- reactiveVal(NULL)
  
  # 2. Use observers to update the tracker
  observeEvent(input$find_degenes, { last_clicked("find_degenes") })
  
  # If the known_genes button is pressed, pop up a box to input the gene list
  observeEvent(input$known_genes, { 
    
    showModal(modalDialog(
      title = "Input Data",
      textInput("raw_input", "Enter items (comma-separated):", 
                placeholder = "e.g. GAPDH, APOE, CD4"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submit_data", "Save & Process", class = "btn-primary")
      )
    ))
    
  })
  
  # If filters are cleared, reset this
  observeEvent(input$clearFilter, { last_clicked(NULL) })
  
  
  observeEvent(input$submit_data, {
    # Validation: Ensure input isn't empty
    req(input$raw_input)
    
    last_clicked(input$raw_input)
    
    # Close the modal
    removeModal()
  })
  
  
  calculate_de_genes <- reactive({  #eventReactive(input$find_degenes, {  # input$known_genes
    
    req(last_clicked()) # Ensure a button has been clicked
    req(rv_anno())
    
    if(last_clicked()=="find_degenes"){
      if(length(rv_sunburst$selected_nodes$foreground)>0)
        
        if(!((length(rv_sunburst$selected_nodes$background)==0)&(input$background_type=="Foreground vs. custom types"))){
          data <- rv_anno()
          
          if(input$background_type=="Trajectory analysis"){
            
            find_trajectory_genes(
              data,
              rv_sunburst$selected_nodes$foreground,
              filter = identical(input$gene_return_mode, "fast")
            )
            
          } else {
            
            find_de_genes(
              data,
              input,
              rv_sunburst$selected_nodes$foreground,
              rv_sunburst$selected_nodes$background,
              filter = identical(input$gene_return_mode, "fast")
            )
            
          }
        }
    } else {
      if(length(rv_sunburst$selected_nodes$foreground)>0)
        
        if(!((length(rv_sunburst$selected_nodes$background)==0)&(input$background_type=="Foreground vs. custom types"))){
          data <- rv_anno()
          
          input_gene_set <- strsplit(last_clicked(), ",")[[1]]
          input_gene_set <- trimws(input_gene_set)
          write("input_gene_set",stderr())
          write(input_gene_set,stderr())
          
          if (input$background_type == "Visualize known genes") {
            return(
              create_known_gene_table(
                data,
                rv_sunburst$selected_nodes$foreground,
                in_genes = input_gene_set
              )
            )
          }
          
          if(input$background_type=="Trajectory analysis"){
            
            find_trajectory_genes(
              data,
              rv_sunburst$selected_nodes$foreground,
              in_genes = input_gene_set,
              filter = TRUE
            )
            
          } else {
            
            find_de_genes(
              data,
              input,
              rv_sunburst$selected_nodes$foreground,
              rv_sunburst$selected_nodes$background,
              in_genes = input_gene_set,
              filter = TRUE
            )
            
          }
        }
    }
  })

  gene_link_lookup <- reactive({
    req(input$select_textbox)
    
    link_name <- table_info[
      table_info$table_name == input$select_textbox,
      "web_urls"
    ]
    
    if (length(link_name) != 1 ||
        is.na(link_name) ||
        !nzchar(link_name)) {
      return(character(0))
    }
    
    link_file <- file.path("links", paste0(link_name, ".csv.gz"))
    
    if (!file.exists(link_file)) {
      warning("Gene link file not found: ", link_file)
      return(character(0))
    }
    
    link_table <- read.csv(
      gzfile(link_file),
      stringsAsFactors = FALSE
    )
    
    if (!all(c("gene", "url") %in% colnames(link_table))) {
      warning("Gene link file must contain 'gene' and 'url' columns: ", link_file)
      return(character(0))
    }
    
    link_table <- link_table[
      !is.na(link_table$gene) &
        !is.na(link_table$url) &
        !duplicated(link_table$gene),
      c("gene", "url")
    ]
    
    setNames(link_table$url, link_table$gene)
  })  

  output$de_table <- renderDataTable({
    req(calculate_de_genes())
    data_df = calculate_de_genes()
    
    if (input$background_type == "Trajectory analysis") {
      
      # Preferred trajectory-table column order
      preferred_order <- c(
        "gene",
        "mean.expression",
        "WLS_Slope",
        "WLS_T_Value",
        "WLS_P_Value",
        "WLS_FDR", 
        "gene_categories________________________________________________________"
      )
      
    } else if (input$background_type == "Visualize known genes") {
      
      # Preferred known-gene-table column order
      preferred_order <- c(
        "gene",
        "mean.expression",
        "gene_categories________________________________________________________"
      )
      
    } else {
      
      # Preferred differential-gene-table column order
      preferred_order <- c(
        "gene", 
        "consensus_score",
        "propMeanScore", 
        "prop_diff", 
        "log2_FC", 
        "gr1_prop", 
        "gr1_mean", 
        "gr2_prop", 
        "gr2_mean", 
        "rank_biserial_corr",
        "overlap_coefficient", 
        "gene_categories________________________________________________________"
      )
      
    }
    
    # keep only columns that actually exist
    preferred_order <- intersect(preferred_order, colnames(data_df))
    
    # reorder dataframe
    data_df <- data_df[, preferred_order, drop = FALSE]
    
    ## Dynamically determine tool tip definitions
    column_definitions <- sapply(colnames(data_df), function(col_name) {
      # Find the matching tooltip from the loaded data
      match <- tooltip_data[tooltip_data$column_names == col_name, "column_definitions"]
      # If a match is not found, provide a default tooltip
      if (length(match) == 0) {
        return(paste("No definition available for", col_name))
      }
      return(match)
    }, USE.NAMES = FALSE)
    
    print(cbind(colnames(data_df),column_definitions))
    
    gene_links <- gene_link_lookup()
    gene_names <- as.character(data_df$gene)
    gene_urls <- unname(gene_links[gene_names])
    
    data_df$gene <- mapply(
      FUN = function(gene, url) {
        if (is.na(url) || !nzchar(url)) {
          return(as.character(htmltools::htmlEscape(gene)))
        }
        
        as.character(
          htmltools::tags$a(
            href = url,
            target = "_blank",
            rel = "noopener noreferrer",
            gene
          )
        )
      },
      gene = gene_names,
      url = gene_urls,
      USE.NAMES = FALSE
    )
    
    datatable(data_df, 
              rownames = FALSE, 
              filter = "top", 
              escape = which(colnames(data_df) != "gene"),
              options = list(
                scrollX = TRUE,
                scrollY = TRUE,
                pageLength = 10,
                lengthMenu = list(c(5, 10, 20), c("5", "10", "20")),
                headerCallback = JS(
                  paste0(
                    "function(thead, data, start, end, display) {",
                    "  var tooltips = ", toJSON(column_definitions), ";",
                    "  $(thead).find('th').each(function(i) {",
                    "    this.setAttribute('title', tooltips[i]);",
                    "  });",
                    "}"
                  )
                )
              )
    )
    
  }, server=TRUE)
  
  output$download_table <- downloadHandler(
    
    filename = function() { paste0(input$background_type,"_results.csv") },
    content = function(file) {
      write.csv(calculate_de_genes()[input$de_table_rows_all,],file)
    }
    
  )
  
  output$download_table_button <- renderUI(
    if(isTruthy(calculate_de_genes())) {
      
      downloadButton("download_table","Download Table", class = "downloads")
      
    }
  )

  
  
  get_genes_dotplot <- reactive({
    req(calculate_de_genes())
    
    cat("de genes dot plot \n")
    
    de_table <- calculate_de_genes()
    
    current_rows <- input$de_table_rows_current
    
    if (is.null(current_rows) || length(current_rows) == 0) {
      current_rows <- seq_len(nrow(de_table))
    }
    
    current_de_table <- de_table[
      current_rows,
      ,
      drop = FALSE
    ]
    # print(current_de_table)
    
    top10_genes <- as.character(
      head(current_de_table[["gene"]], 100)
    )
    
    top10_genes <- top10_genes[
      !is.na(top10_genes) & nzchar(top10_genes)
    ]
    
    top10_genes
    
  })
  
  
  # Group dot plots
  # NOTE:  THIS ALSO IS THE SAME FUNCTION CALL FOR THE TRAJECTORY PLOT!!!
  output$dotplot <- renderPlot({
    
    req(rv_anno())
    
    data <- rv_anno()
    
    if (input$background_type == "Visualize known genes") {
      
      cat("Making known-gene dot plot \n")
      
      known_gene_table <- calculate_de_genes()
      
      shiny::validate(
        shiny::need(
          !is.null(known_gene_table) &&
            is.data.frame(known_gene_table) &&
            nrow(known_gene_table) > 0,
          "None of the submitted genes are available in this data set."
        )
      )
      
      known_genes <- as.character(
        known_gene_table[["gene"]]
      )
      
      known_genes <- known_genes[
        !is.na(known_genes) &
          nzchar(known_genes)
      ]
      
      shiny::validate(
        shiny::need(
          length(known_genes) > 0,
          "None of the submitted genes are available in this data set."
        )
      )
      
      generate_known_gene_dot_plot(
        data,
        rv_sunburst$selected_nodes$foreground,
        known_genes
      )
      
    } else {
      
      req(get_genes_dotplot())
      
      top10_genes <- get_genes_dotplot()
      
      if (input$background_type == "Trajectory analysis") {
        
        cat("Making trajectory plot \n")
        
        generate_trajectory_plot(
          data,
          rv_sunburst$selected_nodes$foreground,
          top10_genes
        )
        
      } else {
        
        cat("Making dot plot \n")
        
        generate_dot_plot(
          input,
          data,
          rv_sunburst$selected_nodes$foreground,
          rv_sunburst$selected_nodes$background,
          top10_genes
        )
      }
    }
  })
  
  output$CHARGE_gene_plot <- downloadHandler(
    filename = "CHARGE_gene_plot.pdf",
    content = function(file) {
      req(get_genes_dotplot())
      req(rv_anno())
      
      top10_genes <- get_genes_dotplot()
      data <- rv_anno()
      
      if(input$background_type=="Trajectory analysis"){
        plot_save <- generate_trajectory_plot(data, rv_sunburst$selected_nodes$foreground, top10_genes)
      } else {
        plot_save <- generate_dot_plot(input, data, rv_sunburst$selected_nodes$foreground, rv_sunburst$selected_nodes$background, top10_genes)
      }
      
      plot_save <- plot_save + theme(text = element_text(size = as.numeric(input$dlf)))

      ggsave(file, 
             plot = plot_save,
             width = as.numeric(input$dlw), 
             height = as.numeric(input$dlh),
             useDingbats = FALSE)
    }
  )
  
  
  
  
  ##################################################
  #####          GENE SET ENRICHMENT           ##### 
  ##################################################
  
  output$gene_analysis_results_text <- renderUI({
    HTML("This section includes an interactive table of genes of potential interest based on the cell types and analysis type above, plots visualizing gene expression in selected cell types for the genes in the table, and buttons for downloading data and images and for performing get set enrichment. <b>We strongly encourage sorting and filtering by the various parameters to optimize gene selection,</b> although reasonable defaults are chosen. Hover over a column name to see a definition of that column. Note: the complete set of <b>filtered</b> genes (not just the shown genes) is used for gene set enrichment analysis.<br>")
  })
  
  
  # Enrichment button
  output$gene_set_enrichment_button <- renderUI(
    if(isTruthy(calculate_de_genes())) {
      
      actionButton("gene_set_enrichment","Calculate Gene Set Enrichment")
      
    }
  )
  
  # Reactive expression to store enrichment results
  enrichment_result <- eventReactive(input$gene_set_enrichment, {

    req(calculate_de_genes())
    req(rv_anno())
    cat("gene set enrichment \n")
    
    # Read in current gene list
    de_table <- calculate_de_genes()
    current_de_table <- de_table[input$de_table_rows_all,]
    genes <- as.character(current_de_table$gene)
    
    # If not human, convert to human gene symbols
    top_col <- apply(orthologs,2,function(x,y) length(intersect(x,y)),genes)
    top_col <- names(sort(-top_col))[1]
    if(top_col!="Human_Symbol"){
      genes <- orthologs[is.element(orthologs[,top_col],genes),"Human_Symbol"]
    }
    
    # Check if any genes remain
    if (length(genes) == 0) {
      showModal(modalDialog(
        title = "No Genes Found",
        "No valid genes found for gene set enrichment."
      ))
      return(NULL)
    }
    
    # Define background for a general enrichment analysis.
    data <- rv_anno()
    background_genes <- rownames(data$counts)
    if(top_col!="Human_Symbol"){
      background_genes <- orthologs[is.element(orthologs[,top_col],background_genes),"Human_Symbol"]
    }
    all_genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"),background_genes)
    
    cat(paste(length(intersect(genes,background_genes)),"genes in gene set \n"))
    
    # Perform Gene Ontology (GO) enrichment analysis
    #   NOTE: `org.Hs.eg.db` uses Entrez IDs. For simplicity, we'll assume the
    #   user's list are gene symbols that match the database.
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
      
      return(go_enrich_results)
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Analysis Error",
        paste("An error occurred during analysis:", e$message),
        footer = modalButton("OK")
      ))
      return(NULL)
    })
  })
  
  # Render the enrichment plot
  output$enrichment_plot <- renderPlot({
    
    # Get the results from the reactive expression
    results <- enrichment_result()
    save(results,file="tmp.RData")
    
    # Check if the enrichment result is a valid object with significant hits
    if (is.null(results) || !inherits(results, "enrichResult") || nrow(results) == 0) {
      # Return a blank plot with a message if no significant results are found
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No significant enrichment results found.",
                 size = 5, color = "grey50") +
        theme_void()
    } else {
      # Use dotplot from clusterProfiler to visualize the results
      dotplot(results, showCategory = 10, title = "Gene Ontology Enrichment Analysis") 
    }
  })
  
  
  output$CHARGE_enrichment_plot <- downloadHandler(
    filename = "CHARGE_enrichment_plot.pdf",
    content = function(file) {
      # Get the results from the reactive expression
      results <- enrichment_result()
      save(results,file="tmp.RData")
      
      # Check if the enrichment result is a valid object with significant hits
      if (is.null(results) || !inherits(results, "enrichResult") || nrow(results) == 0) {
        # Return a blank plot with a message if no significant results are found
        plot_save <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No significant enrichment results found.",
                   size = 5, color = "grey50") +
          theme_void()
      } else {
        # Use dotplot from clusterProfiler to visualize the results
        plot_save <- dotplot(results, showCategory = 10, title = "Gene Ontology Enrichment Analysis") 
      }
      
      plot_save <- plot_save + theme(text = element_text(size = as.numeric(input$enrichment_dlf)))
      
      ggsave(file, 
             plot = plot_save,
             width = as.numeric(input$enrichment_dlw), 
             height = as.numeric(input$enrichment_dlh),
             useDingbats = FALSE)
    }
  )
  
  

}





