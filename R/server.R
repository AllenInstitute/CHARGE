suppressPackageStartupMessages({
  # NOTE: ALL LIBRARIES LOADED FROM ui.R
})
options(stringsAsFactors = F)
options(shiny.maxRequestSize = 1 * 1024^3)  # 1GB max For uploading files

source("initialization.R")
source("sunburst.R")
source("de_genes_functions.R")
source("group_dot_plots.R")

#enableBookmarking("url")  # It was "server", but it doesn't seem to work either way

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
    
    # default values
    # defined in the default_vals list before the server() call, above.
    vals <- default_vals
    
    # URL values
    # defined in the URL
    # These supercede defaults
    if(length(session$clientData$url_search) > 0) {
      
      query <- as.list(parseQueryString(session$clientData$url_search))
      
      for(val in names(query)) {
        vals[[val]] <- query[[val]]
      }
    }
    
    init$vals <- vals
    
  })
  
  
  # Direct link based on input parsing
  # This can be used to provide a direct URL to users that can be bookmarked.
  output$url <- renderUI({
    req(input)
    url <- build_url(session, input)
    a("Direct Link", href = url)
  })
  
  
  output$show_url <- renderText("")
  observeEvent(input$bookmark_url, {
    req(input)
    url <- build_url(session, input)
    output$show_url <- renderText(url)
  })
  
 
  updateSelectInput(session, inputId = "select_category", label = "Choose a category:", choices = names(table_name)) # "Enter your own location"
  
  observeEvent(input$select_category, {
    # Choose a value from the default table, if selected
    # This is updated to be a list of lists
    if(length(input$select_category)>0) category = input$select_category
    updateSelectInput(session, inputId = "select_textbox", label = "Choose an existing data set:", choices = table_name[[category]])
    
  })
  
  
  #########################
  ## General UI Elements ##
  #########################
  
  # Database selection textbox and dropdown.
  # Users provide the network path to the dataset
  # This is in the server.R section so that the default value can be
  #   set using the init$vals reactive values based on defaults, 
  #   a drop-down menu, and URL parsing
  #
  # output$database_textbox - Textbox UI Object
  #
  # input$db - character object
  # 
  
   output$select_category <- renderUI({
     req(init$vals)
     
     id <- "select_category"
     #write("SELECT CATEGORY",stderr())
   
     initial <- NULL
     if(!is.null(init$vals[[id]]))
       initial <- init$vals[[id]]
     #write(initial,stderr())
   
     selectInput("select_category", "choose a category", choices = names(table_name), selected=initial)
     
   })
  
  output$select_textbox <- renderUI({
    req(init$vals)
    
    id <- "select_textbox"
    #write("SELECT TABLE",stderr())
    
    initial <- NULL
    if(!is.null(init$vals[[id]]))
      initial <- init$vals[[id]]
    
    #write(initial,stderr())
    selectizeInput("select_textbox", "select a table", choices = table_name, selected=initial)
    
  })
  

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
    
    text_desc = "README: Select a category and a data set from the boxes above -OR- to compare your own annotation data, choose 'Enter your own location' from the 'Select annotation category' and enter the location of your data file or upload a file yourself. Once a data set is selected, wait for the panels below to refresh."
    
    print("input$select_textbox")
    print(input$select_textbox)
    
    if (length(input$select_textbox)>0){
    
      # If a stored db exists, pull the value from init$vals
      if(length(init$vals[["select_textbox"]]) > 0){
        text_desc <- init$vals[["select_textbox"]]
      } else {
        
        if (input$select_textbox == 'Enter your own location') {
          text_desc = "User-provided data set file, created using the 'chargeTaxonomy' R function (see GitHub page for details)."
        } else if (input$select_textbox == 'Select comparison table...') {
          # Do nothing... text_desc should remain as initialized above
        } else {
          text_desc = table_info[table_info$table_name==input$select_textbox,"description"]
        }
      }
    }
    
    # Return the description text
    div(style = "font-size:14px;", strong("Dataset description"),br(),text_desc)
    
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
     file = rv_path()
     write("Reading file.", stderr())
     
     # THIS IS WHERE THE DATA GETS READ IN.  THERE SHOULD PROBABLY BE MORE CHECKS OF PROPER FORMAT.
     if(substr(file,1,2)=="s3"){
       ## READ FROM s3 bucket
       file2    = substr(file,6,10000)
       file2    = strsplit(file2,"/")[[1]]
       bucket   = file2[1]
       filename = paste(file2[2:length(file2)],collapse="/")

       write("filename:", stderr())
       write(filename, stderr())
       write("bucket:", stderr())
       write(bucket, stderr())
       
       objIn = objects()
       a = try({s3load(object = filename,bucket = bucket)})
       if(class(a)=="try-error"){
         write(paste("s3",file,"does not exist or cannot be accessed."))
         return(NULL)
       }
       objOut = objects()
       objs = setdiff(objOut,objIn)
       
     } else {
       ## READ LOCALLY... THIS MIGHT NOT WORK
       if(file.exists(file)){
         objs <- load(file)
       } else {
         write(paste("Local",file,"does not exist."))
         return(NULL)
       }
     }
     eval(parse(text=paste0("data=list(",paste(objs,collapse=","),")")))  
     names(data) <- objs
     return(data)
  }) # end rv_anno()
   
  
  # Check for valid input
  output$checkInput <- renderUI({
    req(rv_anno)
    if(is.null(rv_anno())){
      p("ENTER VALID DATA SET FILE.")
    } else {
      p(" ")
    }
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
    
    rv_anno()
    data = rv_anno()
    
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
  
  
  ######################################################################
  ##      Constellation plots, Sunburst plots, and plot selection     ##
  ######################################################################
  
  
  
  output$plot_type_selection <- renderUI({

    radioButtons(
      inputId = "plot_selection",
      label = "Choose plot selection type:",
      choices = list(
        "Sunburst" = "Sunburst",
        "Constellation" = "Constellation"
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
      data <- rv_anno()
      sunburst_hierarchy = data$hierarchy[length(data$hierarchy):1]    
      level = which(sunburst_hierarchy==input$hierarchy_level)
      if(length(level)==1)
        sunburst_hierarchy = sunburst_hierarchy[1:level]
      
      sunburstDF <- as.sunburstDF(data$cluster_info, sunburst_hierarchy,rootname="all")
      sunburstDF$key <- 1:dim(sunburstDF)[1]
      
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
      
      # Output 
      write(sunburstDF$label,"label.txt")
      
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
    p
    
  })
  
  # output$clickInfo <- renderPrint({
  #   d <- event_data("plotly_click")
  #   
  #   if (is.null(d)) {
  #     "Click on plot." 
  #     } else {
  #       label <- scan("label.txt",what="character",sep="\n")
  #       label[d$pointNumber+1]
  #     } 
  # })
  
  

  # This function sets the selected nodes
  observeEvent(event_data("plotly_click"), {
    
    req(rv_anno())
    data <- rv_anno()
    constellation <- data$constellation
    
    # Register the event
    d <- event_data("plotly_click")
    
    if(input$plot_selection=="Sunburst"){ 
      # Workaround because keys are not working properly in plotly
      label <- scan("label.txt",what="character",sep="\n")
      clicked_node_id <- label[d$pointNumber+1]
    } else {
      dat  = constellation[[input$hierarchy_level]]$x$layoutAttrs[[1]]$annotations
      xval = as.numeric(lapply(dat,function(x) x$x))
      yval = as.numeric(lapply(dat,function(x) x$y))
      kp   = which((xval==d$x)&(yval==d$y))
      clicked_node_id = as.character(lapply(dat,function(x) x$text))[kp[1]]
    }
    write(clicked_node_id,stderr())
    
    # Determine foreground or background
    which_list     <- "foreground"    # UPDATE THIS
    if(input$background_type=="Foreground vs. custom types")
      if(input$list_selection=="Comparison")
        which_list <- "background"
        
    selected_nodes <- rv_sunburst$selected_nodes[[which_list]]
    
    # Get the current set of selected nodes
    current_selected <- selected_nodes
      
    # Toggle the clicked node's ID in the filter
    if (clicked_node_id %in% current_selected) {
      # If already selected, remove it
      selected_nodes <- setdiff(current_selected, clicked_node_id)
    } else {
      # If not selected, add it
      selected_nodes <- unique(c(current_selected, clicked_node_id))
    }
    
    # Subset selected nodes to only include selections in the current hierarchy
    
    level = input$hierarchy_level
    if(length(level)==0) level = data$hierarchy[1]
    all_types = unique(data$cluster_info[,paste0(level,"_label")])
    selected_nodes <- selected_nodes[selected_nodes %in% all_types]
    
    # Determine foreground or background
    rv_sunburst$selected_nodes[[which_list]] <- selected_nodes
    
  })
  
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
  
  ## BACKGROUND FILTERS
  
  
  output$conditional_background_title <- renderUI({
    
    if(input$background_type=="Foreground vs. custom types"){
      h4("Comparison cell types:")
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
  ##################################################
  
  
  calculate_de_genes <- eventReactive(input$find_degenes, {
    
    req(rv_anno())
    
    if(length(rv_sunburst$selected_nodes$foreground)>0)
      
      if(!((length(rv_sunburst$selected_nodes$background)==0)&(input$background_type=="Foreground vs. custom types"))){
        data <- rv_anno()
        
        if(input$background_type=="Trajectory analysis"){
          
          find_trajectory_genes(data, rv_sunburst$selected_nodes$foreground)
          
        } else {
          
          find_de_genes(data, input, rv_sunburst$selected_nodes$foreground, rv_sunburst$selected_nodes$background)
          
        }
      }
  })
  

  output$de_table <- renderDataTable({
    req(calculate_de_genes())
    data_df = calculate_de_genes()
    
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
    
    datatable(data_df, filter = "top", 
              options = list(scrollX = TRUE, 
                             scrollY = TRUE, 
                             pageLength = 10, 
                             lengthMenu = list(c(5, 10, 20), c('5', '10', '20')),
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
    
  })
  
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

  output$downloadPlot <- downloadHandler(    
    
    filename = "sifter_heatmap.pdf",
    content = function(file) {
      heatmap_plot <- buildplot(pfontsize=as.numeric(input$dlf), showclick = FALSE) + theme(line=element_line(size=0.4))
      legend_plot <- build_legend_plot(pfontsize=as.numeric(input$dlf)) + theme(line=element_line(size=0.4),plot.margin = unit(c(0.1,0.35,0.1,0.35),"npc"))
      
      plot_list <- list(heatmap_plot,legend_plot)
      out_h <- as.numeric(input$dlh)
      out_w <- as.numeric(input$dlw)
      plot <- arrangeGrob(grobs = plot_list,
                          heights = c(out_h/8*7,out_h/8))
      #device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = plot,
             width=as.numeric(input$dlw), 
             height=as.numeric(input$dlh))
    }
  )
  
  
  # Expose input values
  output$show_inputs <- renderTable({
    
    invals <- reactiveValuesToList(input)
    
    df_inputs <- data.frame(variable = rep("",length(invals)),
                            type = rep("",length(invals)),
                            length = rep("",length(invals)),
                            current_value = rep("",length(invals)))
    
    for(i in 1:length(invals)) {
      
      var_name = names(invals)[i]
      df_inputs$variable[i] <- paste0("input$",var_name)
      
      df_inputs$type[i] <- typeof(invals[[i]])
      df_inputs$length[i] <- length(invals[[i]])
      
      if(length(invals[[i]]) == 1) {
        
        df_inputs$current_value[i] <- invals[[i]]
        
      } else if(is.vector(invals[[i]])) {
        
        df_inputs$current_value[i] <- paste0("c(",paste(invals[[i]],collapse=","),")")
        
      } else {
        
        df_inputs$current_value[i] <- typeof(invals[[i]])
        
      }
      
      
    }
    
    df_inputs <- df_inputs %>%
      arrange(variable)
    
    df_inputs
    
  })
  
  
  
  get_genes_dotplot <- reactive({
    req(calculate_de_genes())
    
    cat("de genes dot plot \n")
    
    de_table <- calculate_de_genes()
    
    current_de_table <- de_table[input$de_table_rows_current,]
    # print(current_de_table)
    
    top10_genes <- head(current_de_table$gene,100)
    
    top10_genes
    
  })
  
  
  # Group dot plots
  # NOTE:  THIS ALSO IS THE SAME FUNCTION CALL FOR THE TRAJECTORY PLOT!!!
  output$dotplot <- renderPlot({
    
    req(get_genes_dotplot())
    req(rv_anno())
    
    top10_genes <- get_genes_dotplot()
    data <- rv_anno()
    
    if(input$background_type=="Trajectory analysis"){
      
      cat("Making trajectory plot \n")
      generate_trajectory_plot(data, rv_sunburst$selected_nodes$foreground, top10_genes)
      
    } else {
      
      cat("Making dot plot \n")
      generate_dot_plot(input, data, rv_sunburst$selected_nodes$foreground, rv_sunburst$selected_nodes$background, top10_genes)
      
    }
    
  })
  
  
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
    
    # Show processing message
    output$processing_message <- renderText("Running enrichment analysis...")
    
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
  

}





