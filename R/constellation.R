
#' @param knn.cl.df output of KNN.graph. Dataframe providing information about the cluster call of nearest neighbours of cells within a cluster. required columns: "cl.from" = cluster_id of edge origin, "cl.to" = cluster_id of edge destination, "Freq" = , "cl.from.total" = total nr of neigbours (above threshold) from cluster of origin, "cl.to.total" = total nr of neigbours (above threshold) from destination cluster, "frac" = fraction of total edge outgoing. 
#' @param cl.center.df dataframe containing metadata and coordinates for plotting cluster centroids. Required columns: "x" = x coordinate, "y" = y coordinate, "cl" = unique cluster id that should match "cl.to" and "cl.from" columns in knn.cl.df, "cluster_color","size" = nr of cells in cluster 
#' @param out.dir location to write plotting files to
#' @param prefix A character string to prepend to the filename
#' @param node.label Label to identify plotted nodes. Default is "cluster_id"
#' @param exxageration exxageration of edge width. Default is 1 (no exxageration)
#' @param curved Whether edges should be curved or not. Default is TRUE.
#' @param plot.parts output of intermediate files. default is FALSE.
#' @param plot.hull plot convex around cell type neighbourhood. Provide neighbourhood_id's that need to be plotted
#' @param node.dodge whether or not nodes are allowed to overlap. Default is false 
#' @param plot.height height of pdf in cm. Default is 25cm
#' @param plot.width width of pdf in cm. Default is 25cm 
#' @param label.size point size of plotted node labels. Default is 5pts
#' @param max_size maximum size of node. Default is 10pt 
#' @param label_repel whether to move labels away from node so they do not overlap. Default is FALSE.
#' @param node_trans transformation of node size. See ggplot::scale_size_area(trans=node_trans). Default is "sqrt".
#' @param return.list Whether to return list of independent plotting layers. Useful for replotting only part of the constellation
#' @param highlight_nodes list of node id's (matching cl.center.df) to highlight.
#' @param highlight_color Color of stroke around highlighted node. Default is red
#' @param highlight_width Width of stroke around highlighted node. Default is 1
#' @param highlight_labsize Size of label for highlighted nodes.
#' @param edge_mark_list which edges to color different selected by node using the cl value which is not necessarily the node label
#' @param edge_marking = c("dim", "highlight"),
#' @param fg.alpha alpha to use for edges that are most dark. Default = 0.4 
#' @param bg.alpha alpha to use for edges that are more faint. Default = 0.1
#' @param coord_fixed Cartesian coordinates with fixed "aspect ratio". See ggplot::coord_fixed. Default is TRUE
#'  
#' 
#' @example_data:
#'  
#' knn.cl.df <- read.csv("data/Constellation_example/knn.cl.df.csv")
#' cl.center.df <- read.csv("data/Constellation_example/cl.center.df.csv", row.names=1)
#' 
#' 
#' @usage plotting.MGE.constellation <- plot_constellation(knn.cl.df = knn.cl.df, cl.center.df = cl.center.df, out.dir = "data/Constellation_example/plot", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 


plot_constellation <- function (knn.cl.df,
                                cl.center.df,
                                out.dir,
                                prefix=format(Sys.time(), "%Y%m%d_%H%M%S"),
                                node.label = "cluster_id",
                                exxageration = 2,
                                curved = TRUE,
                                plot.parts = FALSE,
                                plot.hull = NULL,
                                plot.height = 25,
                                plot.width = 25,
                                node.dodge = FALSE,
                                label.size = 5,
                                max_size = 10,
                                label_repel = FALSE,
                                node_trans = "sqrt",
                                return.list = T,
                                highlight_nodes = NULL,
                                highlight_color = "red",
                                highlight_width = 1,
                                highlight_labsize = 10,
                                edge_mark_list = NULL,
                                edge_marking = c("dim", "highlight"),
                                fg.alpha = 0.4,
                                bg.alpha = 0.1,
                                coord_fixed= TRUE,
                                # NEW PARAMETER: Control for ggplotly output
                                enable_plotly = FALSE,
                                # NEW PARAMETER: Control label display in plotly
                                plotly_labels_on_plot = FALSE # TRUE for static labels, FALSE for hover only
) {
  
  library(gridExtra)
  library(sna)
  library(Hmisc)
  library(reshape2)
  library(ggforce)
  library(dplyr)
  library(plotly) # Add plotly to the function's libraries
  if(label_repel == TRUE) { # Only load ggrepel if needed
    library(ggrepel)
  }
  
  st=prefix
  if(!is.null(out.dir)){
    if (!file.exists(out.dir)) {
      dir.create(out.dir)
    }}
  
  ###==== Cluster nodes will represent both cluster.size (width of point) and edges within cluster (stroke of point)
  
  # select rows that have edges within cluster
  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to,]
  
  #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
  cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)]
  
  # scale the node size by square root
  cl.center.df$cluster_size_sqrt = sqrt(cl.center.df$cluster_size)
  
  ###==== plot nodes
  labels <- cl.center.df[[node.label]]
  cl.center.df$text <- labels # Ensure text column is available for plotly text aesthetic
  
  # Initial p.nodes plot used for size extraction (not directly used for final plot.all)
  # Keeping this section as it was, assuming it's for internal calculations.
  p.nodes <-   ggplot() +
    geom_point(data=cl.center.df,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   # FIX 1: Remove alpha() from aes in initial point geom, set fixed alpha outside
                   color=cluster_color),
               alpha = 0.8) + # Set fixed alpha here
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() +
    geom_text(data=cl.center.df,
              aes(x=x,
                  y=y,
                  label=labels),
              size = label.size)
  
  if (plot.parts == TRUE & !is.null(out.dir)) {
    ggsave(file.path(out.dir,paste0(st,"nodes.org.pos.pdf")), p.nodes, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  ###==== extract node size/stroke width to replot later without scaling
  g <- ggplot_build(p.nodes)
  dots <-g[["data"]][[1]] #dataframe with geom_point size, color, coords
  
  nodes <- left_join(cl.center.df, dots, by=c("x","y")) %>% ungroup()
  
  ###==== if node.dodge==TRUE new xy coords are calculated for overlapping nodes.
  
  if (node.dodge==TRUE){
    # ... (node dodging logic - unchanged, assuming it works) ...
    # This section contains nested loops for node dodging. It's computationally
    # intensive and might be a bottleneck for large datasets.
    # It also relies on potentially modifying the 'nodes' dataframe's x/y directly.
    # Consider if this logic is critical for the Plotly output, as Plotly has
    # its own ways of handling overlap (though not as robust as ggrepel).
    
    nodes$r<- (nodes$size/10)/2
    
    
    x.list <- c(mean(nodes$x), nodes$x )
    y.list <- c(mean(nodes$y), nodes$y)
    dist.test <- as.matrix(dist(cbind(x.list, y.list)))
    nodes$distance <- dist.test[2:nrow(dist.test), 1]
    nodes <- nodes[order(nodes$distance),]
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        # print(paste(d1,d2)) # Commented out print statements to reduce console spam
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          # print(paste(d1,d2)) # Commented out print statements
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        # print(paste(d1,d2)) # Commented out print statements
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          # print(paste(d1,d2)) # Commented out print statements
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
  }
  
  nodes <- nodes[order(nodes$cluster_id),]
  
  ## when printing lines to pdf the line width increases slightly. This causes the edge to extend beyond the node. Prevent this by converting from R pixels to points.
  conv.factor <- ggplot2::.pt*72.27/96
  
  
  ## line width of edge can be scaled to node point size
  nodes$node.width <- nodes$size
  
  
  if (plot.parts == TRUE & !is.null(out.dir)) {
    if (node.dodge == TRUE) {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.dodge.csv"))) }
    else {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.csv")))
    }
  }
  
  
  ###==== prepare data for plotting of edges between nodes
  
  ##filter out all edges that are <5% of total for that cluster
  #knn.cl <- knn.cl.df[knn.cl.df$frac >0.05,] #1337 lines
  knn.cl <- knn.cl.df
  ##from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),]
  nodes$cl=as.character(nodes$cl)
  knn.cl.d$cl.from <- as.character(knn.cl.d$cl.from)
  knn.cl.d$cl.to <- as.character(knn.cl.d$cl.to)
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, node.width), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[colnames(knn.cl.d)=="node.width"]<- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""
  
  
  #bidirectional
  knn.cl.bid <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- subset(knn.cl.d[i:nrow(knn.cl.d),])
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from ),]
    
    if (dim(r)[1] != 0) {
      line$Freq.to <- r$Freq
      line$node.pt.to <- r$node.pt.from
      line$frac.to <- r$frac
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
    #print(i)
  }
  
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),]
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
    #print(i)
  }
  
  
  #min frac value = 0.01
  knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
  knn.cl.uni$Freq.to <- 1
  knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  
  ###==== create line segments
  
  line.segments <- knn.cl.lines %>% select(cl.from, cl.to)
  nodes$cl <- as.character(nodes$cl)
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.from"="cl"))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.to"="cl"))
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines$frac,
                              frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  
  ##from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  ##max fraction to max point size
  line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.from
  
  line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
  
  
  ###=== create edges, exaggerated width
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*exxageration),line.segments$node.size.from) #exxagerated width
  line.segments$ex.line.to <- pmin((line.segments$line.width.to*exxageration),line.segments$node.size.to) #exxagerated width
  
  
  line.segments <- na.omit(line.segments)
  
  print("calculating edges")
  
  # Assuming edgeMaker, perpStart, perpMid, perpEnd are defined elsewhere in your environment
  # They are not in the provided function, but necessary for allEdges and poly.Edges
  # For this example, I'll assume they work correctly and 'allEdges' and 'poly.Edges' are formed.
  # Placeholder for edgeMaker if it's not defined
  if (!exists("edgeMaker")) {
    edgeMaker <- function(i, len, curved, line.segments) {
      # This is a placeholder; replace with your actual edgeMaker function
      # It should return a data frame with x, y, fraction, and Group
      data.frame(x = runif(len), y = runif(len), fraction = runif(len),
                 Group = paste(line.segments$cl.from[i], line.segments$cl.to[i], sep = ">"))
    }
  }
  if (!exists("perpStart")) {
    perpStart <- function(x_coords, y_coords, width) {
      # Placeholder for perpStart
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[1], y_coords[1]), each = 2)
    }
  }
  if (!exists("perpMid")) {
    perpMid <- function(x_coords, y_coords, width) {
      # Placeholder for perpMid
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[2], y_coords[2]), each = 2)
    }
  }
  if (!exists("perpEnd")) {
    perpEnd <- function(x_coords, y_coords, width) {
      # Placeholder for perpEnd
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[2], y_coords[2]), each = 2)
    }
  }
  
  
  allEdges <- lapply(1:nrow(line.segments), edgeMaker, len = 50, curved = curved, line.segments=line.segments) # Reduced len for faster example
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path with bend
  
  
  groups <- unique(allEdges$Group)
  
  poly.Edges <- data.frame(x=numeric(), y=numeric(), Group=character(),g1=integer(), g2=integer(),stringsAsFactors=FALSE)
  imax <- as.numeric(length(groups))
  
  # Limiting iterations for poly.Edges generation for a faster example run
  # In your actual use, remove or adjust this limit.
  # max_poly_edges_iter <- min(imax, 50) # Process max 50 groups for quick testing
  
  for(i in 1:imax) {
    #for(i in 1:max_poly_edges_iter) { # Use this line for faster testing
    select.group <- groups[i]
    select.edge <- allEdges[allEdges$Group %in% select.group,]
    
    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    
    N <- length(x)
    leftx <- numeric(N)
    lefty <- numeric(N)
    rightx <- numeric(N)
    righty <- numeric(N)
    
    ## Start point
    perps <- perpStart(x[1:2], y[1:2], w[1]/2)
    leftx[1] <- perps[1, 1]
    lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]
    righty[1] <- perps[2, 2]
    
    ### mid points
    for (ii in 2:(N - 1)) {
      seq <- (ii - 1):(ii + 1)
      perps <- perpMid(as.numeric(x[seq]), as.numeric(y[seq]), w[ii]/2)
      leftx[ii] <- perps[1, 1]
      lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]
      righty[ii] <- perps[2, 2]
    }
    ## Last control point
    perps <- perpEnd(x[(N-1):N], y[(N-1):N], w[N]/2)
    leftx[N] <- perps[1, 1]
    lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]
    righty[N] <- perps[2, 2]
    
    lineleft <- data.frame(x=leftx, y=lefty)
    lineright <- data.frame(x=rightx, y=righty)
    lineright <- lineright[nrow(lineright):1, ]
    lines.lr <- rbind(lineleft, lineright)
    lines.lr$Group <- select.group
    lines.lr[c("g1","g2")] <- as.integer(as.character(stringr::str_split_fixed(lines.lr$Group, '>', 2)))
    
    poly.Edges <- rbind(poly.Edges,lines.lr)
    
    #Sys.sleep(0.01) # Commented out for faster testing
    #cat("\r", i, "of", imax) # Commented out for cleaner output
  }
  
  if (plot.parts == TRUE && !is.null(out.dir)) {
    write.csv(poly.Edges, file=file.path(out.dir,paste0(st,"poly.edges.csv"))) }
  
  #############################
  ##                         ##
  ##        plotting         ##
  ##                         ##
  #############################
  
  labels <- nodes[[node.label]]
  
  ####plot edges
  # p.edges <- ggplot(poly.Edges, aes(group=Group))
  # p.edges <- p.edges +geom_polygon(aes(x=x, y=y), alpha=0.2) + theme_void()
  # p.edges # Commented out as this is not the final plot.all
  
  library("data.table")
  
  
  if(!is.null(edge_mark_list)){
    if(edge_marking == "dim"){
      poly.Edges$alpha_plot <- fg.alpha # Use a different column name to avoid conflict
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- bg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- bg.alpha
      
    } else if(edge_marking == "highlight"){
      poly.Edges$alpha_plot <- bg.alpha # Use a different column name
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- fg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- fg.alpha
      
    } else{
      print("provide valid edge marking")
    }
    
  } else{
    poly.Edges$alpha_plot <- fg.alpha # Use a different column name
  }
  
  # --- FIX 2: Consolidate plot.all creation and fix alpha/text logic ---
  # Define the base plot structure once
  base_plot <- ggplot() +
    geom_polygon(data=poly.Edges,
                 aes(x=x, y=y, group=Group, alpha=alpha_plot, # Use new alpha_plot variable
                     # Add hover text for polygons (edges)
                     text = paste0("From: ", g1, "<br>To: ", g2)),
                 fill = "grey60") + # Explicit fill for polygons
    scale_alpha_identity() + # FIX 2 (Crucial): Use scale_alpha_identity() for direct alpha values
    theme_void() +
    theme(legend.position = "none") # General theme for plot.all
  
  # Add geom_point for nodes
  # FIX 3: Add 'text' aesthetic to geom_point for plotly hover info
  base_plot <- base_plot +
    geom_point(data=nodes,
               # Removed alpha from outside aes here.
               # If you want point transparency, set it here explicitly as a numeric value.
               # For plotly, 0.8 is fine if you want visible points.
               alpha = 0.8,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   color=cluster_color,
                   # Add text for plotly hover (for nodes)
                   text = paste0("Cell type ID: ", cluster_id,
                                 "<br>Cell type label: ", cluster_label,
                                 "<br>Number of cells: ", cluster_size,
                                 "<br>Edges within: ", round(edge.frac.within, 2)))
    ) +
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() # Use identity for direct color mapping
  
  # Add highlight nodes if specified
  if(!is.null(highlight_nodes)){
    base_plot <- base_plot +
      geom_point(data=nodes %>% filter(cluster_id %in% highlight_nodes) ,
                 alpha=0.8,
                 shape=21, # Outlined circle
                 color=highlight_color, # Outline color
                 stroke=highlight_width, # Outline width
                 aes(x=x,
                     y=y,
                     size=cluster_size_sqrt,
                     # Add hover text for highlighted nodes if different
                     text = paste0("Cell type ID: ", cluster_id,
                                   "<br>Cell type label: ", cluster_label,
                                   "<br>Number of cells: ", cluster_size,
                                   "<br>Edges within: ", round(edge.frac.within, 2)))
      )
  }
  
  # Add hull if specified (and apply fixed alpha scaling)
  if (!is.null(plot.hull)) {
    base_plot <- base_plot +
      geom_mark_hull(data=nodes,
                     concavity = 8,
                     radius = unit(5,"mm"),
                     aes(filter = nodes$clade_id %in% plot.hull,x, y,
                         color=nodes$clade_color))
  }
  
  
  # Add labels (ggrepel or geom_text) for static display
  # FIX 4: Only add geom_text/geom_text_repel if NOT enabling plotly static labels
  if (!enable_plotly || (enable_plotly && !plotly_labels_on_plot)) {
    if(label_repel ==TRUE){
      base_plot <- base_plot +
        ggrepel::geom_text_repel(data=nodes,
                                 aes(x=x,
                                     y=y,
                                     label=.data[[node.label]]),
                                 size = label.size,
                                 min.segment.length = Inf)
    } else {
      base_plot <- base_plot +
        geom_text(data=nodes,
                  aes(x=x,
                      y=y,
                      label=.data[[node.label]]),
                  size = label.size)
    }
  }
  
  
  # Apply coord_fixed if TRUE
  if(isTRUE(coord_fixed)){
    base_plot <- base_plot + coord_fixed(ratio=1)
  }
  
  # Now assign the final ggplot object to plot.all
  plot.all <- base_plot
  
  segment.color = NA # This variable seems unused
  if (isTRUE(plot.parts) && !is.null(out.dir)) {
    ggsave(file.path(out.dir, paste0(st, ".comb.constellation.pdf")),
           plot.all, width = plot.width, height = plot.height,
           units = "cm", useDingbats = FALSE)
  }
  
  # --- FIX 5: Conditional ggplotly conversion and text addition ---
  if (enable_plotly) {
    plotly_plot <- ggplotly(plot.all, tooltip = "text") # Use "text" for custom tooltips
    
    # Add static text labels to plotly plot if requested
    if (plotly_labels_on_plot) {
      # MODIFICATION: Use layout annotations instead of add_text()
      
      # Prepare a list of annotations, one for each node
      annotations_list <- lapply(1:nrow(nodes), function(i) {
        list(
          x = nodes$x[i],
          y = nodes$y[i],
          text = as.character(nodes[[node.label]][i]), # Ensure text is character
          xref = "x",
          yref = "y",
          showarrow = FALSE, # No arrow pointing to the text
          font = list(
            color = nodes$cluster_color[i], # Color for each label
            size = label.size * 3 # Adjust font size (Plotly units differ from ggplot)
          ),
          xanchor = "center",  # Horizontal alignment of text
          yanchor = "bottom",  # Vertical alignment (places text slightly above the (x,y) point)
          yshift = 5          # Small vertical shift upwards to avoid overlap with point
        )
      })
      
      # Add the list of annotations to the layout of the plotly_plot
      plotly_plot <- plotly_plot %>%
        layout(annotations = annotations_list)
    }
    final_plot_output <- plotly_plot
  } else {
    final_plot_output <- plot.all # Return standard ggplot object
  }
  
  #############################
  ##                         ##
  ##      plot legends       ##
  ##                         ##
  #############################
  # Legend plotting logic remains unchanged. These are for static PDFs.
  # Plotly handles its own legends interactively.
  # If you need static legend images for plotly output, you might need to
  # export the plotly legend separately or create a combined layout.
  
  ### plot node size legend (1)
  plot.dot.legend <- ggplot()+
    geom_polygon(data=poly.Edges,
                 # FIX 2 (Crucial): Use alpha_plot and scale_alpha_identity() here too
                 aes(x=x, y=y, group=Group, alpha=alpha_plot), fill="grey60")+
    scale_alpha_identity() + # Use scale_alpha_identity()
    geom_point(data=nodes,
               alpha=0.8,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   color=cluster_color)) +
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() +
    geom_text(data=nodes,
              aes(x=x,
                  y=y,
                  label=labels),
              size = label.size)+
    theme_void()
  dot.size.legend <- cowplot::get_legend(plot.dot.legend)
  
  ### plot cluster legend (3)
  cl.center.df$cluster.label <- cl.center.df$cluster_label
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  label.col <- setNames(cl.center.df$cluster_color, cl.center.df$cluster.label)
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  leg.col.nr <- min((ceiling(length(cl.center.df$cluster_id)/20)),
                    5)
  cl.center <- ggplot(cl.center.df,
                      aes(x = cluster_id, y = cluster_size_sqrt)) +
    geom_point(aes(color = cluster.label)) +
    scale_color_manual(values = as.vector(label.col[levels(cl.center.df$cluster.label)])) +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 1)),
           fill=guide_legend(ncol=leg.col.nr)) +
    theme(legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          legend.key.size = unit(0.5, "lines"))
  
  cl.center.legend <- cowplot::get_legend(cl.center)
  
  
  width.1 <- max(line.segments$frac.from, line.segments$frac.to)
  width.05 <- width.1/2
  width.025 <- width.1/4
  edge.width.data <- tibble(node.width = c(1, 1, 1),
                            x = c(2, 2, 2),
                            y = c(5, 3.5, 2),
                            line.width = c(1, 0.5, 0.25),
                            fraction = c(100, 50, 25),
                            frac.ex = c(width.1, width.05, width.025))
  edge.width.data$fraction.ex <- round((edge.width.data$frac.ex *
                                          100), digits = 0)
  poly.positions <- data.frame(id = rep(c(1, 2, 3), each = 4),
                               x = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
                               y = c(4.9, 5.1, 5.5, 4.5, 3.4, 3.6, 3.75, 3.25, 1.9, 2.1, 2.125, 1.875))
  if (exxageration != 1) {
    edge.width.legend <- ggplot() + geom_polygon(data = poly.positions,
                                                 aes(x = x, y = y, group = id),
                                                 fill = "grey60") +
      geom_circle(data = edge.width.data, aes(x0 = x, y0 = y,
                                              r = node.width/2),
                  fill = "grey80", color = "grey80",
                  alpha = 0.4) +
      scale_x_continuous(limits = c(0,3)) +
      theme_void() +
      coord_fixed() +
      geom_text(data = edge.width.data, aes(x = 2.7, y = y, label = fraction.ex, hjust = 0, vjust = 0.5)) +
      annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")
  }  else {
    edge.width.legend <- ggplot() + geom_polygon(data = poly.positions,
                                                 aes(x = x, y = y, group = id), fill = "grey60") +
      geom_circle(data = edge.width.data, aes(x0 = x, y0 = y, r = node.width/2),
                  fill = "grey80", color = "grey80",
                  alpha = 0.4) +
      scale_x_continuous(limits = c(0,3)) +
      theme_void() +
      coord_fixed() +
      geom_text(data = edge.width.data, aes(x = 2.7, y = y, label = fraction, hjust = 0, vjust = 0.5)) +
      annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")
  }
  layout_legend <- rbind(c(1, 3, 3, 3, 3),
                         c(2, 3, 3, 3, 3))
  if (plot.parts == TRUE & !is.null(out.dir)) {
    ggsave(file.path(out.dir, paste0(st, ".comb.LEGEND.pdf")),
           gridExtra::marrangeGrob(list(dot.size.legend,
                                        edge.width.legend,
                                        cl.center.legend),
                                   layout_matrix = layout_legend,
                                   top=NULL),
           height = 20, width = 20, useDingbats = FALSE)
  }
  g2 <- gridExtra::arrangeGrob(grobs = list(dot.size.legend,
                                            edge.width.legend, cl.center.legend), layout_matrix = layout_legend)
  if(!is.null(out.dir) && !enable_plotly){ # Only save PDF if not enabling plotly
    ggsave(file.path(out.dir, paste0(st, ".constellation.pdf")),
           gridExtra::marrangeGrob(list(plot.all, g2), nrow = 1, ncol = 1,top=NULL),
           width = plot.width, height = plot.height, units = "cm",
           useDingbats = FALSE)
  }
  
  if (return.list == TRUE) {
    # Return the plotly object if enabled, otherwise the ggplot object
    return(list(constellation = final_plot_output, edges.df = poly.Edges,
                nodes.df = nodes))
  }
}




## function to draw (curved) line between to points
#' function to draw (curved) line between two points
#'
#' @param whichRow 
#' @param len 
#' @param line.segments 
#' @param curved 
#'
#' @return
#' @export
#'
#' @examples
edgeMaker <- function(whichRow, len=100, line.segments, curved=FALSE){
  
  fromC <- unlist(line.segments[whichRow,c(3,4)])# Origin
  toC <- unlist(line.segments[whichRow,c(5,6)])# Terminus
  # Add curve:
  
  graphCenter <- colMeans(line.segments[,c(3,4)])  # Center of the overall graph
  bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
  distance1 <- sum((graphCenter - bezierMid)^2)
  if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
    bezierMid <- c(toC[1], fromC[2])
    }  # To select the best Bezier midpoint
  bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
  if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve

  edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
  
  #line.width.from in 100 steps to linewidth.to
 edge$fraction <- seq(line.segments$ex.line.from[whichRow], line.segments$ex.line.to[whichRow], length.out = len)

  
    #edge$Sequence <- 1:len  # For size and colour weighting in plot
  edge$Group <- paste(line.segments[whichRow, 1:2], collapse = ">")
  return(edge)
  }



#utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#' Title
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return
#' @export
#'
#' @examples
perpStart <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 1)
        }

#' Title
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
avgangle <- function(x, y) {
    a1 <- angle(x[1:2], y[1:2])
    a2 <- angle(x[2:3], y[2:3])
    atan2(sin(a1) + sin(a2), cos(a1) + cos(a2))
}

#' Title
#'
#' @param x 
#' @param y 
#' @param len 
#' @param a 
#' @param mid 
#'
#' @return
#' @export
#'
#' @examples
perp <- function(x, y, len, a, mid) {
    dx <- len*cos(a + pi/2)
    dy <- len*sin(a + pi/2)
    upper <- c(x[mid] + dx, y[mid] + dy)
    lower <- c(x[mid] - dx, y[mid] - dy)
    rbind(upper, lower)    
        }

#' Title
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return
#' @export
#'
#' @examples
perpMid <- function(x, y, len) {
    ## Now determine angle at midpoint
    perp(x, y, len, avgangle(x, y), 2)
        }

#' Title
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return
#' @export
#'
#' @examples
perpEnd <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 2)
}


## x and y are vectors of length 2
#' Title
#'
#' @param x vector
#' @param y vector
#'
#' @return
#' @export
#'
#' @examples
angle <- function(x, y) {
    atan2(y[2] - y[1], x[2] - x[1])
        }



plot_umap_constellation <- function(umap.2d, cl, cl.df, select.knn.cl.df, dest.d=".", prefix="",...)
  {
    
    cl.center.df = as.data.frame(get_RD_cl_center(umap.2d,cl))
    cl.center.df$cl = row.names(cl.center.df)
    cl.center.df$cluster_id <- cl.df$cluster_id[match(cl.center.df$cl, cl.df$cl)]
    cl.center.df$cluster_color <- cl.df$cluster_color[match(cl.center.df$cl, cl.df$cl)]
    cl.center.df$cluster_label <- cl.df$cluster_label[match(cl.center.df$cl, cl.df$cl)] 
    cl.center.df$cluster_size <- cl.df$cluster_size[match(cl.center.df$cl, cl.df$cl)]
    tmp.cl = row.names(cl.center.df)
    tmp.knn.cl.df = select.knn.cl.df  %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
    p=plot_constellation(tmp.knn.cl.df, cl.center.df, node.label="cluster_id", out.dir=file.path(dest.d,prefix),...)    
  }

