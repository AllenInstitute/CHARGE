
# Function to generate the dot plot
generate_dot_plot <- function(input, data, g1_ids, g2_ids, genes){
  
  ## Read in some data
  cluster_info <- data$cluster_info
  level        <- input$hierarchy_level
  hierarchy    <- data$hierarchy
  
  ## Update g2_ids depending on the requested analysis
  if(input$background_type=="Foreground vs. all other types"){
    g2_ids <- unique(cluster_info[,paste0(level,"_label")])
  }
  if(input$background_type=="Foreground vs. local types"){
    if(level==hierarchy[length(hierarchy)]){ 
      g2_ids <- unique(cluster_info[,paste0(level,"_label")])
    } else{
      level2 <- hierarchy[which(hierarchy==level)+1]
      keep_level2 <- cluster_info[,paste0(level,"_label")] %in% g1_ids
      keep_level2 <- cluster_info[keep_level2,paste0(level2,"_label")]
      g2_ids <- unique(cluster_info[cluster_info[,paste0(level2,"_label")] %in% keep_level2,paste0(level,"_label")])
    }
  }
  
  # Filter g2 to remove any overlap with g1
  g2_ids <- setdiff(g2_ids, g1_ids)
  
  ## Read in the rest of the data
  data_props   <- data$props[,c(g1_ids,g2_ids)]
  data_means   <- data$means[,c(g1_ids,g2_ids)]
  clusters     <- colnames(data_means)
  
  ## Convert matrices to long format
  df_color_long <- as.data.frame(data_means[genes,]) %>%
    tibble::rownames_to_column(var = "Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Cluster", values_to = "ColorValue")
  
  df_size_long <- as.data.frame(data_props[genes,]) %>%
    tibble::rownames_to_column(var = "Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Cluster", values_to = "SizeValue")
  
  ## Combine into a single data frame
  # Ensure Gene and Cluster columns are factors with desired order
  plot_data <- df_color_long %>%
    left_join(df_size_long, by = c("Gene", "Cluster")) %>%
    mutate(
      Gene = factor(Gene, levels = rev(genes)), # Reverse genes for plotting from bottom up
      Cluster = factor(Cluster, levels = clusters)
    )
  
  ## Create the dot plot
  g <- ggplot(plot_data, aes(x = Cluster, y = Gene)) +
    geom_point(aes(color = ColorValue, size = SizeValue)) +
    # --- Color Scale ---
    # Use a diverging color scale if ColorValue has a natural midpoint (e.g., 0 for logFC)
    # Or a sequential scale if values range from low to high (e.g., expression)
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = "Color Value") + # Or use viridis scales like scale_color_viridis_c()
    # --- Size Scale ---
    # Define the range for the dot sizes. Adjust 'range' as needed.
    scale_size_continuous(range = c(1, 10), # Smallest dot size is 1, largest is 10
                          name = "Size Value",
                          limits = c(0, 1)) + # If your proportion is 0-1, set limits
    # --- Labels and Theme ---
    labs(
      title = "",
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      axis.ticks = element_blank(),       # Remove axis ticks
      legend.position = "right"
    )
  
  g
  
}


# Function to generate the trajectory plot
generate_trajectory_plot <- function(data, g1_ids, genes){
  
  scale_data=TRUE
  
  ## Libraries
  library(ggplot2)
  library(dplyr)   # For data manipulation (e.g., bind_rows, mutate)
  library(tidyr)   # For data reshaping (e.g., pivot_longer)
  library(forcats) 
  
  ## Define variables
  means       <- data$means[genes, g1_ids]
  sds         <- data$sds[genes, g1_ids]
  count_n     <- data$count_n[g1_ids]
  num_runs    <- length(genes)
  num_cls     <- length(g1_ids)
  gene_labels <- genes
  
  ## Scale the data if requested
  if(scale_data){
    maxMeans  <- apply(means,1,max)
    means     <- means / maxMeans
    gene_labels <- paste0(genes," (max=",signif(maxMeans,3),")")
  }
  
  # Create a list to store data frames for each gene
  all_series_data <- list()
  
  for (i in 1:num_runs) {
    df_i <- data.frame(
      Cell_Type = 1:num_cls,
      Mean_Value = means[i,],
      SD_Value = sds[i,],
      Sample_Size = count_n,
      Series = gene_labels[i]
    )
    
    # Calculate Weights for each series
    # Variance of the mean = (SD_Value)^2 / Sample_Size
    df_i$Weight <- 1 / ((df_i$SD_Value^2) / df_i$Sample_Size)
    
    all_series_data[[i]] <- df_i
  }

  # Combine all series data frames into one large data frame
  data_df_long <- bind_rows(all_series_data)
  
  
  # Reorder the 'Series' factor from highest to lowest overall mean
  data_df_long <- data_df_long %>%
    group_by(Series) %>%
    mutate(Overall_Mean_For_Series = mean(Mean_Value)) %>% # Calculate mean for each series
    ungroup() %>%
    mutate(Series = fct_reorder(Series, Overall_Mean_For_Series, .desc = TRUE)) # Order factor by overall mean, highest first
  
  
  # --- 2. Plot the means with error bars for multiple series ---
  
  # Remove geom_point to remove dots
  # Map 'color' aesthetic to the 'Series' column
  g <- ggplot(data_df_long, aes(x = Cell_Type, y = Mean_Value, color = Series)) +
    geom_errorbar(aes(ymin = Mean_Value - SD_Value / sqrt(Sample_Size),
                      ymax = Mean_Value + SD_Value / sqrt(Sample_Size)),
                  width = 0.5, alpha = 0.5, linewidth = 1, 
                  # Ensure error bars inherit color from Series
                  position = position_dodge(width = 0.1) # Add a slight dodge if error bars overlap
    ) +
    # Add geom_line to connect the mean points for each series (as dots are removed)
    geom_line(aes(group = Series), linewidth = 0.8) + # Group by Series for distinct lines
    # Use geom_smooth for WLS lines, coloring by Series and weighting
    #geom_smooth(method = "lm", formula = y ~ x, aes(weight = Weight, group = Series), se = TRUE, linetype = "solid") + # WLS lines for each series
    labs(title = "Trajectories of shown genes for foreground cell types",
         x = "",
         y = "Gene Expression") +
    theme_minimal() +
    # Use a color palette suitable for discrete variables
    scale_color_viridis_d(option = "D", direction=1) + # Viridis is colorblind-friendly
    
    scale_x_continuous(
      breaks = unique(data_df_long$Cell_Type), # Ensure a tick for each 'Cell_Type'
      labels = g1_ids # Add the labels
    ) +
    # Optional: Adjust text angle if labels overlap
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    
    theme(legend.position = "right") # Ensure legend is visible
  
  g 
  
}
