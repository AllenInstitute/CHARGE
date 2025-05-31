# This function defines and returns genes and associated statistics for genes differentially expressed between a foreground and background set of cell types.
find_de_genes <- function(data, input, g1_ids, g2_ids) {

  ## Define variables
  counts       <- data$counts
  count_n      <- data$count_n
  sums         <- data$sums
  hierarchy    <- data$hierarchy
  cluster_info <- data$cluster_info
  level        <- input$hierarchy_level
	
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

  # Total number of cells per group
  g1_n <- sum(count_n[names(count_n) %in% g1_ids])
  g2_n <- sum(count_n[names(count_n) %in% g2_ids])
  
  # Subset of count matrix  per group
  g1_data <- counts[, g1_ids]
  g2_data <- counts[, g2_ids]
  
  # Proportions of cells in each group expressing each gene
  if(length(g1_ids)>1){
    g1_counts <- rowSums(g1_data)
  } else {
    g1_counts = g1_data
  }
  if(length(g2_ids)>1){
    g2_counts <- rowSums(g2_data)
  } else {
    g2_counts = g2_data
  }
  g1_props  <- g1_counts/g1_n
  g2_props  <- g2_counts/g2_n

  # Gene names
  all_genes <- rownames(counts)

	# Calculate the log-normalized group sums and means
  if(length(g1_ids)>1){
    g1_sums1 <- rowSums(sums[, g1_ids])
  } else {
    g1_sums1 <- sums[, g1_ids]
  }
  if(length(g2_ids)>1){
    g2_sums1 <- rowSums(sums[, g2_ids])
  } else {
    g2_sums1 <- sums[, g2_ids]
  }
	g1_means <- log2(g1_sums1/g1_n+1)
  g2_means <- log2(g2_sums1/g2_n+1)
	
	# Calculate differnetial gene score based on mean and proportions from Sten Linnarsson's group
	# E[i,j] = ((f[i,j] + epsilon_1)/(f[i,j_hat] + epsilon_1))*((mu[i,j] + epsilon_2)/(mu[i,j_hat] + epsilon_2))
    #where f(i,j) is the fraction of non-zero expression values in the cluster and f(i,j_hat) is the fraction of non-zero expression values for cells not in the cluster. Similarly, mu(i,j) is the mean expression in the cluster and mu(i,j_hat) is the mean expression for cells not in the cluster. Small constants are added to prevent the enrichment score from going to infinity as the mean or non-zero fractions go to zero (we use epsilon_1 = 0.1 and epsilon_2 = 0.01). This formula captures both enrichment in terms of levels (mu) and in terms of fraction expressing cells (f). There is no naturall cutoff, so we usually simply look at the top ten genes, but of course you could compute a null distribution using shuffled data and get a P value.
	epsilon_1 = 0.1
	epsilon_2 = 0.01
	meanAndPropScore <- log2(((g1_props + epsilon_1)/(g2_props + epsilon_1))*
	  ((g1_means + epsilon_2)/(g2_means + epsilon_2)))

	
	# Choose top DEX genes based on difference in proportion
	output <- data.frame(gene = all_genes, 
	                     prop_diff     = round(g1_props - g2_props,5), 
	                     log2_FC       = round(g1_means - g2_means,3),
	                     propMeanScore = round(meanAndPropScore,3),
	                     gr1_prop = round(g1_props,3),
	                     gr1_mean = round(g1_means,3),
	                     gr2_prop = round(g2_props,3), 
	                     gr2_mean = round(g2_means,3),
	                     stringsAsFactors = F)
	
	# Hard-coded filters (could be added as input later)
	meanSum       = 1
	absPropDiff   = 0.2
	
	# Define the output table
	output <- output %>%
	  filter(abs(prop_diff) > absPropDiff) %>%
	  filter(gr1_mean + gr2_mean > meanSum) %>%
	  arrange(-prop_diff)
	
	## Read gene categories (from function in separate file)
	source("read_gene_lists.r", local=TRUE)
	output$ABC_atlas = "Future link-out to ABC Atlas"

	# Return the table
	output
	
}



# This function defines and returns genes and associated statistics for genes showing a trajectory pattern in a single ordered set of cell types.
find_trajectory_genes <- function(data, g1_ids) {
  
  ## Define variables
  means   <- data$means[, g1_ids]
  sds     <- data$sds[, g1_ids]
  count_n <- data$count_n[g1_ids]
  num_runs<- dim(means)[1]
  
  # Create the base data frame that's constant for all runs
  base_data_df <- data.frame(
    Day = 1:length(g1_ids),
    Sample_Size = count_n
  )
  
  # Store means and sds in lists, where each element is a vector for one run
  list_of_actual_means <- asplit(means,MARGIN=1)
  list_of_actual_sds <- asplit(sds,MARGIN=1)
  
  # Calculate Weighted Least Squares (WLS) of the above data for each gene using lapply
  wls_results_lapply <- lapply(1:num_runs, function(i) {
    # Combine with the specific Mean_Value and SD_Value for this run
    current_data_df <- base_data_df
    current_data_df$Mean_Value <- list_of_actual_means[[i]]
    current_data_df$SD_Value <- list_of_actual_sds[[i]]
    
    # Checks for issues and return values that are not significant if there are issues
    if (any(is.na(current_data_df$Mean_Value)) ||
        any(is.infinite(current_data_df$Mean_Value)) ||
        any(current_data_df$SD_Value <= 0) || # SD cannot be zero or negative for variance
        any(is.na(current_data_df$SD_Value)) ||
        any(is.infinite(current_data_df$SD_Value)) ||
        any(current_data_df$Sample_Size <= 0) || # Sample size cannot be zero or negative
        any(is.na(current_data_df$Sample_Size)) ||
        any(is.infinite(current_data_df$Sample_Size))) {
      # If any problematic data, return NA for this run
      return(list(slope = 0, t_value = 0, p_value = 1)) 
    }
    
    # Calculate Weights
    current_data_df$Weight <- 1 / ((current_data_df$SD_Value^2) / current_data_df$Sample_Size)
    
    # Perform WLS
    wls_model <- lm(Mean_Value ~ Day, data = current_data_df, weights = Weight)
    
    # Extract results
    model_summary  <- summary(wls_model)
    slope_estimate <- coef(model_summary)["Day", "Estimate"]
    slope_t_value  <- coef(model_summary)["Day", "t value"]
    slope_p_value  <- coef(model_summary)["Day", "Pr(>|t|)"]
    
    # Return as a named list for easy conversion to data frame later
    list(slope = slope_estimate, t_value = slope_t_value, p_value = slope_p_value)
  })
  
  # Convert the list of results to a data frame
  results_df_lapply <- do.call(rbind, lapply(wls_results_lapply, as.data.frame))
  colnames(results_df_lapply) <- c("WLS_Slope", "WLS_T_Value", "WLS_P_Value")
  rownames(results_df_lapply) <- rownames(means)
  
  # Add FDR values
  output <- results_df_lapply
  output$WLS_FDR <- p.adjust(output[,"WLS_P_Value"], method = "fdr")
  
  # Hard-coded filters (could be added as input later)
  pvalCutoff = 0.1
  
  # Define the output table
  output <- output %>%
    filter(WLS_P_Value < pvalCutoff) %>%
    arrange(-WLS_T_Value)
  
  ## Read gene categories (from function in separate file)
  output <- data.frame(gene=rownames(output),output)
  rownames(output) <- NULL
  source("read_gene_lists.r", local=TRUE)
  output$ABC_atlas = "Future link-out to ABC Atlas"
  
  # Return the table
  output
  
  
  
}

