# This function defines and returns genes and associated statistics for genes differentially expressed between a foreground and background set of cell types.
find_de_genes <- function(data, input, g1_ids, g2_ids, in_genes = NULL, filter=TRUE) {

  ## Define variables
  counts       <- data$counts
  count_n      <- data$count_n
  sums         <- data$sums
  hierarchy    <- data$hierarchy
  local_context_level <- input$local_context_level
  cluster_info <- data$cluster_info
  level        <- input$hierarchy_level
  means        <- data$means
  props        <- data$props
	
	## Update g2_ids depending on the requested analysis
	if(input$background_type=="Foreground vs. all other types"){
	  g2_ids <- unique(cluster_info[,paste0(level,"_label")])
	}
	if(input$background_type=="Foreground vs. local types"){
	  if(level==hierarchy[length(hierarchy)]){ 
	    g2_ids <- unique(cluster_info[,paste0(level,"_label")])
	  } else{
	    level2 <- local_context_level # hierarchy[which(hierarchy==level)+1]
	    keep_level2 <- cluster_info[,paste0(level,"_label")] %in% g1_ids
	    keep_level2 <- cluster_info[keep_level2,paste0(level2,"_label")]
	    g2_ids <- unique(cluster_info[cluster_info[,paste0(level2,"_label")] %in% keep_level2,paste0(level,"_label")])
	    
	    # Deal with edge case where all cell types in a given background are selected (and therefore there are 0 background types)
	    if(length(g2_ids)==length(g1_ids)){
	      if(length(hierarchy)>=(which(hierarchy==level)+2)){
	        showNotification("Warning: no background types one level above. Setting background as two levels above.", type = "warning")
	        level2 <- hierarchy[which(hierarchy==level)+2]
	        keep_level2 <- cluster_info[,paste0(level,"_label")] %in% g1_ids
	        keep_level2 <- cluster_info[keep_level2,paste0(level2,"_label")]
	        g2_ids <- unique(cluster_info[cluster_info[,paste0(level2,"_label")] %in% keep_level2,paste0(level,"_label")])
	      } else{
	        showNotification("Error: no background types available. Please select different options.", type = "warning")
	      }
	    }
	    
	  }
	}
	
	# Filter g2 to remove any overlap with g1
  g2_ids <- setdiff(g2_ids, g1_ids)
  
  ## error checking
  g1_ids <- intersect(
    g1_ids,
    colnames(counts)
  )
  
  g2_ids <- intersect(
    g2_ids,
    colnames(counts)
  )
  
  if (length(g1_ids) == 0) {
    showNotification(
      "Error: no valid foreground cell types were selected.",
      type = "error"
    )
    return(data.frame())
  }
  
  if (length(g2_ids) == 0) {
    showNotification(
      "Error: no comparison cell types are available for this selection.",
      type = "error"
    )
    return(data.frame())
  }
  
  write(paste("Number of g2_ids:",length(g2_ids)),stderr())

  # Total number of cells per group
  g1_n <- sum(count_n[names(count_n) %in% g1_ids])
  g2_n <- sum(count_n[names(count_n) %in% g2_ids])

  # Gene names
  all_genes <- rownames(counts)
  
  #########################################
  ### NEW - FOR USER-PROVIDED GENE SETS ###
    
  if(!is.null(in_genes)){
    use_genes = intersect(all_genes,in_genes)
    if (length(use_genes) < 1) {
      showNotification(
        "Warning: none of the supplied genes are available in this data set.",
        type = "warning"
      )
      return(data.frame())
    } else {
      missing_genes = setdiff(in_genes,use_genes)
      if(length(missing_genes)>0){
        missing_genes <- paste(missing_genes,collapse=", ")
        showNotification(paste("Warning:",missing_genes,"are not valid genes in this data set."), type = "warning")
      }
      rownames(counts) <- rownames(sums) <- rownames(props) <- rownames(means) <- all_genes
      counts <- counts[use_genes, , drop = FALSE]
      sums   <- sums[use_genes, , drop = FALSE]
      props  <- props[use_genes, , drop = FALSE]
      means  <- means[use_genes, , drop = FALSE]
    }
  } else {
    use_genes = all_genes
  }
  ### End new section
  #########################################
  
  # Subset of count matrix  per group
  g1_data <- counts[, g1_ids, drop = FALSE]
  g2_data <- counts[, g2_ids, drop = FALSE]
  
  # Proportions of cells in each group expressing each gene
  g1_counts <- rowSums(
    g1_data,
    na.rm = TRUE
  )
  g2_counts <- rowSums(
    g2_data,
    na.rm = TRUE
  )
  g1_props  <- g1_counts/g1_n
  g2_props  <- g2_counts/g2_n

	# Calculate the log-normalized group sums and means
  g1_sums1 <- rowSums(
    sums[, g1_ids, drop = FALSE],
    na.rm = TRUE
  )
  g2_sums1 <- rowSums(
    sums[, g2_ids, drop = FALSE],
    na.rm = TRUE
  )
	g1_means <- log2(g1_sums1/g1_n+1)
  g2_means <- log2(g2_sums1/g2_n+1)
	
  # Calculate differential gene score based on mean and proportions from Sten Linnarsson's group
  # E[i,j] = ((f[i,j] + epsilon_1)/(f[i,j_hat] + epsilon_1))*((mu[i,j] + epsilon_2)/(mu[i,j_hat] + epsilon_2))
  #where f(i,j) is the fraction of non-zero expression values in the cluster and f(i,j_hat) is the fraction of non-zero expression values for cells not in the cluster. Similarly, mu(i,j) is the mean expression in the cluster and mu(i,j_hat) is the mean expression for cells not in the cluster. Small constants are added to prevent the enrichment score from going to infinity as the mean or non-zero fractions go to zero (we use epsilon_1 = 0.1 and epsilon_2 = 0.01). This formula captures both enrichment in terms of levels (mu) and in terms of fraction expressing cells (f). There is no natural cutoff, so we usually simply look at the top ten genes, but of course you could compute a null distribution using shuffled data and get a P value.
  epsilon_1 = 0.1
  epsilon_2 = 0.01
  propMeanScore <- log2(((g1_props + epsilon_1)/(g2_props + epsilon_1))*
                                ((g1_means + epsilon_2)/(g2_means + epsilon_2)))

  stopifnot(
    length(g1_props) == length(use_genes),
    length(g2_props) == length(use_genes),
    length(g1_means) == length(use_genes),
    length(g2_means) == length(use_genes)
  )
  
	# Choose top DEX genes based on difference in proportion
	output <- data.frame(gene = use_genes, 
	                     prop_diff     = round(g1_props - g2_props,5), 
	                     log2_FC       = round(g1_means - g2_means,3),
	                     propMeanScore = round(propMeanScore,3),
	                     gr1_prop = round(g1_props,3),
	                     gr1_mean = round(g1_means,3),
	                     gr2_prop = round(g2_props,3), 
	                     gr2_mean = round(g2_means,3),
	                     stringsAsFactors = F)
	output$consensus_score = round(output$prop_diff * output$log2_FC * output$propMeanScore,5)
	
	# Hard-coded filters (could be added as input later)
	meanSum       = 1
	absPropDiff   = 0.1
	
	# Define the output table
	if(is.null(in_genes)){
	  if(filter){
	    output <- output %>%
	      filter(abs(prop_diff) > absPropDiff) %>%
	      filter(gr1_mean + gr2_mean > meanSum) %>%
	      arrange(-consensus_score)
	  } else {
	    output <- output %>% arrange(-consensus_score)
	  }
	  
	}
	
	##############################
	## Add additional statistics for the subset of genes still included
	
	genesUse = output$gene
	

	# Calculate the ranked biserial correlation (as a metric for specificity)
	# To get a score that reflects both the purity of the groups and the direction of the sorting (e.g., A's before B's), we use the Rank Biserial Correlation Coefficient (r_b), which is a non-parametric measure of effect size for a two-group ranking. It quantifies how well a binary classification (like being in group A or B) predicts the rank order of the items.
	mean_data <- means[
	  genesUse,
	  c(g1_ids, g2_ids),
	  drop = FALSE
	]
	
	prop_data <- props[
	  genesUse,
	  c(g1_ids, g2_ids),
	  drop = FALSE
	]
	
	calculate_row_rbc <- function(data_matrix) {
	  
	  vapply(
	    seq_len(nrow(data_matrix)),
	    FUN = function(i) {
	      
	      group_a <- as.numeric(
	        data_matrix[i, seq_along(g1_ids), drop = TRUE]
	      )
	      
	      group_b <- as.numeric(
	        data_matrix[
	          i,
	          length(g1_ids) + seq_along(g2_ids),
	          drop = TRUE
	        ]
	      )
	      
	      group_a <- group_a[is.finite(group_a)]
	      group_b <- group_b[is.finite(group_b)]
	      
	      if (length(group_a) == 0 || length(group_b) == 0) {
	        return(NA_real_)
	      }
	      
	      mean(
	        sign(
	          outer(group_a, group_b, FUN = "-")
	        )
	      )
	    },
	    FUN.VALUE = numeric(1)
	  )
	}
	
	mean_rank_biserial_corr <- calculate_row_rbc(mean_data)
	prop_rank_biserial_corr <- calculate_row_rbc(prop_data)
	
	rank_biserial_corr <- rowMeans(
	  cbind(
	    mean_rank_biserial_corr,
	    prop_rank_biserial_corr
	  ),
	  na.rm = TRUE
	)
	
	rank_biserial_corr[
	  !is.finite(rank_biserial_corr)
	] <- NA_real_
	
	rank_biserial_corr <- signif(
	  rank_biserial_corr,
	  5
	)
	
	datIn <- cbind(mean_data, prop_data)
	
	# Calculate the overlap coefficient. A formal statistical metric for quantifying the overlap between two distributions is the overlapping coefficient (OVL). This measures the area of intersection between the probability density functions of two distributions. A low OVL value indicates a high degree of separation between the groups, while a high value means they are largely indistinguishable. The value of OVL ranges from 0 (no overlap) to 1 (complete overlap). This is a general measure of separation
	# In this case we'll take the average value when running this test on means and proportions
	if(min(length(g1_ids),length(g2_ids))>1){
	  overlap_coefficient <- apply(datIn,1,overlap_coefficient_wrapper,
	                               c(rep("A",length(g1_ids)),rep("B",length(g2_ids))))
	} else {
	  overlap_coefficient <- rep(
	    NA_real_,
	    length(genesUse)
	  )
	}
	
	## Add the new statistics and reorder so they show up earlier
	output = cbind(output, rank_biserial_corr, overlap_coefficient)
	#output = output[,c(1:4,7:10,5:6)]
	
	## Read gene categories (from function in separate file)
	source("read_gene_lists.r", local=TRUE)
	rownames(output) = NULL
	
	# Return the table
	output
	
}



# This function defines and returns genes and associated statistics for genes showing a trajectory pattern in a single ordered set of cell types.
find_trajectory_genes <- function(data, g1_ids, in_genes = NULL, filter=TRUE) {
  
  # Deal with edge case where only one cell type is selected
  if(length(g1_ids)<=1){
    showNotification("Error: At least two cell types are required to define a trajectory.", type = "warning")
    return(data.frame())
  }
  
  ## Define variables
  means   <- data$means[, g1_ids, drop = FALSE]
  sds     <- data$sds[, g1_ids, drop = FALSE]
  count_n <- data$count_n[g1_ids]
  num_runs<- dim(means)[1]
  
  #########################################
  ### NEW - FOR USER-PROVIDED GENE SETS ###
  
  if(!is.null(in_genes)){
    use_genes = intersect(rownames(means),in_genes)
    if (length(use_genes) < 1) {
      showNotification(
        "Warning: none of the supplied genes are available in this data set.",
        type = "warning"
      )
      return(data.frame())
    } else {
      missing_genes = setdiff(in_genes,use_genes)
      if(length(missing_genes)>0){
        missing_genes <- paste(missing_genes,collapse=", ")
        showNotification(paste("Warning:",missing_genes,"are not valid genes in this data set."), type = "warning")
      }
      sds      <- sds[use_genes, , drop = FALSE]
      means    <- means[use_genes, , drop = FALSE]
      num_runs <- dim(means)[1]
    }
  } 
  
  ### End new section
  #########################################
  
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
  
  # Add means
  output$mean.expression = rowMeans(means)
  
  # Hard-coded filters (could be added as input later)
  
  pvalCutoff = max(0.1,sort(output$WLS_P_Value)[100])
  
  # Define the output table
  if(is.null(in_genes)){
    if(filter){
      output <- output %>%
        filter(WLS_P_Value <= pvalCutoff) %>%
        arrange(-WLS_T_Value)
    } else {
      output <- output %>% arrange(-WLS_T_Value)
    }
  }
  
  # Round to N significant digits
  output <- signif(output,4)
  
  ## Read gene categories (from function in separate file)
  output <- data.frame(gene=rownames(output),output)
  rownames(output) <- NULL
  source("read_gene_lists.r", local=TRUE)
  rownames(output) = NULL
  
  # Return the table
  output
  
}


# This function returns the gene name, mean expression, and gene sets for input genes
create_known_gene_table <- function(data, g1_ids, in_genes = NULL) {
  
  # Deal with edge case where only one cell type is selected
  if(length(g1_ids)<1){
    showNotification("Error: At least one cell type is required to display genes.", type = "warning")
    return(data.frame())
  }
  
  ## Define variables
  means   <- data$means[, g1_ids, drop = FALSE]
  
  use_genes = intersect(rownames(means),in_genes)
  
  if (length(use_genes) == 0) {
    return(
      data.frame(
        gene = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  if(!is.null(in_genes)){
    if(length(use_genes)<1){
      showNotification("Warning: fewer than one gene included.", type = "warning")
      in_genes = NULL
    } else {
      missing_genes = setdiff(in_genes,use_genes)
      if(length(missing_genes)>0){
        missing_genes <- paste(missing_genes,collapse=", ")
        showNotification(paste("Warning:",missing_genes,"are not valid genes in this data set."), type = "warning")
      }
      means    <- means[use_genes, , drop = FALSE]
    }
  } 
  
  # create data frame
  output = data.frame(gene=rownames(means), mean.expression = round(rowMeans(means),3), stringsAsFactors = FALSE)
  rownames(output) <- NULL
  
  ## Read gene categories (from function in separate file)
  source("read_gene_lists.r", local=TRUE)
  rownames(output) = NULL
  
  # Return the table
  output
  
}



########################################################
########################################################
########################################################
# HELPER FUNCTIONS


overlap_coefficient_wrapper <- function(x, group) {
  
  n_clusters <- length(group)
  
  mean_values <- x[
    seq_len(n_clusters)
  ]
  
  prop_values <- x[
    n_clusters + seq_len(n_clusters)
  ]
  
  mean_overlap <- calculate_overlap_coefficient(
    mean_values[group == "A"],
    mean_values[group == "B"]
  )
  
  prop_overlap <- calculate_overlap_coefficient(
    prop_values[group == "A"],
    prop_values[group == "B"]
  )
  
  valid_overlaps <- c(
    mean_overlap,
    prop_overlap
  )
  
  valid_overlaps <- valid_overlaps[
    is.finite(valid_overlaps)
  ]
  
  if (length(valid_overlaps) == 0) {
    return(NA_real_)
  }
  
  signif(
    mean(valid_overlaps),
    5
  )
}



calculate_overlap_coefficient <- function(x1, x2, n = 512) {
  
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  if (length(x1) < 2 || length(x2) < 2) {
    return(NA_real_)
  }
  
  combined_data <- c(x1, x2)
  
  min_val <- min(combined_data)
  max_val <- max(combined_data)
  data_range <- max_val - min_val
  
  # Both distributions are identical constants
  if (data_range == 0) {
    return(1)
  }
  
  safe_bandwidth <- function(x, fallback_range) {
    
    bandwidth <- suppressWarnings(
      stats::bw.nrd0(x)
    )
    
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      bandwidth <- fallback_range / 100
    }
    
    max(
      bandwidth,
      sqrt(.Machine$double.eps)
    )
  }
  
  bandwidth1 <- safe_bandwidth(x1, data_range)
  bandwidth2 <- safe_bandwidth(x2, data_range)
  
  # Include the effective tails of both kernel densities
  padding <- 4 * max(bandwidth1, bandwidth2)
  
  lower_bound <- min_val - padding
  upper_bound <- max_val + padding
  
  density1 <- stats::density(
    x1,
    bw = bandwidth1,
    from = lower_bound,
    to = upper_bound,
    n = n
  )
  
  density2 <- stats::density(
    x2,
    bw = bandwidth2,
    from = lower_bound,
    to = upper_bound,
    n = n
  )
  
  overlapping_density <- pmin(
    density1$y,
    density2$y
  )
  
  # Trapezoidal numerical integration
  interval_width <- density1$x[2] - density1$x[1]
  
  overlap <- sum(
    (
      overlapping_density[-length(overlapping_density)] +
        overlapping_density[-1]
    ) / 2
  ) * interval_width
  
  # Protect against minor numerical excursions
  max(
    0,
    min(1, overlap)
  )
}


