# This function defines and returns genes and associated statistics for genes differentially expressed between a foreground and background set of cell types.
find_de_genes <- function(data, input, g1_ids, g2_ids, in_genes = NULL) {

  ## Define variables
  counts       <- data$counts
  count_n      <- data$count_n
  sums         <- data$sums
  hierarchy    <- data$hierarchy
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
	    level2 <- hierarchy[which(hierarchy==level)+1]
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
    if(length(use_genes)<2){
      showNotification("Warning: fewer than two genes included. Defaulting to normal differential expression calculation.", type = "warning")
      in_genes = NULL
    } else {
      missing_genes = setdiff(in_genes,use_genes)
      if(length(missing_genes)>0){
        missing_genes <- paste(missing_genes,collapse=", ")
        showNotification(paste("Warning:",missing_genes,"are not valid genes in this data set."), type = "warning")
      }
      rownames(counts) <- rownames(sums) <- rownames(props) <- rownames(means) <- all_genes
      counts <- counts[use_genes,]
      sums   <- sums[use_genes,]
      props  <- props[use_genes,]
      means  <- means[use_genes,]
    }
  } else {
    use_genes = all_genes
  }
  ### End new section
  #########################################
  
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
	
  # Calculate differential gene score based on mean and proportions from Sten Linnarsson's group
  # E[i,j] = ((f[i,j] + epsilon_1)/(f[i,j_hat] + epsilon_1))*((mu[i,j] + epsilon_2)/(mu[i,j_hat] + epsilon_2))
  #where f(i,j) is the fraction of non-zero expression values in the cluster and f(i,j_hat) is the fraction of non-zero expression values for cells not in the cluster. Similarly, mu(i,j) is the mean expression in the cluster and mu(i,j_hat) is the mean expression for cells not in the cluster. Small constants are added to prevent the enrichment score from going to infinity as the mean or non-zero fractions go to zero (we use epsilon_1 = 0.1 and epsilon_2 = 0.01). This formula captures both enrichment in terms of levels (mu) and in terms of fraction expressing cells (f). There is no natural cutoff, so we usually simply look at the top ten genes, but of course you could compute a null distribution using shuffled data and get a P value.
  epsilon_1 = 0.1
  epsilon_2 = 0.01
  propMeanScore <- log2(((g1_props + epsilon_1)/(g2_props + epsilon_1))*
                                ((g1_means + epsilon_2)/(g2_means + epsilon_2)))

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
	
	# Hard-coded filters (could be added as input later)
	meanSum       = 1
	absPropDiff   = 0.1
	
	# Define the output table
	if(is.null(in_genes)){
	  output <- output %>%
	    filter(abs(prop_diff) > absPropDiff) %>%
	    filter(gr1_mean + gr2_mean > meanSum) %>%
	    arrange(-prop_diff)
	}
	
	##############################
	## Add additional statistics for the subset of genes still included
	
	genesUse = output$gene
	

	# Calculate the ranked biserial correlation (as a metric for specificity)
	# To get a score that reflects both the purity of the groups and the direction of the sorting (e.g., A's before B's), we use the Rank Biserial Correlation Coefficient (r_b), which is a non-parametric measure of effect size for a two-group ranking. It quantifies how well a binary classification (like being in group A or B) predicts the rank order of the items.
	datIn <- cbind(means[genesUse, c(g1_ids,g2_ids)],props[genesUse, c(g1_ids,g2_ids)])
	rank_biserial_corr <- apply(datIn,1,score_from_ranks_wrapper,
	                            c(rep("A",length(g1_ids)),rep("B",length(g2_ids))))
	
	# Calculate the overlap coefficient. A formal statistical metric for quantifying the overlap between two distributions is the overlapping coefficient (OVL). This measures the area of intersection between the probability density functions of two distributions. A low OVL value indicates a high degree of separation between the groups, while a high value means they are largely indistinguishable. The value of OVL ranges from 0 (no overlap) to 1 (complete overlap). This is a general measure of separation
	# In this case we'll take the average value when running this test on means and proportions
	if(min(length(g1_ids),length(g2_ids))>1){
	  overlap_coefficient <- apply(datIn,1,overlap_coefficient_wrapper,
	                               c(rep("A",length(g1_ids)),rep("B",length(g2_ids))))
	} else {
	  overlap_coefficient <- 0
	}
	
	## Add the new statistics and reorder so they show up earlier
	output = cbind(output, rank_biserial_corr, overlap_coefficient)
	output = output[,c(1:4,7:10,5:6)]
	
	## Read gene categories (from function in separate file)
	source("read_gene_lists.r", local=TRUE)
	output$ABC_atlas___ = "Coming soon!"
	rownames(output) = NULL
	
	# Return the table
	output
	
}



# This function defines and returns genes and associated statistics for genes showing a trajectory pattern in a single ordered set of cell types.
find_trajectory_genes <- function(data, g1_ids, in_genes = NULL) {
  
  # Deal with edge case where only one cell type is selected
  if(length(g1_ids)<=1){
    showNotification("Error: At least two cell types are required to define a trajectory.", type = "warning")
    return(data.frame())
  }
  
  ## Define variables
  means   <- data$means[, g1_ids]
  sds     <- data$sds[, g1_ids]
  count_n <- data$count_n[g1_ids]
  num_runs<- dim(means)[1]
  
  #########################################
  ### NEW - FOR USER-PROVIDED GENE SETS ###
  
  if(!is.null(in_genes)){
    use_genes = intersect(rownames(means),in_genes)
    if(length(use_genes)<2){
      showNotification("Warning: fewer than two genes included. Defaulting to normal differential expression calculation.", type = "warning")
      in_genes = NULL
    } else {
      missing_genes = setdiff(in_genes,use_genes)
      if(length(missing_genes)>0){
        missing_genes <- paste(missing_genes,collapse=", ")
        showNotification(paste("Warning:",missing_genes,"are not valid genes in this data set."), type = "warning")
      }
      sds      <- sds[use_genes,]
      means    <- means[use_genes,]
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
    output <- output %>%
      filter(WLS_P_Value <= pvalCutoff) %>%
      arrange(-WLS_T_Value)
  }
  
  # Round to N significant digits
  output <- signif(output,4)
  
  ## Read gene categories (from function in separate file)
  output <- data.frame(gene=rownames(output),output)
  rownames(output) <- NULL
  source("read_gene_lists.r", local=TRUE)
  output$ABC_atlas___ = "Coming soon!"
  rownames(output) = NULL
  
  # Return the table
  output
  
}




########################################################
########################################################
########################################################
# HELPER FUNCTIONS

score_from_ranks_wrapper <- function(x,group){
  len <- length(x)/2
  mns <- c(rank(x[1:len]),rank(x[(len+1):(2*len)])+0.5)
  ord <- order(mns)
  group   <- c(group,group)[ord]
  signif(score_from_ranks(which(group=="A"),which(group=="B")),5)
}

score_from_ranks <- function(A_ranks, B_ranks) {
  # --- Core Calculations ---
  n_A <- length(A_ranks)
  n_B <- length(B_ranks)
  n <- n_A + n_B
  
  # Calculate the actual sums of ranks
  W_A <- sum(A_ranks)
  W_B <- sum(B_ranks)
  actual_diff <- W_B - W_A
  
  # --- Determine the Correct Normalization Factor ---
  
  # For a perfect A-B sort, ranks of A are 1:n_A, B are (n_A+1):n
  max_W_B_ab <- sum((n_A + 1):n)
  min_W_A_ab <- sum(1:n_A)
  max_diff_ab <- max_W_B_ab - min_W_A_ab
  
  # For a perfect B-A sort, ranks of B are 1:n_B, A are (n_B+1):n
  max_W_A_ba <- sum((n_B + 1):n)
  min_W_B_ba <- sum(1:n_B)
  max_diff_ba <- max_W_A_ba - min_W_B_ba
  
  # Use the appropriate normalization factor based on the sign of the actual difference
  if (actual_diff >= 0) {
    normalization_factor <- max_diff_ab
  } else {
    normalization_factor <- max_diff_ba
  }
  
  # --- Normalize to get the final score ---
  score <- actual_diff / normalization_factor
  
  return(score)
}


overlap_coefficient_wrapper <- function(x,group){
  len <- length(x)/2
  mns <- x[1:len]
  prp <- x[(len+1):(2*len)]
  signif(0.5*(
      calculate_overlap_coefficient(mns[group=="A"],mns[group=="B"])+
      calculate_overlap_coefficient(prp[group=="A"],prp[group=="B"])
    ),
  5)
}


calculate_overlap_coefficient <- function(x1, x2, n = 512) {
  

  # Combine the data to find the overall range
  combined_data <- c(x1, x2)
  min_val <- min(combined_data)
  max_val <- max(combined_data)
  data_range <- max_val - min_val
  
  # Deal with all the same values
  if(data_range==0) return(1)
  
  # Define 'a' and 'b' automatically based on your formula
  a <- min_val - 0.25 * data_range
  b <- max_val + 0.25 * data_range
  
  # Estimate the density for each group, forcing them to use the same x-points
  density1 <- density(x1, from = a, to = b, n = n)
  density2 <- density(x2, from = a, to = b, n = n)
  
  # The 'density' function's output vectors are already aligned,
  # so no interpolation is needed.
  f1 <- density1$y
  f2 <- density2$y
  
  # Calculate the intersection of the two density curves
  intersection_area <- sum(pmin(f1, f2)) * (density1$x[2] - density1$x[1])
  
  return(intersection_area)
}


