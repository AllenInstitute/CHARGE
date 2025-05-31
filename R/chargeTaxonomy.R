# This function is tested on scrattch docker 1.2

#' Save an RData object for use with CHARGE app
#' 
#' Write an RData object with the following variables
#'  - counts_n - vector of total number of cells per cluster
#'  - counts - cell x cluster matrix where each value represents the number of cells in a cluster with at least one read for each gene 
#'  - props - cell x cluster matrix where each value represents the fraction of cells in a cluster with at least one read for each gene
#'  - sums - cell x cluster matrix where each value represents the total number of counts of that gene for that cluster 
#'  - means - cell x cluster matrix where each value represents the average counts of that gene for that cluster 
#'  - sds - cell x cluster matrix where each value represents the standard deviation of counts of that gene for that cluster
#'  - hierarchy - cell type hierarchy stored in uns$hierarchy, but converted to a corrected ordered character vector
#'  - cluster_info
#'  NOTE: proportions are not provided because they can be calculated trivaially from the above variables
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param check.taxonomy Should the function check if the input variable is a valid AIT taxonomy (default = TRUE)
#' @param subsample How much subsampling should be done before running statistics. Default (recommended!) is none, but subsampling can be done to speed up calculations.
#' @param charge.file.name File name (and path) to write charge file to
#' 
#' @import scrattch.taxonomy
#' @import rsvd
#' 
#' @examples
#' \dontrun{
#'   library(scrattch.taxonomy)
#'   AIT.file    <- "HMBA_BG_Macaque_AIT.h5ad" #"Macaque_HMBA_basalganglia_AIT_pre-print.h5ad"
#'   AIT.anndata <- loadTaxonomy(AIT.file)
#'   #raw = AIT.anndata$raw$X
#'   #raw@x <- round(2^(raw@x)-1)
#'   #AIT.anndata$raw$X@x <- raw@x  # This does not work
#'   chargeTaxonomy(AIT.anndata)
#' }
#'
#' @export
chargeTaxonomy <- function(AIT.anndata,
                           check.taxonomy = TRUE,
                           subsample      = 100000000,
                           charge.file.name = paste0(AIT.anndata$title,"_CHARGE.RData")
){
  #############################################################
  print("===== Checking taxonomy, if requested =====")
  if(check.taxonomy){
    check = checkTaxonomy(AIT.anndata)
    if(!check) warning("Taxonomy does not pass taxonomy checks, potentially causing chargeTaxonomy to fail. Check the log file if needed.") 
  } else {
    warning("Taxonomy is not being checked. If chargeTaxonomy fails, we recommend running checkTaxonomy to ensure valid AIT file.")
  }
  
  cell_counts <- AIT.anndata$raw$X
  if(is.null(cell_counts)){
    error("chargeTaxonomy required raw counts in the raw$X slot to run.")
  }
  
  norm_counts <- AIT.anndata$X
  if(is.null(norm_counts)){
    print("Creating normalized count matrix from cell_counts.")
    norm_counts = scrattch.taxonomy::log2CPM_byRow(cell_counts)
  }

  #############################################################
  print("===== Collecting required cluster information =====")
  # hierarchy
  hierarchy = names(AIT.anndata$uns$hierarchy)[order(-as.numeric(AIT.anndata$uns$hierarchy))]
  
  # cluster_info
  cluster_vector = AIT.anndata$obs[,hierarchy[1]]
  all_clusters   = unique(cluster_vector)
  if(is.factor(all_clusters)) 
    all_clusters = levels(all_clusters)
  if(!is.null(AIT.anndata$uns$cluster_info)){
    cluster_info = AIT.anndata$uns$cluster_info
  } else {
    cluster_info = AIT.anndata$obs
  }
  cluster_info = cluster_info[match(all_clusters,cluster_info[,hierarchy[1]]),hierarchy]
  
  # Rename anything called "_id" or "_label" to avoid breaking shiny
  hierarchy <- gsub("_id","",hierarchy)
  hierarchy <- gsub("_label","",hierarchy)
  colnames(cluster_info) <- hierarchy

  # Subsample vector
  keep_sample <- subsampleCells(cluster_vector,subsample)
  if(mean(keep_sample)<1 ){
    cell_counts    <- cell_counts[,keep_sample]
    cluster_vector <- cluster_vector[keep_sample]
  }
  
  print("===== Building statistics =====")

  # Building sums, counts, and props from the COUNT MATRIX
  print("... transpose count matrix")
  t_cell_counts  <- as(Matrix::t(cell_counts),"dgCMatrix")
  cluster_factor <- factor(cluster_vector,levels=all_clusters)
  rownames(t_cell_counts) <- colnames(AIT.anndata)
  colnames(t_cell_counts) <- names(cluster_factor) <- rownames(AIT.anndata)[keep_sample]

  count_n <- as.numeric(table(cluster_factor))
  names(count_n) <- all_clusters
  print("... sums")
  sums    <- round(get_cl_sums(t_cell_counts,cluster_factor))
  print("... counts")
  counts  <- round(get_cl_sums(t_cell_counts>0,cluster_factor))
  print("... props")
  props   <- t(t(counts)/as.numeric(count_n))
  
  # Building means and sds from LOG2CPM
  print("... transpose logCPM matrix")
  t_norm_counts  <- as(Matrix::t(norm_counts),"dgCMatrix")
  rownames(t_norm_counts) <- colnames(AIT.anndata)
  colnames(t_norm_counts) <- names(cluster_factor) <- rownames(AIT.anndata)[keep_sample]
  
  print("... means")
  means   <- get_cl_means(t_norm_counts,cluster_factor)
  print("... sds")
  sds     <- sqrt(get_cl_vars(t_norm_counts,cluster_factor,means))


  #############################################################
  print("===== Add cluster ids and colors =====")
  ## Add ids and colors to the cluster information
  ## NOTE: I'LL NEED TO UPDATE THIS PART ONCE I SORT OUT HOW TO EMBED CLUSTER COLORS
  ## (AND CLUSTER ORDER?) WITHIN THE CLUSTER_INFO DATA FRAME
  cluster_info$sample_name = paste0("i",1:dim(cluster_info)[1])
  for (i in 1:dim(cluster_info)[2]){
    cluster_info[,i] <- factor(cluster_info[,i],levels = unique(cluster_info[,i]))
  }    
  cluster_info <- auto_annotate(cluster_info)
  
  
  #############################################################
  print("===== Building statistics for the rest of the hierachy =====")
  # NOTE: THIS SUMMARIZES EVERYTHING BY CLUSTER, TREATING EACH CLUSTER WITH THE SAME WEIGHT
  # IF WE WANT TO TREAT EACH CELL WITH THE SAME WEIGHT, WE NEED TO REDO THE "get_cl_XXXX"
  # CALLS ABOVE REPLACING "cluster_factor" FOR EACH LEVEL OF THE HIERARCHY
  sums   <- addHierarchyToStat(sums,hierarchy,cluster_info,"sum")
  counts <- addHierarchyToStat(counts,hierarchy,cluster_info,"sum")
  means  <- addHierarchyToStat(means,hierarchy,cluster_info)
  sds    <- addHierarchyToStat(sds,hierarchy,cluster_info)
  props  <- addHierarchyToStat(props,hierarchy,cluster_info)
  
  count_n2 <- addHierarchyToStat(rbind(count_n,count_n),hierarchy,cluster_info,"sum")
  count_n  <- setNames(as.numeric(count_n2[1,]),colnames(count_n2))
  
 
  #############################################################
  print("===== Create input files for constellation plots =====")
  
  # Create a reduced dimensionality matrix (essentially creating principal components)
  # ---- We are starting from the normalized data here
  # ---- We are also subsampling to random 100 cells per cluster to do this more quickly
  # ---- Finally, we take the top 1000 highly variable genes
  
  ## CHOOSE A SUBSET OF CELLS FOR EVERYTHING
  keep_sample_pca <- subsampleCells(cluster_vector,100)

  ## Get the top variable genes
  if(is.null(AIT.anndata$var$highly_variable_genes_standard)){
    betaScore <- getBetaScore_fast(props[rowMaxs(props)>0.5,1:length(all_clusters)],returnScore=FALSE)
    betaScore <- sort(betaScore)
    top.genes <- names(betaScore)[1:min(1000,length(betaScore))]
  } else {
    top.genes <- AIT.anndata$var$highly_variable_genes_standard
  }
  
  
  ## GET PCS (e.g., reduced dimension data frame)
  print("... calculate principal components")
  
  max_k_to_estimate   <- min(c(nrow(norm_counts),ncol(norm_counts),200))
  pca_result          <- rpca(norm_counts[keep_sample_pca,top.genes], k = max_k_to_estimate, center = TRUE, scale = TRUE)
  exp_variance_ratio  <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(exp_variance_ratio)
  variance_threshold  <- 0.9
  n_components_chosen <- which.max(cumulative_variance >= variance_threshold)
  rd.dat              <- pca_result$x[,1:n_components_chosen]
  colnames(rd.dat)    <- paste0("PC",1:dim(rd.dat)[2])  
  
  
  ## GET THE UMAP
  print("... read in or calculate UMAP")
  embedding <- NULL
  if(length(AIT.anndata$obsm)>0){
    embedding <- AIT.anndata$uns$default_embedding[[1]]
    if(is.null(embedding)) embedding = names(AIT.anndata$obsm)[1]
  }
  if(length(embedding)==1){
    umap.df <- AIT.anndata$obsm[embedding][keep_sample_pca,]
  } else {
    umap.df <- umap(rd.dat)$layout
    umap.df <- as.data.frame(umap.df)
  }
  rownames(umap.df) <- rownames(rd.dat)
  
  
  ##  NOW BUILD THE CONSTELLATION
  constellation <- list()
  for (level in hierarchy){
    print(paste("... creating constellation for",level))
    cl.cl  = AIT.anndata$obs[keep_sample_pca,level]
    names(cl.cl) = rownames(rd.dat)
    result = get_knn_graph(rd.dat, cl=cl.cl, k =50) 
    
    ## Select robust edges for plotting
    knn.cl.df = result$knn.cl.df 
    knn.cl.df = knn.cl.df %>% group_by(cl.from) %>% mutate(cl.from.rank = rank(-Freq))
    knn.cl.df = knn.cl.df %>% group_by(cl.to) %>% mutate(cl.to.rank = rank(-Freq))
    select.knn.cl.df = with(knn.cl.df, knn.cl.df[odds > 1 & pval.log < log(1/100000) & (frac > 0.1 | frac > 0.03 & Freq > 100) & (cl.from.rank < 4| cl.to.rank < 4),])
    
    ## Reorganize the cluster_info matrix for the information required for plotConstellation
    ## --- THIS SECTION NEEDS TO BE EDITED!
    prefix = level
    
    cl.center.df = as.data.frame(get_RD_cl_center(umap.df,cl.cl)) 
    
    ## Define cl.df
    types <- rownames(cl.center.df)
    cl.df <- cluster_info[,paste0(level,c("_id","_label","_color"))]
    cl.df <- cl.df[match(unique(cl.df[,1]),cl.df[,1]),]
    rownames(cl.df) <- cl.df[,2]
    cl.df <- cl.df[types,]
    cl.df$cluster_size <- count_n[types] #as.numeric(table(AIT.anndata$obs[keep_sample_pca,hierarchy[1]])[cl.df[,2]]) #count_n[types]
    
    ## Define cl.center.df (centroids for clusters)
    cl.center.df$cluster_id    <- cl.df[,1]  # id column
    cl.center.df$cluster_label <- rownames(cl.df)
    cl.center.df$cluster_color <- cl.df[,3]  # color column
    cl.center.df$cluster_size  <- cl.df$cluster_size
    rownames(cl.center.df)     <- rownames(cl.df)
    
    # set cl as cluster_id since that was used to summarise the edges
    #cl.center.df$cl = as.integer(as.character(row.names(cl.center.df) ))
    cl.center.df$cl <- cl.df[,1]  # The _id column
    
    # Convert labels to ids 
    select.knn.cl.df$cl.from <- cl.df[,1][match(select.knn.cl.df$cl.from,cl.df[,2])]
    select.knn.cl.df$cl.to   <- cl.df[,1][match(select.knn.cl.df$cl.to,cl.df[,2])]
    
    # Define the edges (I think that is what tmp.knn.cl.df does???)
    tmp.cl = cl.center.df$cluster_id
    tmp.knn.cl.df = select.knn.cl.df %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
    
    # Create the plot
    c.plot=try(plot_constellation(tmp.knn.cl.df, 
                              cl.center.df=cl.center.df, 
                              out.dir=NULL,
                              prefix=prefix,
                              node.label="cluster_label",
                              exxageration=2,
                              plot.parts=FALSE,
                              return.list = T,
                              node.dodge = F,
                              label_repel = TRUE,
                              label.size = 3,
                              plot.height = 15,
                              plot.width = 15,
                              max_size = 5,
                              enable_plotly = TRUE,
                              plotly_labels_on_plot = TRUE))
    if(class(c.plot)[1]=="try-error") {
      constellation[[level]] <- plot_ly() %>%
        add_annotations(
          text = paste("No constellation diagram for",level),
          x = 0.5, y = 0.5,          # Center coordinates
          xref = "paper", yref = "paper", # Relative to plot area
          showarrow = FALSE,
          font = list(size = 36, color = "black") # Basic font styling
        ) %>%
        layout(
          xaxis = list(visible = FALSE), # Hide X-axis
          yaxis = list(visible = FALSE), # Hide Y-axis
          # Optional: make background transparent if embedding or don't want default gray
          plot_bgcolor = 'rgba(0,0,0,0)',
          paper_bgcolor = 'white'
        )
    } else {
      constellation[[level]] <- c.plot$constellation
    }
  } 
  
  #############################################################
  print("===== Save file =====")
  save(
    hierarchy,
    cluster_info,
    count_n,
    counts,
    sums,
    means,
    sds,
    props,
    constellation,
    file = charge.file.name
  )
  
}



### NEED TO UPDATE THE CODE BELOW TO DEAL WITH COLORS, BUT IT ALREADY IS CLOSE

refactorize_annotations <- function(anno, metadata){
  # Do nothing if cell type names are not provided
  if(sum(colnames(metadata)=="cell_type")==0)
    return(anno)
  
  # If provided, look for cell type names and order annotations accordingly by default.
  cns <- colnames(anno)
  cns <- gsub("_label","",cns[grepl("_label",cns)])
  for (cn in cns){
    intersecting_cell_types <- intersect(metadata$cell_type,as.character(anno[,paste0(cn,"_label")]))
    if(length(intersecting_cell_types)>0){
      new_levels <- c(intersecting_cell_types,setdiff(as.character(anno[,paste0(cn,"_label")]),metadata$cell_type))
      anno[,paste0(cn,"_id")]  <- as.numeric(factor(anno[,paste0(cn,"_label")],levels=new_levels))
    }
    
    ## If a column is called "color" in the metadata file, then look for categorical variable colors in the "cell_type" column
    
    if((sum(colnames(metadata)=="cell_type")==1)&(length(intersecting_cell_types)>0)){
      colors <- anno[,paste0(cn,"_color")][match(new_levels,anno[,paste0(cn,"_label")])]
      colors[match(intersecting_cell_types,new_levels)] <- metadata$color[match(intersecting_cell_types,metadata$cell_type)]
      anno[,paste0(cn,"_color")] <- colors[match(anno[,paste0(cn,"_label")],new_levels)]
      new_levels <- c(intersecting_cell_types,setdiff(as.character(anno[,paste0(cn,"_label")]),metadata$cell_type))
      anno[,paste0(cn,"_id")]  <- as.numeric(factor(anno[,paste0(cn,"_label")],levels=new_levels))
    }
    
  }
  
  return(anno)
}




#################################################################################
# UPDATED FUNCTION WITH BUG FIX FOR LARGE CSV FILES AND TO OMIT DIRECTION COLUMNS

auto_annotate <- function (anno, scale_num = "predicted", na_val_num = 0, colorset_num = c("darkblue", 
                                                                                           "white", "red"), sort_label_cat = TRUE, na_val_cat = "ZZ_Missing", 
                           colorset_cat = "varibow", color_order_cat = "sort") 
{
  # Define and properly format a sample name
  anno_out <- as.data.frame(anno)
  cn <- colnames(anno_out)
  if (!is.element("sample_name", cn)) {
    colnames(anno_out) <- gsub("sample_id", "sample_name", cn)
  }
  if (!is.element("sample_name", colnames(anno_out))) {
    anno_out <- cbind(anno_out, paste0("sn_",1:dim(anno_out)[1]))
    colnames(anno_out) <- c(cn,"sample_name")
  }
  anno_out <- anno_out[,c("sample_name",setdiff(colnames(anno_out),"sample_name"))]
  
  # Annotate any columns missing annotations
  cn <- colnames(anno_out)
  convertColumns <- cn[(!grepl("_label", cn)) & (!grepl("_id", 
                                                        cn)) & (!grepl("_color", cn)) & (!grepl("_direction", cn))]  # UPDATE TO OMIT DIRECTION COLUMNS
  convertColumns <- setdiff(convertColumns, "sample_name")
  convertColumns <- setdiff(convertColumns, gsub("_label", 
                                                 "", cn[grepl("_label", cn)]))
  
  # Return input annotation file if there is nothing to convert
  if(length(convertColumns)==0){
    anno_out  <- group_annotations(anno_out)
    return(anno_out)
  }
  
  anno_list <- list()
  for (cc in convertColumns) {
    value <- anno_out[, cc]
    if (sum(!is.na(value)) == 0) 
      value = rep("N/A", length(value))
    if (is.numeric(value)) {
      if (length(table(value)) == 1) 
        value = jitter(value, 1e-06)
      val2 <- value[!is.na(value)]
      if (is.element(scale_num, c("linear", "log10", "log2", 
                                  "zscore"))) {
        out <- annotate_num(df = anno_out[,c("sample_name",cc)], col = cc, 
                            scale = scale_num, na_val = na_val_num, colorset = colorset_num)
      }
      else {
        scalePred <- ifelse(min(val2) < 0, "linear", 
                            "log10")
        if ((max(val2 + 1)/min(val2 + 1)) < 100) {
          scalePred <- "linear"
        }
        if (mean((val2 - min(val2))/diff(range(val2))) < 
            0.01) {
          scalePred <- "log10"
        }
        out <- annotate_num(df = anno_out[,c("sample_name",cc)], col = cc, 
                            scale = scalePred, na_val = na_val_num, colorset = colorset_num)
      }
    }
    else {
      if (is.factor(value)) {
        out <- annotate_factor(df = anno_out[,c("sample_name",cc)], col = cc, 
                               base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                               color_order = color_order_cat)
      }
      else {
        out <- annotate_cat(df = anno_out[,c("sample_name",cc)], col = cc, 
                            base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                            color_order = color_order_cat, sort_label = sort_label_cat)
      }
    }
    anno_list[[cc]] <- out[,colnames(out)!="sample_name"]
  }
  
  # Format the annotations as an appropriate data frame
  for(cc in convertColumns)
    anno_list[[cc]] <- anno_list[[cc]][,colnames(anno_list[[cc]])!="sample_name"]
  anno_out2 <- bind_cols(anno_list)
  anno_out  <- cbind(anno_out[,c(1,which(!is.element(colnames(anno_out),convertColumns)))],anno_out2)
  anno_out  <- group_annotations(anno_out[,c(1,3:dim(anno_out)[2])])
  anno_out
}





#' Returns the summary gene expression value across a group
#' 
#' @param datExpr a matrix of data (rows=genes, columns=samples)
#' @param groupVector character vector corresponding to the group (e.g., cell type)
#' @param fn Summary function to use (default = mean)
#' 
#' @return Summary matrix (genes x groups)
#'
#' @export
findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[1],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[,groupVector==groups[i]]
    if (is.null(dim(datIn)[1])) { 
      datMeans[,i] = as.numeric(datIn)
    } else { 
      datMeans[,i] = as.numeric(apply(datIn,1,fn)) }
  };    
  colnames(datMeans) = groups;
  rownames(datMeans) = rownames(datExpr)
  return(datMeans)
}


#' Adds to the statistics matrix
#' 
#' @param stat a statitics matrix (rows=genes, columns=names of hierarchy[1])
#' @param hierarchy a hierarchy character vector
#' @param cluster_info A cluster_info matrix
#' 
#' @return a statistics matrix with additional columns for other hierarchy levels
#'
#' @export
addHierarchyToStat <- function(stat, hierarchy, cluster_info, fn="mean"){
  if(length(hierarchy)==1)  return(stat)
  stat     <- stat[,cluster_info[,paste0(hierarchy[1],"_label")]]
  stat_out <- stat
  for (i in 2:length(hierarchy)){
    print(paste("........",i,"of",length(hierarchy)))
    vector   <- cluster_info[,paste0(hierarchy[i],"_label")]
    vector   <- factor(vector,levels=unique(vector[order(cluster_info[,paste0(hierarchy[i],"_id")])]))
    stat_out <- cbind(stat_out,findFromGroups(stat,vector,fn=fn))
  }
  stat_out
}




#' Generate colors and ids for categorical annotations that are factors
#'
#' @param df data frame to annotate
#' @param col name of the factor column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'   
#' @export
annotate_factor <- function(df,
                            col = NULL, base = NULL,
                            na_val = "ZZ_Missing",
                            colorset = "varibow", color_order = "sort") {
  
  # library(dplyr)
  # library(lazyeval)
  # library(viridisLite)
  
  if (class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if (class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }
  
  if (class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if (class(base) == "NULL") {
    base <- col
  }
  
  if (!is.factor(df[[col]])) {
    df[[col]] <- as.factor(df[[col]])
  }
  
  # Convert NA values and add NA level to the end
  if (sum(is.na(df[[col]])) > 0) {
    lev <- c(levels(df[[col]]), na_val)
    levels(df[[col]]) <- lev
    df[[col]][is.na(df[[col]])] <- na_val
  }
  
  x <- df[[col]]
  
  annotations <- data.frame(label = as.character(levels(x)), stringsAsFactors = F)
  
  annotations <- annotations %>%
    dplyr::mutate(id = 1:n())
  
  if (colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if (colorset == "rainbow") {
    colors <- sub("FF$", "", grDevices::rainbow(nrow(annotations)))
  } else if (colorset == "viridis") {
    colors <- sub("FF$", "", viridisLite::viridis(nrow(annotations)))
  } else if (colorset == "magma") {
    colors <- sub("FF$", "", viridisLite::magma(nrow(annotations)))
  } else if (colorset == "inferno") {
    colors <- sub("FF$", "", viridisLite::inferno(nrow(annotations)))
  } else if (colorset == "plasma") {
    colors <- sub("FF$", "", viridisLite::plasma(nrow(annotations)))
  } else if (colorset == "terrain") {
    colors <- sub("FF$", "", grDevices::terrain.colors(nrow(annotations)))
  } else if (is.character(colorset)) {
    colors <- grDevices::colorRampPalette(colorset)(nrow(annotations))
  }
  
  if (color_order == "random") {
    colors <- sample(colors, length(colors))
  }
  
  annotations <- dplyr::mutate(annotations, color = colors)
  
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))
  
  names(df)[names(df) == col] <- paste0(base, "_label")
  
  df[[paste0(col,"_label")]] <- as.character(df[[paste0(col,"_label")]]) # convert the factor to a character in the anno
  
  df <- dplyr::left_join(df, annotations, by = paste0(base, "_label"))
  
  df
}




#'Group annotation columns
#'
#'@param df the annotation dataframe to arrange
#'@param sample_col the column with unique sample ids. Default is "cell_id".
#'@param keep_order a logical value. If FALSE, will sort the annotations alphanumerically by base.
#'
#'
#'@return an annotation data frame with reordered columns
#'
#' @export
group_annotations <- function(df, sample_col = "cell_id", keep_order = TRUE) {
  labels <- names(df)[grepl("_label",names(df))]
  if(!keep_order) {
    labels <- labels[order(labels)]
  }
  bases <- sub("_label","",labels)
  
  anno_cols <- c(paste0(rep(bases,each=3),c("_id","_label","_color")))
  extras <- setdiff(names(df),anno_cols)
  
  anno <- select(df,one_of(c(sample_col,anno_cols,extras)))
  
}




#' get knn graph
#' 
#' This function gets a knn graph, which is required for creating a constellation diagram
#'
#' @param rd.dat Reduced dimension data frame where rows are cells and columns are reduced dimensional data (e.g., 50-100 principal components)
#' @param cl Vector of cluster assignments (can be numeric, character, or factor)
#' @param k Number of nearest neighbors (default=50 to align with example)
#' @param ... Other variables that I'm not sure what they are and should probably be left alone
#'
#' @return
#' @export
#'
#' @examples
get_knn_graph <- function(rd.dat, 
                          cl, 
                          k=50, 
                          ref.cells=row.names(rd.dat),
                          method="Annoy.Cosine", 
                          knn.outlier.th=2, 
                          outlier.frac.th=0.5,
                          clean.cells=row.names(rd.dat), 
                          knn.result=NULL,
                          mc.cores=10)
{
  if(is.null(knn.result)){
    ref.rd.dat = rd.dat[ref.cells,]
    knn.result = get_knn_batch(rd.dat, ref.rd.dat, k, method=method, batch.size=10000, mc.cores=mc.cores, transposed=FALSE, return.distance=TRUE)
  }
  knn  = knn.result[[1]]
  knn.dist = knn.result[[2]]
  colnames(knn) = colnames(knn.dist)=1:ncol(knn)
  knn.dist = as.data.frame(as.table(knn.dist),stringsAsFactors=FALSE)
  knn.id = as.data.frame(as.table(knn),stringsAsFactors=FALSE)
  knn.df = cbind(knn.id, knn.dist[,3])
  colnames(knn.df)=c("sample_id","k","ref_id","dist")
  knn.df$cl = cl[knn.df$sample_id]  
  knn.df$knn.cl = cl[ref.cells[knn.df$ref_id]]
  knn.df = knn.df %>% filter(k!=1)
  
  cl.knn.dist.stats = knn.df %>%  group_by(cl) %>% summarise(med=median(dist),mad=mad(dist))
  cl.knn.dist.stats =   cl.knn.dist.stats %>% mutate(th=med + knn.outlier.th * mad)
  th.med = median(cl.knn.dist.stats$th)
  cl.knn.dist.stats =   cl.knn.dist.stats %>% mutate(th=pmax(th, th.med))
  
  outlier.df=knn.df %>% left_join(cl.knn.dist.stats[,c("cl","th")]) %>% group_by(sample_id) %>% summarise(outlier = sum(dist > th))
  
  outlier = outlier.df %>% filter(outlier/(k-1)>outlier.frac.th) %>% pull(sample_id)
  knn.df = knn.df %>% filter(!sample_id %in% outlier)
  knn.cl.df = knn.df %>% group_by(cl, knn.cl) %>% summarise(Freq=n())
  colnames(knn.cl.df)[1:2]=c("cl.from","cl.to")  
  from.size = knn.cl.df %>% group_by(cl.from) %>% summarise(from.total=sum(Freq))
  to.size = knn.cl.df %>% group_by(cl.to) %>% summarise(to.total=sum(Freq))
  total = sum(knn.cl.df$Freq)
  knn.cl.df = knn.cl.df %>% left_join(from.size) %>% left_join(to.size)
  knn.cl.df = knn.cl.df %>% mutate(odds = Freq/(from.total*as.numeric(to.total)/total))
  knn.cl.df = knn.cl.df %>% mutate(pval.log = phyper(q=Freq-1, m=to.total, n=total - to.total, k=from.total, lower.tail=FALSE, log.p=TRUE))
  knn.cl.df$frac = knn.cl.df$Freq/knn.cl.df$from.total
  return(list(knn.result=knn.result, knn.cl.df=knn.cl.df,outlier=outlier))
}



#' Get binary (aka beta) score
#'
#' Returns a beta score which indicates the binaryness of a gene across clusters.  High scores
#'   (near 1) indicate that a gene is either on or off in nearly all cells of every cluster.
#'   Scores near 0 indicate a cells is non-binary (e.g., not expressed, ubiquitous, or
#'   randomly expressed).  This value is used for gene filtering prior to defining clustering.
#'
#' @param propExpr a matrix of proportions of cells (rows) in a given cluster (columns) with
#'   CPM/FPKM > 1 (or 0, HCT uses 1)
#' @param returnScore if TRUE returns the score, if FALSE returns the ranks
#' @param spec.exp scaling factor (recommended to leave as default)
#'
#' @return returns a numeric vector of beta score (or ranks)
#'
#' @export
getBetaScore_fast <- function(propExpr, returnScore = TRUE, spec.exp = 2) {
  
  # Ensure propExpr is a matrix for consistent apply behavior
  propExpr <- as.matrix(propExpr)
  n_cols <- ncol(propExpr) # N in the formula above
  
  # Define a vectorized and optimized calculation for a single row
  calc_beta_optimized <- function(y_row, n_cols, spec.exp) {
    eps1 <- 1e-10
    
    # Numerator: sum of squared differences (spec.exp = 2)
    # sum((yi - yj)^2 for i<j) = N * sum(yk^2) - (sum(yk))^2
    if (spec.exp == 2) {
      sum_sq_diffs <- n_cols * sum(y_row^2) - sum(y_row)^2
    } else {
      # If spec.exp is not 2, we must compute pairwise differences.
      y_diffs <- as.vector(outer(y_row, y_row, FUN = "-"))
      sum_sq_diffs <- sum(y_diffs^spec.exp) / 2 
    }
    
    if (n_cols > 1) {
      sum_abs_diffs <- sum(abs(outer(y_row, y_row, FUN = "-"))) / 2 # Divide by 2 for unique pairs
    } else { # Handle single-column case
      sum_abs_diffs = 0 # No differences if only one element
    }
    
    score1 <- sum_sq_diffs / (sum_abs_diffs + eps1)
    return(score1)
  }
  
  # Apply the optimized calculation across rows
  betaScore <- apply(propExpr, 1, calc_beta_optimized, n_cols = n_cols, spec.exp = spec.exp)
  
  # Handle NA values
  betaScore[is.na(betaScore)] <- 0
  
  if (returnScore) {
    return(betaScore)
  }
  
  # If not returning score, return rank
  scoreRank <- rank(-betaScore, ties.method = "average") # Using average for tie-breaking
  return(scoreRank)
}
