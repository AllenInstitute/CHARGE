# This script compiles all the files in this folder for the subset of requested genes

###################################################################################################
### This code block adds gene categories to the output
wd   <- "gene_lists/"
fnIn <- dir(wd)
fnIn <- fnIn[grepl(".csv",fnIn)|grepl(".txt",fnIn)]
gene_lists <- NULL
if(length(fnIn)>0){
  
  for (i in 1:length(fnIn)) {
    ext = substr(fnIn[i], nchar(fnIn[i]) - 2, nchar(fnIn[i]))
    if (ext == "csv") {
      datIn = read.csv(paste0(wd,fnIn[i]))
      if(dim(datIn)[2]>1){
        datIn = datIn[, 1:2]
      } else {
        datIn = NULL
      }
    } else {
      datIn = scan(paste0(wd,fnIn[i]), what = "character", sep = "\n")
      datIn = cbind(datIn[2:length(datIn)], datIn[1])
    }
    if(!is.null(datIn)){
      colnames(datIn) = c("Gene", "Category")
      gene_lists = rbind(gene_lists, datIn)
    }
  }
  gene_lists$Gene <- toupper(iconv(gene_lists$Gene,"WINDOWS-1252","UTF-8"))
  gene_lists <- gene_lists[rowSums(is.na(gene_lists))==0,]
  
  gene <- toupper(output$gene)
  gene_categories <- rep("",length(gene))
  names(gene_categories) <- gene
  for (gn in intersect(gene_lists[,1],gene)){ 
    cats <- unique(as.character(gene_lists[gene_lists[,1]==gn,2]))
    gene_categories[gn] <- paste(cats,collapse=" | ")
  }
  gene_categories <- as.character(gene_categories)
  output$gene_categories________________________________________________________ = gene_categories
}
###################################################################################################




### COPY OF check_genes from scrattch.vis

#' Check input genes against a vector of gene names
#' 
#' @param genes An input character string (or character vector) of gene symbols
#' @param gene_reference A reference character vector of gene symbols
#' @param result Whether to return a vector of matched or unmatched genes, or a list of both. Options are "matched" (default), "unmatched", or "both".
#' 
check_genes <- function(genes, 
                        gene_reference,
                        result = "matched") {
  
  # Add some input QC
  if (class(genes) != "character")
    stop("Your genes input is not a character string or vector.")
  if (class(gene_reference) != "character")
    stop("Your reference input is not a character vector.")
  
  
  if (length(genes) == 1) {
    raw_genes <- unique(split_text(genes))
  } else {
    raw_genes <- unique(genes)
  }
  
  
  # Convert gene_reference to lowercase and remove "-" and " " for matching
  match_genes <- tolower(gsub("[- .]+","_",gene_reference))
  
  # Remove leading X from Riken genes
  match_genes[grepl("x[0-9]+.+rik",match_genes)] <- sub("^x","",match_genes[grepl("x[0-9]+.+rik",match_genes)])
  
  # for loop will retain the order of the genes.
  good_genes <- character()
  bad_genes <- character()
  
  for (x in raw_genes) {
    this_gene <- tolower(gsub("[- .]+","_", x))
    
    if (this_gene %in% match_genes) {
      good_genes <- c(good_genes, gene_reference[match_genes == this_gene][1])
    } else {
      bad_genes <- c(bad_genes, x)
    }
  }
  
  if (result == "matched") {
    unique(good_genes)
  } else if (result == "not_matched") {
    unique(bad_genes)
  } else if (result == "both") {
    
    list(matched = unique(good_genes),
         not_matched = unique(bad_genes))
    
  }
}
