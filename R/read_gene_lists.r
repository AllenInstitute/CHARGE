# source("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/common/reference_files/gene_lists/read_gene_lists.r", local=TRUE)
# This script compiles all the files in this folder for the subset of requested genes


###################################################################################################
### This code block adds gene categories to the output
wd   <- "gene_lists/"
#wd   <- ("\\\\allen\\programs\\celltypes\\workgroups\\rnaseqanalysis\\shiny\\common\\reference_files\\gene_lists\\")
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

