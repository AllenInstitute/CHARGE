

#########################################################################
###                                                                   ###
### DO NOT ERASE THE FILE IN THIS FOLDER CALLED read_gene_lists.r !!! ###
###                                                                   ###
#########################################################################


User-inputted files can be deposited IN THIS FOLDER in one of two formats:

1) Text files (must end in ".txt") with one list per file, where the first line is the list descriptor and the remaining lines are gene names corresponding to that list, with one gene per line. For example,
  Ribosome
  RPS4
  RPS8
  ...

2) Gene / category files (must be csv files), where the first line is the column headers corresponding to Genes and Lists, and the remaining lines correspond to the genes in each list, for any number of genes and lists. For example: 
  Gene, Category 
  RPS4, Ribosome 
  RPS8, Ribosome 
  ... 
  NDUF1, Mitohcondria 
  NDUF3, Mitochondria 
  ... 
  MAPT, AlzheimersDisease 
  PSEN1, AlzheimersDisease 
  PSEN2, AlzheimersDisease 
  ...