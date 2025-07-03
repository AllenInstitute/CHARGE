# For each table...

#... what is is called in the main page?
table_names <- c(
  "Macaque Basal Ganglia - preprint"
)

#... what category each table is included in on the main page?
categories <- factor(c(
  "Cell type classification"
  ),
# This is the order they will show up on in the list. **MAKE SURE THERE ARE NO TYPOS!**
levels = c("Cell type classification")) 


#... where is the cell x annotation table?
# Note that this can be a csv file, a feather directory, or a scrattch.taxonomy h5ad file.  It can also be local or on the web.
table_locations <- c(
  "s3://pga-genomics-wg-802451596237-us-west-2/CHARGE_data/test_CHARGE.RData" #"../sandbox/test_CHARGE.RData"  # 
)


descriptions   <- c(

  "The Human and Mammalian Brain Atlas (HMBA) consortia has created a unified taxonomy of the mammalian basal ganglia, with >2 million cells collected from humans, macaques, and marmosets, and linkages to existing mouse brain data. More details are available here: https://alleninstitute.github.io/HMBA_BasalGanglia_Consensus_Taxonomy/.  The table here includes data about cell type assignments (cluster, group, subclass, class, and neighborhood) for the macaque taxonomy."
  
)


############################################
## DO NOT EDIT ANYTHING BELOW THIS POINT! ##
############################################

categories = factor(c(as.character(categories)),levels = unique(c(levels(categories))))

# Convert above into a data frame
table_info <- data.frame(table_name   = table_names,
                         table_loc    = table_locations,
                         description  = descriptions
)

# Convert table names into a nested list
table_name <- list()
for (cat in levels(categories)){
  table_name[[cat]] <- c("Select comparison table...", table_names[categories==cat])
}
table_name[["Enter your own location"]] = c("Enter your own location")
category = "Enter your own location"