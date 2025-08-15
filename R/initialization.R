# For each table...

#... what is is called in the main page?
table_names <- c(
  "Middle temporal gyrus (2022; NEW): SMART-seq",
  "Middle temporal gyrus (2018): SMART-seq",
  "Primary motor cortex (2020): SMART-seq",
  "Multiple neocortical areas (2019): SMART-seq",
  "Comparative LGN (2018): SMART-seq",
  "Whole cortex and hippocampus (2021): 10X seq",
  "Whole cortex and hippocampus (2021): SMART-seq",
  "V1 & ALM (2018): SMART-seq",
  "Comparative LGN (2018): SMART-seq",
  "Comparative LGN (2018): SMART-seq"
)

#... what category each table is included in on the main page?
categories <- factor(c(
  "Human brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Mouse brain cell types",
  "Mouse brain cell types",
  "Mouse brain cell types",
  "Mouse brain cell types",
  "Other mammalian brain cell types"
  ),
# This is the order they will show up on in the list. **MAKE SURE THERE ARE NO TYPOS!**
levels = c(
  "Human brain cell types",
  "Mouse brain cell types",
  "Other mammalian brain cell types"
  )
) 


#... where is the CHARGE.RData file?
S3_folder = "s3://pga-genomics-wg-802451596237-us-west-2/CHARGE_data/test_CHARGE.RData"
table_locations <- c(
  paste0(S3_folder,"Human_MTG_SEAAD_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_MTG_SMART_seq_08082025_CHARGE.RData"),
  paste0(S3_folder,"Human_M1_10X_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_neocortex_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_LGN_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_cortex_hippocampus_10X_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_cortex_hippocampus_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_VISp_ALM_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_LGN_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Macaque_LGN_SMART_seq_04042025_CHARGE.RData")
)


descriptions   <- c(

  "Explore cell types in human middle temporal gyrus (MTG) circa 2025 as described in Gabitto, Travaglini, et al 2024 (Nature Neuroscience; https://doi.org/10.1038/s41593-024-01774-5)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad. Note that the data included here is from five adult healthy donors and does NOT include data from the 84 aged donors as part of the SEA-AD study (described at SEA-AD.org).",
  
  "Explore cell types in human middle temporal gyrus (MTG) circa 2018 as described in Hodge, Bakken, et al 2019 (Nature Neuroscience; https://doi.org/10.1038/s41586-019-1506-7)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq.",
  
  "Explore cell types in human primary motor cortex (M1) circa 2021 as described in Bakken et al 2021 (Nature; https://doi.org/10.1038/s41586-021-03465-8)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x, and are also included as part of the Cell Type Knowledge Explorer (https://knowledge.brain-map.org/celltypes).",
  
  "Explore DRAFT cell types in multiple human neocortical areas circa 2023. Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq. These cell types represent an OLD VERSION of the cell types described in Jorstad et al 2023 (Science; https://doi.org/10.1126/science.adf6812).",
  
  "Explore cell types in human lateral geniculate nucleus (LGN) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn.",
  
  "Explore cell types from multiple cortical areas and the hippocampal formation in mouse circa 2021 as described in Yao et al 2021 (Cell; https://doi.org/10.1016/j.cell.2021.04.021)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x. These data include ~1.1 million cells collected using 10X Genomics droplet-based sequencing and are aligned to the same taxonomy as 'Whole Cortex and Hippocampus (2021): SMART seq'.",
  
  "Explore cell types from multiple cortical areas and the hippocampal formation in mouse circa 2021 as described in Yao et al 2021 (Cell; https://doi.org/10.1016/j.cell.2021.04.021)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq. These data include ~77 thousand cells collected using SMART sequencing and are aligned to the same taxonomy as 'Whole Cortex and Hippocampus (2021): 10x seq'.",
  
  "Explore cell types in mouse primary visual cortex (VISp) and anterior lateral motor cortex (ALM) circa 2018 as described in Tasic et al 2018 (Nature; https://doi.org/10.1038/s41586-018-0654-5)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.",
  
  "Explore cell types in mouse dorsolateral geniculate complex (LGd) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn.",
  
  "Explore cell types in macaque lateral geniculate nucleus (LGN) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn."
  
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