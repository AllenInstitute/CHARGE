# For each table...

#... what is is called in the main page?
table_names <- c(
  "Consensus Basal Ganglia (2025): Human, 10x seq",
  "Consensus Basal Ganglia (2025): Macaque, 10x seq",
  "Consensus Basal Ganglia (2025): Marmoset, 10x seq",
  "Human Middle temporal gyrus (2022): 10x seq",
  "Human Middle temporal gyrus (2018): SMART-seq",
  "Human Primary motor cortex (2020): SMART-seq",
  "Human Multiple neocortical areas (2019): SMART-seq",
  "Mouse Whole cortex and hippocampus (2021): 10X seq",
  "Mouse Whole cortex and hippocampus (2021): SMART-seq",
  "Mouse V1 & ALM (2018): , SMART-seq",
  "Comparative LGN (2018): Human, SMART-seq",
  "Comparative LGN (2018): Mouse, SMART-seq",
  "Comparative LGN (2018): Macaque, SMART-seq"
)

#... what category each table is included in on the main page?
categories <- factor(c(
  "Human brain cell types",
  "Other mammalian brain cell types",
  "Other mammalian brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Human brain cell types",
  "Mouse brain cell types",
  "Mouse brain cell types",
  "Mouse brain cell types",
  "Human brain cell types",
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
S3_folder = "s3://pga-genomics-wg-802451596237-us-west-2/CHARGE_data/"
table_locations <- c(
  paste0(S3_folder,"Human_HMBA_basalganglia_AIT_pre-print_CHARGE.RData"),
  paste0(S3_folder,"Macaque_HMBA_basalganglia_AIT_pre-print_CHARGE.RData"),
  paste0(S3_folder,"Marmoset_HMBA_basalganglia_AIT_pre-print_CHARGE.RData"),
  paste0(S3_folder,"Human_MTG_SEAAD_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_MTG_SMART_seq_08082025_CHARGE.RData"),
  paste0(S3_folder,"Human_M1_10X_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_neocortex_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_cortex_hippocampus_10X_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_cortex_hippocampus_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_VISp_ALM_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Human_LGN_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Mouse_LGN_SMART_seq_04042025_CHARGE.RData"),
  paste0(S3_folder,"Macaque_LGN_SMART_seq_04042025_CHARGE.RData")
)


descriptions   <- c(

  "Explore cell types in HUMAN basal ganglia. These data are part of a cross-species consensus taxonomy released 11/13/2025 on Allen Brain Map. Download the data and explore the comprehensive, highly granular cell atlas spanning adult human, macaque, and marmoset brains that links brain structure, function and cellular architecture at https://brain-map.org/consortia/hmba/hmba-release-basal-ganglia. Note that the colors in these plots do not currently match reported cell type colors, and will be updated soon.",
  
  "Explore cell types in MACAQUE basal ganglia. These data are part of a cross-species consensus taxonomy released 11/13/2025 on Allen Brain Map. Download the data and explore the comprehensive, highly granular cell atlas spanning adult human, macaque, and marmoset brains that links brain structure, function and cellular architecture at https://brain-map.org/consortia/hmba/hmba-release-basal-ganglia. Note that the colors in these plots do not currently match reported cell type colors, and will be updated soon.",
  
  "Explore cell types in MARMOSET basal ganglia. These data are part of a cross-species consensus taxonomy released 11/13/2025 on Allen Brain Map. Download the data and explore the comprehensive, highly granular cell atlas spanning adult human, macaque, and marmoset brains that links brain structure, function and cellular architecture at https://brain-map.org/consortia/hmba/hmba-release-basal-ganglia. Note that the colors in these plots do not currently match reported cell type colors, and will be updated soon.",

  "Explore cell types in human middle temporal gyrus (MTG) circa 2025 as described in Gabitto, Travaglini, et al 2024 (Nature Neuroscience; https://doi.org/10.1038/s41593-024-01774-5)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad. Note that the data included here is from five adult healthy donors and does NOT include data from the 84 aged donors as part of the SEA-AD study (described at SEA-AD.org).",
  
  "Explore cell types in human middle temporal gyrus (MTG) circa 2018 as described in Hodge, Bakken, et al 2019 (Nature Neuroscience; https://doi.org/10.1038/s41586-019-1506-7)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq.",
  
  "Explore cell types in human primary motor cortex (M1) circa 2021 as described in Bakken et al 2021 (Nature; https://doi.org/10.1038/s41586-021-03465-8)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x, and are also included as part of the Cell Type Knowledge Explorer (https://knowledge.brain-map.org/celltypes).",
  
  "Explore DRAFT cell types in multiple human neocortical areas circa 2023. Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq. These cell types represent an OLD VERSION of the cell types described in Jorstad et al 2023 (Science; https://doi.org/10.1126/science.adf6812).",
  
  "Explore cell types from multiple cortical areas and the hippocampal formation in mouse circa 2021 as described in Yao et al 2021 (Cell; https://doi.org/10.1016/j.cell.2021.04.021)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x. These data include ~1.1 million cells collected using 10X Genomics droplet-based sequencing and are aligned to the same taxonomy as 'Whole Cortex and Hippocampus (2021): SMART seq'.",
  
  "Explore cell types from multiple cortical areas and the hippocampal formation in mouse circa 2021 as described in Yao et al 2021 (Cell; https://doi.org/10.1016/j.cell.2021.04.021)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq. These data include ~77 thousand cells collected using SMART sequencing and are aligned to the same taxonomy as 'Whole Cortex and Hippocampus (2021): 10x seq'.",
  
  "Explore cell types in mouse primary visual cortex (VISp) and anterior lateral motor cortex (ALM) circa 2018 as described in Tasic et al 2018 (Nature; https://doi.org/10.1038/s41586-018-0654-5)! Underlying data and additonal visualizations are available at https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.",
  
  "Explore cell types in human lateral geniculate nucleus (LGN) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn.",
  
  "Explore cell types in mouse dorsolateral geniculate complex (LGd) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn.",
  
  "Explore cell types in macaque lateral geniculate nucleus (LGN) circa 2021 as described in Bakken, van Velthoven, Menon, et al, 2021 (eLife; https://doi.org/10.7554/eLife.64875)! Underlying data are available at https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn."
  
)

urls <- c(
  "https://knowledge.brain-map.org/abcatlas#AQQBM0NQSVAyODA0U0JEWDFEWDFSVwACUzVLUVRGSENSS0hWWlZSUjVIUgADBQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAEdsdXQgU2VybyBEb3BhAFN1YnBhbGxpdW0gR0FCQQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAAATNYQk9GQlIzMDJFWUMwVVFENTEAAgAAAVNRNldKTzBHTlVKR0ZVOExIVlYAAgAAAUhLMTRPMUFTQU5QWUhPTkZMTFEAAgAAATZFV0JTRUo0WUZWR0tJMkFQQkwAAgAABAEAAoVsLWWGxlUGA4dkj1mGLt%2FTBFdNUDhOTUo1RDhYTzFBT0kxTDkABYHYoeuEthGZhNqUXYUChoYGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISFpFWVhTUU9FRE5EMlE2TTk3TQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNlQ5NjI2MDEyNlRWNEtMR0RYRwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFWTTRFNko5RFBZWUQ3RTBCSE5CAAJLR1pMQVdJUENSNFc5OEY1TEQ2AAMFAUUzSEdCWUhSQzk4R0E2Rlc5RjYAAgROb25uZXVyb24AR2x1dCBTZXJvIERvcGEAU3VicGFsbGl1bSBHQUJBAFN1YnBhbGxpdW0gR0FCQS1HbHV0AAABM1hCT0ZCUjMwMkVZQzBVUUQ1MQACAAABU1E2V0pPMEdOVUpHRlU4TEhWVgACAAABSEsxNE8xQVNBTlBZSE9ORkxMUQACAAABNkVXQlNFSjRZRlZHS0kyQVBCTAACAAAEAQAChTi2eIXEYQYDhrJSsoa4BiYEV01QOE5NSjVEOFhPMUFPSTFMOQAFhBTEFoOMxVKEYmm%2FhLTnSwYABwAABQEBVE1FTTEzMkMAAAYBAQJISzE0TzFBU0FOUFlIT05GTExRAAN%2BAAAABAAACDlBVVRFOVMxMEtBRUMyVVpPMzEACVBPWjJIQ1BCVDYwRFNESjhVQTcACgALAVkzUDFKTTFFTkNHV0VXOFVTVFYAAkNSRUhSQ0syWVhSV0pLMkowQUwAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAMAAABTlhKOU5aQVFBNU1FSFFOT1dKUQACVEdCUTJGRVYyN04wQ1BMS1ZBWgADAQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAFN1YnBhbGxpdW0gR0FCQQBHbHV0IFNlcm8gRG9wYQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAABAEAAoQg%2BpWEuiVQA4aoVQiFduVlBFdNUDhOTUo1RDhYTzFBT0kxTDkABYNsbzGCNuJBg8AJxIM4AlIGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISUE2WTZTTjdRWlNJRUowRlowVQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNzQxNkI1WEpEU1kwSERGOTNSTwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFFSVVBTFBZTkdaNVZVN1kyWTAyAAJRWlEyUUJNVTRUNUk0N0QxSzRGAAMABAEBAoEi0luBCLU%2BA4RS41GD6nvJAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAIRlRMQVVPVjZTTEIwNzFTVk5VOQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAIIAwAA",
  "https://knowledge.brain-map.org/abcatlas#AQQBM0NQSVAyODA0U0JEWDFEWDFSVwACUzVLUVRGSENSS0hWWlZSUjVIUgADBQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAEdsdXQgU2VybyBEb3BhAFN1YnBhbGxpdW0gR0FCQQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAAATNYQk9GQlIzMDJFWUMwVVFENTEAAgAAAVNRNldKTzBHTlVKR0ZVOExIVlYAAgAAAUhLMTRPMUFTQU5QWUhPTkZMTFEAAgAAATZFV0JTRUo0WUZWR0tJMkFQQkwAAgAABAEAAoVsLWWGxlUGA4dkj1mGLt%2FTBFdNUDhOTUo1RDhYTzFBT0kxTDkABYHYoeuEthGZhNqUXYUChoYGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISFpFWVhTUU9FRE5EMlE2TTk3TQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNlQ5NjI2MDEyNlRWNEtMR0RYRwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFWTTRFNko5RFBZWUQ3RTBCSE5CAAJLR1pMQVdJUENSNFc5OEY1TEQ2AAMFAUUzSEdCWUhSQzk4R0E2Rlc5RjYAAgROb25uZXVyb24AR2x1dCBTZXJvIERvcGEAU3VicGFsbGl1bSBHQUJBAFN1YnBhbGxpdW0gR0FCQS1HbHV0AAABM1hCT0ZCUjMwMkVZQzBVUUQ1MQACAAABU1E2V0pPMEdOVUpHRlU4TEhWVgACAAABSEsxNE8xQVNBTlBZSE9ORkxMUQACAAABNkVXQlNFSjRZRlZHS0kyQVBCTAACAAAEAQAChTi2eIXEYQYDhrJSsoa4BiYEV01QOE5NSjVEOFhPMUFPSTFMOQAFhBTEFoOMxVKEYmm%2FhLTnSwYABwAABQEBVE1FTTEzMkMAAAYBAQJISzE0TzFBU0FOUFlIT05GTExRAAN%2BAAAABAAACDlBVVRFOVMxMEtBRUMyVVpPMzEACVBPWjJIQ1BCVDYwRFNESjhVQTcACgALAVkzUDFKTTFFTkNHV0VXOFVTVFYAAkNSRUhSQ0syWVhSV0pLMkowQUwAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAMAAABTlhKOU5aQVFBNU1FSFFOT1dKUQACVEdCUTJGRVYyN04wQ1BMS1ZBWgADAQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAFN1YnBhbGxpdW0gR0FCQQBHbHV0IFNlcm8gRG9wYQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAABAEAAoQg%2BpWEuiVQA4aoVQiFduVlBFdNUDhOTUo1RDhYTzFBT0kxTDkABYNsbzGCNuJBg8AJxIM4AlIGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISUE2WTZTTjdRWlNJRUowRlowVQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNzQxNkI1WEpEU1kwSERGOTNSTwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFFSVVBTFBZTkdaNVZVN1kyWTAyAAJRWlEyUUJNVTRUNUk0N0QxSzRGAAMABAEBAoEi0luBCLU%2BA4RS41GD6nvJAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAIRlRMQVVPVjZTTEIwNzFTVk5VOQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAIIAwAA",
  "https://knowledge.brain-map.org/abcatlas#AQQBM0NQSVAyODA0U0JEWDFEWDFSVwACUzVLUVRGSENSS0hWWlZSUjVIUgADBQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAEdsdXQgU2VybyBEb3BhAFN1YnBhbGxpdW0gR0FCQQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAAATNYQk9GQlIzMDJFWUMwVVFENTEAAgAAAVNRNldKTzBHTlVKR0ZVOExIVlYAAgAAAUhLMTRPMUFTQU5QWUhPTkZMTFEAAgAAATZFV0JTRUo0WUZWR0tJMkFQQkwAAgAABAEAAoVsLWWGxlUGA4dkj1mGLt%2FTBFdNUDhOTUo1RDhYTzFBT0kxTDkABYHYoeuEthGZhNqUXYUChoYGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISFpFWVhTUU9FRE5EMlE2TTk3TQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNlQ5NjI2MDEyNlRWNEtMR0RYRwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFWTTRFNko5RFBZWUQ3RTBCSE5CAAJLR1pMQVdJUENSNFc5OEY1TEQ2AAMFAUUzSEdCWUhSQzk4R0E2Rlc5RjYAAgROb25uZXVyb24AR2x1dCBTZXJvIERvcGEAU3VicGFsbGl1bSBHQUJBAFN1YnBhbGxpdW0gR0FCQS1HbHV0AAABM1hCT0ZCUjMwMkVZQzBVUUQ1MQACAAABU1E2V0pPMEdOVUpHRlU4TEhWVgACAAABSEsxNE8xQVNBTlBZSE9ORkxMUQACAAABNkVXQlNFSjRZRlZHS0kyQVBCTAACAAAEAQAChTi2eIXEYQYDhrJSsoa4BiYEV01QOE5NSjVEOFhPMUFPSTFMOQAFhBTEFoOMxVKEYmm%2FhLTnSwYABwAABQEBVE1FTTEzMkMAAAYBAQJISzE0TzFBU0FOUFlIT05GTExRAAN%2BAAAABAAACDlBVVRFOVMxMEtBRUMyVVpPMzEACVBPWjJIQ1BCVDYwRFNESjhVQTcACgALAVkzUDFKTTFFTkNHV0VXOFVTVFYAAkNSRUhSQ0syWVhSV0pLMkowQUwAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAMAAABTlhKOU5aQVFBNU1FSFFOT1dKUQACVEdCUTJGRVYyN04wQ1BMS1ZBWgADAQFFM0hHQllIUkM5OEdBNkZXOUY2AAIETm9ubmV1cm9uAFN1YnBhbGxpdW0gR0FCQQBHbHV0IFNlcm8gRG9wYQBTdWJwYWxsaXVtIEdBQkEtR2x1dAAABAEAAoQg%2BpWEuiVQA4aoVQiFduVlBFdNUDhOTUo1RDhYTzFBT0kxTDkABYNsbzGCNuJBg8AJxIM4AlIGAAcAAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAISUE2WTZTTjdRWlNJRUowRlowVQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBNzQxNkI1WEpEU1kwSERGOTNSTwACQ1JFSFJDSzJZWFJXSksySjBBTAADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAFFSVVBTFBZTkdaNVZVN1kyWTAyAAJRWlEyUUJNVTRUNUk0N0QxSzRGAAMABAEBAoEi0luBCLU%2BA4RS41GD6nvJAAUABgEBAkhLMTRPMUFTQU5QWUhPTkZMTFEAA34AAAAEAAAIRlRMQVVPVjZTTEIwNzFTVk5VOQAJUE9aMkhDUEJUNjBEU0RKOFVBNwAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAIIAwAA",
  "https://knowledge.brain-map.org/abcatlas?_gl=1*1r4e8w7*_gcl_au*ODk2OTAwNjYwLjE3NzA2NjU3MDY.*_ga*MTI5NjA0MjMyMC4xNzUyNzk0NTk3*_ga_MFZ8RZNMY4*czE3Nzc2NTM3MzckbzE1MCRnMSR0MTc3NzY1NDA2NyRqNTMkbDAkaDQ2MzM3NTA4OQ..#AQIBUU5aTUdPWlhINEVQQ0lVM1YzQQACUlJSNkZOM1dMUjJMOUlWV1hRVAADAAQBAAKMxIiIjQ4RvgONfAb1jrg92AQyTlFUSUU3VEFNUDhQUUFITzRQAAWJT5zAhgPp14lOnMCLmAGzBgAABQAGAQECMFRURjFQUEZSTVhTVlZXTzcyMgADfgAAAAQAAAgyQ0tPUUlLMk45SUU3VDRXSThGAAlGSFVVSzE5QTVYSjdGQzlZWVQ0AAoAAAFMNzNYVkoxR0dURFJTMVJTNUs0AAJHTE1FSjBEODlLODNDQ0FSMUhPAAMABAEBAn0XMUt%2FSDOAA4SSiXOFPDFeAAUABgEBAjBUVEYxUFBGUk1YU1ZWV083MjIAA34AAAAEAAAIOThKRTNaMUlMU0RDSUVNQTZMUQAJVU1TVlhURElBWlRBRktHRTQzVAAKAAACAwA%3D",
  "https://knowledge.brain-map.org/abcatlas?_gl=1*1r4e8w7*_gcl_au*ODk2OTAwNjYwLjE3NzA2NjU3MDY.*_ga*MTI5NjA0MjMyMC4xNzUyNzk0NTk3*_ga_MFZ8RZNMY4*czE3Nzc2NTM3MzckbzE1MCRnMSR0MTc3NzY1NDA2NyRqNTMkbDAkaDQ2MzM3NTA4OQ..#AQIBUU5aTUdPWlhINEVQQ0lVM1YzQQACUlJSNkZOM1dMUjJMOUlWV1hRVAADAAQBAAKMxIiIjQ4RvgONfAb1jrg92AQyTlFUSUU3VEFNUDhQUUFITzRQAAWJT5zAhgPp14lOnMCLmAGzBgAABQAGAQECMFRURjFQUEZSTVhTVlZXTzcyMgADfgAAAAQAAAgyQ0tPUUlLMk45SUU3VDRXSThGAAlGSFVVSzE5QTVYSjdGQzlZWVQ0AAoAAAFMNzNYVkoxR0dURFJTMVJTNUs0AAJHTE1FSjBEODlLODNDQ0FSMUhPAAMABAEBAn0XMUt%2FSDOAA4SSiXOFPDFeAAUABgEBAjBUVEYxUFBGUk1YU1ZWV083MjIAA34AAAAEAAAIOThKRTNaMUlMU0RDSUVNQTZMUQAJVU1TVlhURElBWlRBRktHRTQzVAAKAAACAwA%3D",
  "https://knowledge.brain-map.org/abcatlas#AQIBOU02NEo3UzBBM1pOSUpBNFZTNgACR1o5NFo4M1JYVVRFQ0hZS1FPNAADAQFKR1VSQkNWUlY4QlhKWjY0MlZPAAIBSHVtYW4gTTFDAAAEAQECgiZ0lYKOGw0DhjrDM4YcsqAABQAGAQECTVhFRk9GTlNUTVNGUkhGTkQxOQADfgAAAAQAAAhVTVpKN1pMRjhHNUZINVM5WTRVAAlDM1JSVkFLMThIRzZRMUpONlpRAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAADAAAATlNNjRKN1MwQTNaTklKQTRWUzYAAlI1RldVVlRFQTJPVzlOSlFNWTkAAwEBSkdVUkJDVlJWOEJYSlo2NDJWTwACAUh1bWFuIE0xQwAABAEBAoGAqDN9g8xlA4VQGvSGvI5SAAUABgEBAlgzSzdLMkVHNlVNWE9OMFQ5QkUAA34AAAAEAAAINVQ0QzNNQzNMTlpXRUJERlVFWgAJQzNSUlZBSzE4SEc2UTFKTjZaUQAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAIDAwAA",
  "https://knowledge.brain-map.org/abcatlas#AQIBOU02NEo3UzBBM1pOSUpBNFZTNgACR1o5NFo4M1JYVVRFQ0hZS1FPNAADAgFKR1VSQkNWUlY4QlhKWjY0MlZPAAIAAAE4QVBaRlA0NEoxSEZNVExZNUxKAAIBQ2VyZWJyYWwgY29ydGV4AAAEAQECgiZ0lYKOGw0DhjrDM4YcsqAABQAGAQECTVhFRk9GTlNUTVNGUkhGTkQxOQADfgAAAAQAAAhVTVpKN1pMRjhHNUZINVM5WTRVAAlDM1JSVkFLMThIRzZRMUpONlpRAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAADAAAATlNNjRKN1MwQTNaTklKQTRWUzYAAlI1RldVVlRFQTJPVzlOSlFNWTkAAwIBSkdVUkJDVlJWOEJYSlo2NDJWTwACAAABOEFQWkZQNDRKMUhGTVRMWTVMSgACAUNlcmVicmFsIGNvcnRleAAABAEBAoGAqDN9g8xlA4VQGvSGvI5SAAUABgEBAlgzSzdLMkVHNlVNWE9OMFQ5QkUAA34AAAAEAAAINVQ0QzNNQzNMTlpXRUJERlVFWgAJQzNSUlZBSzE4SEc2UTFKTjZaUQAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwAAAIDAwAA",
  "https://knowledge.brain-map.org/abcatlas?_gl=1*1tvzynf*_gcl_au*ODk2OTAwNjYwLjE3NzA2NjU3MDY.*_ga*MTI5NjA0MjMyMC4xNzUyNzk0NTk3*_ga_MFZ8RZNMY4*czE3Nzc2NTM3MzckbzE1MCRnMSR0MTc3NzY1MzczNyRqNjAkbDAkaDQ2MzM3NTA4OQ..#AQIBSzlKTjIzUDI0S1FDR0s5VTc1QQACSFNZWlBaVzE2NjlVODIxQldZUAADAAQBAAKDgDx7g46YHgOEuBCEhSrCAwQyTlFUSUU3VEFNUDhQUUFITzRQAAWBr6ZKgemsDoGggUeAktXoBgAHAAAFAAYBAQJGUzAwRFhWMFQ5UjFYOUZKNFFFAAN%2BAAAABAAACFZGT0ZZUEZRR1JLVURRVVozRkYACUxWREJKQVc4Qkk1WVNTMVFVQkcACgALAVRMT0tXQ0w5NVJVMDNEOVBFVEcAAjczR1ZURFhERUdFMjdNMlhKTVQAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAAAUFQOEpOTjVMWUFCR1ZNR0tZMUIAAlExTkNXV1BHNkZaMEROSVhKQlEAAwAEAQECgazlhIG8aWIDhMwEVIT%2BMccABQAGAQECRlMwMERYVjBUOVIxWDlGSjRRRQADfgAAAAQAAAhHNEk0R0ZKWEpCOUFUWjNQVFgxAAlMVkRCSkFXOEJJNVlTUzFRVUJHAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAAAAIDAA%3D%3D",
  "https://knowledge.brain-map.org/abcatlas?_gl=1*1tvzynf*_gcl_au*ODk2OTAwNjYwLjE3NzA2NjU3MDY.*_ga*MTI5NjA0MjMyMC4xNzUyNzk0NTk3*_ga_MFZ8RZNMY4*czE3Nzc2NTM3MzckbzE1MCRnMSR0MTc3NzY1MzczNyRqNjAkbDAkaDQ2MzM3NTA4OQ..#AQIBSzlKTjIzUDI0S1FDR0s5VTc1QQACSFNZWlBaVzE2NjlVODIxQldZUAADAAQBAAKDgDx7g46YHgOEuBCEhSrCAwQyTlFUSUU3VEFNUDhQUUFITzRQAAWBr6ZKgemsDoGggUeAktXoBgAHAAAFAAYBAQJGUzAwRFhWMFQ5UjFYOUZKNFFFAAN%2BAAAABAAACFZGT0ZZUEZRR1JLVURRVVozRkYACUxWREJKQVc4Qkk1WVNTMVFVQkcACgALAVRMT0tXQ0w5NVJVMDNEOVBFVEcAAjczR1ZURFhERUdFMjdNMlhKTVQAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAAAUFQOEpOTjVMWUFCR1ZNR0tZMUIAAlExTkNXV1BHNkZaMEROSVhKQlEAAwAEAQECgazlhIG8aWIDhMwEVIT%2BMccABQAGAQECRlMwMERYVjBUOVIxWDlGSjRRRQADfgAAAAQAAAhHNEk0R0ZKWEpCOUFUWjNQVFgxAAlMVkRCSkFXOEJJNVlTUzFRVUJHAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAAAAIDAA%3D%3D",
  "https://knowledge.brain-map.org/abcatlas?_gl=1*1tvzynf*_gcl_au*ODk2OTAwNjYwLjE3NzA2NjU3MDY.*_ga*MTI5NjA0MjMyMC4xNzUyNzk0NTk3*_ga_MFZ8RZNMY4*czE3Nzc2NTM3MzckbzE1MCRnMSR0MTc3NzY1MzczNyRqNjAkbDAkaDQ2MzM3NTA4OQ..#AQIBSzlKTjIzUDI0S1FDR0s5VTc1QQACSFNZWlBaVzE2NjlVODIxQldZUAADAAQBAAKDgDx7g46YHgOEuBCEhSrCAwQyTlFUSUU3VEFNUDhQUUFITzRQAAWBr6ZKgemsDoGggUeAktXoBgAHAAAFAAYBAQJGUzAwRFhWMFQ5UjFYOUZKNFFFAAN%2BAAAABAAACFZGT0ZZUEZRR1JLVURRVVozRkYACUxWREJKQVc4Qkk1WVNTMVFVQkcACgALAVRMT0tXQ0w5NVJVMDNEOVBFVEcAAjczR1ZURFhERUdFMjdNMlhKTVQAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAAAUFQOEpOTjVMWUFCR1ZNR0tZMUIAAlExTkNXV1BHNkZaMEROSVhKQlEAAwAEAQECgazlhIG8aWIDhMwEVIT%2BMccABQAGAQECRlMwMERYVjBUOVIxWDlGSjRRRQADfgAAAAQAAAhHNEk0R0ZKWEpCOUFUWjNQVFgxAAlMVkRCSkFXOEJJNVlTUzFRVUJHAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAAAAIDAA%3D%3D",
  "Comparative LGN (2018): Human, SMART-seq",
  "https://knowledge.brain-map.org/abcatlas#AQIBSzlKTjIzUDI0S1FDR0s5VTc1QQACSFNZWlBaVzE2NjlVODIxQldZUAADAwE3M0dWVERYREVHRTI3TTJYSk1UAAIAAAFJOUxOUDBPMVJOOEs0U04yR1dZAAIBTEdkAAABVEZRRkxORVAzVjIyMlk4OEM2NAACAAAEAQACg7bFJYOU53MDgL5kWoAMANUEMk5RVElFN1RBTVA4UFFBSE80UAAFga%2BmSoHprA6BoIFHgJLV6AYABwAABQAGAQECRlMwMERYVjBUOVIxWDlGSjRRRQADfgAAAAQAAAhWRk9GWVBGUUdSS1VEUVVaM0ZGAAlMVkRCSkFXOEJJNVlTUzFRVUJHAAoACwFUTE9LV0NMOTVSVTAzRDlQRVRHAAI3M0dWVERYREVHRTI3TTJYSk1UAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAADAAAAUFQOEpOTjVMWUFCR1ZNR0tZMUIAAlExTkNXV1BHNkZaMEROSVhKQlEAAwEBWTkzN0NWVVNWWkM3S1lPSFdWTwACAVRIAAAEAQECgazlhIG8aWIDhJIyUIT%2BMccABQAGAQECRlMwMERYVjBUOVIxWDlGSjRRRQADfgAAAAQAAAhHNEk0R0ZKWEpCOUFUWjNQVFgxAAlMVkRCSkFXOEJJNVlTUzFRVUJHAAoACwFub25lAAJub25lAAMBBAEAAiMwMDAwMDAAA8gBAAUBAQIjMDAwMDAwAAPIAQAADAAAAgMDAAA%3D",
  "Comparative LGN (2018): Macaque, SMART-seq"
)


############################################
## DO NOT EDIT ANYTHING BELOW THIS POINT! ##
############################################

categories = factor(c(as.character(categories)),levels = unique(c(levels(categories))))

# Convert above into a data frame
table_info <- data.frame(table_name   = table_names,
                         table_loc    = table_locations,
                         description  = descriptions,
                         url          = urls
)

# Convert table names into a nested list
table_name <- list()
for (cat in levels(categories)){
  table_name[[cat]] <- c("Select comparison table...", table_names[categories==cat])
}
table_name[["Enter your own location"]] = c("Enter your own location")
category = "Enter your own location"
