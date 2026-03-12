# Assign cell types based on RCTD weights
assign_celltypes <- function(seurat_md){
  if("Fibroblast" %in% colnames(seurat_md)){
    seurat_md <- seurat_md %>%
      mutate(assignment1 = case_when(B > 0.3 ~ "B",
                                     `Tcell` > 0.4 ~ "T",
                                     Myeloid > 0.3 ~ "Myeloid",
                                     Tumor > 0.2 ~ "Tumor",
                                     Fibroblast > 0.4 ~ "Fibroblast", 
                                     Endothelial > 0.7 ~ "Endothelial",
                                     NK > 0.3 ~ "NK"))
  } else{
    seurat_md <- seurat_md %>%
      mutate(assignment1 = case_when(B > 0.3 ~ "B",
                                     `Tcell` > 0.4 ~ "T",
                                     Myeloid > 0.3 ~ "Myeloid",
                                     Tumor > 0.2 ~ "Tumor",
                                     Endothelial > 0.7 ~ "Endothelial",
                                     NK > 0.3 ~ "NK"))
  }
  return(seurat_md)
}

order_celltypes <- function(seurat_md){
  if("Fibroblast" %in% colnames(seurat_md)){
    seurat_md <- seurat_md %>%
      mutate(assignment1_ord = case_when(assignment1 == "Endothelial" ~ "1Endothelial",
                                         assignment1 == "Fibroblast" ~ "2Fibroblast",
                                         assignment1 == "Tumor" ~ "3Tumor",
                                         assignment1 == "Myeloid" ~ "4Myeloid",
                                         assignment1 == "NK" ~ "5NK",
                                         assignment1 == "B" ~ "6B",
                                         assignment1 == "T" ~ "7T"))
      
  } else{
    seurat_md <- seurat_md %>%
      mutate(assignment1_ord = case_when(assignment1 == "Endothelial" ~ "1Endothelial",
                                         assignment1 == "Tumor" ~ "2Tumor",
                                         assignment1 == "Myeloid" ~ "3Myeloid",
                                         assignment1 == "NK" ~ "4NK",
                                         assignment1 == "B" ~ "5B",
                                         assignment1 == "T" ~ "6T"))
  }
  return(seurat_md)
}
