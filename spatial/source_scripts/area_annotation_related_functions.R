#### Description: Source script holding re-usable chunks of code for analysis
# related to area annotations (from ImageScope) ####
#### Author: Chloe Tu ####


############################################################################
######## CODE RELATED TO IMAGESCOPE ANALYSIS ----

## Function description: Need anndata object to convert ImageScope XML to PNG in 
# Mehdi's python script, so convert seurat object to anndata and saves the ad obj.
## NOTES: the functions SaveH5Seurat() and Convert() as provided by mojaveazure 
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html do NOT work
# Instead, convertFormat by sceasy works. convertFormat is only compatible with 
# Seurat v3, so I converted my v5 object to v3. ConvertFormat will error in the 
# function Seurat::GetAssayData() with the warning "assay must be one of "TCR" or 
# "Spatial3", not "RNA"", and since I can't seem to change the default variables
# inside that function, I simply renamed my "SpatialV3" assay to "RNA"

save_seurat_as_anndata <- function(seurat_obj, out_path, puck_name){
  # Convert obj to V3 
  to_anndata <- seurat_obj
  to_anndata[["Spatial3"]] <- as(object = to_anndata[["Spatial"]], Class = "Assay")
  DefaultAssay(to_anndata) <- "Spatial3"
  to_anndata[["Spatial"]] <- NULL
  # Rename my Spatial assay
  # ! `assay` must be one of "TCR" or "Spatial3", not "RNA".
  to_anndata <- RenameAssays(object = to_anndata, Spatial3 = 'RNA')
  
  # Save the object
  sceasy::convertFormat(to_anndata, from="seurat", to="anndata",
                        outFile=paste0(out_path, "/", puck_name, ".h5ad"))
  
}

## Function description: After obtaining PNG versions of the masks, this code will
# convert the masks into area annotations using the package imager to line up pixels
# NOTE: SORRY- i can't figure out why the coordinates change between anndata+Seurat.
# Just play with adjustment factor and make sure things line up. The adjustment is often between 0-200.

add_masks_to_obj <- function(rna_obj, mask_dir, adjust_x, adjust_y){
  files <- list.files(mask_dir)
  
  mask_paths <- paste0(mask_dir, files)
  
  plot_coords <- rna_obj@images$image@coordinates %>%
    mutate(barcode = rownames(.))
  
  for(i in c(1:length(mask_paths))){
    mask <- load.image(mask_paths[i])
    mask.df <- as.data.frame(mask) %>%
      dplyr::rename("new_x" = "y", 
             "new_y" = "x") %>%
      select(new_x, new_y, value) %>%
      dplyr::rename("x" = "new_x",
             "y" = "new_y") %>%
      mutate(x = x + adjust_x,
             y = y + adjust_y,
             value = value*i)
    
    plot_coords <- plot_coords %>%
      mutate(x = as.integer(x),
             y = as.integer(y))
    
    beads <- inner_join(mask.df, plot_coords, by = c("x", "y")) %>%
      select(value, barcode) %>%
      column_to_rownames("barcode") 
    
    colnames(beads)[1] <- paste0("mask_", i)
    
    rna_obj <- AddMetaData(rna_obj, beads)
  }
  
  return(rna_obj)
}


## Function description: Removes any traces of TCRs in the areas annotated as empty
#  Extract beads in "empty" area, update the TCR assay so that those beads have 0 TCRs, and metadata columns so that those beads have 0 TCRs
remove_tcrs_in_empty_area <- function(seurat_obj){
  # Get empty beads
  empty_beads <- WhichCells(seurat_obj, expression = area == "Empty")
  # Get TCR assay matrix
  obj_TCR_assay <- as.matrix(seurat_obj@assays$TCR$counts)
  
  # Remove empty beads from this matrix (columns)
  obj_TCR_assay <- obj_TCR_assay[,!(colnames(obj_TCR_assay) %in% empty_beads)]
  
  # Add back empty beads to this matrix w/ 0's in all rows
  TCR_assay_to_add <- matrix(data = 0, nrow = nrow(obj_TCR_assay), ncol = length(empty_beads))
  colnames(TCR_assay_to_add) <- empty_beads
  rownames(TCR_assay_to_add) <- rownames(obj_TCR_assay)
  
  obj_TCR_assay <- cbind(obj_TCR_assay, TCR_assay_to_add)
  
  # Remove old TCR assay and add new TCR assay
  seurat_obj[["TCR"]] <- NULL
  seurat_obj[["TCR"]] <- CreateAssay5Object(counts = obj_TCR_assay) 
  
  # Update metadata by replacing TCR-related columns with NA for empty beads
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(chains = case_when(barcode %in% empty_beads ~ NA,
                              !(barcode %in% empty_beads) ~ chains),
           identity = case_when(barcode %in% empty_beads ~ NA,
                                !(barcode %in% empty_beads) ~ identity),
           all_CDR3 = case_when(barcode %in% empty_beads ~ NA,
                                !(barcode %in% empty_beads) ~ all_CDR3),
           sctcr_clone = case_when(barcode %in% empty_beads ~ NA,
                                   !(barcode %in% empty_beads) ~ sctcr_clone),
           tumor_specific = case_when(barcode %in% empty_beads ~ NA,
                                      !(barcode %in% empty_beads) ~ tumor_specific),
           tumor_reactive = case_when(barcode %in% empty_beads ~ NA,
                                      !(barcode %in% empty_beads) ~ tumor_reactive),
           virus_reactive = case_when(barcode %in% empty_beads ~ NA,
                                      !(barcode %in% empty_beads) ~ virus_reactive), 
           antigen = case_when(barcode %in% empty_beads ~ NA,
                               !(barcode %in% empty_beads) ~ antigen),
           virus_antigen = case_when(barcode %in% empty_beads ~ NA,
                                     !(barcode %in% empty_beads) ~ virus_antigen),
           tumor_antigen_specific = case_when(barcode %in% empty_beads ~ NA,
                                              !(barcode %in% empty_beads) ~ tumor_antigen_specific),
           PrimaryCluster = case_when(barcode %in% empty_beads ~ NA,
                                      !(barcode %in% empty_beads) ~ PrimaryCluster),
           clone_id = case_when(barcode %in% empty_beads ~ NA,
                                !(barcode %in% empty_beads) ~ clone_id))
  
  # Recalculate log feature and counts
  seurat_obj$log_nFeature_TCR <- log(seurat_obj$nFeature_TCR)
  seurat_obj$log_nCount_TCR <- log(seurat_obj$nCount_TCR)
  
  return(seurat_obj)
}

################## end- CODE TO ADD AREA MASKS VIA IMAGESCOPE -----------------


