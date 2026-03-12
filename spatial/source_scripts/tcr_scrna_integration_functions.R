#### Description: Source script holding re-usable chunks of code to integrate
# TCR data from nanoranger and scRNA data from slideseq ####
#### Author: Chloe Tu ####

############################################################
##################### GENERAL ------------------------------

# Function to read in Slide-seq data
# gene_exp_dir: directory where barcodes.tsv, matrix.tsv, genes.tsv are held
# coords_dir: file where barcode_xy.txt is held
# project_name: desired project name
load_slideseq <- function(gene_exp_dir, coords_dir, project_name = "SlideSeq"){
  # Load gene exp
  gene_exp <- Read10X(gene_exp_dir)
  # Create Seurat obj using constructor
  slideseq_obj <- CreateSeuratObject(gene_exp, project = project_name, assay = "Spatial")
  # Load co-ordinates associated with barcodes
  barcode_xy <- read.table(coords_dir, row.names = 1)
  colnames(barcode_xy) <- c("x", "y")
  # Reorder barcodes to match order of scRNA object
  barcode_xy <- barcode_xy[match(rownames(slideseq_obj@meta.data), rownames(barcode_xy)),]
  # Add coordinates to seurat obj
  slideseq_obj[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = barcode_xy
  )
  # # Remove beads outside of radius
  # slideseq_obj <- FilterSlideSeq(slideseq_obj)
  return(slideseq_obj)
}

###################################### end- GENERAL -------------

############################################################
######## CODE TO INTEGRATE TCR (nanoranger) AND RNA (slideseq) OUTPUT ----

## Function description: helper function to map TCRs that do not have a 1-1 
# barcode match in RNA. Mapping is done using XY locations.
## Input: 
## tcr_unique_barcode must be a nx3 dataframe holding barcodes unique to 
# the TCR object, thus not seen in RNA. The columns must be barcode, x and y.
## rna_barcode must be a nx3 dataframe holding all barcodes in the RNA object.
# The columns must be barcode, x and y.

map_tcr_to_closest_rna <- function(tcr_unique_xy, rna_xy){
  library(RANN)
  tcr_unique <- tcr_unique_xy %>%
    select(TCR_x, TCR_y)
  rna <- rna_xy %>%
    select(RNA_x, RNA_y)
  # Use the nearest neighbor algorithm to obtain the closest RNA barcode 
  # to every TCR barcode
  nearest_rna <- nn2(tcr_unique, data = rna, k = 1)
  reorder_rna <- rna_xy[nearest_rna$nn.idx,]
  mapped <- cbind(tcr_unique_xy, reorder_rna)
  
  # Calculate the distance between the TCR and RNA barcode
  mapped <- mapped %>%
    mutate(delta_x = (RNA_x-TCR_x)^2,
           delta_y = (RNA_y-TCR_y)^2,
           distance = sqrt(delta_x+delta_y))
  
  return(mapped)
}

# Function description: Having RNA data processed via Slideseq, and TCR data
# processed via Nanoranger, assign all TCR barcodes to their closest RNA barcode, 
# using their exact barcode or XY locations.
## Input: 
## tcr_unique_barcode must be a nx3 dataframe holding barcodes unique to 
# the TCR object, thus not seen in RNA. The columns must be barcode, x and y.
## rna_barcode must be a nx3 dataframe holding all barcodes in the RNA object.
# Columns must be TCR_barcode, TCR_x and TCR_y, OR RNA_barcode, RNA_x and RNA_y
map_tcr_and_rna_barcodes <- function(tcr_xy, rna_xy){
  ## Handling barcodes shared between TCR and RNA
  # Format so that both dataframes look the same
  shared_xy <- tcr_xy %>%
    inner_join(rna_xy, by = c("TCR_barcode" = "RNA_barcode")) %>%
    # Calculating distance between barcodes for sanity checks
    mutate(RNA_barcode = TCR_barcode,
           delta_x = (RNA_x-TCR_x)^2,
           delta_y = (RNA_y-TCR_y)^2,
           distance = sqrt(delta_x+delta_y))
  
  # Sanity check: Same barcode in RNA-TCR should be at the same XY location
  
  print(paste0("Do matching barcodes share the same location? ", 
               paste0(shared_xy %>% count(delta_x == 0 & delta_y == 0), collapse = " ")))
  
  ## Handling barcodes seen in TCR, and not seen in RNA
  tcr_unique_xy <- tcr_xy %>%
    filter(!(TCR_barcode %in% rna_xy$RNA_barcode))
  # Finds nearest RNA barcode to every TCR barcode, and calculates distance
  # between barcodes
  mapped_xy <- map_tcr_to_closest_rna(tcr_unique_xy, rna_xy)
  
  # Merges shared barcodes and TCR-unique barcodes together
  final_xy <- rbind(shared_xy, mapped_xy)
  
  # Sanity checks
  print(paste0("Number of rows: ", nrow(final_xy)))
  print(paste0("Number of distinct TCR barcodes: ", final_xy %>% 
                 distinct(TCR_barcode) %>% 
                 nrow()))
  print(paste0("Number of distinct RNA barcodes: ", final_xy %>% 
                 distinct(RNA_barcode) %>% 
                 nrow()))
  return(final_xy)
}

## Function description: Update the TCR count matrix after mapping the TCR 
# barcodes to the correct RNA barcodes. Then add the new TCR count matrix
# to the RNA seurat object.
## Input: 
## tcr_anndata must be the anndata object holding the Nanoranger TCR data
## final_xy must be the output from map_tcr_and_rna_barcodes, specifically a df
# containing all TCR barcodes assigned to its closest RNA barcode
## rna_obj must be a Seurat object holding Slide-seq RNA data
update_and_add_tcr_assay <- function(tcr_anndata, final_xy, rna_obj){
  # Convert sparse counts matrix to more friendly matrix
  tcr_matrix <- as.matrix(tcr_anndata$X) 
  rownames(tcr_matrix) <- paste0(rownames(tcr_matrix), "-1")
  
  # Sanity check: Make sure there's the same TCRs between the xy coordinate DF and counts matrix
  print(paste0("Sanity: Number of same TCRs between xy coordinate DF and counts matrix: ", 
               final_xy %>% filter(TCR_barcode %in% rownames(tcr_matrix)) %>% count()))
  
  # Find all RNA barcodes with a TCR associated to it
  T_rna_barcodes <- unique(final_xy$RNA_barcode)
  
  # Find every RNA barcode without a TCR associated to it
  all_rna_barcodes <- Cells(rna_obj)
  non_T_rna_barcodes <- all_rna_barcodes[!(all_rna_barcodes %in% T_rna_barcodes)]
  
  # Create new TCR matrix for TCR reads
  unified_tcr_matrix <- matrix(nrow = 0, ncol = length(colnames(tcr_matrix)))
  colnames(unified_tcr_matrix) <- colnames(tcr_matrix)
  
  # Create a TCR matrix for barcodes without TCRs associated to it
  empty_tcr_matrix <- matrix(nrow = length(non_T_rna_barcodes), ncol = length(colnames(tcr_matrix)), data = 0)
  colnames(empty_tcr_matrix) <- colnames(tcr_matrix)
  rownames(empty_tcr_matrix) <- non_T_rna_barcodes
  
  # Loop over every RNA barcode with a TCR associated to it
  for(rna_barcode in T_rna_barcodes){
    # Find all TCR barcodes associated with a given RNA barcode
    associated_TCRs <- final_xy %>%
      filter(RNA_barcode == rna_barcode) %>%
      pull(TCR_barcode)
    
    # Pull these TCR barcodes from the TCR matrix
    associated_matrix <- tcr_matrix[rownames(tcr_matrix) %in% associated_TCRs,]
    
    # Sum up counts matrix by columns
    if(!is.null(nrow(associated_matrix))){
      sum_matrix <- t(as.matrix(colSums(associated_matrix)))
    } else {
      sum_matrix <- t(as.matrix(associated_matrix))
    }
    
    # Replace the barcode with the unifying RNA barcode
    rownames(sum_matrix) <- rna_barcode
    
    # Add this row to new TCR matrix
    unified_tcr_matrix <- rbind(unified_tcr_matrix, sum_matrix)
  }
  
  # Sanity check: Make sure there's the same number of counts between the xy coordinate DF and counts matrix
  print(paste0("Sanity: Total TCR counts from original TCR matrix: ", 
               sum(colSums(tcr_anndata$X))))
  print(paste0("Sanity: Total TCR counts from new TCR matrix: ", 
               sum(colSums(unified_tcr_matrix))))
  
  # Append both TCR matricies together- 1. RNA barcodes w/ TCR, 2. RNA barcodes w/o
  unified_tcr_matrix <- rbind(unified_tcr_matrix, empty_tcr_matrix)
  
  rna_obj[["TCR"]] <- CreateAssayObject(counts = t(unified_tcr_matrix))
  
  # Add CDR3 metadata
  
  return(rna_obj)
}

## Function description: Add the associated CDR3 to each bead from the TCR assay inside the seurat object
# in the RNA metadata
add_CDR3s_to_RNA_md <- function(rna_obj){
  cdr3_counts <- as.data.frame(t(as.matrix(rna_obj@assays$TCR@counts)))
  nonzero_cdr3_counts <- cdr3_counts[rowSums(cdr3_counts) != 0,]
  
  # Change non-zero entries into the CDR3 seq itself
  for (i in c(1:nrow(nonzero_cdr3_counts))){
    for (j in c(1:ncol(nonzero_cdr3_counts))){
      if (nonzero_cdr3_counts[i,j] != 0){
        nonzero_cdr3_counts[i,j] <- colnames(nonzero_cdr3_counts)[j]
      }
    }
  }
  
  # Summarize the CDR3 sequences for each bead (barcode)
  nonzero_cdr3_counts[nonzero_cdr3_counts == 0] <- NA
  
  # Create new df to hold the summarized CDR3 for each barcode
  summarized_cdr3 <- as.data.frame(matrix(nrow = nrow(nonzero_cdr3_counts), ncol = 1))
  rownames(summarized_cdr3) <- rownames(nonzero_cdr3_counts)
  colnames(summarized_cdr3) <- "all_CDR3"
  summarized_cdr3 <- summarized_cdr3
  
  # Add all visible CDR3 to the bead
  for (i in c(1:nrow(nonzero_cdr3_counts))){
    # Remove all CDR3 not seen in the bead
    temp_cdr3 <- nonzero_cdr3_counts[i,] %>%
      select_if(~ !any(is.na(.)))
    # Paste the CDR3 together, separated by a comma if there are multiple
    summarized_cdr3[i, 1] <- paste0(temp_cdr3, collapse = ",")
  }
  
  # Sanity check: make sure my summarizing code didn't drop any CDR3s
  number_summarized <- paste0(summarized_cdr3$all_CDR3, collapse = ",") %>%
    # Counting comma's, so +1 for the last CDR3 that doesn't have an associated comma
    str_count(",") + 1
  number_nfeature <- sum(rna_obj$nFeature_TCR)
  
  print(paste0("Sanity check: Summarized CDR3s properly? ", number_summarized == number_nfeature))
  
  # Add to summarized CDR3s metadata
  rna_obj <- AddMetaData(rna_obj, summarized_cdr3)
  
  # Add barcode
  rna_obj@meta.data <- rna_obj@meta.data %>%
    mutate(barcode = rownames(.))
  
  return(rna_obj)
}

## Function description: Get all scTCR clones seen in this specific puck, and
# keep the most expanded scTCR clone across clones that share the same CDR3
get_puck_specific_sctcr_clones <- function(tcr_seq, sctcr_clones){
  # Get all CDR3 sequences seen in the given puck
  tcr_seq <- tcr_seq %>%
    select(chains, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore, CDR3)
  
  # Subset sc clones for the CDR3 seen in the given puck
  puck_specific_sc_clones <- sctcr_clones %>%
    filter(CDR3 %in% tcr_seq$CDR3)
  
  # Join CDR3 and puck specific clones defined in single cell analysis
  merged_tcr_seq <- tcr_seq %>%
    left_join(puck_specific_sc_clones, by = c("chains", "CDR3"))
  
  ## Handling CDR3 that map to multiple scTCR clones
  ## Conversely, it's ok if a scTCR clone has multiple CDR3 sequences
  # Find the CDR3 that map to multiple scTCR clones from the same patient
  cdr3_mapping_multiple_sc_clones <- merged_tcr_seq %>%
    group_by(CDR3) %>%
    distinct(CDR3, sctcr_clone, chains) %>%
    count() %>%
    filter(n > 1)
  
  # Print the multiple scTCR clones that these CDR3 map to
  print("Print the single cell clones with the same CDR3- the most highly expanded scTCR clone for each CDR3 will be kept")
  print(merged_tcr_seq %>% 
          filter(CDR3 %in% cdr3_mapping_multiple_sc_clones$CDR3) %>%
          arrange(CDR3))
  
  # It's most likely that the CDR3 found in the spatial data came from the most expanded scTCR clone, 
  # so find the most expanded scTCR clone (in sc data) for all clones
  top_sc_clones <- merged_tcr_seq %>% 
    filter(CDR3 %in% cdr3_mapping_multiple_sc_clones$CDR3) %>%
    group_by(CDR3) %>%
    distinct(.keep_all = TRUE) %>%
    arrange(desc(total_count)) %>%
    dplyr::slice(1)
  
  # Remove the scTCR clones that less expanded
  sc_clones_to_remove <- merged_tcr_seq %>% 
    filter(CDR3 %in% cdr3_mapping_multiple_sc_clones$CDR3) %>%
    group_by(CDR3) %>%
    distinct(.keep_all = TRUE) %>%
    arrange(desc(total_count)) %>%
    dplyr::slice(-1)
  
  # Remove the sc clones from the list of puck-specific single cell clones, 
  # using the clone id AND the CDR3 sequence. ie, both chains from a clone are
  # seen in the data, but one of those chains belongs to another more expanded clone.
  # we could remove the duplicated chain, but keep the unique chain.
  puck_specific_sc_clones <- puck_specific_sc_clones %>%
    anti_join(sc_clones_to_remove, by = c("sctcr_clone", "CDR3"))
  
  # Join CDR3 and puck specific clones defined in single cell analysis
  puck_specific_sc_clones <- tcr_seq %>%
    left_join(puck_specific_sc_clones, by = c("chains", "CDR3")) %>%
    distinct(.keep_all = TRUE)
  
  return(puck_specific_sc_clones)
}

## Helper function: Pull phenotypes, specificities, chain information,
# CDR3 in the puck for multiple CDR3s
summarize_beads_multiple_cdr3 <- function(puck_specific_sc_clones, rna_obj){
  # Find beads with multiple CDR3s
  beads_with_multiple_cdr3 <- rna_obj@meta.data %>%
    filter(nFeature_TCR > 1)
  
  if(nrow(beads_with_multiple_cdr3) == 0){
    return()
  }
  
  summarized_df <- data.frame()
  for(i in c(1:nrow(beads_with_multiple_cdr3))){
    # Get all CDR3 associated with one bead
    CDR3_i <- str_split_1(beads_with_multiple_cdr3[i,"all_CDR3"], ",")
    
    # If all clones are found in the scTCR data, grab all clones associated with those CDR3
    # Create comma sep list for chain, clone ID, CDR3, primary cluster, identity,
    # tumor mutation and virus antigen
    if(all(CDR3_i %in% puck_specific_sc_clones$CDR3)){
      pasteable <- puck_specific_sc_clones %>%
        filter(CDR3 %in% CDR3_i) %>%
        reframe(chains = paste0(unique(chains), collapse = ","),
                CDR3 = paste0(unique(CDR3), collapse = ","),
                PrimaryCluster = paste0(unique(PrimaryCluster), collapse = ","),
                identity = paste0(unique(identity), collapse = ","),
                sctcr_clone = paste0(unique(sctcr_clone), collapse = ","),
                antigen = paste0(unique(antigen), collapse = ","),
                virus_antigen = paste0(unique(virus_antigen), collapse = ",")) %>%
        distinct()
      
      # Merge one value for tumor antigen specific, tumor specific, virus specific 
      bools <- puck_specific_sc_clones %>%
        filter(CDR3 %in% CDR3_i) %>%
        mutate(tumor_reactive = any(tumor_reactive),
               tumor_specific = any(tumor_specific),
               tumor_antigen_specific = any(tumor_antigen_specific), 
               virus_reactive = any(virus_reactive)) %>%
        distinct(tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive)
      
    } else {
      # If not all of the clones are found in the scTCR data, create dummy row(s) to summarize
      # Get number of/ the CDR3s to create dummies for
      count_dummy <- sum(CDR3_i %in% puck_specific_sc_clones$CDR3 == FALSE)
      unique_CDR3_i <- CDR3_i[!(CDR3_i %in% puck_specific_sc_clones$CDR3)]
      
      dummy_pasteable <- puck_specific_sc_clones %>% 
        select(chains, CDR3, PrimaryCluster, identity, sctcr_clone, tumor_mutation, 
               virus_antigen) %>%
        dplyr::slice(1:count_dummy) %>%
        replace(!is.na(.), NA) %>%
        mutate(CDR3 = unique_CDR3_i)
      
      pasteable <- puck_specific_sc_clones %>%
        filter(CDR3 %in% CDR3_i) %>%
        select(chains, CDR3, PrimaryCluster, identity, sctcr_clone, tumor_mutation, 
               virus_antigen) %>%
        rbind(dummy_pasteable) %>%
        reframe(chains = paste0(unique(chains), collapse = ","),
                CDR3 = paste0(unique(CDR3), collapse = ","),
                PrimaryCluster = paste0(unique(PrimaryCluster), collapse = ","),
                identity = paste0(unique(identity), collapse = ","),
                sctcr_clone = paste0(unique(sctcr_clone), collapse = ","),
                tumor_mutation = paste0(unique(tumor_mutation), collapse = ","),
                virus_antigen = paste0(unique(virus_antigen), collapse = ",")) %>%
        distinct()
      
      # Merge one value for tumor antigen specific, tumor specific, virus specific 
      
      dummy_bools <- puck_specific_sc_clones %>% 
        select(tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive) %>%
        dplyr::slice(1:count_dummy) %>%
        replace(!is.na(.), NA)
      
      bools <- puck_specific_sc_clones %>%
        filter(CDR3 %in% CDR3_i) %>%
        select(tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive) %>%
        rbind(dummy_bools) %>%
        mutate(tumor_reactive = any(tumor_reactive),
               tumor_specific = any(tumor_specific),
               tumor_antigen_specific = any(tumor_antigen_specific), 
               virus_reactive = any(virus_reactive)) %>%
        distinct(tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive)
    }
    # Get barcode
    barcode <- beads_with_multiple_cdr3[i, "barcode"]
    
    # Push summarized data into a new dataframe
    summarized_bead <- cbind(pasteable, bools, barcode)
    summarized_df <- rbind(summarized_df, summarized_bead)
  }
  
  return(summarized_df)
}

## Helper function: Pull phenotypes, specificities, chain information for
# beads in the puck that have 1 CDR3
summarize_beads_one_cdr3 <- function(puck_specific_sc_clones, rna_obj){
  # Find beads with one CDR3
  beads_with_one_cdr3 <- rna_obj@meta.data %>%
    filter(nFeature_TCR == 1)
  
  summarized_df <- data.frame()
  for(i in c(1:nrow(beads_with_one_cdr3))){
    # Get all CDR3 associated with one bead
    CDR3_i <- beads_with_one_cdr3[i, "all_CDR3"]
    
    # If the clone is found in the scRNA data, then grab the clone associated 
    # with that CDR3 and grab additional data
    if(CDR3_i %in% puck_specific_sc_clones$CDR3){
      pasteable <- puck_specific_sc_clones %>%
        filter(CDR3 %in% CDR3_i) %>%
        select(chains, CDR3, PrimaryCluster, identity, sctcr_clone, antigen, 
               virus_antigen, tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive)
    } else { 
      # If the clone is not found in the scRNA data, make a dummy row by keeping 
      # the CDR3 but assigning NA for other values
      pasteable <- puck_specific_sc_clones %>% 
        select(chains, CDR3, PrimaryCluster, identity, sctcr_clone, antigen, 
               virus_antigen, tumor_reactive, tumor_specific, tumor_antigen_specific, virus_reactive) %>%
        dplyr::slice(1) %>%
        replace(!is.na(.), NA) %>%
        mutate(CDR3 = CDR3_i)
    }
    
    # Get barcode
    barcode <- beads_with_one_cdr3[i, "barcode"]
    
    if(nrow(pasteable) > 1){
      stop("wtff")
    }
    
    # Push summarized data into a new dataframe
    summarized_bead <- cbind(pasteable, barcode)
    summarized_df <- rbind(summarized_df, summarized_bead)
  }
  return(summarized_df)
}

## Helper function description: Clean up phenotypes, specificities, chain information
# if all values are NA
clean_up_md <- function(summarized){
  cleaned <- summarized %>%
    mutate(across(c(PrimaryCluster, identity, sctcr_clone, antigen, virus_antigen), ~ case_when((. %in% c("NA", "NA,NA", "NA,NA,NA")) ~ NA,
                                                                                                !(. %in% c("NA", "NA,NA", "NA,NA,NA")) ~ .)))
  
  return(cleaned)
}

## Function description: Attach phenotypes, specificities, chain information,
# CDR3 in the puck
attach_phenotype_specificity_to_RNA_md <- function(puck_specific_sc_clones, rna_obj){
  # Handle clones with multiple CDR3
  summarized_multiple <- summarize_beads_multiple_cdr3(puck_specific_sc_clones, rna_obj)
  # Handle clones with 1 CDR3
  summarized_one <- summarize_beads_one_cdr3(puck_specific_sc_clones, rna_obj)
  # Append together
  summarized <- rbind(summarized_multiple, summarized_one)
  # Sanity check: should be true
  tcr_barcodes <- rna_obj@meta.data %>%
    filter(nFeature_TCR > 0) %>%
    pull(barcode)
  
  print(paste0("Were all RNA barcodes with a TCR assigned phenotype/specificity metadata? ", 
               all(tcr_barcodes %in% summarized$barcode)))
  
  # Create dummy df for barcodes without TCRs
  non_tcr_barcodes <- rna_obj@meta.data %>%
    filter(nFeature_TCR == 0) %>%
    pull(barcode)
  
  non_tcr_summarized <- as.data.frame(matrix(nrow = length(non_tcr_barcodes), ncol = 12))
  colnames(non_tcr_summarized) <- colnames(summarized)
  non_tcr_summarized <- non_tcr_summarized %>%
    mutate(barcode = non_tcr_barcodes)
  
  # Append TCR-containing beads and non-TCR containing beads
  summarized <- rbind(summarized, non_tcr_summarized)
  
  # Cleaning up- remove CDR3 column which had previously been added, and
  # adding back NAs instead of NA,NA etc
  summarized <- summarized %>%
    select(-CDR3)
  summarized <- clean_up_md(summarized)
  
  # Make the barcodes the row names- necessary for correct usage of AddMetaData
  rownames(summarized) <- summarized$barcode
  
  # Append metadata to rna obj
  rna_obj <- AddMetaData(rna_obj, summarized)
  
  return(rna_obj)
}

create_final_clone_ids <- function(rna_obj, sctcr_clones){
  tcr_barcodes <- rna_obj@meta.data %>%
    filter(nFeature_TCR > 0)
  
  puck_all_cdr3 <- str_split_1(paste0(tcr_barcodes$all_CDR3, collapse = ","), pattern = fixed(","))
  puck_all_cdr3 <- puck_all_cdr3[puck_all_cdr3!="NA"]
  
  # spatial unique cdr3 for each barcode
  puck_spatial_unique_cdr3 <- puck_all_cdr3[!(puck_all_cdr3 %in% sctcr_clones$CDR3)]
  
  clone_id.df <- data.frame() 
  for(i in c(1:nrow(tcr_barcodes))){
    # Extract the CDR3 in to a particular barcode
    barcode_n_clones <- str_split_1(tcr_barcodes[i,"all_CDR3"], fixed(","))
    
    clone_id <- c() 
    # Extract the scTCR clones
    if(any(barcode_n_clones %in% sctcr_clones$CDR3)){
      sctcr_clone <- str_split_1(paste0(tcr_barcodes[i, "sctcr_clone"], collapse = ","), pattern = fixed(","))
      sctcr_clone <- sctcr_clone[sctcr_clone!="NA"]
      
      # Push clone ID into vector
      clone_id <- c(clone_id, sctcr_clone)
    }
    
    # Extract the CDR3 which are unique to spatial
    if(any(barcode_n_clones %in% puck_spatial_unique_cdr3)){
      # Find the CDR3 that's unique to spatial
      spatial_cdr3 <- barcode_n_clones[barcode_n_clones %in% puck_spatial_unique_cdr3]
      # Push the CDR3 into vector
      clone_id <- c(clone_id, spatial_cdr3)
    }
    
    # Make unique
    clone_id <- unique(clone_id)
    
    # Paste clones/CDR3 into one line
    clone_id <- paste0(clone_id, collapse = ",")
    
    # Create new row w/ associated barcode
    new_row <- c(tcr_barcodes[i,"barcode"], clone_id)
    
    # Push into DF
    clone_id.df <- rbind(clone_id.df, new_row)
  }
  colnames(clone_id.df) <- c("barcode", "clone_id")
  
  rna_obj@meta.data <- rna_obj@meta.data %>%
    left_join(clone_id.df, by = "barcode")
  
  rownames(rna_obj@meta.data) <- rna_obj@meta.data$barcode
  return(rna_obj)  
}

add_total_clone_count_to_puck <- function(puck, total_clone_count){
  clone_counts.df <- data.frame()
  # Remove non-TCR containing beads
  puck_T <- puck@meta.data %>%
    filter(nFeature_TCR > 0)
  # Beads can have multiple TCRs- assign the largest count
  for(i in c(1:nrow(puck_T))){
    # Get the clones/CDR3 associated with each barcode
    associated_clones <- puck_T[i,"clone_id"] %>%
      str_split_1(",")
    # Find these clones in the total counts df and pick whichever clone is the largest
    associated_counts <- total_clone_count %>%
      filter(clone_id %in% associated_clones) %>%
      filter(total_expansion == max(total_expansion)) %>%
      slice_head(n = 1) %>%
      dplyr::rename("largest_clone_id" = "clone_id") %>%
      select(largest_clone_id, total_expansion) %>%
      mutate(log_total_expansion = log10(total_expansion))
    # Push this barcode, it's associated clone and counts to a dataframe
    new_row <- cbind(associated_counts, barcode = puck_T[i,"barcode"])
    clone_counts.df <- rbind(clone_counts.df, new_row)
  }
  
  # Update metadata of the puck
  puck@meta.data <- puck@meta.data %>%
    full_join(clone_counts.df, by = "barcode") %>%
    mutate(total_expansion = case_when(is.na(total_expansion) ~ 0,
                                       !is.na(total_expansion) ~ total_expansion),
           log_total_expansion = case_when(is.na(log_total_expansion) ~ 0,
                                           !is.na(log_total_expansion) ~ log_total_expansion))
  
  rownames(puck@meta.data) <- puck@meta.data$barcode
  
  return(puck)
}

######## end- CODE TO INTEGRATE TCR (nanoranger) AND RNA (slideseq) OUTPUT ----

