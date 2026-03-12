#### Description: Source script holding re-usable chunks of code for analysis
# related spatial distribution of TCRs ####
#### Author: Chloe Tu ####


############################################################################
######## CODE RELATED TO SPATIAL DISTRIBUTION ANALYSIS ----


## Function description: From a vector of strings, normally from object@meta.data$clone_id
get_unique_productive_clones <- function(vec){
  # Collapse all clones from puck into one string
  collapsed <- paste0(vec, collapse = ",")
  
  # Reseparate clones from single string
  separated <- strsplit(collapsed, ",")
  
  # Get all unique clones/CDR3, and remove the clones with an * in them 
  uniq <- unique(separated[[1]])
  uniq <- uniq[uniq != "NA"]
  uniq <- uniq[!str_detect(uniq, fixed("*"))]
  return(uniq)
}

calculate_clonotype_enrichment <- function(others, current, id, given_area){
  # Clonotype enrichment calculation depending on which area to calculate enrichment FOR
  other_in_given_area <- as.integer(others[others$area == given_area, "count"])
  other_in_all <- sum(others$count)
  
  current_in_given_area <- as.integer(current[current$area == given_area, "count"])
  current_in_all <- sum(current$count)
  
  # Handle when there aren't any clones in given area
  if(is.na(current_in_given_area)){current_in_given_area = 0}
  
  # Calculate the enrichment of a given clone in a given area
  enrichment_in_area <- (current_in_given_area/current_in_all) * (other_in_all/other_in_given_area) - 1
  clonotype_enrichment <- c(id, enrichment_in_area, current_in_given_area, current_in_all, other_in_given_area, other_in_all)
  
  return(clonotype_enrichment)
}

#### CODE FOR CLONES FOUND IN SINGLE CELL AND SPATIAL DATA
## Function description: From a spatial Seurat object or dataframe that has been annotated 
# (in ImageScope etc) with specific physical areas, count the number of clones 
# in all annotated areas of the puck. Flexible- works with n+ no. of areas
## Output: a dataframe containing the counts of scTCR and spatial clones in the puck, grouped
# by the annotated areas
## NOTE: STR_DETECT DOESN'T WORK WITH CLONES THAT HAVE A * IN THEM.
## SECOND NOTE: DO NOT USE THIS FOR FISHERS TEST SINCE A BEAD WITH TWO TCRS WILL BE COUNTED TWICE
count_clones_by_annotated_area <- function(dat, clones){
  tcr_area_counts <- as.data.frame(matrix(nrow = 0, ncol = 3))
  colnames(tcr_area_counts) <- c("clone_id", "area", "count")
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  for(clone in clones){
    # Count expansion of a specific clone in all areas
    temp_count <- dat %>%
      mutate(tcr_of_interest = str_detect(clone_id, paste0(clone, "($|,)"))) %>%
      group_by(tcr_of_interest) %>%
      dplyr::count(area) %>%
      filter(tcr_of_interest == TRUE) %>%
      mutate(clone_id = clone, .before = area) %>%
      dplyr::rename("count" = "n") %>%
      ungroup %>%
      select(clone_id, area, count)
    
    tcr_area_counts <- rbind(tcr_area_counts, temp_count)
  }
  
  colnames(tcr_area_counts) <- c("clone_id", "area", "count")
  
  return(tcr_area_counts)
}

## Function: Creates a contingency table for fisher's exact test
create_contingency_table <- function(others, current, given_area){
  # Set focus on a specific area- the annotations of other areas do not matter 
  # anymore, so we add up counts from other areas
  others <- others %>%
    mutate(area = case_when(area == given_area ~ given_area,
                            !(area == given_area) ~ "Other")) %>%
    group_by(area) %>%
    summarize(count = sum(count))  %>%
    mutate(id = "others_id", .before = "area")
  
  current <- current %>%
    mutate(area = case_when(area == given_area ~ given_area,
                            !(area == given_area) ~ "Other")) %>%
    group_by(area) %>%
    summarize(count = sum(count)) %>%
    mutate(id = "current_id", .before = "area")
  
  # Create a table to run Fisher's exact test
  fishers_table <- rbind(others, current)
  
  # Cast to wide table
  fishers_table_wide <- acast(fishers_table, area ~ id, value.var = "count")
  
  # Replace NAs with 0
  fishers_table_wide[is.na(fishers_table_wide)] <- 0
  
  return(fishers_table_wide)
}

## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas, create a 2x2 table suitable for Fisher's exact test,
# then run a two-sided Fisher's exact test
## Example of the table: 
#                Other clones    Specific clone
# Other areas         10            1
# Specific area       2             10
## Output: a dataframe containing the scTCR clone and the associated, newly calculated
# p-value that is NOT FDR adjusted.
spatial_fishers_clones <- function(dat, clones, given_area){
  fishers_pval.list <- list()
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }

  dat <- dat %>%
    filter(nFeature_TCR > 0)

  for(clone in clones){
    # Calculate counts of all other clones, in all annotated areas
    others <- dat %>% 
      filter(!str_detect(clone_id, paste0(clone, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Calculate counts of the given clone, in all annotated areas
    current <- dat %>%
      filter(str_detect(clone_id, paste0(clone, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    contigency_table <- create_contingency_table(others, current, given_area)
    
    # Run Fisher's exact test
    fishers_res <- fisher.test(contigency_table)
    
    # Extract pval
    fishers_pval.list <- append(fishers_pval.list, list(fishers_res[[1]]))
  }
  
  pval.df <- fishers_pval.list %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    cbind(clones) %>%
    dplyr::rename("pval" = "V1",
                  "clone_id" = "clones")
  
  rownames(pval.df) <- NULL
  
  return(pval.df)
}


## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas, create a 2x2 table suitable for Fisher's exact test,
# then run a two-sided Fisher's exact test. Traits currently include the 
# columns "antigen", "virus_antigen", "tumor_reactive", "tumor_specific", 
# "tumor_antigen_specific", "virus_reactive"
## Example of the table: 
#                Number of beads with other trait    Number of beads with specific trait
# Other areas         10            1
# Specific area       2             10
## Output: a dataframe containing the clone and the associated, newly calculated
# p-value that is NOT FDR adjusted.
spatial_fishers_specificity <- function(dat, given_area, 
                                        traits = c("tumor_reactive", "tumor_specific", 
                                                   "tumor_antigen_specific", "virus_reactive")){
  fishers_pval.list <- list()
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Loop over traits 
  for(trait in traits){
    dat.temp <- dat[,c("area", trait)]
    colnames(dat.temp) <- c("area", "specificity")
    # Calculate counts of beads without the given trait, in all annotated areas
    others <- dat.temp %>% 
      filter(is.na(specificity)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Calculate counts of beads with the given trait, in all annotated areas
    current <- dat.temp %>%
      filter(specificity == TRUE) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Create contingency table
    contigency_table <- create_contingency_table(others, current, given_area)
    
    # Run Fisher's exact test
    fishers_res <- fisher.test(contigency_table, alternative = "two.sided")
    
    # Extract pval
    fishers_pval.list <- append(fishers_pval.list, list(fishers_res[[1]]))
  }
  
  pval.df <- fishers_pval.list %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    cbind(traits) %>%
    dplyr::rename("pval" = "V1")
  
  rownames(pval.df) <- NULL
  
  return(pval.df)
}


## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas, create a 2x2 table suitable for Fisher's exact test,
# then run a two-sided Fisher's exact test. "Antigen" and "virus_antigen"
## Example of the table: 
#                Number of beads with other/no antigen    Number of beads with other/no antigen
# Other areas         10            1
# Specific area       2             10
## Output: a dataframe containing the clone and the associated, newly calculated
# p-value that is NOT FDR adjusted.
spatial_fishers_antigen <- function(dat, given_area){
  fishers_pval.list <- list()
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  # Only care about beads with TCR
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Get all tumor/pool antigens
  t_antigens <- get_unique_productive_clones(dat$antigen)
  
  # Get all virus antigens 
  v_antigens <- get_unique_productive_clones(dat$virus_antigen)
  
  # Loop over tumor/pool antigens 
  for(ant in t_antigens){
    # Some antigens have parenthesis in them. Reformat for grep
    if(str_detect(ant, "\\(")){
      ant <- str_replace_all(ant, fixed("("), "\\(")
      ant <- str_replace_all(ant, fixed(")"), "\\)")
    }
    
    # Calculate counts of all other clones, in all annotated areas
    others <- dat %>% 
      filter(!str_detect(antigen, paste0(ant, "($|,)")) | is.na(antigen)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Calculate counts of the given clone, in all annotated areas
    current <- dat %>%
      filter(str_detect(antigen, paste0(ant, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Create contingency table
    contigency_table <- create_contingency_table(others, current, given_area)
    
    # Run Fisher's exact test
    fishers_res <- fisher.test(contigency_table, alternative = "two.sided")
    
    # Extract pval
    fishers_pval.list <- append(fishers_pval.list, list(fishers_res[[1]]))
  }
  
  # Loop over virus antigens 
  for(ant in v_antigens){
    # Calculate counts of all other clones, in all annotated areas
    others <- dat %>% 
      filter(!str_detect(virus_antigen, paste0(ant, "($|,)")) | is.na(virus_antigen)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Calculate counts of the given clone, in all annotated areas
    current <- dat %>%
      filter(str_detect(virus_antigen, paste0(ant, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Create contingency table
    contigency_table <- create_contingency_table(others, current, given_area)
    
    # Run Fisher's exact test
    fishers_res <- fisher.test(contigency_table, alternative = "two.sided")
    
    # Extract pval
    fishers_pval.list <- append(fishers_pval.list, list(fishers_res[[1]]))
  }
  
  pval.df <- fishers_pval.list %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    cbind(c(t_antigens, v_antigens)) %>%
    dplyr::rename("antigens" = "c(t_antigens, v_antigens)",
                  "pval" = "V1")
  
  rownames(pval.df) <- NULL
  
  return(pval.df)
  
}

## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas, create a 2x2 table suitable for Fisher's exact test,
# then run a two-sided Fisher's exact test.
## Example of the table: 
#                Number of beads with other phenotype    Number of beads with specific phenotype
# Other areas         10            1
# Specific area       2             10
## Output: a dataframe containing the clone and the associated, newly calculated
# p-value that is NOT FDR adjusted.
spatial_fishers_phenotype <- function(dat, given_area, phenotypes){
  fishers_pval.list <- list()
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Loop over traits 
  for(phenotype in phenotypes){
    # Calculate counts of beads without the given phenotype, in all annotated areas
    others <- dat %>% 
      filter(!str_detect(PrimaryCluster, paste0(phenotype, "($|,)")) | is.na(PrimaryCluster)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Calculate counts of beads with the given phenotype, in all annotated areas
    current <- dat %>%
      filter(str_detect(PrimaryCluster, paste0(phenotype, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    # Create contingency table
    contigency_table <- create_contingency_table(others, current, given_area)
    
    # Run Fisher's exact test
    fishers_res <- fisher.test(contigency_table, alternative = "two.sided")
    
    # Extract pval
    fishers_pval.list <- append(fishers_pval.list, list(fishers_res[[1]]))
  }
  
  pval.df <- fishers_pval.list %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    cbind(phenotypes) %>%
    dplyr::rename("pval" = "V1")
  
  rownames(pval.df) <- NULL
  
  return(pval.df)
}

## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas and a dataframe of previously calculated p-values,
# calculate the spatial enrichment of clones in a given area, and correct for FDR
# using the BH method
# This formula follows the method in "Spatially mapping T cell receptors and 
# transcriptomes reveals distinct immune niches and interactions underlying 
# the adaptive immune response", where clonotype enrichment in each compartment 
# was calculated in the formula
# Ec = (Cc/Call) * (Aall/Ac) - 1, where 
# Ec = clonotype enrichment for the compartment 
# Cc = number of beads for a clonotype in the compartment 
# Call = total number of beads for a clonotype (across all compartments);
# Ac = number of beads for all other clonotypes in the compartment
# Aall = total number of beads for all other clonotypes (across all compartments),
## Output: a dataframe
spatial_enrichment_clones <- function(dat, clones_to_keep, given_area){
  clone_enrichment <- as.data.frame(matrix(nrow = 0, ncol = 6))
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  # Only care about beads with TCR
  dat <- dat %>%
    filter(!is.na(clone_id))
  
  # Calculate clonotype enrichment in each compartment
  for(clone in clones_to_keep){
    # Calculate counts of a specific clone, and all other clones, in all areas
    others <- dat %>% 
      filter((!str_detect(clone_id, paste0(clone, "($|,)")))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "others", .before = "area")
    
    current <- dat %>%
      filter(str_detect(clone_id, paste0(clone, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "current", .before = "area")
    
    # Calculate clonotype enrichment per sophias paper
    enrichment_row <- calculate_clonotype_enrichment(others, current, clone, given_area)
    
    clone_enrichment <- rbind(clone_enrichment, enrichment_row)
  }
  
  colnames(clone_enrichment) <- c("clone_id", "enrichment_in_area", "current_in_given_area", 
                                  "current_in_all", "other_in_given_area", "other_in_all")
  
  return(clone_enrichment)
}

## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas and a dataframe of previously calculated p-values,
# calculate the spatial enrichment of clones in a given area, and correct for FDR
# using the BH method
# This formula follows the method in "Spatially mapping T cell receptors and 
# transcriptomes reveals distinct immune niches and interactions underlying 
# the adaptive immune response", where clonotype enrichment in each compartment 
# was calculated in the formula
# Ec = (Cc/Call) * (Aall/Ac) - 1, where 
# Ec = clonotype enrichment for the compartment 
# Cc = number of beads for a clonotype in the compartment 
# Call = total number of beads for a clonotype (across all compartments);
# Ac = number of beads for all other clonotypes in the compartment
# Aall = total number of beads for all other clonotypes (across all compartments),
## Output: a dataframe
spatial_enrichment_specificity <- function(dat, given_area, 
                                           traits = c("tumor_reactive", "tumor_specific", 
                                                      "tumor_antigen_specific", "virus_reactive")){
  trait_enrichment <- as.data.frame(matrix(nrow = 0, ncol = 6))
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  # Only care about beads with TCR
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Calculate clonotype enrichment in each compartment
  for(trait in traits){
    # Calculate counts of a specific clone, and all other clones, in all areas
    dat.temp <- dat[,c("area", trait)]
    colnames(dat.temp) <- c("area", "specificity")
    # Calculate counts of beads without the given trait, in all annotated areas
    others <- dat.temp %>% 
      filter(is.na(specificity)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "others", .before = "area")
    
    current <- dat.temp %>%
      filter(!is.na(specificity)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "current", .before = "area")
    
    # calculate clonotype enrichment per sophias paper
    enrichment_row <- calculate_clonotype_enrichment(others, current, trait, given_area)
    
    trait_enrichment <- rbind(trait_enrichment, enrichment_row)
  }
  
  colnames(trait_enrichment) <- c("traits", "enrichment_in_area", "current_in_given_area", 
                                  "current_in_all", "other_in_given_area", "other_in_all")
  return(trait_enrichment)
}


## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas and a dataframe of previously calculated p-values,
# calculate the spatial enrichment of clones in a given area, and correct for FDR
# using the BH method
# This formula follows the method in "Spatially mapping T cell receptors and 
# transcriptomes reveals distinct immune niches and interactions underlying 
# the adaptive immune response", where clonotype enrichment in each compartment 
# was calculated in the formula
# Ec = (Cc/Call) * (Aall/Ac) - 1, where 
# Ec = clonotype enrichment for the compartment 
# Cc = number of beads for a clonotype in the compartment 
# Call = total number of beads for a clonotype (across all compartments);
# Ac = number of beads for all other clonotypes in the compartment
# Aall = total number of beads for all other clonotypes (across all compartments),
## Output: a dataframe
spatial_enrichment_antigen <- function(dat, given_area){
  trait_enrichment <- as.data.frame(matrix(nrow = 0, ncol = 6))
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  # Only care about beads with TCR
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Get all tumor/pool antigens
  t_antigens <- get_unique_productive_clones(dat$antigen)
  
  # Get all virus antigens 
  v_antigens <- get_unique_productive_clones(dat$virus_antigen)

  # Calculate clonotype enrichment in each compartmen. loop over tumor/pool antigens 
  for(ant in t_antigens){
    # Some antigens have parenthesis in them. Reformat for grep
    if(str_detect(ant, fixed("("))){
      ant_temp <- str_replace_all(ant, fixed("("), "\\(")
      ant_temp <- str_replace_all(ant_temp, fixed(")"), "\\)")
    } else {
      ant_temp <- ant
    }
    
    # Calculate counts of all other clones, in all annotated areas
    others <- dat %>% 
      filter(!str_detect(antigen, paste0(ant_temp, "($|,)")) | is.na(antigen)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n")
    
    current <- dat %>%
      filter(str_detect(antigen, paste0(ant_temp, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "current", .before = "area")
    
    # calculate clonotype enrichment per sophias paper
    enrichment_row <- calculate_clonotype_enrichment(others, current, ant, given_area)
    
    trait_enrichment <- rbind(trait_enrichment, enrichment_row)
  }
  
  # Calculate clonotype enrichment in each compartmen. loop over virus antigens 
  for(ant in v_antigens){
    # Calculate counts of beads with a tumor specific antigen, and all other beads, in all areas
    others <- dat %>% 
      filter(!str_detect(virus_antigen, paste0(ant, "($|,)")) | is.na(virus_antigen)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
    mutate(clone_id = "others", .before = "area")
    
    current <- dat %>%
      filter(str_detect(virus_antigen, paste0(ant, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "current", .before = "area")
    
    # calculate clonotype enrichment per sophias paper
    enrichment_row <- calculate_clonotype_enrichment(others, current, ant, given_area)
    
    trait_enrichment <- rbind(trait_enrichment, enrichment_row)
  }
  
  colnames(trait_enrichment) <- c("antigens", "enrichment_in_area", "current_in_given_area", 
                                  "current_in_all", "other_in_given_area", "other_in_all")
  return(trait_enrichment)
}

## Function description: From a spatial Seurat object that has been annotated 
# with specific physical areas and a dataframe of previously calculated p-values,
# calculate the spatial enrichment of beads with a phenotype in a given area, and correct for FDR
# using the BH method
# This formula follows the method in "Spatially mapping T cell receptors and 
# transcriptomes reveals distinct immune niches and interactions underlying 
# the adaptive immune response", where clonotype enrichment in each compartment 
# was calculated in the formula
# Ec = (Cc/Call) * (Aall/Ac) - 1, where 
# Ec = clonotype enrichment for the compartment 
# Cc = number of beads for a clonotype in the compartment 
# Call = total number of beads for a clonotype (across all compartments);
# Ac = number of beads for all other clonotypes in the compartment
# Aall = total number of beads for all other clonotypes (across all compartments),
## Output: a dataframe
spatial_enrichment_phenotypes <- function(dat, phenotypes, given_area){
  clone_enrichment <- as.data.frame(matrix(nrow = 0, ncol = 6))
  
  if(class(dat)[1]== "Seurat"){
    dat <- dat@meta.data
  }
  
  # Only care about beads with TCR
  dat <- dat %>%
    filter(nFeature_TCR > 0)
  
  # Calculate clonotype enrichment in each compartment
  for(phenotype in phenotypes){
    # Calculate counts of a specific clone, and all other clones, in all areas
    others <- dat %>% 
      filter((!str_detect(PrimaryCluster, paste0(phenotype, "($|,)"))) | is.na(PrimaryCluster)) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "others", .before = "area")
    
    current <- dat %>%
      filter(str_detect(PrimaryCluster, paste0(phenotype, "($|,)"))) %>%
      group_by(area) %>%
      dplyr::count() %>%
      dplyr::rename("count" = "n") %>%
      mutate(clone_id = "current", .before = "area")
    
    # Calculate clonotype enrichment per sophias paper
    enrichment_row <- calculate_clonotype_enrichment(others, current, phenotype, given_area)
    
    clone_enrichment <- rbind(clone_enrichment, enrichment_row)
  }
  
  colnames(clone_enrichment) <- c("phenotypes", "enrichment_in_area", "current_in_given_area", 
                                  "current_in_all", "other_in_given_area", "other_in_all")
  
  return(clone_enrichment)
}


## Function description: Calculate the distance of each barcode to a given point
# puck_metadata must have columns for the XY coordinates, the barcode, and clone ids
calculate_distance_of_beads_to_a_point <- function(puck_metadata, given_point, area_name){
  # Calculate euclidean distance between the given point and all beads provided
  distance.df <- nn2(data = given_point, query = puck_metadata[c("x", "y")], k = 1)[[2]]
  # Add barcodes and clone ids back to the dataframe
  distance.df <- cbind(puck_metadata[c("barcode", "clone_id")], distance.df)
  colnames(distance.df)[3] <- paste0("distance_to_", area_name)
  return(distance.df)
}


## Function description: Plot the distribution of distances between given clones
# and the point of interest
plot_distance_distribution_using_clone_id <- function(clone_vector, distance.df, distance_column_name, color_pal){
  for(clone in clone_vector){
    ### Fixed background line- better for comparison between clones
    distance_clone <- distance.df %>%
      mutate(clone = case_when(str_detect(clone_id, paste0(clone, "($|,)")) ~ clone,
                               !(str_detect(clone_id, paste0(clone, "($|,)"))) ~ "others")) %>%
      filter(clone != "others")
    
    plot <- ggplot() +
      geom_density(data = distance.df, aes(x = !!sym(distance_column_name))) +
      geom_density(data = distance_clone, aes(x = !!sym(distance_column_name), color = clone)) +
      ggtitle("Fixed background line") +
      scale_fill_manual(values = color_pal)
    
    print(plot)
    
    plot <- ggplot() +
      geom_histogram(data = distance.df, aes(x = !!sym(distance_column_name))) +
      geom_histogram(data = distance_clone, aes(x = !!sym(distance_column_name), fill = clone))  +
      ggtitle("Fixed background line") +
      scale_fill_manual(values = color_pal)
    
    print(plot)
    # 
    # ### Changing background line- better to compare a particular clone against the distribution of all other clones
    # distance_clone <- distance.df %>%
    #   mutate(clone = case_when(str_detect(clone_id, paste0(clone, "($|,)")) ~ clone,
    #                          !(str_detect(clone_id, paste0(clone, "($|,)"))) ~ "others"))
    # 
    # plot <- ggplot() +
    #   geom_density(data = distance_clone, aes(x = distance_to_tls, color = clone, fill = clone), alpha = 0.4) +
    #   ggtitle("Per-clone background line") 
    #   
    # print(plot)
    
  }
}

## Function description: Plot the distribution of distances between given clones
# and the point of interest, grouped by whether they are seen in scTCR or unique to spatial
plot_distance_distribution_by_uniqueness <- function(distance.df, distance_column_name, color_pal){
  # Mutate ordering column
  distance_by_type <- distance.df %>%
    mutate(clone = case_when(str_detect(clone_id, fixed("p1")) & str_detect(clone_id, "C") ~ "both",
                             str_detect(clone_id, "C") & !str_detect(clone_id, fixed("p1")) ~ "spatial-unique",
                             !str_detect(clone_id, "C") & str_detect(clone_id, fixed("p1")) ~ "scTCR"))
  # Plot density
  p <- ggplot(distance_by_type, aes(x = !!sym(distance_column_name), fill = clone)) +
    geom_density(position = "identity", alpha = 0.4) +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  # Plot histogram
  p <- ggplot(distance_by_type, aes(x = !!sym(distance_column_name), fill = clone)) +
    geom_histogram(position = "identity", alpha = 0.4) +
    scale_fill_manual(values = color_pal)
  
  print(p)
}


# Plot distribution of distance of clones grouped by phenotypes
plot_distance_distribution_by_phenotype <- function(distance.df, clone_metadata, distance_column_name, color_pal){
  # Add phenotype labels to beads with calculated distances
  distance_phenotype <- distance.df %>%
    left_join(clone_metadata[,c("barcode", "PrimaryCluster", "identity", "tumor_mutation", "virus_antigen", "tumor_specific", "tumor_antigen_specific", "virus_specific")], by = "barcode")
  
  # Plot distance grouped by T cell phenotype
  p <- distance_phenotype %>%
    filter(!is.na(PrimaryCluster)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = PrimaryCluster)) +
    geom_density(alpha = 0.4) +
    ggtitle("Distribution of T cell phenotype") +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  p <- distance_phenotype %>%
    filter(!is.na(PrimaryCluster)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = PrimaryCluster)) +
    geom_histogram(position = "identity", alpha = 0.4) +
    ggtitle("Distribution of T cell phenotype") +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  # Plot distance grouped by CD4/CD8 identity
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = identity)) +
    geom_density(alpha = 0.4) +
    ggtitle("Distribution of CD4/CD8 identity") +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = identity)) +
    geom_histogram(position = "identity",alpha = 0.4) +
    ggtitle("Distribution of CD4/CD8 identity") +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  # Plot distance grouped by whether they are neoantigen reactive
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = tumor_antigen_specific)) +
    geom_density(alpha = 0.4) +
    ggtitle("Distribution of tumor antigen specific clones") +
    scale_fill_manual(values = color_pal)
  
  print(p)
  
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = tumor_antigen_specific)) +
    geom_histogram(position = "identity",alpha = 0.4) +
    ggtitle("Distribution of tumor antigen specific clones") +
    scale_fill_manual(values = color_pal)
  print(p)
  
  # Plot distance grouped by whether they are tumor specific
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = tumor_specific)) +
    geom_density(alpha = 0.4) +
    ggtitle("Distribution of tumor specific clones") +
    scale_fill_manual(values = color_pal)
  print(p)
  
  p <- distance_phenotype %>%
    filter(!is.na(identity)) %>%
    ggplot(aes(x = !!sym(distance_column_name), fill = tumor_specific)) +
    geom_histogram(position = "identity",alpha = 0.4) +
    ggtitle("Distribution of tumor specific clones") +
    scale_fill_manual(values = color_pal)
  print(p)
}

## Function description: plot clones given
plot_clones_spatially <- function(rna_obj, clones_to_plot, area_fill_env){
  for(i in c(1:length(clones_to_plot))){
    beads_to_highlight <- rna_obj@meta.data %>%
      filter(str_detect(clone_id, paste0(paste0(clones_to_plot[i], "($|,)"), collapse = "|"))) %>%
      mutate(trim = str_extract(clone_id, paste0(clones_to_plot[i], collapse = "|")))
    
    p <- SpatialDimPlot(rna_obj, group.by = "area", alpha = 0.3, stroke = 0) +
      # Add beads of interest
      geom_point(data = beads_to_highlight, aes(x = plot_x, y = plot_y, color = trim), 
                 size = 1) +
      # Add outline to beads of interest
      geom_point(data = beads_to_highlight, aes(x = plot_x, y = plot_y), 
                 shape = 1, size = 1.5, alpha = 0.8, stroke = 0.5) +
      area_fill +
      ggtitle(paste0("Distribution of ",clones_to_plot[i]))
    
    print(p)
  }
}

# Use calculate_distance_of_beads_to_a_point() to calculate distance_df please.
# The 3rd column needs to refer to the distance calculations!!
plot_clone_distribution <- function(distance_df, clones_to_plot, area_name, adjust = 0.5){
  colnames(distance_df)[3] <- "distance"
  for(i in c(1:length(clones_to_plot))){
    distance_subset <- distance_df %>%
      filter(str_detect(clone_id, paste0(clones_to_plot[i], "($|,)")))
    
    p <- ggplot() +
      geom_density(data = distance_df, aes(x = distance), linetype = "dashed", adjust = adjust) + 
      geom_density(data = distance_subset, aes(x = distance), adjust = adjust) +
      xlab(paste0("Distance to ", area_name)) +
      ggtitle(paste0("Distr. of ", clones_to_plot[i]))
    
    print(p)
  }
}

######## end -- CODE RELATED TO SPATIAL DISTRIBUTION ANALYSIS
