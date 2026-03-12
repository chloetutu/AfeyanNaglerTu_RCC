# Functions

get_beads_near_specific_TCR_phenotype <- function(obj, phenotype = c("Tex", "Tmem"), lower_radius = lower_radius_near, upper_radius = upper_radius_near){
  # Ensure we're taking a phenotype
  phenotype <- match.arg(phenotype)
  
  # Get Tex or Tmem TCRs
  tcrs <- obj@meta.data %>%
    filter(str_detect(PrimaryCluster, regex(phenotype, ignore_case = FALSE))) %>%
    select(x, y)
  
  all_beads <- obj@meta.data %>%
    select(x, y)
  
  # Get all beads of the lower pixel radius around TCRs
  near_beads_l <- nn2(all_beads, tcrs, radius = lower_radius, k = 5000, searchtype = "radius")
  
  # Get all beads of the upper pixel radius around TCRs
  near_beads_u <- nn2(all_beads, tcrs, radius = upper_radius, k = 5000, searchtype = "radius")
  
  # Collapse matrix of row indexes into a vector of unique values
  near_beads_l <- unique(as.vector(near_beads_l$nn.idx))
  near_beads_u <- unique(as.vector(near_beads_u$nn.idx))
  
  # Remove the beads in the lower radius from the upper radius
  if(lower_radius > 0){
    near_beads <- near_beads_u[!(near_beads_u %in% near_beads_l)]
    near_beads_md <- obj@meta.data[near_beads,]
  } else {
    near_beads_md <- obj@meta.data[near_beads_u,]
  }
  
  return(list(TCRs = tcrs, near_beads = near_beads_md, lower_radius = lower_radius, upper_radius = upper_radius))
}

# When Tex and Tmem TCRs are close to each other, beads may be called within range of both types of TCRs. This function helps reassign those beads to the closer TCR type.
reassign_overlapping_beads <- function(obj, tex_beads, tmem_beads){
  
  # Get beads near both exhausted and memory TCRs
  overlapped_beads <- intersect(tex_beads[[2]]$barcode, tmem_beads[[2]]$barcode)
  
  overlapped_xy <- obj@meta.data %>%
    filter(barcode %in% overlapped_beads) %>%
    select(x, y)
  
  # Fork to make code compatible with both Tex vs Tmem and neoantigen vs neoantigen analysis
  # Fork for Tex vs Tmem analysis
  if(!is.null(tex_beads[["lower_radius"]])){ 
    # Ensure lower and upper radius are the same
    if(tex_beads$lower_radius != tmem_beads$lower_radius | tex_beads$upper_radius != tmem_beads$upper_radius){
      print("Radii aren't matching up. Is this on purpose?")
    }
    
    # Calculate the distance of overlapping beads to exhausted and memory TCRs
    overlap_to_tex <- nn2(tex_beads[[1]], overlapped_xy, radius = tex_beads$upper_radius, searchtype = "radius")
    overlap_to_tmem <- nn2(tmem_beads[[1]], overlapped_xy, radius = tmem_beads$upper_radius, searchtype = "radius")
  } else if (!is.null(tex_beads[["radius"]])){ # Fork for neoantigen clone vs neoantigen clone analysis
    # Ensure the radius are the same
    if(tex_beads$radius != tmem_beads$radius){
      print("Radii aren't matching up. Is this on purpose?")
    }
    
    # Calculate the distance of overlapping beads to exhausted and memory TCRs
    overlap_to_tex <- nn2(tex_beads[[1]], overlapped_xy, radius = tex_beads$radius, searchtype = "radius")
    overlap_to_tmem <- nn2(tmem_beads[[1]], overlapped_xy, radius = tmem_beads$radius, searchtype = "radius")
    
  } else {print("Where's the radius?")}
  
  
  # For each overlapping bead, see whether it's closer to the exhausted or memory TCR and push barcode name into corresponding vector
  closer_to_tex <- c()
  closer_to_tmem <- c()
  for(i in c(1:length(overlapped_beads))){
    if(overlap_to_tex$nn.dists[[i]] <= overlap_to_tmem$nn.dists[[i]]){
      closer_to_tex <- c(closer_to_tex, rownames(overlapped_xy[i,]))
    } else {
      closer_to_tmem <- c(closer_to_tmem, rownames(overlapped_xy[i,]))
    }
  }
  
  # Remove beads closer to Tex from original Tmem object, and vice versa
  tex_beads$near_beads <- tex_beads$near_beads %>%
    filter(!(barcode %in% closer_to_tmem))
  tmem_beads$near_beads <- tmem_beads$near_beads %>%
    filter(!(barcode %in% closer_to_tex))
  
  return(list(reassigned_tex_beads = tex_beads, reassigned_tmem_beads = tmem_beads))
}

# Get beads near a specific neoantigen TCR inside or outside of the TLS
get_beads_near_specific_neoantigen_TCR <- function(obj, tls_stat = c("inTLS", "outTLS"), neoant = c("p108_271", "p108_2", "p108_215", "p108_266", "p108_1", "p108_208", "p108_205"), radius = rad){
  # Ensure we're taking a neoantigen TCR
  neoant <- match.arg(neoant)
  # Ensure we're taking neoantigen TCRs inside or outside the TLS
  tls_stat <- match.arg(tls_stat)
  
  # Get specific neoantigen TCRs either inside or outside the TLS
  neoant_tcrs <- obj@meta.data %>%
    filter(tls_status == tls_stat,
           str_detect(neoantigen_id, paste0(neoant, "($|,)"))) %>%
    select(x, y)
  
  if(nrow(neoant_tcrs) > 0){
    # Get beads near the neoantigen TCR
    near_beads <- nn2(obj@meta.data[,c("x", "y")], neoant_tcrs[,c("x", "y")], radius = radius, k = 10, searchtype = "radius")
    
    # Collapse matrix of row indexes into a vector of unique values
    near_beads <- unique(as.vector(near_beads$nn.idx))
    
    # Collect metadata of nearby beads
    near_beads_md <- obj@meta.data[near_beads,]
  } else{
    near_beads_md <- neoant_tcrs
  }
  
  return(list(TCRs = neoant_tcrs, near_beads = near_beads_md, radius = radius))
}

# Get beads near a specific neoantigen TCR inside or outside of the TLS
get_beads_near_specific_neoantigen_TCR2 <- function(obj, tls_stat = c("inTLS", "outTLS"), neoant = c("p108_271", "p108_2", "p108_215", "p108_266", "p108_1", "p108_208", "p108_205"), radius = rad){
  # Ensure we're taking a neoantigen TCR
  neoant <- match.arg(neoant)
  # Ensure we're taking neoantigen TCRs inside or outside the TLS
  tls_stat <- match.arg(tls_stat)
  
  # Get specific neoantigen TCRs either inside or outside the TLS
  neoant_tcrs <- obj@meta.data %>%
    filter(tls_status == tls_stat,
           str_detect(neoantigen_id, paste0(neoant, "($|,)"))) %>%
    select(x, y)
  
  if(nrow(neoant_tcrs) > 0){
    near_beads_md <- data.frame()
    for(i in c(1:nrow(neoant_tcrs))){
      # Get beads near the neoantigen TCR
      near_beads <- nn2(obj@meta.data[,c("x", "y")], neoant_tcrs[i, c("x", "y")], radius = radius, k = 10, searchtype = "radius")
      # near_beads <- nn2(obj@meta.data[,c("x", "y")], neoant_tcrs["AAGGGGCGAGCTTA-1", c("x", "y")], radius = radius, k = 10, searchtype = "radius")
      
      # Collapse matrix of row indexes into a vector of unique values
      near_beads <- unique(as.vector(near_beads$nn.idx))
      
      # Collect metadata of nearby beads
      near_beads_md_i <- obj@meta.data[near_beads,] %>%
        mutate(partition = rownames(neoant_tcrs[i,]))
      
      near_beads_md <- rbind(near_beads_md, near_beads_md_i)
    }
    # Rarely a single bead is in the radius of two neoantigens. Remove the second instance of that bead for simplicity
    near_beads_md <- near_beads_md %>%
      distinct(barcode, .keep_all = TRUE)
  } else{
    near_beads_md <- neoant_tcrs
  }
  

  
  return(list(TCRs = neoant_tcrs, near_beads = near_beads_md, radius = radius))
}

get_beads_near_specific_TCR_reactivity <- function(obj, reactivity = c("TS", "VR"), lower_radius = lower_radius_near, upper_radius = upper_radius_near){
  # Ensure we're taking a reactivity
  reactivity <- match.arg(reactivity)

  # Get TS or VR TCRs
  if(reactivity == "TS"){
    tcrs <- obj@meta.data %>%
      filter(tumor_specific == TRUE) %>%
      select(x, y)
  } else if(reactivity == "VR"){
    tcrs <- obj@meta.data %>%
      filter(virus_reactive == TRUE) %>%
      select(x, y)
  }

  all_beads <- obj@meta.data %>%
    select(x, y)

  # Get all beads of the lower pixel radius around TCRs
  near_beads_l <- nn2(all_beads, tcrs, radius = lower_radius, k = 1000, searchtype = "radius")

  # Get all beads of the upper pixel radius around TCRs
  near_beads_u <- nn2(all_beads, tcrs, radius = upper_radius, k = 1000, searchtype = "radius")

  # Collapse matrix of row indexes into a vector of unique values
  near_beads_l <- unique(as.vector(near_beads_l$nn.idx))
  near_beads_u <- unique(as.vector(near_beads_u$nn.idx))

  # Remove the beads in the lower radius from the upper radius, unless the radius is 0
  if(lower_radius > 0){
    near_beads <- near_beads_u[!(near_beads_u %in% near_beads_l)]
    near_beads_md <- obj@meta.data[near_beads,]
  } else {
    near_beads_md <- obj@meta.data[near_beads_u,]
  }

  return(list(TCRs = tcrs, near_beads = near_beads_md, lower_radius = lower_radius, upper_radius = upper_radius))
}

# Collapse outputs from get_beads_near_specific_TCR
collapse_list <- function(clone_list, name){
  # Extract first element from each sublist
  tcr_beads <- lapply(clone_list, "[[", 1)
  # Extract second element from each sublist
  near_beads <- lapply(clone_list, "[[", 2)
  
  # Combine the dataframes using rbind
  tcr_beads <- do.call(rbind, tcr_beads)
  near_beads <- do.call(rbind, near_beads)
  
  # Remove duplicated rows
  tcr_beads <- tcr_beads %>% distinct(.keep_all = TRUE)
  near_beads <- near_beads %>% distinct(.keep_all = TRUE)
  
  # Create new list
  new_list <- list(list(TCRs = tcr_beads, 
                   near_beads = near_beads, 
                   lower_radius = clone_list[[1]]$lower_radius, 
                   upper_radius = clone_list[[1]]$upper_radius))

  names(new_list) <- eval(parse(text = "name"))
  
  return(new_list)
}


get_beads_near_specific_TCR <- function(obj, cln, lower_radius = lower_radius_near, upper_radius = upper_radius_near){
  # Get specific TS or VR TCRs
  tcrs <- obj@meta.data %>%
    filter(str_detect(clone_id, paste0(cln, "($|,)"))) %>%
    select(x, y)

  all_beads <- obj@meta.data %>%
    select(x, y)

  if(nrow(tcrs) > 0){
    # Get all beads of the lower pixel radius around TCRs
    near_beads_l <- nn2(all_beads, tcrs, radius = lower_radius, k = 1000, searchtype = "radius")
  
    # Get all beads of the upper pixel radius around TCRs
    near_beads_u <- nn2(all_beads, tcrs, radius = upper_radius, k = 1000, searchtype = "radius")
  
    # Collapse matrix of row indexes into a vector of unique values
    near_beads_l <- unique(as.vector(near_beads_l$nn.idx))
    near_beads_u <- unique(as.vector(near_beads_u$nn.idx))
  
    # Remove the beads in the lower radius from the upper radius, unless the radius is 0
    if(lower_radius > 0){
      near_beads <- near_beads_u[!(near_beads_u %in% near_beads_l)]
      near_beads_md <- obj@meta.data[near_beads,]
    } else {
      near_beads_md <- obj@meta.data[near_beads_u,]
    }
  } else {
    near_beads_md <- tcrs
  }
  
  return(list(TCRs = tcrs, near_beads = near_beads_md, lower_radius = lower_radius, upper_radius = upper_radius))
}

# Downsample beads based on the column "type"
downsample_beads <- function(obj, downsample){
  ds <- obj@meta.data %>%
    select(barcode, type) %>%
    group_by(type) %>%
    slice_sample(n = downsample, replace = FALSE) %>%
    pull(barcode)
  
  new_obj <- subset(x = obj, cells = ds)
  
  return(new_obj)
}