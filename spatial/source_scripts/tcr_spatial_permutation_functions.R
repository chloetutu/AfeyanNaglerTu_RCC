# ------------------------------------------------------------------------------
# Assuming that distance.df is an NxM dataframe where N is the number of beads 
# with a clonotype, and that the column "median_dist" in the dataframe represents 
# the distance of each bead to the centroid of a chosen TLS.
# p-value is calculated as a one sided test.
distance_permutation_test <- function(distance.df, n_perm, tcr_name) {
  # Step 1: Get the observed test statistic - the median distance of TCR_X 
  # containing beads from the centroid of TLSX, minus the median distance of 
  # non-TCR_X containing beads from the centroid of TLSX.
  obs_temp.df <- distance.df %>%
    # The order of group_by() is important
    group_by(label) %>%
    # Get the median distance to TLSX of (1) beads with TCR_X, and (2) beads with all other clones
    summarize(median_dist = median(distance_to))
  
  median_dist.obs <- obs_temp.df %>%
    # Calculate the median distance to TLSX of beads with TCR_X minus beads with all other clones
    summarize(median_diff = diff(median_dist)) %>% 
    pull()
  
  # Step 2: Create n replicates of our dataset, and permute the TCR_x label within 
  # each puck and replicate- this handles confounding variable of puck identity
  distance_perm.df <- distance.df %>%
    rep_sample_n(size = nrow(distance.df), replace = FALSE, reps = n_perm) %>%
    group_by(replicate, Puck) %>%
    mutate(perm_label = sample(label, size = n(), replace = FALSE))
  
  # Make sure there's the same number of TCR containing and non TCR containing labels 
  # within each puck and replicate
  original_check <- distance.df %>%
    group_by(Puck) %>%
    count(label) %>%
    group_by(Puck, label) %>%
    # Collapse replicates
    summarise(n = unique(n), .groups = "keep")
  
  perm_check <- distance_perm.df %>%
    group_by(replicate, Puck) %>%
    count(perm_label) %>%
    group_by(Puck, perm_label) %>%
    # Collapse replicates
    summarise(n = unique(n), .groups = "keep")
  
  print(table(perm_check == original_check))
  
  # Step 3: For each permutated replicate, get the test statistic - the median distance of TCR_X 
  # containing beads from the centroid of TLSX, minus the median distance of 
  # non-TCR_X containing beads from the centroid of TLSX.
  perm_temp.df <- distance_perm.df %>% 
    group_by(replicate, perm_label) %>% 
    summarise(perm_median_dist = median(distance_to), .groups = "drop_last")
  
  median_dist.perm <- perm_temp.df %>%
    group_by(replicate) %>%
    # Calculate the median distance to TLSX of beads with TCR_X minus beads with all other clones
    summarize(perm_stat = diff(perm_median_dist))
  
  # Step 4: Compare the permuted test statistics with the observed test statistic
  # If the alternative hypothesis is true, then the median distance of TCR_X 
  # is smaller than the median distance of other clones. So a smaller test statistic
  # means that TCR_X is closer to the TLS than other clones
  median_dist.perm <- median_dist.perm %>% 
    mutate(observed_stat = median_dist.obs,
           verdict = case_when(perm_stat <= observed_stat ~ "Closer to TLS",
                                 perm_stat > observed_stat ~ "Farther from TLS"))
  
  plot <- ggplot(median_dist.perm, aes(x = perm_stat, fill = verdict)) +
    geom_histogram(bins = 35, color = "white")    +  
    scale_fill_manual(values = c("Farther from TLS" = "grey", "Closer to TLS" = "black")) +
    theme_classic() +
    theme(legend.position = "bottom")             + 
    xlab("Test statistic") +
    ylab("Number of permutations") +
    ggtitle(paste0("Test statistic distribution of ", tcr_name))         +
    geom_vline(xintercept = median_dist.obs, color = "blue") 
  
  print(plot)
  
  pval <- sum(median_dist.perm[,"verdict"] == "Closer to TLS", na.rm = TRUE)/n_perm
  
  # Collect output dataframes in a list
  results.list <- list(obs_median_df = obs_temp.df, obs_test_stat = median_dist.obs, 
                    perm_median_df = perm_temp.df, perm_test_stat = median_dist.perm, pval = pval)
  
  return(results.list)
}

# ------------------------------------------------------------------------------
# Get p-val nested in my lists
get_pval <- function(lis){
  return(lis[["pval"]])
}

# ------------------------------------------------------------------------------
# Get enrichment nested in my lists
get_enrichment <- function(lis){
  return(lis[["obs_test_stat"]])
}

# ------------------------------------------------------------------------------
# Get bead count of the clone nested in my lists
get_count <- function(lis){
  return(lis[["current_in_area"]])
}

# ------------------------------------------------------------------------------
# Get median distances nested in my lists
get_dist <- function(lis){
  dist <- lis[["obs_median_df"]] %>%
    pivot_wider(names_from = "label", values_from = "median_dist")
  return(dist)
}

# ------------------------------------------------------------------------------
# Get median distances nested in my lists
get_dist_df <- function(lis){
  dist.df <- lapply(lis, get_dist) %>%
    do.call(rbind, .)
  return(dist.df)
}

# ------------------------------------------------------------------------------
# Get p-val, enrichment and the bead count of the clone in the given area
# nested in my lists and format into a dataframe
get_vals_df <- function(lis){
  pval.df <- lapply(lis, get_pval) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    rownames_to_column("clone") %>%
    dplyr::rename("pval" = "V1")
  enrichment.df <- lapply(lis, get_enrichment) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    rownames_to_column("clone") %>%
    dplyr::rename("enrichment" = "V1")
  count.df <- lapply(lis, get_count) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    rownames_to_column("clone") %>%
    dplyr::rename("count_in_area" = "V1")
  
  df <- pval.df %>%
    full_join(enrichment.df, by = "clone") %>%
    full_join(count.df, by = "clone")
  
  return(df)
}

# ------------------------------------------------------------------------------
# Assuming that md is an NxM dataframe where N is the number of beads 
# with a clonotype, and that the column "area_label" in the dataframe represents 
# whether the bead is in the chosen biological structure, and the column 
# "tcr_label" in the dataframe represents whether the bead contains TCR_X
# p-value is calculated as a two sided test
enrichment_permutation_test <- function(md, n_perm, clonotype){
  # Step 1: Get the observed test statistic, the clonotyp enrichment of TCR_X in the structure
  # Clonotype enrichment = Cc/Call * Aall/Ac - 1
  # where Cc = number of beads for clonotype X in the structure 
  #       Call = total number of beads for clonotype X across all structures
  #       Ac = number of beads for all other clonotypes in the structure
  #       Aall = total number of beads for all other clonotypes across all structures
  # In other words, Enrichment = (fraction of TCR_X in the structure) / (fraction of all other clones in the structure) - 1
  enrichment_obs <- md %>%
    group_by(tcr_label, area_label) %>%
    dplyr::count() %>%
    pivot_wider(names_from = area_label, values_from = "n") %>%
    # 0's are translated to NAs in dplyr::count(), so convert back NAs back to 0.
    replace(is.na(.), 0) %>%
    # Calculate Call and Aall
    mutate(all_areas = sum(area_of_interest, other_area),
    # Calculate the fractions: Cc/Call and Ac/All
           fraction = area_of_interest/all_areas) %>%
    select(tcr_label, fraction) %>%
    # Calculate clonotype enrichment
    pivot_wider(names_from = tcr_label, values_from = fraction) %>%
    mutate(clonotype_enrichment = (TCR_of_interest/other_tcr) - 1) %>%
    pull(clonotype_enrichment)
  
  # Get the number of TCR_x beads in the structure
  num_beads_in_area <- md %>%
    filter(tcr_label == "TCR_of_interest",
           area_label == "area_of_interest") %>%
    dplyr::count() %>%
    pull(n)
  
  # Step 2: Create n replicates of our dataset, and permute the TCR_x label within 
  # each puck and replicate- this handles confounding variable of puck identity
  md_perm <- md %>%
    rep_sample_n(size = nrow(md), replace = FALSE, reps = n_perm) %>%
    group_by(replicate, orig.ident) %>%
    mutate(perm_label = sample(tcr_label, size = n(), replace = FALSE))
  
  # Make sure there's the same number of TCR containing and non TCR containing labels 
  # within each puck and replicate. i.e. There should be two values of n per puck
  original_check <- md %>%
    group_by(orig.ident, tcr_label) %>%
    count(tcr_label) %>%
    summarise(n = unique(n), .groups = "keep")
  
  perm_check <- md_perm %>%
    group_by(replicate, orig.ident, perm_label) %>%
    count(perm_label) %>%
    group_by(orig.ident, perm_label) %>%
    # Collapse replicates
    summarise(n = unique(n))
  
  table(perm_check == original_check)
  
  # Step 3: For each permutated replicate, recalculate the test statistic (clonotype enrichment score)
  enrichment_perm <- md_perm %>% 
    group_by(replicate, perm_label, area_label) %>% 
    dplyr::count() %>%
    pivot_wider(names_from = area_label, values_from = "n") %>%
    # For less highly expanded clonotypes, the number of TCR_X in the structure or outside the structure
    # is 0 in some permutations due to chance. 0's are translated to NA using dplyr::count(), so convert these NAs back to 0.
    replace(is.na(.), 0) %>%
    # Calculate Call and Aall
    mutate(all_areas = sum(area_of_interest, other_area),
           # Calculate enrichment = Cc/Call and Ac/All
           fraction = area_of_interest/all_areas) %>%
    select(replicate, perm_label, fraction) %>%
    # Calculate clonotype enrichment
    pivot_wider(names_from = perm_label, values_from = fraction) %>%
    mutate(perm_clonotype_enrichment = (TCR_of_interest/other_tcr) - 1) %>%
    select(replicate, perm_clonotype_enrichment)
  
  # Step 4: Compare the permuted test statistics with the observed test statistic
  # A larger test statistic means that TCR_X is more enriched structure than other clones
  enrichment_perm <- enrichment_perm %>% 
    mutate(abs_obs_enrichment = abs(enrichment_obs),
           abs_perm_enrichment = abs(perm_clonotype_enrichment),
           verdict = case_when(abs_perm_enrichment >= abs_obs_enrichment ~ "More extreme",
                               abs_perm_enrichment < abs_obs_enrichment ~ "Less extreme"))
  
  # Plot test statistic distribution
  plot <- ggplot(enrichment_perm, aes(x = perm_clonotype_enrichment, fill = verdict)) +
    geom_histogram(bins = 35, color = "white")    +  
    scale_fill_manual(values = c("Less extreme" = "grey", "More extreme" = "black")) +
    theme_classic() +
    theme(legend.position = "bottom")             + 
    xlab("Test statistic") +
    ylab("Number of permutations") +
    ggtitle(paste0("Test statistic distribution of ", clonotype))         +
    geom_vline(xintercept = enrichment_obs, color = "blue") +
    geom_vline(xintercept = -enrichment_obs, color = "blue")
  
  print(plot)
  
  # Calculate pval by summing up the number of permutations in which the test 
  # statistic was more extreme, and then dividing it by the number of permutations made
  pval <- sum(enrichment_perm[,"verdict"] == "More extreme")/n_perm
  pval
  
  # Collect output dataframes in a list
  data_list <- list(obs_test_stat = enrichment_obs, perm_test_stat = enrichment_perm, 
                    pval = pval, current_in_area = num_beads_in_area)
  
  return(data_list)
}


# ------------------------------------------------------------------------------
# To interrogate whether our enrichment analysis is robust, remove a layer of cells
# X pixels away from a given annotated area. Those beads get cast to "Other" which
# is compatible with the current way enrichment_permutation_test is calculated. This
# may not be compatible for other applications.
# The number of neaby beads tested may need to increase if the pixel distance is increased
remove_perimeter_beads_from_given_area <- function(obj, area_to_trim, pixel_distance = 15){
  if (pixel_distance > 15){
    print("WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          5 BEADS MAY NOT BE ENOUGH TO CAPTURE ALL BEADS AT THAT DISTANCE")
  }
  obj_current_area <- obj@meta.data %>%
    filter(area == area_to_trim) %>%
    select(x, y)
  
  obj_other_area <- obj@meta.data %>%
    filter(area != area_to_trim) %>%
    select(x, y)
  
  # For each bead in the area, calculate the distance to the 5 closest beads outside of the area
  # i.e. find 5 beads in the area closest to the perimeter of the area.
  bead_distances <- nn2(data = obj_other_area, query = obj_current_area, k = 5)[[2]]
  bead_distances <- cbind(obj_current_area, bead_distances)
  
  # For the layer of closest beads, find those closer than the given pixel distance
  beads1 <- bead_distances %>%
    select(x, y, `1`) %>%
    filter(`1` < pixel_distance) %>%
    rownames()
  # For the second closest layer of beads, find those closer than the given pixel distance
  beads2 <- bead_distances %>%
    select(x, y, `2`) %>%
    filter(`2` < pixel_distance) %>%
    rownames()
  # etc
  beads3 <- bead_distances %>%
    select(x, y, `3`) %>%
    filter(`3` < pixel_distance) %>%
    rownames()
  # etc
  beads4 <- bead_distances %>%
    select(x, y, `4`) %>%
    filter(`4` < pixel_distance) %>%
    rownames()
  # etc
  beads5 <- bead_distances %>%
    select(x, y, `5`) %>%
    filter(`5` < pixel_distance) %>%
    rownames()
  
  beads_to_trim <- c(beads1, beads2, beads3, beads4, beads5)
  
  # Plot
  print(SpatialDimPlot(subset(obj, cells = beads_to_trim)) + 
    coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000), ratio = 1) + ggtitle("Beads to trim"))
  print(SpatialDimPlot(obj, group.by = "area") + 
      coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("Original data"))
  
  # Remove the layer of beads by casting their area annotation to other
  obj@meta.data <- obj@meta.data %>%
    mutate(area = case_when(barcode %in% beads_to_trim ~ "Other",
                            T ~ area))
  
  print(SpatialDimPlot(obj, group.by = "area") + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("New data"))
  
  cat("Beads in area: ", nrow(obj_current_area), "\nNo. of beads trimmed: ", length(beads_to_trim),
        "\nPct of beads trimmed: ", (length(beads_to_trim)/nrow(obj_current_area))*100)
  
  return(obj)
}

# ------------------------------------------------------------------------------
# To interrogate whether our enrichment analysis is robust, add a layer of cells
# X pixels away from a given annotated area. Beads in "other" areas get cast to the given which
# is compatible with the current way enrichment_permutation_test is calculated. This
# may not be compatible for other applications.
# The number of neaby beads tested may need to increase if the pixel distance is increased
add_perimeter_beads_to_given_area <- function(obj, area_to_add, pixel_distance = 15){
  if (pixel_distance > 15){
    print("WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          5 BEADS MAY NOT BE ENOUGH TO CAPTURE ALL BEADS AT THAT DISTANCE")
  }
  obj_current_area <- obj@meta.data %>%
    filter(area == area_to_add) %>%
    select(x, y)
  
  obj_other_area <- obj@meta.data %>%
    filter(area != area_to_add) %>%
    select(x, y)
  
  # For each bead outside the area, calculate the distance to the 5 closest beads inside the area
  # i.e. find 5 beads outside the area closest to the perimeter of the area.
  bead_distances <- nn2(data = obj_current_area, query = obj_other_area, k = 5)[[2]]
  bead_distances <- cbind(obj_other_area, bead_distances)
  
  # For the layer of closest beads, find those closer than the given pixel distance
  beads1 <- bead_distances %>%
    select(x, y, `1`) %>%
    filter(`1` < pixel_distance) %>%
    rownames()
  # For the second closest layer of beads, find those closer than the given pixel distance
  beads2 <- bead_distances %>%
    select(x, y, `2`) %>%
    filter(`2` < pixel_distance) %>%
    rownames()
  # etc
  beads3 <- bead_distances %>%
    select(x, y, `3`) %>%
    filter(`3` < pixel_distance) %>%
    rownames()
  # etc
  beads4 <- bead_distances %>%
    select(x, y, `4`) %>%
    filter(`4` < pixel_distance) %>%
    rownames()
  # etc
  beads5 <- bead_distances %>%
    select(x, y, `5`) %>%
    filter(`5` < pixel_distance) %>%
    rownames()
  
  beads_to_add <- c(beads1, beads2, beads3, beads4, beads5)
  
  # Plot
  print(SpatialDimPlot(subset(obj, cells = beads_to_add)) + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000), ratio = 1) + ggtitle("Beads to add"))
  print(SpatialDimPlot(obj, group.by = "area") + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("Original data"))
  
  # Add the layer of beads by casting their area annotation to the area
  obj@meta.data <- obj@meta.data %>%
    mutate(area = case_when(barcode %in% beads_to_add ~ area_to_add,
                            T ~ area))
  
  print(SpatialDimPlot(obj, group.by = "area") + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("New data"))
  
  cat("Beads in area: ", nrow(obj_current_area), "\nNo. of beads added: ", length(beads_to_add),
      "\nPct of beads added: ", (length(beads_to_add)/nrow(obj_current_area))*100)
  
  return(obj)
}



# ------------------------------------------------------------------------------
# To interrogate whether our enrichment analysis is robust, remove a layer of cells
# X pixels away in both directions from a given annotated area.
# The number of neaby beads tested may need to increase if the pixel distance is increased
remove_perimeter_beads_completely <- function(obj, area_to_trim, pixel_distance = 15){
  if (pixel_distance > 15){
    print("WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          5 BEADS MAY NOT BE ENOUGH TO CAPTURE ALL BEADS AT THAT DISTANCE")
  }
  obj_current_area <- obj@meta.data %>%
    filter(area == area_to_trim) %>%
    select(x, y)
  
  obj_other_area <- obj@meta.data %>%
    filter(area != area_to_trim) %>%
    select(x, y)
  
  # For each bead in the area, calculate the distance to the 5 closest beads outside of the area
  # i.e. find 5 beads in the area closest to the perimeter of the area.
  bead_distances <- nn2(data = obj_other_area, query = obj_current_area, k = 5)[[2]]
  bead_distances <- cbind(obj_current_area, bead_distances)
  
  # For the layer of closest beads, find those closer than the given pixel distance
  beads1 <- bead_distances %>%
    select(x, y, `1`) %>%
    filter(`1` < pixel_distance) %>%
    rownames()
  # For the second closest layer of beads, find those closer than the given pixel distance
  beads2 <- bead_distances %>%
    select(x, y, `2`) %>%
    filter(`2` < pixel_distance) %>%
    rownames()
  # etc
  beads3 <- bead_distances %>%
    select(x, y, `3`) %>%
    filter(`3` < pixel_distance) %>%
    rownames()
  # etc
  beads4 <- bead_distances %>%
    select(x, y, `4`) %>%
    filter(`4` < pixel_distance) %>%
    rownames()
  # etc
  beads5 <- bead_distances %>%
    select(x, y, `5`) %>%
    filter(`5` < pixel_distance) %>%
    rownames()
  
  inside_beads_to_trim <- c(beads1, beads2, beads3, beads4, beads5)
  
  # For each bead outside the area, calculate the distance to the 5 closest beads inside the area
  # i.e. find 5 beads outside the area closest to the perimeter of the area.
  bead_distances <- nn2(data = obj_current_area, query = obj_other_area, k = 5)[[2]]
  bead_distances <- cbind(obj_other_area, bead_distances)
  
  # For the layer of closest beads, find those closer than the given pixel distance
  beads1 <- bead_distances %>%
    select(x, y, `1`) %>%
    filter(`1` < pixel_distance) %>%
    rownames()
  # For the second closest layer of beads, find those closer than the given pixel distance
  beads2 <- bead_distances %>%
    select(x, y, `2`) %>%
    filter(`2` < pixel_distance) %>%
    rownames()
  # etc
  beads3 <- bead_distances %>%
    select(x, y, `3`) %>%
    filter(`3` < pixel_distance) %>%
    rownames()
  # etc
  beads4 <- bead_distances %>%
    select(x, y, `4`) %>%
    filter(`4` < pixel_distance) %>%
    rownames()
  # etc
  beads5 <- bead_distances %>%
    select(x, y, `5`) %>%
    filter(`5` < pixel_distance) %>%
    rownames()
  
  outside_beads_to_trim <- c(beads1, beads2, beads3, beads4, beads5)
  
  # Merge beads to trim
  beads_to_trim <- c(inside_beads_to_trim, outside_beads_to_trim)
  
  # Plot
  print(SpatialDimPlot(subset(obj, cells = beads_to_trim)) + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000), ratio = 1) + ggtitle("Beads to trim"))
  print(SpatialDimPlot(obj, group.by = "area") + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("Original data"))
  
  # Remove the two layers of beads COMPLETELY
  obj <- subset(obj, subset = barcode %in% beads_to_trim, invert = TRUE)
  
  print(SpatialDimPlot(obj, group.by = "area") + 
          coord_fixed(xlim = c(0, 5000), ylim = c(0, 5000)) + ggtitle("New data"))
  
  cat("Beads in area: ", nrow(obj_current_area), "\nNo. of beads trimmed: ", length(beads_to_trim),
      "\nPct of beads trimmed: ", (length(beads_to_trim)/nrow(obj_current_area))*100)
  
  return(obj)
}
