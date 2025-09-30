# ---- Functions for reading and handling data --------------------------

read_metadata <- function(path, sort_table = FALSE){
  md <- read.csv(path, row.names = 1)
  if(isTRUE(sort_table)){
    md <- md[order(row.names(md)), ] # sort my row names (sample names)
  }
  return(md)
}

# Read Hitchhikers guide style export feature table
read_ft <- function(path, sort_by_names = FALSE, p_sep = ","){
  ft <- read.csv2(path, header = TRUE, row.names = 1, sep = p_sep, dec = ".") #read csv table
  if(isTRUE(sort_by_names)){
    ft <- ft[order(row.names(ft)), ] # sort my row names (sample names)
  }
  
  rownames(ft) <- gsub("\\.mzML$", "", rownames(ft))
  #col_names <- colnames(ft)
  #ft <- sapply(ft, as.numeric)
  #ft <- as.data.frame(ft)
  #colnames(ft) <- col_names
  return(t(ft))
}

remove_feature_by_prefix <- function(df, patterns) {
  # Create a single regex pattern that matches any of the species names at the start
  combined_pattern <- paste0("^(", paste(patterns, collapse = "|"), ")")
  
  # Filter the dataframe: keep rows that do NOT match the pattern
  df_filtered <- df[!grepl(combined_pattern, rownames(df)), ]
  
  return(df_filtered)
}

filter_features_by_col_counts <- function(feature_table, min_count, col_number){
  if (ncol(feature_table) > 1) {
    return(feature_table[which(rowSums(feature_table >= min_count) >= col_number), ])
  }
  else if(ncol(feature_table) == 1){
    ft <- feature_table[feature_table >= min_count, ,drop=FALSE]
    return(ft)
  }
  else{
    print("Dataframe has no columns")
  }
}

# Sort otu table in barcodes numeration
sort_nanopore_table_by_barcodes <- function(df, new_names = NULL){
  cn <- colnames(df) # store column names
  sorted_names <- cn[order(nchar(cn), cn)] # order columns names
  df_sorted <- df[, sorted_names] # order data frame using colnames
  if (!is.null(new_names) && ncol(df_sorted) == length(new_names)) {
    colnames(df_sorted) <- new_names
  }
  return(df_sorted)
}

# Takes an feature table (OTUs) and removes the strain information from the species NOT in the passed vector of species
merge_non_target_strains <- function(df, target_species) {
  # Extract species names (first two words) from rownames
  species_names <- sapply(strsplit(rownames(df), " "), function(x) paste(x[1:2], collapse = " "))
  #print(species_names)
  # Identify which rows belong to target or non-target species
  is_target <- species_names %in% target_species
  #print(is_target)
  # Separate target and non-target
  target_df <- df[is_target, , drop = FALSE]
  #print(target_df)
  non_target_df <- df[!is_target, , drop = FALSE]
  #print(non_target_df)
  non_target_species <- species_names[!is_target]
  #print(non_target_species)
  # 
  # # Aggregate non-target strains by species
  if (nrow(non_target_df) > 0) {
    aggregated <- aggregate(non_target_df, by = list(Species = non_target_species), FUN = sum)
    # Set the species name as rownames and remove the Group column
    # >>> THIS IS WHERE YOU ADD " 1" TO THE SPECIES NAMES <<<
    rownames(aggregated) <- paste(aggregated$Species, "1")
    #rownames(aggregated) <- aggregated$Species
    aggregated$Species <- NULL
  } else {
    aggregated <- NULL
  }
  print(aggregated)
  # Combine target and aggregated non-target dataframes
  result <- rbind(target_df, aggregated)
  
  return(result)
}

##### Function to convert a OTU table to a strain-level-table
# It takes the otu table at species level and a second dataframe including strain-level data.
# First dataframe should be a dataframe containing species level data.
# Second dataframe should be a dataframe containing
# the strain inoculation data in the following format:
merge_abundance_by_strain <- function(df1, df2) {
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  # Extract species names from df1
  species_names_df1 <- rownames(df1)
  #print(species_names_df1)
  
  # Extract species names and strain names from df2
  strain_names_df2 <- df2[, 1]  # Full strain names (including strain number)
  #print(strain_names)
  
  # Create an empty matrix to store the new abundance data
  new_abundance_matrix <- matrix(0, nrow = nrow(df2), ncol = ncol(df1))
  rownames(new_abundance_matrix) <- strain_names_df2
  colnames(new_abundance_matrix) <- colnames(df1)
  
  #print(head(new_abundance_matrix))
  #print(nrow(new_abundance_matrix))
  #print(ncol(new_abundance_matrix))
  
  samples <- colnames(new_abundance_matrix)
  # Iterate over each sample of the new DF
  for (i in seq_along(samples)) {
    #print(samples[i])
    # Get the SC to which it belongs to.
    current_sc <- strsplit(x = samples[i], split = "_")[[1]][1]
    print(current_sc)
    # get the list of strains inoculated in that sample.
    inoc_strains_per_sample <- get_inoculated_strains(df2 = df2, sample_name = current_sc)
    print(inoc_strains_per_sample)
    
    for (x in seq_along(inoc_strains_per_sample)) {
      strain_name <- inoc_strains_per_sample[x]
      print(strain_name)
      # get the index where the data is going to be inserted. The index is the same row as in the df2
      index_strain_df2 <- which(strain_names_df2 == strain_name) # this is also the same in the new df
      #print(index_strain_df2)
      # get the name of the species.
      species_name <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", strain_name)  # Remove strain number, keeping species
      print(species_name)
      if (species_name %in% species_names_df1) {
        index_species_df1 <- which(species_names_df1 == species_name)
        #print(index_species_df1)
        #print(species_names_df1[index_species_df1])
        # get the actual data, that corresponds to the species in df1
        current_abundance <- df1[index_species_df1, i]
        #print(current_abundance)
        # paste the data
        new_abundance_matrix[index_strain_df2, i] <- current_abundance
      }
    }
  }
  return(as.data.frame(new_abundance_matrix))
}

get_inoculated_strains <- function(df2, sample_name) {
  # Select the column corresponding to the sample
  sample_column <- df2[[sample_name]]
  
  # Get row indices where the value is 1 (inoculated strains)
  inoculated_indices <- which(sample_column == 1)
  
  # Extract the strain names based on the indices
  inoculated_strains <- df2[inoculated_indices, 1]  # First column contains strain names
  
  return(inoculated_strains)
}

# Function to set selected species/sample combinations to zero
zero_out_species_in_samples <- function(df, species_name, sample_names) {
  # Safety check: does the species exist?
  if (!(species_name %in% rownames(df))) {
    stop(paste("Species", species_name, "not found in rownames"))
  }
  
  # Safety check: do all samples exist?
  if (!all(sample_names %in% colnames(df))) {
    missing_samples <- sample_names[!sample_names %in% colnames(df)]
    stop(paste("Samples not found in dataframe:", paste(missing_samples, collapse = ", ")))
  }
  
  # Set the selected cells to zero
  df[species_name, sample_names] <- 0
  
  return(df)
}

# ---- Cluster Barplots --------------------------

cluster_barplot_panels <- function(abundance_df, k = NULL, colour_palette = NULL) {
  require(cluster)
  require(ggplot2)
  require(reshape2)
  
  # Convert to matrix
  mat <- as.matrix(abundance_df)
  
  # Relative abundance
  mat_rel <- sweep(mat, 2, colSums(mat), "/")
  
  #print(mat_rel)
  
  mat_input <- mat_rel
  
  # Distance and clustering
  dist_mat <- dist(t(mat_input), method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Number of clusters
  if (is.null(k)) {
    sil_widths <- sapply(2:10, function(k_try) {
      pam_fit <- pam(dist_mat, k_try)
      pam_fit$silinfo$avg.width
    })
    best_k <- which.max(sil_widths) + 1
  } else {
    best_k <- k
  }
  
  # Assign clusters
  pam_fit <- pam(dist_mat, best_k)
  clusters <- pam_fit$clustering
  cluster_df <- data.frame(Sample = names(clusters), Cluster = as.factor(clusters))
  
  # Reorder samples by dendrogram order
  sample_order <- hc$labels[hc$order]
  mat_rel_ordered <- mat_rel[, sample_order]
  
  # Melt for ggplot
  df_long <- melt(mat_rel_ordered)
  colnames(df_long) <- c("Bacteria", "Sample", "Abundance")
  df_long <- merge(df_long, cluster_df, by = "Sample")
  
  # Barplot
  barplot <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Bacteria)) +
    geom_bar(stat = "identity") +
    facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab("Relative Abundance") +
    ggtitle(paste("Stacked Barplot with Clusters (k =", best_k, ")"))
  
  if (!is.null(colour_palette)) {
    # Expecting a named vector: names must match rownames(abundance_df) (Bacteria)
    barplot <- barplot + scale_fill_manual(values = colour_palette, drop = FALSE)
  }
  
  print(barplot)
  
  return(list(
    clusters = cluster_df,
    best_k = best_k
  ))
}

#' Add a cluster column to a metadata dataframe by mapping keys across dataframes
#'
#' @param meta_df data.frame containing the destination column to receive clusters
#' @param clusters_df data.frame with the key->cluster mapping
#' @param meta_key_col character; column in meta_df used to look up the cluster (e.g., "ATTRIBUTE_SynCom")
#' @param cluster_key_col character; key column in clusters_df (e.g., "Sample")
#' @param cluster_value_col character; value column in clusters_df to copy over (e.g., "Cluster")
#' @param new_col_name character; name of the new column to create in meta_df (default "ATTRIBUTE_Cluster")
#' @param warn_missing logical; warn if some meta keys aren’t found in clusters_df (default TRUE)
#' @return meta_df with an added column `new_col_name`
add_cluster_column <- function(meta_df,
                               clusters_df,
                               meta_key_col,
                               cluster_key_col,
                               cluster_value_col,
                               new_col_name,
                               warn_missing = TRUE) {
  # Basic checks
  for (nm in c(meta_key_col)) {
    if (!nm %in% names(meta_df)) stop(sprintf("Column '%s' not found in meta_df.", nm))
  }
  for (nm in c(cluster_key_col, cluster_value_col)) {
    if (!nm %in% names(clusters_df)) stop(sprintf("Column '%s' not found in clusters_df.", nm))
  }
  
  # Optional: detect duplicate keys in the cluster table
  dups <- duplicated(clusters_df[[cluster_key_col]])
  if (any(dups)) {
    dup_vals <- unique(clusters_df[[cluster_key_col]][dups])
    warning(sprintf(
      "Duplicate keys in clusters_df$%s for: %s. Using the last occurrence.",
      cluster_key_col, paste(head(dup_vals, 10), collapse = ", ")
    ))
  }
  
  # Build a named lookup vector: names = keys, values = clusters
  lookup <- setNames(clusters_df[[cluster_value_col]],
                     as.character(clusters_df[[cluster_key_col]]))
  
  # Map onto meta_df using the key column
  keys <- as.character(meta_df[[meta_key_col]])
  mapped <- unname(lookup[keys])
  
  # Warn if some keys were not found
  if (warn_missing) {
    missing_idx <- which(is.na(mapped) & !is.na(keys))
    if (length(missing_idx) > 0) {
      missing_vals <- unique(keys[missing_idx])
      warning(sprintf(
        "No cluster found for %d key(s) in meta_df$%s. Examples: %s",
        length(missing_vals), meta_key_col,
        paste(head(missing_vals, 10), collapse = ", ")
      ))
    }
  }
  
  # Add the new column
  meta_df[[new_col_name]] <- mapped
  meta_df
}


get_palette <- function(nColors = 60, replace_cols = FALSE){
  colors_vec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2",
                  "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3",
                  "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
                  "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
                  "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
                  "lavenderblush4", "lawngreen", "indianred1", "lightblue1", "honeydew4",
                  "hotpink", "#e3ae78", "#a23f3f", "#290f76", "#ce7e00",
                  "#386857", "#738564", "#e89d56", "#cd541d", "#1a3a46",
                  "#9C4A1A", "#ffe599", "#583E26", "#A78B71", "#F7C815",
                  "#EC9704", "#4B1E19", "firebrick2", "#C8D2D1", "#14471E",
                  "#6279B8", "#DA6A00", "#C0587E", "#FC8B5E", "#FEF4C0",
                  "#EA592A", "khaki3", "lavenderblush3", "indianred4", "lightblue",
                  "honeydew1", "hotpink4", "ivory3", "#49516F", "#502F4C",
                  "#A8C686", "#669BBC", "#29335C", "#E4572E", "#F3A712",
                  "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD")
  
  #set.seed(1)
  
  return(colors_vec[sample(1:length(colors_vec), size = nColors, replace = replace_cols)])
}


align_samples_attr <- function(metab_df, metadata_df, sample_col = NULL) {
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df
  
  # Figure out where sample IDs live in metadata
  if (is.null(sample_col)) {
    if ("Sample" %in% colnames(md)) {
      sample_col <- "Sample"
    } else if (!is.null(rownames(md))) {
      md$Sample <- rownames(md)
      sample_col <- "Sample"
    } else stop("metadata_df must have a 'Sample' column or rownames = sample IDs.")
  }
  
  # Ensure required columns exist; fill from IDs if missing
  req <- c("ATTRIBUTE_Cluster","ATTRIBUTE_Time","ATTRIBUTE_SynCom")
  missing <- setdiff(req, colnames(md))
  if (length(missing) > 0) {
    if ("ATTRIBUTE_Cluster" %in% missing) md$ATTRIBUTE_Cluster <- NA
    if ("ATTRIBUTE_Time" %in% missing) {
      md$ATTRIBUTE_Time <- suppressWarnings(
        as.integer(sub(".*_T(\\d+)_.*", "\\1", md[[sample_col]]))
      )
    }
    if ("ATTRIBUTE_SynCom" %in% missing) {
      md$ATTRIBUTE_SynCom <- sub("_.*", "", md[[sample_col]])  # "SC##"
    }
  }
  #print(md$ATTRIBUTE_Time)
  # Intersect & align
  common <- intersect(colnames(metab_df), md[[sample_col]])
  if (length(common) == 0) stop("No overlapping sample IDs between metabolome and metadata.")
  X  <- as.matrix(metab_df[, common, drop = FALSE])
  md <- md[match(common, md[[sample_col]]), , drop = FALSE]
  rownames(md) <- md[[sample_col]]
  
  # Types
  md$ATTRIBUTE_Cluster <- factor(md$ATTRIBUTE_Cluster)
  #md$ATTRIBUTE_Time    <- suppressWarnings(as.integer(md$ATTRIBUTE_Time))
  md$ATTRIBUTE_Time    <- factor(md$ATTRIBUTE_Time)
  md$ATTRIBUTE_SynCom  <- as.character(md$ATTRIBUTE_SynCom)
  
  print("Data alligned succesfully")
  
  list(X = X, meta = md)
}

pcoa_flex <- function(
    metab_df, metadata_df,
    sample_id_col = NULL,              # "Sample" if IDs live in a column; NULL -> use rownames(metadata_df)
    color_var,                          # metadata column for point color (factor)
    shape_var = NULL,                   # metadata column for point shape (factor) or NULL (Optional)
    ellipse_var = NULL,                 # metadata column for ellipses (factor) or NULL (Optional)
    distance = c("bray","euclidean"),
    preprocess = c("none","hellinger","sqrt"),
    k_axes = 2,
    ellipse = TRUE,
    min_n_for_ellipse = 3,
    label_points = FALSE,
    points_palette  = NULL,             # named vector for levels(color_var) (optional)
    ellipse_palette = NULL,             # named vector for levels(ellipse_var) (optional)
    color_var_leg_columns = 1,
    # PERMANOVA options (single variable test)
    permanova_var = NULL,               # metadata column to test (factor or numeric)
    strata_var = NULL,                  # optional blocking factor (e.g., SynCom)
    permutations = 999
) {
  distance   <- match.arg(distance)
  preprocess <- match.arg(preprocess)
  
  # Package loading
  pkgs <- c("ggplot2","vegan")
  to_load <- pkgs[!pkgs %in% (.packages())]
  if (length(to_load)) suppressPackageStartupMessages(sapply(to_load, require, character.only = TRUE))
  
  # ---------- Align samples ----------
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df
  
  # find sample IDs in metadata (column or rownames)
  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(md)) stop("sample_id_col not found in metadata.")
    ids <- as.character(md[[sample_id_col]])
  } else {
    if (is.null(rownames(md))) stop("metadata_df must have rownames or provide sample_id_col.")
    ids <- rownames(md)
  }
  common <- intersect(colnames(metab_df), ids)
  if (!length(common)) stop("No overlapping sample IDs between matrix columns and metadata IDs.")
  
  X  <- as.matrix(metab_df[, common, drop = FALSE]) # convert to matrix
  
  # Put md’s rows samples in the same order as X column of feature table
  md <- if (!is.null(sample_id_col)) md[match(common, ids), , drop = FALSE]
  else                         md[match(common, rownames(md)), , drop = FALSE]
  rownames(md) <- common
  
  # Ensure aesthetics columns are in metadata
  if (!color_var %in% colnames(md))   stop("color_var not found in metadata.")
  if (!is.null(shape_var)   && !shape_var   %in% colnames(md)) stop("shape_var not found in metadata.")
  if (!is.null(ellipse_var) && !ellipse_var %in% colnames(md)) stop("ellipse_var not found in metadata.")
  # Factor de aesthetic variables
  md[[color_var]] <- droplevels(factor(md[[color_var]]))
  if (!is.null(shape_var))   md[[shape_var]]   <- droplevels(factor(md[[shape_var]]))
  if (!is.null(ellipse_var)) md[[ellipse_var]] <- droplevels(factor(md[[ellipse_var]]))
  
  # ---------- Build samples x features and preprocess ----------
  S <- t(X)  # samples x metabolites
  
  # Check requirements for the different prepropcessing and distance methods.
  if (distance == "bray") {
    if (any(S < 0, na.rm = TRUE))
      stop("Bray–Curtis requires non-negative data. Use a non-negative matrix or switch distance='euclidean'.")
    if (preprocess == "hellinger")      S <- vegan::decostand(S, method = "hellinger")
    else if (preprocess == "sqrt")      S <- sqrt(S)
    dd <- vegan::vegdist(S, method = "bray")
  } else { # euclidean
    if (preprocess %in% c("hellinger","sqrt")) {
      if (any(S < 0, na.rm = TRUE)) stop("Selected preprocessing requires non-negative data.")
      if (preprocess == "hellinger") S <- vegan::decostand(S, method = "hellinger")
      if (preprocess == "sqrt")      S <- sqrt(S)
    }
    dd <- vegan::vegdist(S, method = "euclidean")
  }
  
  # ---------- PCoA ----------
  pc <- cmdscale(dd, eig = TRUE, k = max(2, k_axes), add = TRUE)
  
  scores <- as.data.frame(pc$points)
  colnames(scores) <- paste0("PCo", seq_len(ncol(scores)))
  scores$Sample <- rownames(scores)
  scores <- cbind(scores, md[rownames(scores), , drop = FALSE])
  
  eig <- pc$eig; pos <- eig[eig > 0]
  explained <- 100 * pos / sum(pos)
  
  # ---------- Plot (no aes_string) ----------
  # build dynamic mapping for points
  pt_map <- if (is.null(shape_var)) {
    ggplot2::aes(color = .data[[color_var]])
  } else {
    ggplot2::aes(color = .data[[color_var]], shape = .data[[shape_var]])
  }
  
  p <- ggplot2::ggplot(scores, ggplot2::aes(PCo1, PCo2)) +
    ggplot2::geom_point(mapping = pt_map, size = 2, alpha = 0.9) +
    ggplot2::labs(
      x = paste0("PCo1 (", round(explained[1], 1), "%)"),
      y = paste0("PCo2 (", round(explained[2], 1), "%)"),
      color = color_var,
      shape = if (!is.null(shape_var)) shape_var else NULL
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw()
  
  # custom palette for point colors (optional)
  if (!is.null(points_palette)) {
    p <- p + ggplot2::scale_color_manual(values = points_palette, name = color_var, guide = guide_legend(ncol = color_var_leg_columns))
  }else p <- p + ggplot2::scale_color_discrete(name = color_var, guide = guide_legend(ncol = color_var_leg_columns))
  
  # Ellipses
  if (isTRUE(ellipse) && !is.null(ellipse_var)) {
    grp_counts <- table(scores[[ellipse_var]])
    ok_groups  <- names(grp_counts)[grp_counts >= min_n_for_ellipse]
    if (length(ok_groups)) {
      dat_ell <- scores[scores[[ellipse_var]] %in% ok_groups, , drop = FALSE]
      if (ellipse_var == color_var) {
        # reuse the same color scale as points
        p <- p + ggplot2::stat_ellipse(
          data = dat_ell,
          ggplot2::aes(x = PCo1, y = PCo2,
                       group = .data[[ellipse_var]],
                       color = .data[[ellipse_var]]),
          inherit.aes = FALSE, level = 0.95, type = "norm",
          linewidth = 0.6, alpha = 0.9
        )
      } else {
        if (!"ggnewscale" %in% (.packages())) suppressPackageStartupMessages(require(ggnewscale))
        p <- p +
          ggnewscale::new_scale_color() +
          ggplot2::stat_ellipse(
            data = dat_ell,
            ggplot2::aes(x = PCo1, y = PCo2,
                         group = .data[[ellipse_var]],
                         color = .data[[ellipse_var]]),
            inherit.aes = FALSE, level = 0.95, type = "norm",
            linewidth = 0.7, alpha = 0.9
          ) +
          {
            if (!is.null(ellipse_palette))
              ggplot2::scale_color_manual(values = ellipse_palette,
                                          name = paste0(ellipse_var, " (ellipse)"))
            else
              ggplot2::scale_color_discrete(name = paste0(ellipse_var, " (ellipse)"))
          }
      }
    } else {
      message("Skipping ellipses: fewer than ", min_n_for_ellipse,
              " samples in all groups of '", ellipse_var, "'.")
    }
  }
  
  # Labels on points (optional)
  if (label_points) {
    if (!"ggrepel" %in% (.packages())) suppressPackageStartupMessages(require(ggrepel))
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 2, max.overlaps = 60)
  }
  
  
  # ---------- PERMANOVA (single variable) ----------
  permanova <- NULL
  if (!is.null(permanova_var)) {
    if (!permanova_var %in% colnames(md)) stop("permanova_var not found in metadata.")
    # keep complete cases for the tested (and strata) variables
    keep <- complete.cases(md[, permanova_var, drop = FALSE]) &
      (if (!is.null(strata_var)) complete.cases(md[, strata_var, drop = FALSE]) else TRUE)
    
    if (!all(keep)) message("PERMANOVA: dropping ", sum(!keep), " sample(s) with NA in selected variables.")
    md_perm <- md[keep, , drop = FALSE]
    
    # subset the distance object to kept samples
    labs <- rownames(md_perm)
    dd_mat <- as.matrix(dd)
    dd_perm <- stats::as.dist(dd_mat[labs, labs])
    
    # build formula dd_perm ~ .__var
    md_perm$.__var <- md_perm[[permanova_var]]
    form <- stats::as.formula("dd_perm ~ .__var")
    
    permanova <- vegan::adonis2(
      formula = form,
      data = md_perm,
      permutations = permutations,
      strata = if (!is.null(strata_var)) md_perm[[strata_var]] else NULL
    )
  }
  
  list(
    pcoa      = pc,
    explained = explained,
    scores    = scores,
    plot      = p,
    dist      = dd,
    dist_mat  = as.matrix(dd),
    permanova = permanova,
    settings  = list(distance = distance, preprocess = preprocess,
                     color_var = color_var, shape_var = shape_var,
                     ellipse_var = ellipse_var, permanova_var = permanova_var,
                     strata_var = strata_var, permutations = permutations)
  )
}

#### Heatmap functions
# ---- Main: limma markers (generalized variable names) --------------------------
# Align metabolite matrix (rows = metabolites, cols = samples) to metadata
.align_to_metadata <- function(metab_df, metadata_df, sample_id_col = NULL) {
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df
  
  if (is.null(sample_id_col)) {
    if ("Sample" %in% colnames(md)) sample_id_col <- "Sample"
    else if (!is.null(rownames(md))) { md$Sample <- rownames(md); sample_id_col <- "Sample" }
    else stop("metadata_df must have a 'Sample' column or rownames = sample IDs.")
  }
  
  common <- intersect(colnames(metab_df), md[[sample_id_col]])
  if (!length(common)) stop("No overlapping sample IDs between metabolome and metadata.")
  X  <- as.matrix(metab_df[, common, drop = FALSE])
  md <- md[match(common, md[[sample_id_col]]), , drop = FALSE]
  rownames(md) <- md[[sample_id_col]]
  
  list(X = X, meta = md)
}


limma_markers_by_cluster_general <- function(
    metab_df, metadata_df,
    sample_id_col = NULL,          # e.g., "Sample" (or set rownames(metadata) accordingly and pass "Sample")
    cluster_var,            # e.g., "ATTRIBUTE_Cluster"
    covariates = NULL,      # e.g., c("ATTRIBUTE_Time")
    block_var = NULL,       # e.g., "ATTRIBUTE_SynCom"  (optional)
    log_transform = TRUE, log_offset = 1,
    do_pairwise = TRUE,
    adjust_method = "BH"
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Install it via install.packages('limma').")
  }
  
  aligned <- .align_to_metadata(metab_df, metadata_df, sample_id_col)
  X  <- aligned$X
  md <- aligned$meta
  
  # Build model pieces
  if (!cluster_var %in% colnames(md)) stop("cluster_var not found in metadata.")
  cluster <- factor(md[[cluster_var]])
  # Make cluster design explicitly so we know its columns
  Z_cluster <- model.matrix(~ 0 + cluster)
  colnames(Z_cluster) <- paste0("cluster", levels(cluster))  # e.g., cluster1, cluster2...
  
  # Covariate design (optional). Character -> factor; numeric left as is.
  cov_mat <- NULL
  if (!is.null(covariates) && length(covariates) > 0) {
    miss <- setdiff(covariates, colnames(md))
    if (length(miss)) stop("Missing covariate(s) in metadata: ", paste(miss, collapse = ", "))
    cov_df <- md[, covariates, drop = TRUE]
    cov_df <- as.data.frame(cov_df, stringsAsFactors = FALSE)
    for (nm in names(cov_df)) {
      if (is.character(cov_df[[nm]])) cov_df[[nm]] <- factor(cov_df[[nm]])
    }
    cov_mat <- model.matrix(~ 0 + ., data = cov_df)  # dummy-code factors, keep numerics
    colnames(cov_mat) <- paste0("cov_", colnames(cov_mat))
  }
  
  design <- if (is.null(cov_mat)) Z_cluster else cbind(Z_cluster, cov_mat)
  
  # Optional log transform (recommended for TIC-normalized non-negative data)
  if (log_transform) {
    minX <- suppressWarnings(min(X, na.rm = TRUE))
    if (is.finite(minX) && minX < 0) X <- X - minX
    X <- log(X + log_offset)
  }
  
  # Duplicate correlation (block by, e.g., SynCom) if provided
  use_block <- !is.null(block_var)
  if (use_block) {
    if (!block_var %in% colnames(md)) stop("block_var not found in metadata.")
    block <- md[[block_var]]
    corfit <- limma::duplicateCorrelation(X, design, block = block)
    fit <- limma::lmFit(X, design, block = block, correlation = corfit$consensus)
  } else {
    corfit <- list(consensus = NA_real_)
    fit <- limma::lmFit(X, design)
  }
  
  # One-vs-rest contrasts over cluster columns
  cl_cols <- colnames(Z_cluster)             # e.g., cluster1, cluster2, ...
  K <- length(cl_cols)
  if (K < 2) stop("Need at least 2 cluster levels.")
  
  L_ovr <- sapply(seq_len(K), function(i) {
    v <- rep(-1/(K - 1), K); v[i] <- 1
    names(v) <- cl_cols
    v
  })
  colnames(L_ovr) <- paste0(sub("^cluster", "", cl_cols), "_vs_rest")
  
  # Embed in full design (zeros elsewhere)
  L_full <- matrix(0, nrow = ncol(design), ncol = ncol(L_ovr),
                   dimnames = list(colnames(design), colnames(L_ovr)))
  L_full[cl_cols, colnames(L_ovr)] <- L_ovr
  
  fit_ovr <- limma::contrasts.fit(fit, L_full)
  fit_ovr <- limma::eBayes(fit_ovr)
  
  markers_ovr <- lapply(seq_len(ncol(L_full)), function(i) {
    limma::topTable(fit_ovr, coef = i, number = Inf, adjust.method = adjust_method) %>%
      rownames_to_column("Metabolite") %>%
      mutate(Contrast = colnames(L_full)[i],
             TargetCluster = sub("_vs_rest$", "", colnames(L_full)[i]))
  }) %>% bind_rows()
  
  # Optional: pairwise contrasts among clusters
  markers_pairwise <- NULL
  if (isTRUE(do_pairwise) && K >= 2) {
    combs <- t(combn(cl_cols, 2))
    L_pw <- sapply(seq_len(nrow(combs)), function(i) {
      v <- setNames(rep(0, K), cl_cols)
      A <- combs[i, 1]; B <- combs[i, 2]
      v[B] <-  1; v[A] <- -1
      v
    })
    colnames(L_pw) <- apply(combs, 1, function(x) {
      paste0(sub("^cluster", "", x[2]), "_vs_", sub("^cluster", "", x[1]))
    })
    
    Lpw_full <- matrix(0, nrow = ncol(design), ncol = ncol(L_pw),
                       dimnames = list(colnames(design), colnames(L_pw)))
    Lpw_full[cl_cols, colnames(L_pw)] <- L_pw
    
    fit_pw <- contrasts.fit(fit, Lpw_full)
    fit_pw <- eBayes(fit_pw)
    
    markers_pairwise <- lapply(seq_len(ncol(Lpw_full)), function(i) {
      topTable(fit_pw, coef = i, number = Inf, adjust.method = adjust_method) %>%
        rownames_to_column("Metabolite") %>%
        mutate(Contrast = colnames(Lpw_full)[i])
    }) %>% bind_rows()
  }
  
  list(
    markers_one_vs_rest = markers_ovr,
    markers_pairwise    = markers_pairwise,
    duplicate_correlation = corfit$consensus,
    design = design,
    cluster_levels = sub("^cluster", "", cl_cols),
    used_args = list(cluster_var = cluster_var, covariates = covariates,
                     block_var = block_var, log_transform = log_transform,
                     log_offset = log_offset, adjust_method = adjust_method)
  )
}

summarize_markers_and_heatmap_with_classes <- function(
    # --- core inputs ---
  metab_df, metadata_df,
  sample_id_col = NULL,             # e.g. "Sample"; if NULL, uses metadata rownames
  cluster_var,                       # e.g. "ATTRIBUTE_Cluster"
  # --- SIRIUS inputs ---
  sirius_df,
  id_col = "row.ID",
  class_cols = c("SIRIUS_ClassyFire.class",
                 "SIRIUS_ClassyFire.most.specific.class",
                 "SIRIUS_ClassyFire.subclass",
                 "SIRIUS_ClassyFire.level.5"),
  id_pattern = "^X(\\d+).*",
  # --- limma selection ---
  limma_res,                         # output of limma_markers_by_cluster_general()
  top_n = 15,
  p_adj_thresh = 0.05,
  min_logFC = 0,
  log_transform = TRUE, log_offset = 1,
  scale_rows = TRUE,
  heatmap_colors = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),
  # --- plotting options ---
  cluster_colors = NULL,             # named vector for sample clusters (optional)
  class_na_label = "Unclassified",
  class_na_color = "#BDBDBD",
  out_file   = NULL,  # e.g. "markers_heatmap.pdf" or "markers_heatmap.png"
  out_width  = 9,     # inches
  out_height = 12,    # inches
  out_dpi    = 300,   # used for raster formats (png/jpg/tiff)
  c_legend_ncol = 2,     # columns for both heatmap & annotation legends,
  r_legend_ncol = 2,
  legend_side = "bottom",   # "bottom", "top", "left", or "right"
  merge_legends = FALSE     # TRUE = put all legends into one combined block
) {
  ## ---------- Step A: align + column annotation ----------
  if (is.null(colnames(metab_df))) stop("metab_df must have sample column names.")
  if (!(cluster_var %in% colnames(metadata_df))) {
    stop("Cluster column '", cluster_var, "' not found in metadata_df.")
  }
  if (!is.null(sample_id_col) && sample_id_col %in% colnames(metadata_df)) {
    md_ids <- as.character(metadata_df[[sample_id_col]]); md_id_src <- paste0("column '", sample_id_col, "'")
  } else if (!is.null(rownames(metadata_df)) && all(nchar(rownames(metadata_df)) > 0)) {
    md_ids <- rownames(metadata_df); md_id_src <- "metadata rownames"; sample_id_col <- "<rownames>"
  } else stop("Could not find sample IDs in metadata_df (need sample_id_col or non-empty rownames).")
  
  common <- intersect(colnames(metab_df), md_ids)
  if (length(common) == 0) stop("No overlapping sample IDs between metab_df and metadata_df (", md_id_src, ").")
  
  X <- as.matrix(metab_df[, common, drop = FALSE])
  if (sample_id_col == "<rownames>") {
    md <- metadata_df[match(common, rownames(metadata_df)), , drop = FALSE]; rownames(md) <- rownames(md)
  } else {
    md <- metadata_df[match(common, metadata_df[[sample_id_col]]), , drop = FALSE]; rownames(md) <- md[[sample_id_col]]
  }
  
  ann_col_all <- data.frame(Cluster = factor(md[[cluster_var]]))
  rownames(ann_col_all) <- rownames(md)
  message("Prepared alignment: ", ncol(X), " samples; clusters = {", paste(levels(ann_col_all$Cluster), collapse = ", "), "}")
  
  ## ---------- Step B: SIRIUS row annotations ----------
  if (!id_col %in% colnames(sirius_df)) {
    stop("SIRIUS id_col '", id_col, "' not found in sirius_df. Available: ", paste(colnames(sirius_df), collapse = ", "))
  }
  class_cols_present <- intersect(class_cols, colnames(sirius_df))
  if (length(class_cols_present) == 0) {
    stop("None of the requested class_cols found in sirius_df. Available: ", paste(colnames(sirius_df), collapse = ", "))
  }
  
  metas <- rownames(X); if (is.null(metas)) stop("metab_df must have metabolite rownames.")
  ids_extracted <- sub(id_pattern, "\\1", metas)
  
  sirius_min <- sirius_df |>
    dplyr::mutate(.id = as.character(.data[[id_col]])) |>
    dplyr::distinct(.id, .keep_all = TRUE) |>
    dplyr::select(.id, dplyr::all_of(class_cols_present))
  
  id_map <- data.frame(Metabolite = metas, .id = ids_extracted, stringsAsFactors = FALSE)
  ann_row_full <- id_map |>
    dplyr::left_join(sirius_min, by = ".id") |>
    as.data.frame()
  rownames(ann_row_full) <- ann_row_full$Metabolite
  ann_row_full$Metabolite <- NULL
  if (".id" %in% colnames(ann_row_full)) ann_row_full$.id <- NULL
  
  if (ncol(ann_row_full) > 0) {
    all_na <- vapply(ann_row_full, function(x) all(is.na(x)), logical(1))
    if (any(all_na)) ann_row_full <- ann_row_full[, !all_na, drop = FALSE]
    for (nm in colnames(ann_row_full)) ann_row_full[[nm]] <- droplevels(factor(ann_row_full[[nm]]))
  }
  match_count <- sum(ids_extracted %in% unique(sirius_min$.id))
  message("SIRIUS matched IDs: ", match_count, " / ", length(ids_extracted))
  
  ## ---------- Step C: select top-N markers per cluster (from limma) ----------
  if (!("markers_one_vs_rest" %in% names(limma_res))) {
    stop("limma_res must contain $markers_one_vs_rest.")
  }
  ovr <- limma_res$markers_one_vs_rest
  req_cols <- c("Metabolite","adj.P.Val","logFC","TargetCluster")
  if (!all(req_cols %in% colnames(ovr))) {
    stop("limma_res$markers_one_vs_rest must have columns: ", paste(req_cols, collapse=", "))
  }
  
  top_by_cluster <- ovr |>
    dplyr::filter(adj.P.Val <= p_adj_thresh, logFC >= min_logFC) |>
    dplyr::group_by(TargetCluster) |>
    dplyr::arrange(adj.P.Val, dplyr::desc(logFC), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
  
  top_metabs <- unique(top_by_cluster$Metabolite)
  if (!length(top_metabs)) stop("No metabolites passed thresholds; relax p_adj_thresh or min_logFC.")
  
  keep <- intersect(rownames(X), top_metabs)
  if (!length(keep)) stop("Selected metabolites not found in rownames(metab_df). Check naming.")
  Xsub <- X[keep, , drop = FALSE]
  
  # optional transform for display (match limma)
  if (isTRUE(log_transform)) {
    minX <- suppressWarnings(min(Xsub, na.rm = TRUE)); if (is.finite(minX) && minX < 0) Xsub <- Xsub - minX
    Xsub <- log(Xsub + log_offset)
  }
  
  # order columns by cluster
  md[[cluster_var]] <- factor(md[[cluster_var]])
  col_order <- order(md[[cluster_var]], decreasing = FALSE)
  Xsub <- Xsub[, col_order, drop = FALSE]
  ann_col <- ann_col_all[colnames(Xsub), , drop = FALSE]
  
  # subset row annotations
  ann_row <- ann_row_full[rownames(Xsub), , drop = FALSE]
  
  # NEW: drop unused levels so the legend shows only present classes
  # if (ncol(ann_row) > 0) {
  #   for (nm in colnames(ann_row)) {
  #     ann_row[[nm]] <- droplevels(factor(ann_row[[nm]]))  # force factor + drop unused
  #   }
  # }
  ### end of NEW
  # New 2
  if (ncol(ann_row) > 0) {
    for (nm in colnames(ann_row)) {
      v <- as.character(ann_row[[nm]])
      v[is.na(v) | trimws(v) == ""] <- class_na_label
      ann_row[[nm]] <- droplevels(factor(v))
    }
  }
  # end of new 2
  # optional row scaling for display
  X_display <- Xsub
  if (isTRUE(scale_rows) && nrow(X_display) > 1) {
    X_display <- t(scale(t(X_display))); X_display[!is.finite(X_display)] <- 0
  }
  
  # ---------- Step D: build annotation colors & plot ----------
  # Column (Cluster) palette
  if (is.null(cluster_colors)) {
    k <- nlevels(ann_col$Cluster)
    cluster_colors <- setNames(brewer.pal(max(3, min(8, k)), "Set2")[seq_len(k)],
                               levels(ann_col$Cluster))
  }
  ann_colors <- list(Cluster = cluster_colors)
  
  # New 2
  # NEW: Row-annotation palettes — discrete, present levels only, add NA/"" color
  if (!is.null(ann_row) && ncol(ann_row) > 0) {
    make_discrete_pal <- function(vals) {
      lev <- levels(vals)
      # Separate placeholder label if present
      has_na_lab <- class_na_label %in% lev
      lev_core   <- setdiff(lev, class_na_label)
      n <- length(lev_core)
      base <- if (n <= 12) RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)]
      else          scales::hue_pal()(n)
      pal <- setNames(base, lev_core)
      if (has_na_lab) pal[[class_na_label]] <- class_na_color
      pal
    }
    row_palettes <- lapply(ann_row, make_discrete_pal)
    ann_colors   <- c(ann_colors, row_palettes)
  }
  # end of new 2
  
  # Now plot
  # --- replace pheatmap block with ComplexHeatmap for multi-column legends ---
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install packages 'ComplexHeatmap' and 'circlize' for multi-column legends.")
  }
  
  # heatmap color function (symmetric around 0; X_display is row-scaled z-scores)
  zlim <- max(abs(range(X_display, finite = TRUE)))
  col_fun <- circlize::colorRamp2(c(-zlim, 0, zlim), c("#2166AC","#F7F7F7","#B2182B"))
  
  # split palettes for annotations
  col_ann_cols <- list(Cluster = cluster_colors)  # column annotation colors
  row_ann_cols <- lapply(ann_row, function(v) {
    lev <- levels(v); n <- length(lev)
    cols <- if (n <= 12) RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)] else scales::hue_pal()(n)
    stats::setNames(cols, lev)
  })
  
  ha_top <- ComplexHeatmap::HeatmapAnnotation(
    df  = ann_col,
    col = col_ann_cols,
    annotation_legend_param = list(ncol = c_legend_ncol)   # <-- set legend columns here
  )
  
  ha_left <- if (ncol(ann_row) > 0) {
    ComplexHeatmap::rowAnnotation(
      df  = ann_row,
      col = row_ann_cols,
      annotation_legend_param = list(ncol = r_legend_ncol) # <-- and here for row-annotation legends
    )
  } else NULL
  
  ht <- ComplexHeatmap::Heatmap(
    X_display,
    name = "z",
    col  = col_fun,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    top_annotation  = ha_top,
    left_annotation = ha_left,
    heatmap_legend_param = list(ncol = 2),     # <-- heatmap (z) legend columns
    column_title = sprintf("Top markers per cluster (top %d, FDR≤%.02f)", top_n, p_adj_thresh)
  )
  
  # when saving
  save_ht <- function(ht, file, width, height, dpi, legend_side, merge_legends) {
    ext <- tolower(tools::file_ext(file))
    if (ext %in% c("pdf"))  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    else if (ext %in% c("svg")) grDevices::svg(file, width = width, height = height)
    else if (ext %in% c("png")) grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi, units = "px")
    else if (ext %in% c("tif","tiff")) grDevices::tiff(file, width = width * dpi, height = height * dpi, res = dpi, units = "px", compression = "lzw")
    else if (ext %in% c("jpg","jpeg")) grDevices::jpeg(file, width = width * dpi, height = height * dpi, res = dpi, units = "px", quality = 95)
    else stop("Unsupported file extension: ", ext)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = legend_side,
      annotation_legend_side = legend_side,
      merge_legend = merge_legends
    )
  }
  
  # call it
  if (!is.null(out_file)) {
    save_ht(ht, out_file, out_width, out_height, out_dpi, legend_side, merge_legends)
  } else {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = legend_side,
      annotation_legend_side = legend_side,
      merge_legend = merge_legends
    )
  }
  
  
  # return pieces + plot object
  list(
    heatmap              = ht,
    X_display            = X_display,
    X_values             = Xsub,
    ann_col              = ann_col,
    ann_row              = ann_row,
    annotation_colors    = ann_colors,
    top_table            = top_by_cluster,
    selected_metabolites = rownames(Xsub)
  )
}




barplots_grid <- function(feature_tables, experiments_names, shared_samples = FALSE, strains = FALSE, plot_title = "",
                          plot_title_size = 14, x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                          y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 0,
                          legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3, legend_key_size = 1, 
                          colour_palette = NULL){
  # Creates a grid of Barplots
  
  ### Step 1. Clean, join and gather the otu tables.
  sample_names = c() # to keep track of the sample names
  for (table in seq(from = 1, to = length(feature_tables), by=1)) { # iterate over all the feature tables
    # copy current feature table to avoid modifying the original table.
    feature_table <- feature_tables[[table]]
    
    #print(head(feature_table2)) # check the working feature table
    
    if (isTRUE(strains)) {
      # Convert table with strain names to a strain-number table
      feature_table <- strain_name2strain_number(feature_table)
    }
    
    # Remove rows with Zero counts
    feature_table <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1)
    
    #print(head(feature_table2))
    
    # save names of species
    species_names <- row.names(feature_table)
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table) > 1) {
      feature_table <- feature_table[, colSums(feature_table != 0) > 0]
    }
    
    sample_names <- c(sample_names, colnames(feature_table))
    
    #print(head(feature_table2))
    
    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table["species"] <- species_names
    #print(feature_table2$species)
    
    # Use dplyr gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table, 1:(ncol(feature_table) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table]) # check experiment name that corresponds to working feature table.
    
    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector
    
    #print(head(feature_table_g))
    
    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data for all tables to make a barplot.
    if (table == 1) {
      plot_df <- feature_table_g
    }else{
      plot_df <- rbind(plot_df, feature_table_g)
    }
  }
  print(sample_names) # check sample_names
  print(head(plot_df)) # check gathered table
  
  ### Step 2. Convert Strain data to a graphing-compatible format.
  # Add strain data column to long dataframe
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      mutate(
        strain = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)  # Remove strain number from species name
      )
  }
  
  print(head(plot_df))
  
  ### Step 3. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    filter(!is.na(abundance) & abundance != 0)
  
  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      filter(!is.na(strain) & strain != 0)
  }
  
  plot_df_filtered$experiment <- factor(plot_df_filtered$experiment, levels = experiments_names)
  
  ### Step 4. Plotting
  # get color palette
  if (is.null(colour_palette)) {
    colour_palette <- get_palette(nColors = length(unique(plot_df$species)))
  }
  
  print(plot_df_filtered) # check final table prevouos to plotting
  
  # Create base plot.
  if (shared_samples) {
    p1 <- ggplot(data = plot_df_filtered, aes(x = experiment, y=abundance)) +
      facet_grid(~sample)
  } else{
    p1 <- ggplot(data = plot_df_filtered, aes(x = sample, y=abundance)) +
      facet_grid(~experiment, scales = "free", space = "free")
  }
  
  # Add elements based on graph type.
  if (isTRUE(strains)) {
    print("strains processing")
    p1 <- p1 + ggpattern::geom_bar_pattern(aes(fill = species2, pattern = strain, pattern_density = strain),
                                           position = "fill",
                                           stat="identity",
                                           show.legend = TRUE,
                                           pattern_color = "white",
                                           pattern_fill = "white",
                                           pattern_angle = 45,
                                           pattern_spacing = 0.025) +
      ggpattern::scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(pattern = guide_legend(override.aes = list(fill = "black")),
             fill = guide_legend(override.aes = list(pattern = "none")))
  } else{
    print("no strains")
    p1 <- p1 + geom_bar(aes(fill = species),
                        position = position_fill(),
                        stat = "identity")
  }
  
  if (!is.null(colour_palette)) {
    p1 <- p1 + ggplot2::scale_fill_manual(values=colour_palette)
  } else{
    print("Colours vec is null, using standard color palette.")
  }
  
  p1 <- p1 +
    theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                   axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                   axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, vjust = 0.5, hjust=1, size = x_axis_text_size),
                   axis.title.y = ggplot2::element_text(size=y_axis_title_size, angle = 90),
                   axis.text.y = ggplot2::element_text(size = x_axis_text_size, angle = y_axis_text_angle),
                   legend.title=ggplot2::element_text(size=legend_title_size),
                   legend.text=ggplot2::element_text(size=legend_text_size),
                   legend.position=legend_pos, legend.key.size = unit(legend_key_size, "cm")) + 
    guides(fill = guide_legend(ncol = legend_cols))
  
  # Show plot
  p1
  
  return(p1)
}


# This function takes a dataframe where the rownames are strain level OTUs/ASVs in the form:
# Genera species strain data. The two first words are used a the Species names that are numbered then as:
# Genera species 1; Genera species 2; Genera species 3
strain_name2strain_number <- function(df){
  # Extract only the "Genus species" part
  species_names <- sub(" \\S+$", "", rownames(df))  
  
  # Create a numeric ID for each strain within the same species
  species_ids <- ave(species_names, species_names, FUN = function(x) seq_along(x))
  
  # Create new rownames with species + strain ID
  new_rownames <- paste(species_names, species_ids)
  
  # Assign new rownames to the dataframe
  rownames(df) <- new_rownames
  
  # Print the updated dataframe
  #print(df)
  return(df)
}
