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


ft_pca_2 <- function(ft, metadata_table, grouping_col = NULL, p_shape = NULL, dist_method = "euclidean", colour_palette = NULL){
  ft_t <- t(ft)
  if (isTRUE(all.equal(rownames(ft_t), rownames(metadata_table)))) {
    print("Sample names in feature table and metadatable are identical :)")
    # Compute distance matrix according to dist_method
    # 1. Calculate Bray-Curtis dissimilarity
    dist_mat <- vegan::vegdist(ft_t, method = dist_method)
    
    # 2. Convert the distance matrix to a matrix object
    dist_mat <- as.matrix(dist_mat)
    
    print(head(dist_mat))
    
    # Perform PCoA
    pcoa_results <- cmdscale(dist_mat, k = 2, eig = TRUE)
    
    pcoa_df <- as.data.frame(pcoa_results$points)
    colnames(pcoa_df) <- c("PC1", "PC2")  # Rename axes
    #print(head(pcoa_df))
    
    pcoa_df$Sample <- rownames(ft_t)  # Add sample names to pcoa_df
    
    metadata_table$Sample <- rownames(metadata_table)  # Add sample names to pcoa_df
    
  }else{
    print("Sample names in feature table and metadatable are not identical")
    #print(all.equal(colnames(ft_t),metadata_table$Sample))
    return()
  }
  
  #print(pcoa_df$Sample == row.names(metadata_table))
  if (identical(pcoa_df$Sample,row.names(metadata_table))) {
    print("Sample names in feature table and metadatable are identical :)")
    
  } else{
    print("Sample names in feature table and metadatable are not identical")
    return()
  }
  
  pcoa_df <- left_join(pcoa_df, metadata_table, by = c("Sample" = "Sample"))
  
  if(is.null(colour_palette)){
    colour_palette <- get_palette(nColors = 60)
  }
  
  #print(colour_palette)
  
  if (is.null(grouping_col) && is.null(p_shape)) {
    ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2")
  }else if(!is.null(grouping_col) && is.null(p_shape)){
    print("Plotting with grouping variable")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }else if(!is.null(grouping_col) && !is.null(p_shape)){
    print("Plotting with two grouping variables")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]], shape = .data[[p_shape]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }
  
}

