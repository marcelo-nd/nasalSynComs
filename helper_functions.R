
remove_feature_by_prefix <- function(df, patterns) {
  # Create a single regex pattern that matches any of the species names at the start
  combined_pattern <- paste0("^(", paste(patterns, collapse = "|"), ")")
  
  # Filter the dataframe: keep rows that do NOT match the pattern
  df_filtered <- df[!grepl(combined_pattern, rownames(df)), ]
  
  return(df_filtered)
}


cluster_barplot_panels <- function(abundance_df, metadata_df = NULL,
                                   k = NULL, colour_palette = NULL) {
  require(cluster)
  require(ggplot2)
  require(reshape2)
  
  # Convert to matrix
  mat <- as.matrix(abundance_df)
  
  # Relative abundance
  mat_rel <- sweep(mat, 2, colSums(mat), "/")
  
  print(mat_rel)
  
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
  
  # Update metadata
  if (!is.null(metadata_df)) {
    if (!"Sample" %in% colnames(metadata_df)) {
      metadata_df$Sample <- rownames(metadata_df)
    }
    metadata_updated <- merge(metadata_df, cluster_df, by = "Sample")
  } else {
    metadata_updated <- NULL
  }
  
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
    metadata = metadata_updated,
    best_k = best_k
  ))
}