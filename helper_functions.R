# ---- Functions for reading data --------------------------

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
