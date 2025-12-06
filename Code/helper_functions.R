########################################################
## Code: Marcelo Navarro Diaz
## Contact: marcelo.n.d@ciencias.unam.mx
########################################################
# Load libraries
suppressPackageStartupMessages({
  library(cluster)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(limma)
  library(vegan)
  library(ComplexHeatmap)
  library(purrr)
})

# ----- Functions for reading and handling data -----
#' Read metadata table from CSV
#'
#' Loads a metadata table from a CSV file, optionally sorting it by row names.
#' The first column of the CSV file is treated as row names (typically sample IDs).
#'
#' @param path Character. Path to the metadata CSV file.
#' @param sort_table Logical. If TRUE, metadata rows are sorted alphabetically
#'   by their row names (default = FALSE).
#'
#' @return A data.frame containing the metadata, with row names assigned from
#'   the first column of the CSV file.
#'
#' @details
#' The function expects the CSV file to have:
#' - A header row  
#' - A first column containing unique sample identifiers  
#'
#' @examples
#' # Read metadata as-is
#' md <- read_metadata("metadata.csv")
#'
#' # Read metadata and sort by sample names
#' md_sorted <- read_metadata("metadata.csv", sort_table = TRUE)
read_metadata <- function(path, sort_table = FALSE) {
  
  # Input checks
  if (!is.character(path) || length(path) != 1) {
    stop("`path` must be a single character string.")
  }
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  # Read CSV
  md <- read.csv(path, row.names = 1, stringsAsFactors = FALSE)
  
  # Optional: sort by row names
  if (isTRUE(sort_table)) {
    md <- md[order(row.names(md)), , drop = FALSE]
  }
  
  return(md)
}


#' Read feature table exported from FBMN/Hitchhiker's Guide workflow
#'
#' Imports a feature table generated following the workflow described in
#' Pakkir Shah et al. (Nature Protocols) for feature-based molecular networking.
#' The function reads the table, optionally sorts rows by feature names, removes
#' trailing `.mzML` from row names, and returns the **transposed** table so that
#' features become columns and samples become rows.
#'
#' @param path Character. Path to the exported feature table (CSV or delimited).
#' @param sort_by_names Logical. If TRUE, row names (usually sample names or
#'   feature IDs) are sorted alphabetically (default = FALSE).
#' @param p_sep Character. Field separator used in the file. Default is `","`.
#'
#' @return A transposed numeric data.frame where:
#'   - rows correspond to samples  
#'   - columns correspond to features  
#'   Feature IDs have `.mzML` removed.
#'
#' @details
#' The function expects an export format compatible with the Hitchhiker's Guide /
#' FBMN workflow (GNPS), typically including:
#' - a header row  
#' - first column containing feature identifiers  
#' - sample intensities across columns  
#'
#' The returned object is transposed because many downstream analyses in this
#' project expect **samples as rows** and **features as columns**.
#'
#' @examples
#' # Read feature table as exported by the FBMN workflow
#' ft <- read_ft("feature_table.csv")
#'
#' # Read and sort by row names
#' ft_sorted <- read_ft("feature_table.csv", sort_by_names = TRUE)
read_ft <- function(path, sort_by_names = FALSE, p_sep = ",") {
  
  # ---- Input checks ----
  if (!is.character(path) || length(path) != 1) {
    stop("`path` must be a single character string.")
  }
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  if (!is.character(p_sep) || length(p_sep) != 1) {
    stop("`p_sep` must be a single character separator string.")
  }
  
  # Read feature table
  ft <- read.csv2(
    file   = path,
    header = TRUE,
    row.names = 1,
    sep    = p_sep,
    dec    = ".",
    stringsAsFactors = FALSE
  )
  
  # Optional: sort by row names
  if (isTRUE(sort_by_names)) {
    ft <- ft[order(row.names(ft)), , drop = FALSE]
  }
  
  # Clean row names (remove trailing .mzML)
  rownames(ft) <- gsub("\\.mzML$", "", rownames(ft))
  
  # Return transposed table
  # Samples become rows, features become columns.
  return(t(ft))
}

#' Filter features by minimum count across samples
#'
#' Filters rows (features) of a count table, keeping only those that have
#' at least `min_count` in at least `col_number` columns.
#'
#' @param feature_table A numeric matrix or data.frame with features in rows
#'   and samples in columns.
#' @param min_count Numeric scalar. Minimum count value for a column to be
#'   considered "present".
#' @param col_number Integer scalar. Minimum number of columns that must meet
#'   `min_count` for a feature to be kept.
#'
#' @return A subset of `feature_table` containing only rows that meet the
#'   filtering criteria. If no features pass the filter, an object with zero
#'   rows (but same columns) is returned.
#'
#' @examples
#' # Keep features present (>= 10 counts) in at least 3 samples
#' filtered <- filter_features_by_col_counts(ft, min_count = 10, col_number = 3)
filter_features_by_col_counts <- function(feature_table, min_count, col_number) {
  
  # Input checks
  if (!is.matrix(feature_table) && !is.data.frame(feature_table)) {
    stop("`feature_table` must be a matrix or data.frame.")
  }
  if (!is.numeric(min_count) || length(min_count) != 1) {
    stop("`min_count` must be a single numeric value.")
  }
  if (!is.numeric(col_number) || length(col_number) != 1) {
    stop("`col_number` must be a single numeric value.")
  }
  
  # Handle "no columns" edge case explicitly
  if (ncol(feature_table) == 0) {
    stop("`feature_table` has no columns.")
  }
  
  # Main filtering logic
  if (ncol(feature_table) > 1) {
    # Count how many columns per row are >= min_count, then filter
    keep <- rowSums(feature_table >= min_count) >= col_number
    return(feature_table[keep, , drop = FALSE])
  } else {
    # Single-column case: only keep rows where that column meets min_count
    ft <- feature_table[feature_table[, 1] >= min_count, , drop = FALSE]
    return(ft)
  }
}

#' Set selected species–sample combinations to zero
#'
#' Sets specified entries of a species-by-sample abundance matrix to zero.
#' This is useful for removing detected species from selected samples
#' (e.g., manual curation or experimental design adjustments).
#'
#' @param df A numeric matrix or data.frame where rows are species
#'   and columns are samples.
#' @param species_name Character. Single species name (must match a row name).
#' @param sample_names Character vector of sample names (must match column names).
#'
#' @return The same object as `df`, with the selected species in the selected
#'   samples set to zero.
#'
#' @details
#' The function performs two safety checks:
#'   - Species must exist in `rownames(df)`.  
#'   - All provided sample names must exist in `colnames(df)`.  
#'
#' If any name is missing, the function stops with an informative error.
#'
#' @examples
#' # Remove species "Staphylococcus aureus" from samples SC7_1 and SC7_2
#' cleaned_df <- zero_out_species_in_samples(
#'   df = abundance_table,
#'   species_name = "Staphylococcus aureus",
#'   sample_names = c("SC7_1", "SC7_2")
#' )
zero_out_species_in_samples <- function(df, species_name, sample_names) {
  
  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }
  if (!is.character(species_name) || length(species_name) != 1) {
    stop("`species_name` must be a single character string.")
  }
  if (!is.character(sample_names)) {
    stop("`sample_names` must be a character vector of column names.")
  }
  
  # Safety check: species must exist
  if (!(species_name %in% rownames(df))) {
    stop("Species '", species_name, "' not found in rownames(df).")
  }
  
  # Safety check: samples must exist
  missing_samples <- sample_names[!sample_names %in% colnames(df)]
  if (length(missing_samples) > 0) {
    stop("Samples not found in df: ",
         paste(missing_samples, collapse = ", "))
  }
  
  # Zero out selected species in selected samples
  df[species_name, sample_names] <- 0
  
  return(df)
}

#' Remove features whose names begin with specified prefixes
#'
#' Filters a feature table by removing any rows (features) whose row names
#' start with one or more given prefixes. Matching is performed using a
#' combined regular expression and is case-sensitive unless patterns
#' include case-insensitive constructs.
#'
#' @param df A matrix or data.frame with features in rows.
#' @param patterns Character vector of prefixes that should trigger removal.
#'   Each element is interpreted as a *literal prefix* used to identify
#'   unwanted features.
#'
#' @return The input `df` with all rows removed whose row names begin with
#'   any of the specified prefixes.
#'
#' @details
#' The function builds a single regex of the form:
#' \code{"^(prefix1|prefix2|...)"}  
#' and removes all rows whose row names match this expression.
#'
#' If no rows match any prefix, the original table is returned unchanged.
#'
#' @examples
#' # Remove all features beginning with "Blank" or "Control"
#' ft_clean <- remove_feature_by_prefix(
#'   df = feature_table,
#'   patterns = c("Blank", "Control")
#' )
remove_feature_by_prefix <- function(df, patterns) {
  
  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("`patterns` must be a non-empty character vector.")
  }
  
  # Construct combined regex pattern
  # Escape regex metacharacters to ensure literal prefix matching
  safe_patterns <- vapply(patterns, utils::glob2rx, character(1))
  combined_pattern <- paste0("^(", paste(safe_patterns, collapse = "|"), ")")
  
  # Filter rows: keep only those NOT matching the combined prefix
  matched <- grepl(combined_pattern, rownames(df))
  df_filtered <- df[!matched, , drop = FALSE]
  
  return(df_filtered)
}

#' Map species-level abundances to strain-level abundances
#'
#' Converts a species-level abundance table into a strain-level table using an
#' inoculation design matrix. For each sample, the abundance of a species in
#' `df1` is assigned to all strains of that species that were inoculated in the
#' corresponding SynCom (as specified in `df2`).
#'
#' @param df1 A numeric matrix or data.frame with **species-level abundances**.
#'   Rows are species, columns are samples (e.g., "SC7_1", "SC7_2", ...).
#' @param df2 A data.frame describing **strain inoculation**. The first column
#'   contains strain names (e.g., "Corynebacterium propinquum 16"), and the
#'   remaining columns are SynCom IDs (e.g., "SC7", "SC12", ...) with values
#'   0/1 indicating whether that strain was inoculated in that SynCom.
#'
#' @return A data.frame of strain-level abundances with:
#'   - rows = strains (from the first column of `df2`)
#'   - columns = samples (same as `df1` column names)
#'
#' @details
#' Species–strain mapping:
#' - Species names in `df1` are expected to match the **first two words** of the
#'   strain names in `df2` (e.g., `"Corynebacterium propinquum"` in `df1`
#'   matches `"Corynebacterium propinquum 16"` in `df2`).
#' - Species names are extracted from strain names using the pattern
#'   `^([A-Za-z]+ [A-Za-z]+).*`.
#'
#' Samples and SynComs:
#' - Sample names in `df1` (columns) are assumed to follow the pattern
#'   `"SCx_rep"` (e.g., `"SC7_1"`, `"SC7_2"`).
#' - The SynCom ID used to query `df2` is obtained by taking the part before
#'   the first underscore (e.g., `"SC7"` from `"SC7_1"`).
#'
#' Strains whose species are not present in `df1` remain at zero.  
#'
#' @examples
#' # df1: species-level OTU/abundance table
#' # df2: inoculation design (first column = strains, other columns = SynComs)
#' strain_abundance <- merge_abundance_by_strain(df1 = species_abundance,
#'                                               df2 = inoculation_design)
merge_abundance_by_strain <- function(df1, df2) {
  
  # Input checks
  if (!is.matrix(df1) && !is.data.frame(df1)) {
    stop("`df1` must be a matrix or data.frame with species-level abundances.")
  }
  if (!is.matrix(df2) && !is.data.frame(df2)) {
    stop("`df2` must be a matrix or data.frame with strain design information.")
  }
  
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  
  if (is.null(rownames(df1))) {
    stop("`df1` must have species names as rownames.")
  }
  if (ncol(df2) < 2) {
    stop("`df2` must have at least two columns: strains and SynCom indicators.")
  }
  
  # Helper: get inoculated strains for a given SynCom
  get_inoculated_strains <- function(df2, sample_name) {
    if (!sample_name %in% colnames(df2)) {
      stop("SynCom column '", sample_name, "' not found in df2.")
    }
    sample_column <- df2[[sample_name]]
    inoculated_indices <- which(sample_column == 1)
    df2[inoculated_indices, 1]  # first column = strain names
  }
  
  # Extract names
  species_names_df1 <- rownames(df1)
  strain_names_df2  <- df2[, 1]
  
  # Prepare empty strain-level abundance matrix
  new_abundance_matrix <- matrix(
    0,
    nrow = nrow(df2),
    ncol = ncol(df1)
  )
  rownames(new_abundance_matrix) <- strain_names_df2
  colnames(new_abundance_matrix) <- colnames(df1)
  
  samples <- colnames(new_abundance_matrix)
  
  # Optional consistency check: required SynCom IDs in df2
  syncom_ids_from_samples <- vapply(strsplit(samples, "_"), `[`, character(1), 1)
  missing_syncoms <- setdiff(unique(syncom_ids_from_samples), colnames(df2))
  if (length(missing_syncoms) > 0) {
    stop("SynCom IDs missing in df2: ",
         paste(missing_syncoms, collapse = ", "))
  }
  
  # Main loop over samples
  for (i in seq_along(samples)) {
    sample_name <- samples[i]
    
    # Extract SynCom ID from sample name (e.g., "SC7" from "SC7_1")
    current_sc <- strsplit(sample_name, "_")[[1]][1]
    
    # Get strains inoculated in this SynCom
    inoc_strains_per_sample <- get_inoculated_strains(df2 = df2,
                                                      sample_name = current_sc)
    
    # Loop over strains inoculated in this SynCom
    for (x in seq_along(inoc_strains_per_sample)) {
      strain_name <- inoc_strains_per_sample[x]
      
      # Row index in df2 / new_abundance_matrix
      index_strain_df2 <- which(strain_names_df2 == strain_name)
      
      # Derive species name from strain name (first two words)
      species_name <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", strain_name)
      
      if (species_name %in% species_names_df1) {
        index_species_df1 <- which(species_names_df1 == species_name)
        
        # Abundance of this species in this sample
        current_abundance <- df1[index_species_df1, sample_name]
        
        # Assign abundance to the corresponding strain row + sample column
        new_abundance_matrix[index_strain_df2, i] <- current_abundance
      }
    }
  }
  
  return(as.data.frame(new_abundance_matrix))
}

#' Merge non-target strains to species-level
#'
#' Takes a strain-level feature (OTU) table and collapses all **non-target**
#' strains to a single pseudo-strain per species, while keeping **target**
#' strains unchanged.
#'
#' @param df A numeric matrix or data.frame with strains in rows and samples
#'   in columns. Row names are expected to contain species and strain
#'   information, where the **first two words** encode the species name
#'   (e.g., "Corynebacterium propinquum 16").
#' @param target_species Character vector of species names (e.g.,
#'   `"Corynebacterium propinquum"`) for which **individual strains should be kept**.
#'
#' @return A data.frame where:
#'   - Rows corresponding to `target_species` remain at strain-level (unchanged).
#'   - All other strains are **aggregated by species**, and each species is
#'     represented by a single row whose name is `"<Genus> <species> 1"`.
#'
#' @details
#' Species names are derived from the row names by taking the first two words:
#' \code{paste(x[1:2], collapse = " ")} after splitting by spaces.
#'
#' For non-target species:
#'   - All strain rows for a given species are summed across rows
#'   - The resulting aggregated row is named `"<Genus> <species> 1"` to preserve
#'     compatibility with functions that expect a strain-like naming convention.
#'
#' If there are no non-target strains, only the target strain rows are returned.
#'
#' @examples
#' # Keep all strains of C. propinquum separate, merge all others by species
#' merged_df <- merge_non_target_strains(
#'   df = strain_table,
#'   target_species = c("Corynebacterium propinquum")
#' )
merge_non_target_strains <- function(df, target_species) {
  
  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame with strains in rows.")
  }
  if (is.null(rownames(df))) {
    stop("`df` must have row names containing species and strain information.")
  }
  if (!is.character(target_species)) {
    stop("`target_species` must be a character vector of species names.")
  }
  
  df <- as.data.frame(df)
  
  # Extract species names from rownames (first two words)
  species_names <- sapply(
    strsplit(rownames(df), " "),
    function(x) paste(x[1:2], collapse = " ")
  )
  
  # Identify target vs non-target rows
  is_target <- species_names %in% target_species
  
  target_df     <- df[is_target, , drop = FALSE]
  non_target_df <- df[!is_target, , drop = FALSE]
  non_target_species <- species_names[!is_target]
  
  # Aggregate non-target strains by species
  if (nrow(non_target_df) > 0) {
    aggregated <- aggregate(
      non_target_df,
      by = list(Species = non_target_species),
      FUN = sum
    )
    
    # Row names become "<Genus> <species> 1"
    rownames(aggregated) <- paste(aggregated$Species, "1")
    aggregated$Species <- NULL
  } else {
    aggregated <- NULL
  }
  
  # Combine target strains with aggregated non-target species
  result <- rbind(target_df, aggregated)
  
  return(result)
}


#' Cluster samples based on relative abundance profiles
#'
#' Performs sample clustering on a feature (OTU) table using relative abundance,
#' hierarchical clustering (Ward.D2), and PAM clustering with silhouette-based
#' estimation of the optimal number of clusters (K) when not provided.
#'
#' @param abundance_df A numeric matrix or data.frame with features (e.g. OTUs,
#'   species) in rows and samples in columns.
#' @param k Optional integer. Number of clusters to use for PAM. If \code{NULL}
#'   (default), K is estimated using the average silhouette width over
#'   K = 2,...,K_max, where K_max = \code{min(10, n_samples - 1)}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{clusters}{A data.frame with columns \code{Sample} and \code{Cluster}
#'     containing the PAM cluster assignment for each sample.}
#'   \item{best_k}{The number of clusters used (either provided via \code{k} or
#'     estimated using the silhouette method).}
#'   \item{rel_abundance_ordered}{A data.frame of relative abundances with
#'     samples ordered according to the hierarchical clustering dendrogram.}
#'   \item{sample_order}{Character vector of sample names in dendrogram order.}
#' }
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Convert the input to a matrix.
#'   \item Convert counts to relative abundance per sample
#'         (each column sums to 1).
#'   \item Compute a distance matrix between samples using Euclidean distance
#'         on transposed relative abundance (\code{dist(t(mat_rel))}).
#'   \item Perform hierarchical clustering with method \code{"ward.D2"}.
#'   \item If \code{k} is \code{NULL}, estimate the best K by fitting PAM
#'         clustering models for K = 2,...,K_max and choosing the K that
#'         maximizes the average silhouette width.
#'   \item Fit the final PAM model with \code{best_k} clusters and extract the
#'         cluster assignments.
#'   \item Order samples according to the hierarchical clustering dendrogram and
#'         reorder the relative abundance matrix accordingly.
#' }
#'
#' Requires the \pkg{cluster} package for \code{cluster::pam()}.
#'
#' @examples
#' # Cluster samples using silhouette-based K selection
#' res <- cluster_samples(abundance_df = otu_table)
#'
#' # Use a fixed number of clusters (e.g., K = 3)
#' res_k3 <- cluster_samples(abundance_df = otu_table, k = 3)
cluster_samples <- function(abundance_df, k = NULL) {
  
  # Input checks
  if (!is.matrix(abundance_df) && !is.data.frame(abundance_df)) {
    stop("`abundance_df` must be a numeric matrix or data.frame.")
  }
  
  mat <- as.matrix(abundance_df)
  
  if (!is.numeric(mat)) {
    stop("`abundance_df` must be numeric (counts or abundances).")
  }
  
  n_samples <- ncol(mat)
  if (n_samples < 2) {
    stop("`abundance_df` must contain at least 2 samples (columns) for clustering.")
  }
  
  # Relative abundance per sample
  col_sums <- colSums(mat, na.rm = TRUE)
  if (any(col_sums == 0)) {
    stop("One or more samples have total abundance of zero; cannot compute relative abundances.")
  }
  
  mat_rel <- sweep(mat, 2, col_sums, "/")
  
  mat_input <- mat_rel
  
  # Distance matrix and hierarchical clustering
  dist_mat <- stats::dist(t(mat_input), method = "euclidean")
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  
  # Determine number of clusters (best_k)
  if (is.null(k)) {
    # Upper bound for K: at most 10 or n_samples - 1
    k_max <- min(10, n_samples - 1)
    
    if (k_max < 2) {
      stop("Not enough samples to estimate clusters (need at least 3 samples).")
    }
    
    ks_to_try <- 2:k_max
    
    sil_widths <- sapply(ks_to_try, function(k_try) {
      pam_fit <- cluster::pam(dist_mat, k_try)
      pam_fit$silinfo$avg.width
    })
    
    best_k <- ks_to_try[which.max(sil_widths)]
  } else {
    if (!is.numeric(k) || length(k) != 1 || k <= 1) {
      stop("`k` must be a single integer greater than 1.")
    }
    if (k >= n_samples) {
      stop("`k` must be less than the number of samples.")
    }
    best_k <- as.integer(k)
  }
  
  # Final PAM clustering with best_k
  pam_fit <- cluster::pam(dist_mat, best_k)
  clusters <- pam_fit$clustering
  
  cluster_df <- data.frame(
    Sample  = names(clusters),
    Cluster = as.factor(clusters),
    row.names = NULL
  )
  
  # Reorder samples by dendrogram order
  sample_order <- hc$labels[hc$order]
  mat_rel_ordered <- mat_rel[, sample_order, drop = FALSE]
  
  # Return results
  return(list(
    clusters = cluster_df,
    best_k = best_k,
    rel_abundance_ordered = as.data.frame(mat_rel_ordered),
    sample_order = sample_order
  ))
}

#' Transform a feature table using scaling or normalization methods
#'
#' Applies one of three common transformations to a feature table:
#' \describe{
#'   \item{\code{"zscale"}}{Z-score standardization of each column using \code{scale()}.}
#'   \item{\code{"min_max"}}{Min–max normalization of each numeric column to the range [0, 1].}
#'   \item{\code{"rel_abundance"}}{Conversion to relative abundance where each column sums to 1.}
#' }
#'
#' @param feature_table A numeric matrix or data.frame with features in rows
#'   and samples in columns.
#' @param transform_method Character string. One of:
#'   \code{"zscale"}, \code{"min_max"}, \code{"rel_abundance"}.
#'
#' @return A transformed data.frame with the same dimensions as \code{feature_table}.
#'
#' @details
#' \itemize{
#'   \item For \code{"zscale"}, columns are centered and scaled to unit variance.
#'   \item For \code{"min_max"}, non-numeric columns are left unchanged.
#'   \item For \code{"rel_abundance"}, columns with sum zero will trigger an error
#'         to avoid division by zero.
#' }
#'
#' @examples
#' scaled  <- transform_feature_table(otu_table, "zscale")
#' normed  <- transform_feature_table(otu_table, "min_max")
#' relab   <- transform_feature_table(otu_table, "rel_abundance")
transform_feature_table <- function(feature_table, transform_method) {
  
  # ---- Input checks ----
  if (!is.matrix(feature_table) && !is.data.frame(feature_table)) {
    stop("`feature_table` must be a matrix or data.frame.")
  }
  
  if (!transform_method %in% c("zscale", "min_max", "rel_abundance")) {
    stop("Invalid `transform_method`. Must be 'zscale', 'min_max', or 'rel_abundance'.")
  }
  
  feature_table <- as.data.frame(feature_table)
  
  # ---- Transformation logic ----
  if (transform_method == "zscale") {
    
    df_transformed <- as.data.frame(scale(feature_table))
    
  } else if (transform_method == "min_max") {
    
    df_transformed <- feature_table
    numeric_cols <- sapply(df_transformed, is.numeric)
    
    normalize <- function(x) {
      rng <- max(x) - min(x)
      if (rng == 0) {
        warning("A column has zero variance; min-max scaling will return zeros.")
        return(rep(0, length(x)))
      }
      (x - min(x)) / rng
    }
    
    df_transformed[numeric_cols] <- lapply(df_transformed[numeric_cols], normalize)
    
  } else if (transform_method == "rel_abundance") {
    
    col_sums <- colSums(feature_table, na.rm = TRUE)
    if (any(col_sums == 0)) {
      stop("One or more columns sum to 0; cannot compute relative abundance normalization.")
    }
    
    df_transformed <- sweep(feature_table, 2, col_sums, FUN = "/")
    df_transformed <- as.data.frame(df_transformed)
  }
  
  return(df_transformed)
}

#' Sort Nanopore feature table columns by barcode order
#'
#' Sorts the columns of a Nanopore feature/OTU table by barcode-like column
#' names, using a mixed strategy that orders first by the character length of
#' the column name and then lexicographically. This is useful when barcodes are
#' named like "BC1", "BC2", ..., "BC10", so that "BC2" comes before "BC10".
#'
#' Optionally, new column names can be assigned after sorting.
#'
#' @param df A matrix or data.frame with features in rows and Nanopore samples
#'   (barcodes) in columns.
#' @param new_names Optional character vector of new column names. If provided,
#'   its length must match the number of columns in \code{df}.
#'
#' @return A data.frame with columns sorted by barcode order. If
#'   \code{new_names} is provided, columns are renamed accordingly.
#'
#' @details
#' Column names are sorted using:
#' \itemize{
#'   \item primary key: \code{nchar(colname)}
#'   \item secondary key: lexicographic order of \code{colname}
#' }
#'
#' This avoids the usual problem where purely lexicographic sorting would place
#' "BC10" before "BC2".
#'
#' @examples
#' # Sort a Nanopore OTU table by barcode-like column names
#' df_sorted <- sort_nanopore_table_by_barcodes(nanopore_otu)
#'
#' # Sort and assign human-readable sample names
#' df_sorted_named <- sort_nanopore_table_by_barcodes(
#'   df = nanopore_otu,
#'   new_names = paste0("Sample_", seq_len(ncol(nanopore_otu)))
#' )
sort_nanopore_table_by_barcodes <- function(df, new_names = NULL) {
  
  # ---- Input checks ----
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }
  
  df <- as.data.frame(df)
  
  if (is.null(colnames(df))) {
    stop("`df` must have column names (barcodes).")
  }
  
  if (!is.null(new_names)) {
    if (!is.character(new_names)) {
      stop("`new_names` must be a character vector.")
    }
    if (length(new_names) != ncol(df)) {
      stop("Length of `new_names` (", length(new_names),
           ") must match the number of columns in `df` (", ncol(df), ").")
    }
  }
  
  # ---- Sort column names by length, then lexicographically ----
  cn <- colnames(df)
  sorted_names <- cn[order(nchar(cn), cn)]
  
  df_sorted <- df[, sorted_names, drop = FALSE]
  
  # ---- Optionally rename columns ----
  if (!is.null(new_names)) {
    colnames(df_sorted) <- new_names
  }
  
  return(df_sorted)
}

# ----- Cluster Barplots -----
#' Generate a Palette of Distinct Colors
#'
#' Returns a vector of distinct colors for plotting. By default, the function
#' samples \code{nColors} colors from a predefined palette of visually distinct
#' color values. Colors may be sampled with or without replacement.
#'
#' @param nColors Integer. Number of colors to return. Default is \code{60}.
#' @param replace_cols Logical. Whether to sample colors with replacement.
#'   Default is \code{FALSE}.
#'
#' @return A character vector of color hex codes or named R colors of length
#'   \code{nColors}.
#'
#' @examples
#' # Generate 10 distinct colors
#' get_palette(10)
#'
#' # Sample colors with replacement
#' get_palette(10, replace_cols = TRUE)
#'
#' @export
get_palette <- function(nColors = 60, replace_cols = FALSE){
  colors_vec <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
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
    "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD"
  )
  
  sample(colors_vec, nColors, replace = replace_cols)
}


#' Create Relative Abundance Barplots by Cluster
#'
#' Creates a stacked relative abundance barplot with one panel per cluster.
#' Samples are grouped by their assigned cluster and displayed as stacked bars
#' of bacterial relative abundance. Optionally, strain-level information can be
#' encoded using fill colors for species and patterns for strains.
#'
#' @param abundance_df A matrix or data frame of relative abundances with
#'   taxa (e.g., species or strain labels) in rows and samples in columns.
#'   Row names are used as the \code{Bacteria} identifiers in the plot.
#' @param cluster_df A data frame containing at least two columns:
#'   \code{Sample} (sample IDs matching the column names of \code{abundance_df})
#'   and \code{Cluster} (cluster assignment for each sample).
#' @param sample_order Optional character vector specifying the order of
#'   samples (column names of \code{abundance_df}) to be used on the x-axis.
#'   If \code{NULL}, the original column order of \code{abundance_df} is used.
#' @param colour_palette Optional named character vector of colors used to
#'   fill taxa. Names must match the \code{Bacteria} values (i.e., row names of
#'   \code{abundance_df}). If \code{NULL}, ggplot2's default color palette is
#'   used.
#' @param strains Logical. If \code{TRUE}, the function assumes that
#'   \code{Bacteria} names encode strain information in the last numeric token
#'   (e.g., \code{"Corynebacterium propinquum 1"}) and will represent species
#'   by fill color and strains by bar patterns. Default is \code{FALSE}.
#' @param best_k Integer or character. The number of clusters (k) to display in
#'   the plot title. This argument must be provided; otherwise, the function
#'   will stop with an error.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Orders the abundance matrix according to \code{sample_order}
#'     if provided.
#'   \item Optionally converts strain names to numeric strain IDs using
#'     \code{strain_name2strain_number()} when \code{strains = TRUE}.
#'   \item Reshapes the abundance matrix to long format and merges it with
#'     \code{cluster_df} via the \code{Sample} column.
#'   \item When \code{strains = FALSE}, creates a standard stacked barplot of
#'     taxa relative abundances per sample, faceted by cluster.
#'   \item When \code{strains = TRUE}, uses \pkg{ggpattern} to encode species
#'     as fill colors and strains as patterns within each stacked bar.
#' }
#'
#' The function relies on \pkg{ggplot2}, \pkg{reshape2}, \pkg{dplyr},
#' \pkg{ggpattern}, and a user-defined helper
#' \code{strain_name2strain_number()} that converts strain names to numeric
#' strain labels.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} object representing the relative abundance
#'     barplot faceted by cluster.}
#'   \item{df_long}{A data frame in long format containing the abundance data,
#'     cluster annotations, and (optionally) species/strain columns used for
#'     plotting.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage (no strain patterns, default colors)
#' res <- cluster_barplot_panels(
#'   abundance_df = abundance_mat,
#'   cluster_df   = sample_clusters,
#'   best_k       = 4
#' )
#' res$plot
#'
#' # With custom sample order and custom colors
#' res <- cluster_barplot_panels(
#'   abundance_df   = abundance_mat,
#'   cluster_df     = sample_clusters,
#'   sample_order   = c("S1", "S2", "S3"),
#'   colour_palette = my_colors,
#'   best_k         = 4
#' )
#'
#' # With strain-level encoding (requires strain_name2strain_number())
#' res <- cluster_barplot_panels(
#'   abundance_df = abundance_mat,
#'   cluster_df   = sample_clusters,
#'   strains      = TRUE,
#'   best_k       = 4
#' )
#' }
#'
#' @export
cluster_barplot_panels <- function(
    abundance_df,
    cluster_df,
    sample_order   = NULL,
    colour_palette = NULL,
    strains        = FALSE,
    best_k         = NULL
) {
  require(cluster)
  require(ggplot2)
  require(reshape2)
  require(dplyr)
  require(ggpattern)
  
  if (is.null(best_k)) {
    stop("Argument 'best_k' must be provided.")
  }
  
  # Base matrix (all samples), optionally reordered
  mat_rel_ordered <- as.matrix(abundance_df)
  if (!is.null(sample_order)) {
    mat_rel_ordered <- mat_rel_ordered[, sample_order, drop = FALSE]
  }
  
  if (isTRUE(strains)) {
    message("Using strain data")
    # Convert table with strain names to a strain-number table
    mat_rel_ordered <- strain_name2strain_number(mat_rel_ordered)
  }
  
  # Melt for ggplot
  df_long <- reshape2::melt(mat_rel_ordered)
  colnames(df_long) <- c("Bacteria", "Sample", "Abundance")
  df_long <- merge(df_long, cluster_df, by = "Sample")
  
  # Add strain data columns if needed
  if (isTRUE(strains)) {
    df_long <- df_long %>%
      dplyr::mutate(
        strain   = paste0("Strain ", sub(".* ", "", Bacteria)),  # last token as strain
        species2 = sub(" \\d+$", "", Bacteria)                   # drop trailing number
      )
  }
  
  if (isFALSE(strains)) {
    p1 <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Bacteria)) +
      geom_bar(stat = "identity") +
      facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ylab("Relative Abundance") +
      ggtitle(paste("Stacked Barplot with Clusters (k =", best_k, ")"))
    message("Created plot without strain data")
  } else if (isTRUE(strains)) {
    # Clean the long-format table
    df_long <- df_long %>%
      dplyr::filter(!is.na(Abundance) & Abundance != 0) %>%
      dplyr::filter(!is.na(strain) & strain != 0)
    
    p1 <- ggplot(data = df_long, aes(x = Sample, y = Abundance)) +
      ggpattern::geom_bar_pattern(
        aes(fill = species2, pattern = strain),
        position        = "fill",
        stat            = "identity",
        show.legend     = TRUE,
        pattern_spacing = unit(2.5, "mm"),
        pattern_density = 0.0050,
        pattern_color   = "white",
        pattern_fill    = "white",
        pattern_angle   = 45
      ) +
      facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
      ggpattern::scale_pattern_manual(
        values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")
      ) +
      ggpattern::scale_pattern_spacing_manual(
        values = c(0, unit(0.025, "mm"), unit(0.025, "mm"))
      ) +
      guides(
        pattern = guide_legend(override.aes = list(fill = "grey")),
        fill    = guide_legend(override.aes = list(pattern = "none"))
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    message("Created plot with strain data")
  }
  
  if (!is.null(colour_palette)) {
    # Expecting a named vector: names must match Bacteria (rownames(abundance_df))
    p1 <- p1 + scale_fill_manual(values = colour_palette, drop = FALSE)
    message("Added custom color scale")
  }
  
  return(list(
    plot   = p1,
    df_long = df_long
  ))
}


#' Convert Strain-Level OTU/ASV Names to Species-Level Numbered Labels
#'
#' Converts row names of a data frame from full strain-level identifiers
#' (e.g., \code{"Genus species strainX"}) to a standardized species + strain ID
#' format (e.g., \code{"Genus species 1"}, \code{"Genus species 2"}).
#'
#' The function assumes that row names contain at least three whitespace-
#' separated components, where the first two correspond to the genus and species,
#' and the final component encodes the strain identity. Rows belonging to the
#' same species receive sequential numeric identifiers in the order they appear.
#'
#' @param df A data frame or matrix whose row names contain strain-level OTU/ASV
#'   identifiers. The row names must follow the structure:
#'   \code{"Genus species strainInfo"}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts the first two words (\code{"Genus species"}) from each row name.
#'   \item Assigns sequential numeric strain IDs within each species.
#'   \item Replaces the original row names with the standardized format
#'     \code{"Genus species <ID>"}.
#' }
#'
#' @return The input data frame with updated row names representing species-level
#'   groups with numeric strain identifiers.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6)
#' rownames(df) <- c("Corynebacterium propinquum A1",
#'                   "Corynebacterium propinquum B7",
#'                   "Staphylococcus aureus T3")
#'
#' strain_name2strain_number(df)
#'
#' @export
strain_name2strain_number <- function(df){
  # Extract only the "Genus species" part
  species_names <- sub(" \\S+$", "", rownames(df))
  
  # Create a numeric ID for each strain within the same species
  species_ids <- ave(species_names, species_names, FUN = function(x) seq_along(x))
  
  # Create new row names with species + strain ID
  new_rownames <- paste(species_names, species_ids)
  
  # Assign the new rownames to the dataframe
  rownames(df) <- new_rownames
  
  return(df)
}

#' Cluster Samples and Compute Mean Abundance of a Given Species
#'
#' Performs hierarchical clustering of samples based on their abundance profiles
#' and computes the mean relative abundance of a specified species within each
#' cluster. Optionally prints the sample names assigned to each cluster.
#'
#' @param df A numeric data frame or matrix where rows represent taxa (e.g.,
#'   species or strains) and columns represent samples.
#' @param species_name Character string giving the row name of the species whose
#'   mean abundance should be computed within each cluster.
#' @param k Integer. Number of clusters to cut the hierarchical clustering tree
#'   into. Default is \code{2}.
#' @param method Character string specifying the distance metric passed to
#'   \code{\link{dist}}. Default is \code{"euclidean"}.
#' @param show_samples Logical. If \code{TRUE}, prints the sample names belonging
#'   to each cluster. Default is \code{FALSE}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Computes a distance matrix between samples (\code{dist(t(df))}).
#'   \item Performs hierarchical clustering with \code{\link{hclust}}.
#'   \item Cuts the tree into \code{k} clusters using \code{\link{cutree}}.
#'   \item Computes the mean abundance of the target species within each cluster.
#' }
#'
#' Results (cluster means and optionally sample lists) are printed to the console.
#'
#' @return This function is called for its side effects (printed output) and does
#'   not return a value.
#'
#' @examples
#' \dontrun{
#' cluster_mean_abundance(df = abundance_matrix,
#'                        species_name = "Corynebacterium propinquum",
#'                        k = 3,
#'                        method = "euclidean",
#'                        show_samples = TRUE)
#' }
#'
#' @export
cluster_mean_abundance <- function(df, species_name, k = 2, method = "euclidean", show_samples = FALSE) {
  # Check species
  if (!(species_name %in% rownames(df))) {
    stop("Species not found in the dataframe.")
  }
  
  # Transpose for clustering samples
  dist_matrix <- dist(t(df), method = method)
  hc <- hclust(dist_matrix)
  
  # Cut tree into k groups
  groups <- cutree(hc, k = k)
  
  # Group sample names by cluster
  cluster_samples <- split(names(groups), groups)
  
  # Print results
  cat("Number of clusters:", length(cluster_samples), "\n\n")
  
  for (i in seq_along(cluster_samples)) {
    samples <- cluster_samples[[i]]
    mean_abund <- mean(as.numeric(df[species_name, samples]))
    
    cat("Cluster", i, "- Mean relative abundance of",
        species_name, ":", round(mean_abund, 5), "\n")
    
    if (show_samples) {
      cat("  Samples in cluster", i, ":\n")
      cat("   ", paste(samples, collapse = ", "), "\n\n")
    }
  }
}


#' Add a Cluster Column to a Metadata Data Frame
#'
#' Maps cluster assignments from a lookup table onto a metadata data frame by
#' matching key columns. The function performs sanity checks, warns about
#' duplicated keys in the cluster table, optionally warns about metadata rows
#' with missing matches, and appends a new column containing the mapped cluster
#' values.
#'
#' @param meta_df A data frame containing metadata. A key column within this
#'   data frame will be used to map cluster values onto it.
#' @param clusters_df A data frame containing cluster assignments. Must include
#'   a key column and a value column.
#' @param meta_key_col Character string. Name of the column in \code{meta_df}
#'   used to match entries against \code{clusters_df}.
#' @param cluster_key_col Character string. Name of the key column in
#'   \code{clusters_df} that corresponds to \code{meta_key_col}.
#' @param cluster_value_col Character string. Name of the column in
#'   \code{clusters_df} containing the cluster labels or values to be mapped.
#' @param new_col_name Character string. Name of the new column to be added to
#'   \code{meta_df} containing the mapped cluster entries.
#' @param warn_missing Logical. If \code{TRUE} (default), warns when some values
#'   in \code{meta_df[[meta_key_col]]} have no corresponding match in
#'   \code{clusters_df}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Checks that all required columns exist in both data frames.
#'   \item Warns if duplicate keys appear in the cluster table; in such cases,
#'         the last occurrence is used for mapping.
#'   \item Constructs a named lookup vector mapping keys to cluster values.
#'   \item Uses this lookup to assign cluster values to the metadata table.
#'   \item Optionally warns about unmatched keys in the metadata.
#' }
#'
#' This is a general-purpose utility for augmenting metadata with additional
#' annotations derived from separate data tables.
#'
#' @return A modified version of \code{meta_df} containing a new column
#'   \code{new_col_name} with mapped cluster values.
#'
#' @examples
#' \dontrun{
#' meta <- data.frame(Sample = c("S1", "S2", "S3"))
#' clusters <- data.frame(SampleID = c("S1", "S2"), Cluster = c(1, 2))
#'
#' add_cluster_column(
#'   meta_df          = meta,
#'   clusters_df      = clusters,
#'   meta_key_col     = "Sample",
#'   cluster_key_col  = "SampleID",
#'   cluster_value_col = "Cluster",
#'   new_col_name     = "Cluster_Assignment"
#' )
#' }
#'
#' @export
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


#' Filter Features by Minimum Counts Across Columns
#'
#' Filters rows (features) in a feature table based on how many columns meet or
#' exceed a specified minimum count. This is a general utility for retaining
#' only features that are sufficiently present across samples or conditions.
#'
#' @param feature_table A numeric matrix or data frame in which rows represent
#'   features (e.g., taxa, metabolites, genes) and columns represent samples or
#'   observations.
#' @param min_count Numeric. Minimum value a feature must reach in a column to
#'   be considered present.
#' @param col_number Integer. Minimum number of columns in which the feature
#'   must meet or exceed \code{min_count} in order to be retained.
#'
#' @details
#' The function behaves differently depending on the number of columns:
#' \itemize{
#'   \item If \code{ncol(feature_table) > 1}, a feature is retained if it has at
#'         least \code{col_number} columns with values \code{>= min_count}.
#'   \item If \code{ncol(feature_table) == 1}, the function simply keeps rows
#'         whose single column value meets the threshold.
#'   \item If the table has zero columns, a message is printed and no value is
#'         returned.
#' }
#'
#' @return A filtered version of \code{feature_table} containing only rows that
#'   meet the specified abundance criteria. For one-column tables, a subset of
#'   the original table is returned with \code{drop = FALSE}.
#'
#' @examples
#' ft <- data.frame(
#'   f1 = c(5, 0, 3),
#'   f2 = c(2, 1, 4)
#' )
#'
#' # Keep features present in at least 2 columns with count >= 3
#' filter_features_by_col_counts(ft, min_count = 3, col_number = 2)
#'
#' # One-column case
#' ft2 <- data.frame(f1 = c(0, 5, 10))
#' filter_features_by_col_counts(ft2, min_count = 5, col_number = 1)
#'
#' @export
filter_features_by_col_counts <- function(feature_table, min_count, col_number){
  if (ncol(feature_table) > 1) {
    return(feature_table[which(rowSums(feature_table >= min_count) >= col_number), ])
  }
  else if(ncol(feature_table) == 1){
    ft <- feature_table[feature_table >= min_count, , drop = FALSE]
    return(ft)
  }
  else{
    message("Dataframe has no columns")
  }
}

#' Order Samples According to Hierarchical Clustering
#'
#' Computes a hierarchical clustering of samples based on their feature
#' profiles and returns the sample names in the order determined by the
#' clustering dendrogram. This is useful for ordering heatmaps, barplots,
#' or other visualizations consistently with sample similarity.
#'
#' @param feature_table A numeric data frame or matrix where rows represent
#'   features (e.g., taxa, genes, metabolites) and columns represent samples.
#'   Row names are treated as feature identifiers.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Converts row names into a column and transposes the data so that
#'         samples become rows for clustering.
#'   \item Computes Euclidean distances between samples.
#'   \item Performs hierarchical clustering using Ward's \code{D2} method.
#'   \item Extracts sample names in the order they appear in the dendrogram.
#' }
#'
#' The resulting sample order is often used for heatmap column ordering or to
#' ensure consistent clustering-based visualization across analyses.
#'
#' @return A character vector of sample names ordered according to the
#'   hierarchical clustering.
#'
#' @examples
#' \dontrun{
#' sample_order <- order_samples_by_clustering(feature_table)
#' heatmap(feature_table[, sample_order])
#' }
#'
#' @export
order_samples_by_clustering <- function(feature_table){
  # Add feature names as a column (Species)
  df_otu <- feature_table %>% rownames_to_column(var = "Species")
  
  # Transpose table so samples become rows (remove Species column first)
  df_t <- as.matrix(t(df_otu[, -1]))
  
  # Hierarchical clustering of samples
  d <- dist(df_t, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  # Extract ordered sample names (skip the Species column name)
  ordered_samples_cluster <- colnames(df_otu)[-1][hc$order]
  
  return(ordered_samples_cluster)
}

#' Create a Grid of Relative Abundance Barplots Across Experiments
#'
#' Builds a faceted grid of stacked relative abundance barplots from multiple
#' feature tables (e.g., OTU/ASV count tables or relative abundance tables),
#' typically corresponding to different experiments or treatments.
#'
#' Each feature table is reshaped to long format, combined into a single data
#' frame, and plotted as stacked bars either:
#' \itemize{
#'   \item per sample, faceted by experiment (default), or
#'   \item per experiment, faceted by sample (when \code{shared_samples = TRUE}).
#' }
#'
#' Optionally, strain-level information can be encoded using patterns
#' (\pkg{ggpattern}) while species are encoded by fill color.
#'
#' @param feature_tables A list of numeric matrices or data frames, where rows
#'   represent features (e.g., taxa/strains) and columns represent samples. Each
#'   element of the list corresponds to one experiment or condition.
#' @param experiments_names Character vector of the same length as
#'   \code{feature_tables}, giving the experiment/treatment name for each
#'   feature table. These names are used to annotate and facet the plots.
#' @param shared_samples Logical. If \code{TRUE}, assumes that the same samples
#'   are shared across experiments and plots \code{experiment} on the x-axis
#'   with \code{sample} as facets. If \code{FALSE} (default), plots
#'   \code{sample} on the x-axis and facets by \code{experiment}.
#' @param strains Logical. If \code{TRUE}, treats row names of each feature
#'   table as strain-level identifiers of the form
#'   \code{"Genus species strainID"} and uses \code{\link{strain_name2strain_number}}
#'   and pattern aesthetics (\pkg{ggpattern}) to represent strains. Default
#'   is \code{FALSE}.
#' @param plot_title Character string. Overall title for the plot. Default is
#'   an empty string.
#' @param plot_title_size Numeric. Text size for the plot title. Default is
#'   \code{14}.
#' @param x_axis_text_size Numeric. Text size for x-axis tick labels. Default
#'   is \code{12}.
#' @param x_axis_title_size Numeric. Text size for the x-axis title. Default
#'   is \code{12}.
#' @param x_axis_text_angle Numeric. Rotation angle (in degrees) for x-axis
#'   tick labels. Default is \code{0}.
#' @param y_axis_title_size Numeric. Text size for the y-axis title. Default
#'   is \code{12}.
#' @param y_axis_text_size Numeric. Text size for y-axis tick labels. Default
#'   is \code{12}.
#' @param y_axis_text_angle Numeric. Rotation angle (in degrees) for y-axis
#'   tick labels. Default is \code{0}.
#' @param legend_pos Character string specifying the legend position (e.g.,
#'   \code{"right"}, \code{"bottom"}, \code{"none"}). Passed to
#'   \code{theme(legend.position = ...)}. Default is \code{"right"}.
#' @param legend_title_size Numeric. Text size for legend title. Default is
#'   \code{12}.
#' @param legend_text_size Numeric. Text size for legend text. Default is
#'   \code{12}.
#' @param legend_cols Integer. Number of columns in the legend for the
#'   \code{fill} aesthetic. Default is \code{3}.
#' @param legend_key_size Numeric. Size (in cm) of the legend key boxes.
#'   Default is \code{1}.
#' @param colour_palette Optional named character vector of colors used for
#'   species (or species groups). Names should match the values in
#'   \code{species} (or \code{species2} when \code{strains = TRUE}). If
#'   \code{NULL}, a palette is generated by \code{\link{get_palette}}.
#'
#' @details
#' For each feature table, the function:
#' \itemize{
#'   \item Optionally renames strain-level rows using
#'         \code{\link{strain_name2strain_number}} when \code{strains = TRUE}.
#'   \item Filters out features with all-zero counts using
#'         \code{\link{filter_features_by_col_counts}}.
#'   \item Drops samples (columns) that contain only zeros.
#'   \item Converts row names to a \code{species} column and reshapes the table
#'         to long format using \pkg{tidyr::gather}.
#'   \item Adds an \code{experiment} column from \code{experiments_names}.
#' }
#'
#' When \code{strains = TRUE}, additional columns \code{strain} and
#' \code{species2} are created:
#' \itemize{
#'   \item \code{strain}: a label like \code{"Strain 1"}, \code{"Strain 2"}, etc.
#'   \item \code{species2}: the species name without the trailing numeric ID.
#' }
#'
#' The final combined long data frame is filtered to remove zero abundances,
#' converted to factors with \code{experiment} levels matching
#' \code{experiments_names}, and plotted using \pkg{ggplot2}:
#' \itemize{
#'   \item If \code{shared_samples = FALSE}, x-axis = \code{sample}, faceted
#'         by \code{experiment}.
#'   \item If \code{shared_samples = TRUE}, x-axis = \code{experiment}, faceted
#'         by \code{sample}.
#'   \item For \code{strains = TRUE}, \pkg{ggpattern} is used to distinguish
#'         strains by patterns while species are distinguished by fill colors.
#'   \item For \code{strains = FALSE}, a standard stacked barplot is created
#'         with \code{species} as fill.
#' }
#'
#' Relative abundances are shown as fractions of 1 via
#' \code{position = "fill"}.
#'
#' @return A \code{ggplot} object representing the faceted grid of barplots.
#'
#' @examples
#' \dontrun{
#' p <- barplots_grid(
#'   feature_tables   = list(exp1_abund, exp2_abund),
#'   experiments_names = c("Experiment 1", "Experiment 2"),
#'   shared_samples   = FALSE,
#'   strains          = FALSE,
#'   plot_title       = "Relative abundance across experiments"
#' )
#' p
#' }
#'
#' @export
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


#' Create a Relative Abundance Barplot from a Feature Table
#'
#' Generates a stacked relative abundance barplot from a single feature table
#' (e.g., OTU/ASV or species-by-sample table), with several options for sorting
#' samples and optionally encoding strain-level information using patterns.
#'
#' @param feature_table A numeric matrix or data frame where rows represent
#'   features (e.g., species or strains) and columns represent samples.
#'   Row names are assumed to be feature identifiers.
#' @param sort_type Character string specifying how to order samples on the
#'   x-axis. One of:
#'   \itemize{
#'     \item \code{"none"} (default): keep the original column order.
#'     \item \code{"feature_value"}: order samples by the relative contribution
#'           of a given feature (see \code{feature_to_sort}).
#'     \item \code{"similarity"}: order samples according to hierarchical
#'           clustering via \code{\link{order_samples_by_clustering}}.
#'   }
#' @param feature_to_sort Character string giving the feature (row name) used
#'   for ordering samples when \code{sort_type = "feature_value"}. Ignored
#'   otherwise.
#' @param strains Logical. If \code{TRUE}, treats row names as strain-level
#'   identifiers of the form \code{"Genus species strainID"} and uses
#'   \code{\link{strain_name2strain_number}} plus pattern aesthetics
#'   (\pkg{ggpattern}) to represent strains; species are encoded as fill colors.
#'   Default is \code{FALSE}.
#' @param plot_title Character. Overall plot title. (Note: currently only the
#'   theme title size is controlled; you can extend to use this string explicitly
#'   if desired.) Default is \code{""}.
#' @param plot_title_size Numeric. Text size for the plot title. Default is
#'   \code{14}.
#' @param x_axis_text_size Numeric. Text size for x-axis tick labels. Default
#'   is \code{12}.
#' @param x_axis_title_size Numeric. Text size for the x-axis title. Default
#'   is \code{12}.
#' @param x_axis_text_angle Numeric. Rotation angle (degrees) for x-axis tick
#'   labels. Default is \code{0}.
#' @param y_axis_title_size Numeric. Text size for the y-axis title. Default
#'   is \code{12}.
#' @param y_axis_text_size Numeric. Text size for y-axis tick labels. Default
#*   is \code{12}.
#' @param y_axis_text_angle Numeric. Rotation angle (degrees) for y-axis tick
#'   labels. Default is \code{90}.
#' @param legend_pos Character string specifying the legend position (e.g.,
#'   \code{"right"}, \code{"bottom"}, \code{"none"}). Default is \code{"right"}.
#' @param legend_title_size Numeric. Text size for the legend title. Default is
#'   \code{12}.
#' @param legend_text_size Numeric. Text size for legend text. Default is
#'   \code{12}.
#' @param legend_cols Integer. Number of columns in the legend for the
#'   \code{fill} aesthetic. Default is \code{3}.
#' @param x_vjust Numeric. Vertical justification for x-axis text. Default is
#'   \code{0.5}.
#' @param x_hjust Numeric. Horizontal justification for x-axis text. Default
#'   is \code{1}.
#' @param transform_table Logical. If \code{TRUE} and \code{sort_type =
#'   "similarity"}, the feature table is transformed using
#'   \code{\link{transform_feature_table}} with method \code{"min_max"} prior
#'   to clustering. Default is \code{TRUE}.
#' @param colour_palette Optional named character vector of colors to use for
#'   features (or species groups). If \code{NULL}, a palette is generated with
#'   \code{\link{get_palette}}. When \code{strains = TRUE}, the palette is
#'   applied to \code{species2}.
#' @param replace_c Logical. Passed to \code{\link{get_palette}} as
#'   \code{replace_cols}, controlling whether colors may be sampled with
#'   replacement when generating a palette. Default is \code{FALSE}.
#'
#' @details
#' The function proceeds in several steps:
#' \itemize{
#'   \item Filters out rows (features) with zero counts using
#'         \code{\link{filter_features_by_col_counts}}.
#'   \item Removes samples (columns) that are all zeros.
#'   \item Optionally converts strain-level row names into a standardized
#'         species + numeric strain ID format with
#'         \code{\link{strain_name2strain_number}}.
#'   \item Depending on \code{sort_type}, orders samples by:
#'         \itemize{
#'           \item original column order (\code{"none"}),
#'           \item decreasing relative contribution of a specific feature
#'                 (\code{"feature_value"}),
#'           \item similarity-based clustering
#'                 (\code{"similarity"}, via \code{\link{order_samples_by_clustering}}).
#'         }
#'   \item Reshapes the data into long format with columns \code{species},
#'         \code{sample}, and \code{abundance}. When \code{strains = TRUE},
#'         also creates \code{strain} and \code{species2} columns.
#'   \item Filters out zero and \code{NA} abundances, and, if applicable,
#'         invalid strain entries.
#'   \item Constructs a stacked barplot with relative abundances
#'         (\code{position = "fill"}), using either:
#'         \itemize{
#'           \item \pkg{ggpattern} to encode strains via patterns and species
#'                 via fill colors, or
#'           \item standard stacked bars colored by species.
#'         }
#' }
#'
#' @return A \code{ggplot} object representing the relative abundance barplot.
#'
#' @examples
#' \dontrun{
#' # Basic usage (no sorting)
#' p <- barplot_from_feature_table(feature_table = abundance_mat)
#' p
#'
#' # Sort samples by the relative contribution of a given species
#' p <- barplot_from_feature_table(
#'   feature_table   = abundance_mat,
#'   sort_type       = "feature_value",
#'   feature_to_sort = "Corynebacterium propinquum"
#' )
#'
#' # Sort samples by similarity (clustering-based order)
#' p <- barplot_from_feature_table(
#'   feature_table = abundance_mat,
#'   sort_type     = "similarity",
#'   transform_table = TRUE
#' )
#'
#' # Strain-level visualization with patterns (requires ggpattern)
#' p <- barplot_from_feature_table(
#'   feature_table = strain_level_mat,
#'   strains       = TRUE
#' )
#' }
#'
#' @export
barplot_from_feature_table <- function(feature_table, sort_type = "none", feature_to_sort = NULL, strains = FALSE,
                                       plot_title = "", plot_title_size = 14,
                                       x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                                       y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 90,
                                       legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3,
                                       x_vjust = 0.5, x_hjust = 1, transform_table = TRUE,
                                       colour_palette = NULL, replace_c = FALSE){
  ### Step 1. Clean feature table
  # Remove empty rows (features)
  feature_table2 <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1)
  
  # Remove columns (samples) with zero count
  if (ncol(feature_table2) > 1) {
    feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
  }
  
  if (isTRUE(strains)) {
    # Convert table with strain names to a strain-number table
    feature_table2 <- strain_name2strain_number(feature_table2)
  }
  
  # Save species names from row names
  species <- row.names(feature_table2)
  print(head(feature_table2))
  
  ### Step 2. If sorting, determine sample order.
  if (sort_type == "feature_value" && !is.null(feature_to_sort)) {
    print("Sort samples by feature_value")
    # Make "species" column with the rownames 
    df1 <- feature_table2 %>% rownames_to_column(var = "species")
    
    total_abundance <- colSums(df1[, -1])
    
    # Filter the row of the species of interest and calculate its proportion with respect to total abundance
    df_proportion <- df1 %>%
      dplyr::filter(species == feature_to_sort) %>%
      dplyr::select(-species)
    
    # calculate species of interest proportion
    df_proportion <- df_proportion[1, ] / total_abundance
    
    # Get sample names sorted by the species of interest proportion
    ordered_samples <- df_proportion %>%
      unlist() %>%
      sort(decreasing = TRUE) %>%
      names()
    
  } else if (sort_type == "similarity") {
    print("Sort samples by similarity")
    
    # transform table
    if (transform_table) {
      df1 <- transform_feature_table(feature_table = feature_table2, transform_method = "min_max")
    } else {
      df1 <- feature_table2
    }
    
    # Get the order of samples based on clustering
    ordered_samples <- order_samples_by_clustering(df1)
    
    df1 <- df1 %>% rownames_to_column(var = "species")
    
  } else if (sort_type == "none") {
    print("No sorting chosen")
    df1 <- feature_table2
    ordered_samples <- colnames(feature_table2)
    # Generate a column with the names of ASVs/OTUs using rownames.
    df1["species"] <- species
  } else {
    print("No valid sorting option chosen")
    return()
  }
  
  ### Step 3. Process features table to plotting table.
  plot_df <- df1 %>%
    tidyr::pivot_longer(-species, names_to = "sample", values_to = "abundance")
  
  # If strain processing has to be done.
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      dplyr::mutate(
        strain   = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)                   # Remove strain number from species name
      )
  }
  
  ### Step 4. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    dplyr::filter(!is.na(abundance) & abundance != 0)
  
  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      dplyr::filter(!is.na(strain) & strain != 0)
  }
  
  # Factor the "sample" variable so the order of samples is as in "ordered_samples"
  plot_df_filtered$sample <- factor(plot_df_filtered$sample, levels = ordered_samples)
  
  print(head(plot_df_filtered))
  
  ### Step 5. Create plot.
  if (is.null(colour_palette)) { # get colour palette
    print("Colour pallette generated")
    nfeatures <- length(unique(plot_df_filtered$species))
    colour_palette <- get_palette(nColors = nfeatures, replace_cols = replace_c)
    print(colour_palette)
  }
  
  # Create base plot.
  ft_barplot <- ggplot2::ggplot(plot_df_filtered, ggplot2::aes(x = sample, y = abundance, fill = species))
  
  if (isTRUE(strains)) {
    print("strains processing")
    ft_barplot <- ft_barplot +
      ggpattern::geom_bar_pattern(
        ggplot2::aes(fill = species2, pattern = strain, pattern_density = strain),
        position        = "fill",
        stat            = "identity",
        show.legend     = TRUE,
        pattern_color   = "white",
        pattern_fill    = "white",
        pattern_angle   = 45,
        pattern_spacing = 0.025
      ) +
      ggpattern::scale_pattern_manual(
        values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")
      ) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(
        pattern = guide_legend(override.aes = list(fill = "black")),
        fill    = guide_legend(override.aes = list(pattern = "none"))
      )
  } else {
    print("no strains")
    ft_barplot <- ft_barplot +
      ggplot2::geom_bar(
        ggplot2::aes(fill = species),
        position = ggplot2::position_fill(),
        stat     = "identity"
      )
  }
  
  # add theme options
  ft_barplot <- ft_barplot +
    ggplot2::theme_void() +
    ggplot2::scale_fill_manual(values = colour_palette) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = plot_title_size, face = "bold"),
      axis.title.x = ggplot2::element_text(size = x_axis_title_size),
      axis.text.x  = ggplot2::element_text(
        angle = x_axis_text_angle,
        vjust = x_vjust,
        hjust = x_hjust,
        size  = x_axis_text_size
      ),
      axis.title.y = ggplot2::element_text(
        size   = y_axis_title_size,
        angle  = y_axis_text_angle,
        margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)
      ),
      axis.text.y  = ggplot2::element_text(size = y_axis_text_size),
      legend.position = legend_pos,
      legend.title    = ggplot2::element_text(size = legend_title_size),
      legend.text     = ggplot2::element_text(size = legend_text_size)
    ) +
    guides(fill = guide_legend(ncol = legend_cols))
  
  ft_barplot
}


# ----- PCoA -----
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
    #ggplot2::coord_equal() +
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

# ----- Targeted metabolomics analyses -----
get_sample_info <- function(df, replicate_regex = "^[^_]+_\\d+$") {
  # 1) Find replicate columns like PREFIX_1, PREFIX_2, ...
  sample_cols <- grep(replicate_regex, colnames(df), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No replicate columns found. Check your column names and 'replicate_regex'.")
  }
  
  # 2) Extract base names (prefix before final underscore)
  base_names <- sub("_[^_]+$", "", sample_cols)
  
  # 3) Unique sample prefixes
  unique_samples <- unique(base_names)
  
  # 4) Return everything useful
  list(
    sample_cols   = sample_cols,
    base_names    = setNames(base_names, sample_cols), # named by column for clarity
    unique_samples = unique_samples
  )
}

build_mats_from_df <- function(df, sample_cols, base_names) {
  # ---- 0) Rownames for metabolites ----
  if (is.null(rownames(df))) {
    if ("Metabolite" %in% names(df)) {
      rn <- as.character(df$Metabolite)
    } else {
      stop("No rownames and no 'Metabolite' column found to set rownames.")
    }
  } else {
    rn <- rownames(df)
  }
  rn[is.na(rn) | rn == ""] <- "NA_metabolite"
  if (anyDuplicated(rn)) rn <- make.unique(rn)
  
  # ---- 1) Pull the replicate columns and coerce to numeric matrix ----
  if (length(sample_cols) == 0) stop("sample_cols is empty.")
  X <- df[, sample_cols, drop = FALSE]
  
  # guard against accidental non-numeric columns
  if (!all(vapply(X, is.numeric, logical(1)))) {
    X <- data.matrix(X)  # safe numeric coercion
  } else {
    X <- as.matrix(X)
  }
  rownames(X) <- rn
  
  # ---- 2) Compute per-sample (prefix) means across replicates ----
  # base_names is a named vector mapping each replicate col -> its prefix
  if (is.null(names(base_names))) {
    # try to align names with sample_cols if missing
    names(base_names) <- sample_cols
  }
  unique_samples <- unique(unname(base_names))
  
  mat_mean <- sapply(unique_samples, function(smpl) {
    cols <- names(base_names)[base_names == smpl]
    rowMeans(X[, cols, drop = FALSE], na.rm = TRUE)
  })
  mat_mean <- as.matrix(mat_mean)
  rownames(mat_mean) <- rownames(X)
  
  list(
    mat_raw = X,
    mat_mean = mat_mean,
    unique_samples = unique_samples
  )
}

compute_lfc_and_stars <- function(mat_raw,
                                  mat_mean,
                                  base_names,             # named vector: replicate_col -> prefix
                                  control_prefix = "CTRL",
                                  alpha = 0.05,           # FDR cutoff (BH)
                                  lfc_gate = 2,           # |log2 FC| threshold (2 = 4x)
                                  pseudocount_test = 1,   # added to replicate intensities before log2 for t-tests
                                  pseudocount_disp = 1e-8 # added to means before building heatmap LFC
) {
  # Checks
  stopifnot(!is.null(rownames(mat_raw)), !is.null(rownames(mat_mean)))
  unique_samples <- colnames(mat_mean)
  if (is.null(unique_samples)) stop("mat_mean must have column names (sample prefixes).")
  if (!control_prefix %in% unique_samples)
    stop(sprintf("control_prefix '%s' not found in mat_mean columns.", control_prefix))
  
  # --- 1) Heatmap LFC vs control (means across replicates, same as you plot) ---
  ctrl_mean <- mat_mean[, control_prefix, drop = TRUE]
  lfc_heatmap_full <- log2((mat_mean + pseudocount_disp) / (ctrl_mean + pseudocount_disp))
  # Drop control column for plotting/stars
  lfc <- lfc_heatmap_full[, setdiff(colnames(lfc_heatmap_full), control_prefix), drop = FALSE]
  
  # --- 2) Prepare stars matrix (same shape as lfc) ---
  stars <- matrix("",
                  nrow = nrow(lfc),
                  ncol = ncol(lfc),
                  dimnames = dimnames(lfc))
  
  # --- 3) T-test on replicate-level log2 values (as before) ---
  rep_cols  <- names(base_names)
  ctrl_cols <- rep_cols[base_names == control_prefix]
  log_ctrl  <- log2(mat_raw[, ctrl_cols, drop = FALSE] + pseudocount_test)
  
  for (smpl in colnames(lfc)) {
    trt_cols <- rep_cols[base_names == smpl]
    if (length(trt_cols) == 0) next
    
    log_trt <- log2(mat_raw[, trt_cols, drop = FALSE] + pseudocount_test)
    
    # p-values per metabolite
    pvals <- vapply(seq_len(nrow(mat_raw)), function(i) {
      x <- log_trt[i, ]
      y <- log_ctrl[i, ]
      if (all(is.na(x)) || all(is.na(y)) ||
          length(na.omit(x)) < 2 || length(na.omit(y)) < 2) return(NA_real_)
      tt <- try(t.test(x, y, var.equal = FALSE), silent = TRUE)
      if (inherits(tt, "try-error")) NA_real_ else tt$p.value
    }, numeric(1))
    
    padj <- p.adjust(pvals, method = "BH")
    
    # --- 4) Gate on the SAME effect that the heatmap shows: |LFC_heatmap| >= lfc_gate ---
    # Use the lfc from the heatmap matrix (already aligns by rownames/colnames)
    lfc_vec <- lfc[, smpl, drop = TRUE]  # log2((mean_trt+eps)/(mean_ctrl+eps))
    
    sig <- !is.na(padj) & (padj < alpha) & (abs(lfc_vec) >= lfc_gate)
    
    stars[, smpl] <- ifelse(sig, "*", "")
  }
  
  list(lfc = lfc, stars = stars)
}

plot_metabolites_lfc_panel <- function(df,
                                       metabolites,
                                       ctrl_prefix = "CTRL",
                                       n_rows = 2,
                                       n_cols = 3,
                                       replicate_regex = "^[^_]+_\\d+$",
                                       tiny_pseudocount = 0,
                                       y_limits = c(-10, 10),
                                       show_guides = TRUE,
                                       palette = NULL,
                                       facet_label_width = 28,
                                       debug = FALSE) {
  # --- identify replicate/sample columns ---
  sample_cols <- grep(replicate_regex, colnames(df), value = TRUE)
  if (!length(sample_cols)) stop("No replicate columns found.")
  unique_samples <- unique(sub("_[^_]+$", "", sample_cols))
  if (!ctrl_prefix %in% unique_samples) {
    stop(sprintf("CTRL prefix '%s' not among samples: %s",
                 ctrl_prefix, paste(unique_samples, collapse = ", ")))
  }
  samp_levels <- c(ctrl_prefix, setdiff(unique_samples, ctrl_prefix))
  
  # --- lift rownames & subset metabolites ---
  df2 <- as.data.frame(df[, sample_cols, drop = FALSE], stringsAsFactors = FALSE)
  df2$Metabolite <- rownames(df)
  metabolites <- as.character(metabolites)
  present <- intersect(metabolites, df2$Metabolite)
  if (!length(present)) stop("None of the requested metabolites found.")
  met_order <- present
  df2 <- dplyr::filter(df2, Metabolite %in% present)
  
  # --- long format ---
  long <- df2 |>
    tidyr::pivot_longer(cols = tidyselect::all_of(sample_cols),
                        names_to = "Replicate", values_to = "Abundance") |>
    dplyr::mutate(
      Sample     = sub("_[^_]+$", "", Replicate),
      Sample     = factor(Sample, levels = samp_levels),
      Metabolite = as.character(Metabolite),
      Abundance  = suppressWarnings(as.numeric(Abundance))
    )
  
  # --- CTRL means ---
  ctrl_means <- long |>
    dplyr::filter(Sample == ctrl_prefix) |>
    dplyr::group_by(Metabolite) |>
    dplyr::summarise(ctrl_mean = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  bad <- ctrl_means |>
    dplyr::filter(!is.finite(ctrl_mean) | (tiny_pseudocount == 0 & ctrl_mean == 0)) |>
    dplyr::pull(Metabolite)
  
  keep_mets <- setdiff(met_order, bad)
  if (!length(keep_mets)) stop("All requested metabolites had invalid CTRL means; nothing to plot.")
  
  # --- normalize to CTRL & compute log2FC ---
  long2 <- long |>
    dplyr::filter(Metabolite %in% keep_mets) |>
    dplyr::left_join(ctrl_means, by = "Metabolite") |>
    dplyr::mutate(
      Relative = if (tiny_pseudocount > 0)
        (Abundance + tiny_pseudocount) / (ctrl_mean + tiny_pseudocount)
      else
        Abundance / ctrl_mean,
      log2FC = log2(Relative)
    ) |>
    dplyr::filter(is.finite(log2FC))
  
  # remove CTRL from plot (keep legend for other samples)
  long2 <- dplyr::filter(long2, Sample != ctrl_prefix)
  samp_levels <- setdiff(samp_levels, ctrl_prefix)
  long2 <- dplyr::mutate(long2, Metabolite = factor(Metabolite, levels = met_order))
  
  # --- PLOT: now SHOW x-axis sample names and keep legend ---
  p <- ggplot2::ggplot(long2, ggplot2::aes(x = Sample, y = log2FC, fill = Sample)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(ggplot2::aes(color = Sample),
                         width = 0.2, size = 1.8, alpha = 0.85, show.legend = FALSE) +
    ggplot2::facet_wrap(
      ~ Metabolite, nrow = n_rows, ncol = n_cols, scales = "fixed", drop = FALSE,
      labeller = ggplot2::labeller(Metabolite = function(x) stringr::str_wrap(x, width = facet_label_width))
    ) +
    ggplot2::ggtitle(paste("log2 Fold-Change vs", ctrl_prefix)) +
    ggplot2::labs(x = NULL, y = expression(log[2]("sample / CTRL"))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      # SHOW x-axis labels (was hidden before)
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks.x = ggplot2::element_line(),
      legend.title = ggplot2::element_text(),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Sample"),
                    color = "none")
  
  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits) +
      ggplot2::scale_y_continuous(breaks = pretty(y_limits, n = 7))
  }
  if (show_guides) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "solid") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed")
  }
  
  # palette (named or unnamed)
  if (!is.null(palette)) {
    if (is.null(names(palette))) {
      if (length(palette) < length(samp_levels)) {
        stop(sprintf("Palette has %d colors but needs at least %d.",
                     length(palette), length(samp_levels)))
      }
      pal_vec <- setNames(palette[seq_along(samp_levels)], samp_levels)
    } else {
      missing_cols <- setdiff(samp_levels, names(palette))
      if (length(missing_cols)) {
        stop(sprintf("Palette missing colors for: %s", paste(missing_cols, collapse = ", ")))
      }
      pal_vec <- palette[samp_levels]
    }
    p <- p + ggplot2::scale_fill_manual(values = pal_vec) +
      ggplot2::scale_color_manual(values = pal_vec)
  }
  
  return(p)
}


# ----- limma markers analysis -----
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




# ----- Analysis of Human Microbiome Porject data -----
# This function takes a "biom_path" to a biom file with an otu_table and a tax_table,
# a string "tax_rank" which indicates the level of analyses, and bool "order_table".
# tax_rank parameter must be a value of greengenes ranks format; if not an error is returned.
# The ASVs/OTUs in the biom file are agglomerated by the "tax_rank" provided
# "order_table" indicates if the table should be ordered by larger to smaller values of rowMeans.
# Generally, ASV/OTU tables from QIIME2 are already ordered by row sums.
# This function returns a dataframe where rows are the ASVs and the columns are samples,
# "rownames" are ASVs taxonomy at the selected rank, and "colnames" are samples names.
# Taxonomy is dereplicated so that no row has the same name (which is not allowed in R dataframes).
# The output format is useful for using in other packages like vegan and to generate plots like barplots and heatmaps.
load_biom_as_table <- function(biom_path, tax_rank = "Species", strain_taxonomy = FALSE, order_table = FALSE){
  
  unite_colNames <- get_colNames_per_rank(tax_rank)
  
  if(strain_taxonomy) {
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_strain)
  }else{
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=phyloseq::parse_taxonomy_greengenes)
  }
  
  extracted_feature_table <- extract_table(biom_object, tax_rank, unite_colNames)
  
  return(clean_table(extracted_feature_table, order_table = order_table))
}

# This function takes a "tax_rank" string that correspond to a taxonomic rank in Greengenes format.
# Returns a list of strings which represent the columns in the tax_table of a biom file 
# that have to be joined to get the taxonomy assignment of each AVS/OTU as a string.
# If a not valid tax_rank is provided it returns an error.
get_colNames_per_rank <- function(tax_rank){
  colNames = NULL
  switch(tax_rank,
         Strain = {
           colNames = c("Genus", "Species", "Strain")
         },
         Species = {
           # Species level
           colNames = c("Genus", "Species")
         },
         Genus = {
           # Genus level
           colNames = c("Genus")
         },
         Family = {
           # Family
           colNames = c("Family")
         },
         Order = {
           # Order
           colNames = c("Order")}
  )
  if (!is.null(colNames)){
    return(colNames)
  }else{
    stop("Please choose a valid taxonomy rank!", call. = FALSE)
  }
}

# This function takes a character vector containing the result of splitting a taxonomy vector in the greenegenes format.
# It returns a named vector where each field is a taxonomic rank for the passed taxonomy entry.
# The taxonomic ranks are the same as in the greengenes taxonomy format but include a "Strain" rank.
# This function is used by phyloseq's "import_biom" function to parse taxonomy.
# import_biom splits taxonomy vectors automatically when they are the in the greengenes format.
parse_taxonomy_strain <- function(char.vec){
  # Remove the greengenes taxonomy rank id prefix.
  named.char.vec <- substring(char.vec, first = 4)
  # Set the names for each rank.
  names(named.char.vec) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  return(named.char.vec)
}

# This function takes a biom object and extracts it's "tax_table" and the "otu table".
# Then cbinds both dataframes to obtain a dataframe where the 1st column is the taxonomy and the 
extract_table <- function(biom_object, tax_rank, col_names){
  # Agglomerate tax_table by the chosen tax_rank
  biom_object <- phyloseq::tax_glom(biom_object, taxrank = tax_rank, NArm=TRUE)
  # cbind tax_table and otu_table
  feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(phyloseq::tax_table(biom_object)),
                                                    taxonomy,
                                                    all_of(col_names),
                                                    sep = "_"), "taxonomy"),
                         data.frame(phyloseq::otu_table(biom_object)))
}

# This function takes a feature table. First it make all the values in the "taxonomy" column unique.
# Then it makes the "taxonomy" column the rownames of the table.
# If "order_table" is TRUE it orders the table by ASVs/OTUs abundance.
clean_table <- function(feature_table, order_table){
  # Get valid (unique) names for all ASVs/OTUs.
  feature_table["taxonomy"] <- make.unique(feature_table$taxonomy, sep = "_")
  # Set taxonomy column as rownames
  feature_table <- tibble::column_to_rownames(tibble::remove_rownames(feature_table), var = "taxonomy")
  if (order_table) {
    # Order by abundances mean, from higher to lower.
    feature_table <- feature_table[order(rowMeans(feature_table), decreasing = TRUE),]
  }else{
    return(feature_table)
  }
}

filter_low_abundance <- function(rel_abundance, threshold = 0.01) {
  # rel_abundance: species x samples matrix/dataframe of relative abundances
  # threshold: minimum relative abundance (e.g., 0.01 = 1%)
  
  # Compute maximum abundance of each species across samples
  species_max <- apply(rel_abundance, 1, max)
  
  # Keep only species with max abundance >= threshold
  filtered_df <- rel_abundance[species_max >= threshold, ]
  
  message(paste("Filtered from", nrow(rel_abundance), 
                "to", nrow(filtered_df), "species"))
  
  return(filtered_df)
}

# ----- Within-timepoint replicate distances ----
# For replicate similarity
prepare_data <- function(abund, meta) {
  # Ensure column/rownames are present
  stopifnot(!is.null(colnames(abund)), !is.null(rownames(meta)))
  
  # intersect samples
  common_samples <- intersect(colnames(abund), rownames(meta))
  if (length(common_samples) == 0) stop("No overlapping sample IDs between abund colnames and meta rownames.")
  
  abund2 <- abund[, common_samples, drop = FALSE]
  meta2  <- meta[common_samples, , drop = FALSE]
  
  # Make a clean sample_id column
  meta2 <- meta2 %>%
    mutate(sample_id = rownames(meta2),
           syncom_id = .data$ATTRIBUTE_SynCom,
           time_raw  = .data$ATTRIBUTE_Time,
           rep_raw   = .data$ATTRIBUTE_Replicate)
  
  # Map times to numeric order 1..4 (T1,T2,T3,TF->4)
  time_map <- c(T1 = 1, T2 = 2, T3 = 3, TF = 4)
  if (!all(meta2$time_raw %in% names(time_map))) {
    bad_levels <- setdiff(unique(meta2$time_raw), names(time_map))
    stop("Unexpected time labels in ATTRIBUTE_Time: ", paste(bad_levels, collapse = ", "))
  }
  
  meta2 <- meta2 %>%
    mutate(
      time_num   = unname(time_map[time_raw]),
      time_label = factor(time_raw, levels = c("T1","T2","T3","TF"))
    )
  
  # Optionally normalize abundances per sample (just to be safe)
  # Transpose to samples x taxa because vegdist expects samples in rows
  X <- t(abund2)
  X[is.na(X)] <- 0
  row_sums <- rowSums(X)
  row_sums[row_sums == 0] <- 1
  X <- sweep(X, 1, row_sums, "/")
  
  list(meta = meta2, X = as.matrix(X))
}
# Compute within-timepoint replicate distances
compute_within_tp_distances <- function(meta, X, method = "bray") {
  dat <- meta %>% select(sample_id, syncom_id, time_num, time_label, rep_raw)
  
  # helper to compute pairwise distances for a block of sample_ids
  compute_block <- function(sample_ids) {
    if (length(sample_ids) < 2) return(tibble())
    d <- vegdist(X[sample_ids, , drop = FALSE], method = method)  # distances among rows (samples)
    dv <- as.numeric(d)
    
    # label pairs as e.g., "R1_vs_R2"
    reps <- meta$rep_raw[match(sample_ids, meta$sample_id)]
    pair_names <- combn(reps, 2, FUN = function(xx) paste0(xx[1], "_vs_", xx[2]))
    tibble(replicate_pair = pair_names, pairwise_dist = dv)
  }
  
  # group by SynCom × Time and compute the 3 pairwise distances
  dist_tbl <- dat %>%
    group_by(syncom_id, time_num, time_label) %>%
    group_modify(~{
      sample_ids <- .x$sample_id
      compute_block(sample_ids)
    }) %>%
    ungroup()
  
  # summarize with mean/sd (retaining individual pairs too)
  dist_tbl %>%
    group_by(syncom_id, time_num, time_label) %>%
    mutate(mean_dist = mean(pairwise_dist),
           sd_dist   = sd(pairwise_dist),
           n_pairs   = dplyr::n()) %>%
    ungroup()
}

# Plot: one panel per SynCom, x = time, y = dissimilarity
plot_replicate_similarity <- function(dist_tbl) {
  dist_tbl <- dist_tbl %>%
    mutate(
      syncom_order = as.numeric(gsub("\\D", "", syncom_id)),
      syncom_id = factor(syncom_id, levels = unique(syncom_id[order(syncom_order)]))
    )
  
  ggplot(dist_tbl, aes(x = time_num, y = pairwise_dist)) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.05, height = 0)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    facet_wrap(~ syncom_id, scales = "free_y") +
    scale_x_continuous(breaks = 1:4, labels = c("T1","T2","T3","TF")) +
    labs(
      x = "Time point",
      y = "Replicate dissimilarity (Bray–Curtis)",
      title = "Within–time point replicate dissimilarity per SynCom"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

compute_distance_to_final <- function(meta, X, method = "bray",
                                      mode = c("centroid", "replicate", "allpair_mean")) {
  mode <- match.arg(mode)
  stopifnot(nrow(meta) == nrow(X))
  
  # Coerce to tibble to avoid rowname weirdness
  meta <- tibble::as_tibble(meta)
  
  # Ensure required columns exist; if not, try to create them from ATTRIBUTE_* columns
  if (!"sample_id" %in% names(meta)) {
    stop("meta must contain a 'sample_id' column. Did you pass prepared$meta?")
  }
  if (!"syncom_id" %in% names(meta)) {
    if (all(c("ATTRIBUTE_SynCom") %in% names(meta))) {
      meta <- meta %>% mutate(syncom_id = .data$ATTRIBUTE_SynCom)
    } else {
      stop("meta lacks 'syncom_id' and 'ATTRIBUTE_SynCom'.")
    }
  }
  if (!"time_label" %in% names(meta) || !"time_num" %in% names(meta)) {
    if ("ATTRIBUTE_Time" %in% names(meta)) {
      time_map <- c(T1 = 1, T2 = 2, T3 = 3, TF = 4)
      meta <- meta %>%
        mutate(time_label = factor(.data$ATTRIBUTE_Time, levels = c("T1","T2","T3","TF")),
               time_num   = unname(time_map[as.character(.data$ATTRIBUTE_Time)]))
    } else {
      stop("meta lacks 'time_label/time_num' and 'ATTRIBUTE_Time'.")
    }
  }
  if (!"rep_raw" %in% names(meta)) {
    if ("ATTRIBUTE_Replicate" %in% names(meta)) {
      meta <- meta %>% mutate(rep_raw = .data$ATTRIBUTE_Replicate)
    } else {
      stop("meta lacks 'rep_raw' and 'ATTRIBUTE_Replicate'.")
    }
  }
  
  # Build the compact 'dat'
  dat <- meta %>%
    dplyr::select(sample_id, syncom_id, time_label, time_num, rep_raw)
  
  # Helper to compute Bray–Curtis between a sample row and a reference vector
  bc_to_vec <- function(sample_row, ref_vec) {
    d <- vegan::vegdist(rbind(sample_row, ref_vec), method = method)
    as.numeric(d)[1]
  }
  
  print(head(dat))
  
  results <- dat %>%
    dplyr::group_by(syncom_id) %>%
    dplyr::group_modify(~{
      block_meta <- .x
      ids <- block_meta$sample_id
      tf_ids <- block_meta$sample_id[block_meta$time_label == "TF"]
      
      if (length(tf_ids) == 0) {
        return(tibble::tibble(sample_id = ids, dist_to_TF = NA_real_))
      }
      
      if (mode == "centroid") {
        tf_centroid <- colMeans(X[tf_ids, , drop = FALSE])
        tibble::tibble(
          sample_id = ids,
          dist_to_TF = apply(X[ids, , drop = FALSE], 1, function(r) {
            d <- vegan::vegdist(rbind(r, tf_centroid), method = method)
            as.numeric(d)[1]
          })
        )
        
      } else if (mode == "replicate") {
        tf_centroid <- colMeans(X[tf_ids, , drop = FALSE])
        tf_rep_map <- tibble::tibble(
          tf_id  = tf_ids,
          tf_rep = block_meta$rep_raw[match(tf_ids, block_meta$sample_id)]
        )
        purrr::map_dfr(ids, function(sid) {
          rep_s <- block_meta$rep_raw[block_meta$sample_id == sid]
          tf_match <- tf_rep_map$tf_id[tf_rep_map$tf_rep == rep_s]
          dval <- if (length(tf_match) >= 1) {
            if (length(tf_match) == 1) {
              as.numeric(vegan::vegdist(rbind(X[sid, ], X[tf_match, ]), method = method))[1]
            } else {
              mean(as.numeric(vegan::vegdist(rbind(X[sid, ], X[tf_match, ]), method = method))[seq_along(tf_match)])
            }
          } else {
            as.numeric(vegan::vegdist(rbind(X[sid, ], tf_centroid), method = method))[1]
          }
          tibble::tibble(sample_id = sid, dist_to_TF = dval)
        })
        
      } else { # mode == "allpair_mean"
        purrr::map_dfr(ids, function(sid) {
          d <- vegan::vegdist(rbind(X[sid, ], X[tf_ids, ]), method = method)
          tibble::tibble(sample_id = sid,
                         dist_to_TF = mean(as.numeric(d)[seq_along(tf_ids)]))
        })
      }
    }) %>%
    dplyr::ungroup() %>%
    # <<< drop the group key that group_modify appended
    dplyr::select(sample_id, dist_to_TF) %>%
    # reattach a single clean copy of metadata
    dplyr::left_join(dat, by = "sample_id")
  
  summary_tbl <- results %>%
    dplyr::group_by(syncom_id, time_label, time_num) %>%
    dplyr::summarize(
      mean_dist_to_TF = mean(dist_to_TF, na.rm = TRUE),
      sd_dist_to_TF   = sd(dist_to_TF, na.rm = TRUE),
      n               = sum(!is.na(dist_to_TF)),
      .groups = "drop"
    )
  
  
  list(per_sample = results, summary = summary_tbl)
}

plot_distance_to_final <- function(per_sample, summary_tbl) {
  # Consistent numeric ordering of SynCom panels
  order_levels <- summary_tbl %>%
    mutate(syncom_order = as.numeric(gsub("\\D", "", syncom_id))) %>%
    arrange(syncom_order) %>%
    pull(syncom_id) %>%
    unique()
  
  per_sample <- per_sample %>%
    mutate(syncom_id = factor(syncom_id, levels = order_levels))
  summary_tbl <- summary_tbl %>%
    mutate(syncom_id = factor(syncom_id, levels = order_levels))
  
  ggplot() +
    geom_point(
      data = per_sample,
      aes(x = time_num, y = dist_to_TF),
      alpha = 0.6, position = position_jitter(width = 0.05, height = 0)
    ) +
    geom_line(
      data = summary_tbl,
      aes(x = time_num, y = mean_dist_to_TF, group = 1),
      linewidth = 0.8
    ) +
    geom_point(
      data = summary_tbl,
      aes(x = time_num, y = mean_dist_to_TF),
      size = 2
    ) +
    facet_wrap(~ syncom_id, scales = "free_y") +
    scale_x_continuous(breaks = 1:4, labels = c("T1","T2","T3","TF")) +
    labs(
      x = "Time point",
      y = "Bray–Curtis distance to final state (TF)",
      title = "Stabilization toward final community composition"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


