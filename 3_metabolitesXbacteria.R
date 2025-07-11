# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_graphing.R")


# Load OTU tables
# Species-level
otu_table <- read.csv(file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/2_20sc_species_ot.csv",
                      row.names = 1, header = TRUE)
# Strain-level
otu_table <- read.csv(file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/3_20sc_strain_ot.csv",
                      row.names = 1, header = TRUE)

# Filter and sort OTU table
otu_table <- remove_feature_by_prefix(otu_table, c("Corynebacterium tuberculostearicum DSM44922"))
# Filter and sort OTU table
otu_table <- otu_table[, order(colnames(otu_table))]

# Load feature tables for metabolomics data
# Feature table - only imputated
feature_table_imputated <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_quantTable_imp2.csv",
                                   sort_by_names = TRUE, p_sep = ";")
# Filter and sort feature table
feature_table_imputated <- feature_table_imputated[, order(colnames(feature_table_imputated))]

# Feature table - imputated and scaled
feature_table_scaled <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_QuantTable_scaled2.csv",
                                sort_by_names = TRUE, p_sep = ";")
# Filter and sort feature table
feature_table_scaled <- feature_table_scaled[, order(colnames(feature_table_scaled))]

# Feature table - random forest-selected metabolites
feature_table_rf <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled2.csv",
                             sort_by_names = TRUE, p_sep = ";")
# Filter and sort feature table
feature_table_rf <- feature_table_rf[, order(colnames(feature_table_rf))]

# Read metadata
metadata <- read_metadata("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/SC100_metadata_noqcs_nosinStrs.csv",
                          sort_table = TRUE)
metadata <- metadata[7:nrow(metadata),]
rownames(metadata) <- gsub("\\.mzML$", "", rownames(metadata))

# Exploratory Heatmap
metadata2 <- flag_samples_by_abundance(feature_df = otu_table, metadata_df = metadata,
                          feature_name = "Staphylococcus aureus USA300",
                          percentage_threshold = 25)



### Correlation Heatmaps - Bacteria and metabolites selected by random forest
colnames(otu_table) == colnames(feature_table_rf)
feature_table_heatmap(ft1 = otu_table, ft2 = feature_table_rf, corr_type = "pearson", sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr", text_angle = 45) # can use also "by" correction

### PCAs metabolites


# PCA metabolome
ft_pca_1(ft = feature_table_imputated, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = FALSE, dist_method = "bray")

# PCA with vegan
ft_pca_2(ft = feature_table_imputated, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", dist_method = "bray")

# Remove highly variable metabolites (intrareplicate variability).
#filtered_ft <- filter_by_error(feature_table = feature_table, metadata_table = metadata, grouping_var = "ATTRIBUTE_Sample", error_threshold = 25)
#rownames(filtered_ft) <- rownames(feature_table)


################################## Sparce Canonical Correspondance Analysis ##################################
install.packages("PMA")  # if not already installed
library(PMA)

X <- as.data.frame(t(otu_table))

Y <- as.data.frame(t(feature_table_scaled))

apply(X, 2, sd)  # Species
apply(Y, 2, sd)  # Metabolites

# Find shared samples
common_samples <- intersect(rownames(X), rownames(Y))

# Subset and reorder both matrices by the common sample names
X <- X[common_samples, , drop = FALSE]
Y <- Y[common_samples, , drop = FALSE]

# Check final alignment
all(rownames(X) == rownames(Y))  # should return TRUE

# Lets scale species too
X_scaled <- scale(X)
apply(X_scaled, 2, sd)  # Species


scca_result <- CCA(x = X_scaled, z = Y, typex = "standard", typez = "standard",
                   penaltyx = 0.4, penaltyz = 0.3, standardize = TRUE, K = 1)

# Canonical correlation
scca_result$cors

# Selected species (non-zero u)
selected_species <- colnames(X_scaled)[scca_result$u[,1] != 0]
selected_species

# Selected metabolites (non-zero v)
selected_metabs <- colnames(Y)[scca_result$v != 0]
selected_metabs

U_scores <- X_scaled %*% scca_result$u
V_scores <- as.matrix(Y) %*% scca_result$v

# Combine into a data frame
df_scores <- data.frame(
  Sample = rownames(X_scaled),
  Species_CV = U_scores,
  Metabolite_CV = V_scores
)

# Scatter plot
ggplot(df_scores, aes(x = Species_CV, y = Metabolite_CV)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, size = 3) +
  labs(x = "Canonical Variable 1 (Species)",
       y = "Canonical Variable 1 (Metabolites)",
       title = "sCCA Sample Projection") +
  theme_minimal()


# Top loadings
top_species <- sort(scca_result$u, decreasing = TRUE)[1:15]
top_metabs  <- sort(scca_result$v, decreasing = TRUE)[1:15]

names(scca_result$u) <- colnames(X_scaled)
names(scca_result$v) <- colnames(Y)

head(scca_result$u)
head(scca_result$v)
# Add names
names(top_species) <- names(scca_result$u)[order(scca_result$u, decreasing = TRUE)][1:15]
names(top_metabs)  <- names(scca_result$v)[order(scca_result$v, decreasing = TRUE)][1:15]

# Now plot
barplot(top_species, las = 2, main = "Top 10 Species Loadings", col = "steelblue")
barplot(top_metabs, las = 2, main = "Top 10 Metabolite Loadings", col = "darkgreen")

######################################################################
######################################################################
########################################## Heatmaps, only metaboloites

plot_heatmap_by_time <- function(feature_df, metadata_df, time_point, syncom_colors = NULL) {
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  
  # 1. Check sample name alignment
  if (!identical(colnames(feature_df), rownames(metadata_df))) {
    stop("Sample names in feature_df and metadata_df do not match.")
  }
  
  # 2. Filter metadata for the given time point
  filtered_meta <- metadata_df %>% filter(Time == time_point)
  
  if (nrow(filtered_meta) == 0) {
    stop("No samples found for the given time point.")
  }
  
  # 3. Subset the feature dataframe
  selected_samples <- rownames(filtered_meta)
  filtered_features <- feature_df[, selected_samples, drop = FALSE]
  
  # 4. Order samples by SynCom
  filtered_meta$SynCom <- as.character(filtered_meta$SynCom)
  ordered_meta <- filtered_meta[order(filtered_meta$SynCom), ]
  ordered_samples <- rownames(ordered_meta)
  ordered_features <- filtered_features[, ordered_samples]
  
  # 5. Scale the data by row (z-score)
  scaled_matrix <- t(scale(t(ordered_features)))  # Row-wise scaling
  
  # 6. Prepare SynCom colors
  syncoms <- ordered_meta$SynCom
  if (is.null(syncom_colors)) {
    unique_syncoms <- unique(syncoms)
    syncom_colors <- setNames(rainbow(length(unique_syncoms)), unique_syncoms)
  }
  
  # 7. Create annotation
  top_annotation <- HeatmapAnnotation(
    SynCom = syncoms,
    col = list(SynCom = syncom_colors),
    annotation_name_side = "left"
  )
  
  # 8. Plot heatmap
  Heatmap(
    scaled_matrix,
    name = "z-score",
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    column_title = paste("Samples at Time =", time_point),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
}
