source("C:/Users/marce/Documents/GitHub/nasalSynComs/helper_functions.R")

# ---------- Screening Results ----------

otu_table_screening <- read.csv("D:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1, sep = ";")


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "SC18", "SC19", "SC20",
                                   "SC21", "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30",
                                   "SC31", "SC32", "SC33", "SC34", "SC35", "SC36", "SC37", "SC38", "SC39", "SC40",
                                   "SC41", "SC42", "SC43", "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")


cluster_barplot_result <- cluster_barplot_panels(ot_scree_filtered, colour_palette = colours_vec)


# ---------- Selected SynComs ----------
metadata <- read_metadata("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/SC100_metadata_noqcs_nosinStrs.csv",
                          sort_table = TRUE)
metadata <- metadata[7:nrow(metadata),]
metadata_or_names <- rownames(metadata)
rownames(metadata) <- gsub("\\.mzML$", "", rownames(metadata))

# Add cluster result to metadata of Selected SynComs

meta_df <- add_cluster_column(
  meta_df = metadata,
  clusters_df = cluster_barplot_result$clusters,
  meta_key_col      = "ATTRIBUTE_SynCom",
  cluster_key_col   = "Sample",
  cluster_value_col = "Cluster",
  new_col_name      = "ATTRIBUTE_Cluster"
)

# Barplot with strain-level information for C. propinquum and D. pigrum

# ---------- Metabolites PCoA  ----------
feature_table_tic <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_quantTable_ticNorm2.csv",
                             sort_by_names = TRUE, p_sep = ";")

# Sort feature table by sample names
feature_table_tic <- feature_table_tic[, order(colnames(feature_table_tic))]

###
res_euc <- pcoa_flex(
  metab_df      = feature_table_tic,
  metadata_df   = meta_df,
  color_var     = "ATTRIBUTE_SynCom",
  shape_var     = "ATTRIBUTE_Cluster",
  ellipse_var   = "ATTRIBUTE_Cluster",
  color_var_leg_columns = 4,
  distance      = "bray",
  preprocess    = "hellinger",
  permanova_var = "ATTRIBUTE_Cluster",
  permutations  = 999
)

print(res_euc$plot)
res_euc$permanova

# ---------- Metabolites + Clusters Markers Heatmap  ----------
an_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/2025-05-10_Merged_Annotations_GNPS_SIRIUS (1).csv", row.names=1)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(limma)
})

res_limma <- limma_markers_by_cluster_general(
  metab_df      = feature_table_tic,   # ~6000 x ~200 (non-negative)
  metadata_df   = meta_df,               # has ATTRIBUTE_* columns + Sample
  #  sample_id_col = "Sample",
  cluster_var   = "ATTRIBUTE_Cluster",
  covariates    = c("ATTRIBUTE_Time"),       # optional; drop or add more if you like
  block_var     = "ATTRIBUTE_SynCom",        # optional; recommended for repeated measures
  log_transform = TRUE, log_offset = 1,
  do_pairwise   = TRUE
)

sum_ht_sirius <- summarize_markers_and_heatmap_with_classes(
  metab_df      = feature_table_tic,
  metadata_df   = meta_df,
  sample_id_col = "Sample",
  cluster_var   = "ATTRIBUTE_Cluster",
  sirius_df     = an_table,
  id_col        = "row.ID",
  class_cols    = c("SIRIUS_ClassyFire.most.specific.class",
                    "SIRIUS_ClassyFire.subclass",
                    "SIRIUS_ClassyFire.level.5"),
  id_pattern    = "^X(\\d+).*",
  limma_res     = res_limma,
  top_n = 25, p_adj_thresh = 0.05, min_logFC = 0,
  log_transform = TRUE, log_offset = 1,
  scale_rows    = TRUE,
  out_file      = file.path("C:/Users/marce/Desktop/markers_heatmap2.pdf"),  # <- save here
  out_width     = 15,
  out_height    = 12,
  class_na_label = "Unclassified",
  class_na_color = "#BDBDBD",
  c_legend_ncol = 2,
  r_legend_ncol = 4,
  legend_side = "bottom"   # "bottom", "top", "left", or "right"
)


# c("SIRIUS_ClassyFire.class",
#  "SIRIUS_ClassyFire.most.specific.class",
#  "SIRIUS_ClassyFire.subclass",
#  "SIRIUS_ClassyFire.level.5")


# ---------- Repetition Experiment (Diversity and Aminoacids)  ----------






#ot_scree_filtered_rel_ab <- transform_feature_table(feature_table = ot_scree_filtered,
#                                                   transform_method = "rel_abundance")
#result <- cluster_barplot_panels(ot_scree_filtered_rel_ab, colour_palette = colours_vec)