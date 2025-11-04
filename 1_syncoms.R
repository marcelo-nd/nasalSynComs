source("C:/Users/marce/Documents/GitHub/nasalSynComs/helper_functions.R")
library(dplyr)

# ---------- Figure 2. Screening Results with Add strain level info----------
# Get data
otu_table_screening <- read.csv("D:/SequencingData/SynCom100/1_Screening/emu_results/otu_table.csv", row.names=1, sep = ";")


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "SC18", "SC19", "SC20",
                                   "SC21", "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30",
                                   "SC31", "SC32", "SC33", "SC34", "SC35", "SC36", "SC37", "SC38", "SC39", "SC40",
                                   "SC41", "SC42", "SC43", "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)

### Strain data processing
strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data/nasal_syncom_strains.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)

strain_data <- tibble::column_to_rownames(strain_data, "Species")

strain_data <- remove_feature_by_prefix(strain_data, species_to_remove)

strain_data <- tibble::rownames_to_column(strain_data, "Species")

# To use strain-level data
strain_ft <- merge_abundance_by_strain(ot_scree_filtered, strain_data)

otu_table <- strain_ft

##### Run only for creating barplots with strain-level data for certain species.
otu_table <- merge_non_target_strains(otu_table, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

strain_level_sel = TRUE

clustering_results <- cluster_samples(ot_scree_filtered)

clusters <- clustering_results$clusters
rel_abundance_ordered <- clustering_results$rel_abundance_ordered
sample_order <- clustering_results$sample_order
k <- clustering_results$best_k

cluster_barplot_result <- cluster_barplot_panels(abundance_df = rel_abundance_ordered,
                                                 cluster_df = clusters,
                                                 sample_order = sample_order,
                                                 best_k = k,
                                                 strains = FALSE,
                                                 colour_palette = colours_vec)

cluster_barplot_result <- cluster_barplot_panels(abundance_df = calculate_relative_abundance(otu_table),
                                                 cluster_df = clusters,
                                                 sample_order = sample_order,
                                                 best_k = k,
                                                 strains = TRUE,
                                                 colour_palette = get_palette(20))

df_test <- cluster_barplot_result$df_long

#barplots1 <- barplots1 + xlab("Time") + # for the x axis label
#  ylab("Relative abundance")

# ---------- Figure 3. Selected SynComs Barplots ----------
# Barplot with strain-level information for C. propinquum and D. pigrum
###### Time-series analyses
otu_table_sctp <- read.csv("D:/SequencingData/SynCom100/2_TheChampions/emu_results/otu_table.csv",
                           row.names=1, sep = ";")

otu_table_sctp_sorted <- sort_nanopore_table_by_barcodes(df = otu_table_sctp,
                                                         new_names = c("SC4_T1_R1", "SC4_T1_R2", "SC4_T1_R3",
                                                                       "SC4_T2_R1", "SC4_T2_R2", "SC4_T2_R3",
                                                                       "SC4_T3_R1", "SC4_T3_R2", "SC4_T3_R3",
                                                                       "SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                                                                       "SC7_T1_R1", "SC7_T1_R2", "SC7_T1_R3",
                                                                       "SC7_T2_R1", "SC7_T2_R2", "SC7_T2_R3",
                                                                       "SC7_T3_R1", "SC7_T3_R2", "SC7_T3_R3",
                                                                       "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                                                                       "SC9_T1_R1", "SC9_T1_R2", "SC9_T1_R3",
                                                                       "SC9_T2_R1", "SC9_T2_R2", "SC9_T2_R3",
                                                                       "SC9_T3_R1", "SC9_T3_R2", "SC9_T3_R3",
                                                                       "SC9_TF_R1", "SC9_TF_R2", "SC9_TF_R3",
                                                                       "SC10_T1_R1", "SC10_T1_R2", "SC10_T1_R3",
                                                                       "SC10_T2_R1", "SC10_T2_R2", "SC10_T2_R3",
                                                                       "SC10_T3_R1", "SC10_T3_R2", "SC10_T3_R3",
                                                                       "SC10_TF_R1", "SC10_TF_R2", "SC10_TF_R3",
                                                                       "SC11_T1_R1", "SC11_T1_R2", "SC11_T1_R3",
                                                                       "SC11_T2_R1", "SC11_T2_R2", "SC11_T2_R3",
                                                                       "SC11_T3_R1", "SC11_T3_R2", "SC11_T3_R3",
                                                                       "SC11_TF_R1", "SC11_TF_R2", "SC11_TF_R3",
                                                                       "SC12_T1_R1", "SC12_T1_R2", "SC12_T1_R3",
                                                                       "SC12_T2_R1", "SC12_T2_R2", "SC12_T2_R3",
                                                                       "SC12_T3_R1", "SC12_T3_R2", "SC12_T3_R3",
                                                                       "SC12_TF_R1", "SC12_TF_R2", "SC12_TF_R3",
                                                                       "SC13_T1_R1", "SC13_T1_R2", "SC13_T1_R3",
                                                                       "SC13_T2_R1", "SC13_T2_R2", "SC13_T2_R3",
                                                                       "SC13_T3_R1", "SC13_T3_R2", "SC13_T3_R3",
                                                                       "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                                                                       "SC14_T1_R1", "SC14_T1_R2", "SC14_T1_R3",
                                                                       "SC14_T2_R1", "SC14_T2_R2", "SC14_T2_R3",
                                                                       "SC14_T3_R1", "SC14_T3_R2", "SC14_T3_R3",
                                                                       "SC14_TF_R1", "SC14_TF_R2", "SC14_TF_R3",
                                                                       "SC19_T1_R1", "SC19_T1_R2", "SC19_T1_R3",
                                                                       "SC19_T2_R1", "SC19_T2_R2", "SC19_T2_R3",
                                                                       "SC19_T3_R1", "SC19_T3_R2", "SC19_T3_R3",
                                                                       "SC19_TF_R1", "SC19_TF_R2", "SC19_TF_R3",
                                                                       "SC22_T1_R1", "SC22_T1_R2", "SC22_T1_R3",
                                                                       "SC22_T2_R1", "SC22_T2_R2", "SC22_T2_R3",
                                                                       "SC22_T3_R1", "SC22_T3_R2", "SC22_T3_R3",
                                                                       "SC22_TF_R1", "SC22_TF_R2", "SC22_TF_R3",
                                                                       "SC23_T1_R1", "SC23_T1_R2", "SC23_T1_R3",
                                                                       "SC23_T2_R1", "SC23_T2_R2", "SC23_T2_R3",
                                                                       "SC23_T3_R1", "SC23_T3_R2", "SC23_T3_R3",
                                                                       "SC23_TF_R1", "SC23_TF_R2", "SC23_TF_R3",
                                                                       "SC24_T1_R1", "SC24_T1_R2", "SC24_T1_R3",
                                                                       "SC24_T2_R1", "SC24_T2_R2", "SC24_T2_R3",
                                                                       "SC24_T3_R1", "SC24_T3_R2", "SC24_T3_R3",
                                                                       "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                                                                       "SC25_T1_R1", "SC25_T1_R2", "SC25_T1_R3",
                                                                       "SC25_T2_R1", "SC25_T2_R2", "SC25_T2_R3",
                                                                       "SC25_T3_R1", "SC25_T3_R2", "SC25_T3_R3",
                                                                       "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                                                                       "SC27_T1_R1", "SC27_T1_R2", "SC27_T1_R3",
                                                                       "SC27_T2_R1", "SC27_T2_R2", "SC27_T2_R3",
                                                                       "SC27_T3_R1", "SC27_T3_R2", "SC27_T3_R3",
                                                                       "SC27_TF_R1", "SC27_TF_R2", "SC27_TF_R3",
                                                                       "SC31_T1_R1", "SC31_T1_R2", "SC31_T1_R3",
                                                                       "SC31_T2_R1", "SC31_T2_R2", "SC31_T2_R3",
                                                                       "SC31_T3_R1", "SC31_T3_R2", "SC31_T3_R3",
                                                                       "SC31_TF_R1", "SC31_TF_R2", "SC31_TF_R3",
                                                                       "SC34_T1_R1", "SC34_T1_R2", "SC34_T1_R3",
                                                                       "SC34_T2_R1", "SC34_T2_R2", "SC34_T2_R3",
                                                                       "SC34_T3_R1", "SC34_T3_R2", "SC34_T3_R3",
                                                                       "SC34_TF_R1", "SC34_TF_R2", "SC34_TF_R3",
                                                                       "SC39_T1_R1", "SC39_T1_R2", "SC39_T1_R3",
                                                                       "SC39_T2_R1", "SC39_T2_R2", "SC39_T2_R3",
                                                                       "SC39_T3_R1", "SC39_T3_R2", "SC39_T3_R3",
                                                                       "SC39_TF_R1", "SC39_TF_R2", "SC39_TF_R3",
                                                                       "SC40_T1_R1", "SC40_T1_R2", "SC40_T1_R3",
                                                                       "SC40_T2_R1", "SC40_T2_R2", "SC40_T2_R3",
                                                                       "SC40_T3_R1", "SC40_T3_R2", "SC40_T3_R3",
                                                                       "SC40_TF_R1", "SC40_TF_R2", "SC40_TF_R3",
                                                                       "SC44_T1_R1", "SC44_T1_R2", "SC44_T1_R3",
                                                                       "SC44_T2_R1", "SC44_T2_R2", "SC44_T2_R3",
                                                                       "SC44_T3_R1", "SC44_T3_R2", "SC44_T3_R3",
                                                                       "SC44_TF_R1", "SC44_TF_R2", "SC44_TF_R3",
                                                                       "SC50_T1_R1", "SC50_T1_R2", "SC50_T1_R3",
                                                                       "SC50_T2_R1", "SC50_T2_R2", "SC50_T2_R3",
                                                                       "SC50_T3_R1", "SC50_T3_R2", "SC50_T3_R3",
                                                                       "SC50_TF_R1", "SC50_TF_R2", "SC50_TF_R3"))

# Remove species with no counts
otu_table_sctp_filt <- filter_features_by_col_counts(otu_table_sctp_sorted,
                                                     min_count = 10,
                                                     col_number = 1)

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

otu_table_sctp_filt <- remove_feature_by_prefix(otu_table_sctp_filt, species_to_remove)

# For inoculum with out strain-level data
inoculum_spp_df <- strain_data %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
  ungroup()

inoc_spps <- inoculum_spp_df$Species

inoculum_spp_df <- select(inoculum_spp_df, -1)

rownames(inoculum_spp_df) <- inoc_spps


### Run to include inoculation to barplot
strain_data2 <- as.data.frame(strain_data)

strain_data2 <- strain_data2[,3:ncol(strain_data2)]

rownames(strain_data2) <- strain_data$Species


# To use strain-level data
strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)

otu_table <- strain_ft

##### Run only for creating barplots with strain-level data for certain species.
otu_table <- merge_non_target_strains(otu_table, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

### If inoculation included and strain-level data for certain species is going to be used.
strain_data2 <- zero_out_species_in_samples(df = strain_data2, species_name = "Staphylococcus aureus 1", sample_names = colnames(strain_data2))

strain_data2 <- merge_non_target_strains(strain_data2, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

time_names <- c("Inoc", "T1", "T2", "T3", "T4")

sc4 <- cbind(strain_data2["SC4"], otu_table[c(2,5,8,11)])
colnames(sc4) <- time_names
sc7 <- cbind(strain_data2["SC7"], otu_table[c(14,17,20,23)])
colnames(sc7) <- time_names
sc9 <- cbind(strain_data2["SC9"], otu_table[c(26,29,32,35)])
colnames(sc9) <- time_names
sc10 <- cbind(strain_data2["SC10"], otu_table[c(38,41,44,47)])
colnames(sc10) <- time_names
sc11 <- cbind(strain_data2["SC11"], otu_table[c(50,53,56,59)])
colnames(sc11) <- time_names
sc12 <- cbind(strain_data2["SC12"], otu_table[c(62,65,68,71)])
colnames(sc12) <- time_names
sc13 <- cbind(strain_data2["SC13"], otu_table[c(74,77,80,83)])
colnames(sc13) <- time_names
sc14 <- cbind(strain_data2["SC14"], otu_table[c(86,89,92,95)])
colnames(sc14) <- time_names
sc19 <- cbind(strain_data2["SC19"], otu_table[c(98,101,104,107)])
colnames(sc19) <- time_names
sc22 <- cbind(strain_data2["SC22"], otu_table[c(110,113,116,119)])
colnames(sc22) <- time_names
sc23 <- cbind(strain_data2["SC23"], otu_table[c(122,125,128,131)])
colnames(sc23) <- time_names
sc24 <- cbind(strain_data2["SC24"], otu_table[c(134,137,140,143)])
colnames(sc24) <- time_names
sc25 <- cbind(strain_data2["SC25"], otu_table[c(146,149,152,155)])
colnames(sc25) <- time_names
sc27 <- cbind(strain_data2["SC27"], otu_table[c(158,161,164,167)])
colnames(sc27) <- time_names
sc31 <- cbind(strain_data2["SC31"], otu_table[c(170,173,176,179)])
colnames(sc31) <- time_names
sc34 <- cbind(strain_data2["SC34"], otu_table[c(182,185,188,191)])
colnames(sc34) <- time_names
sc39 <- cbind(strain_data2["SC39"], otu_table[c(194,197,200,203)])
colnames(sc39) <- time_names
sc40 <- cbind(strain_data2["SC40"], otu_table[c(206,209,212,215)])
colnames(sc40) <- time_names
sc44 <- cbind(strain_data2["SC44"], otu_table[c(218,221,224,227)])
colnames(sc44) <- time_names
sc50 <- cbind(strain_data2["SC50"], otu_table[c(230,233,236,239)])
colnames(sc50) <- time_names

### Barplots

strain_level_sel = FALSE
strain_level_sel = TRUE

barplots1 <- barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                                 sc12, sc13, sc14,sc19,sc22),
                           strains = strain_level_sel, shared_samples = FALSE,
                           experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11",
                                                 "SC12", "SC13", "SC14", "SC19","SC22"),
                           x_axis_title_size = 12, x_axis_text_size = 12,
                           y_axis_title_size = 12, y_axis_text_size = 12,
                           legend_pos = "none", legend_cols = 2,
                           legend_title_size = 12, legend_text_size = 12,
                           legend_key_size = 0.3, colour_palette = colours_vec)

barplots1 <- barplots1 + xlab("Time") + # for the x axis label
  ylab("Relative abundance")

barplots1

barplots2 <- barplots_grid(feature_tables = list(sc23, sc24, sc25, sc27, sc31,
                                                 sc34, sc39, sc40, sc44, sc50),
                           strains = strain_level_sel, shared_samples = FALSE,
                           experiments_names = c("SC23", "SC24", "SC25", "SC27", "SC31",
                                                 "SC34", "SC39", "SC40", "SC44","SC50"),
                           x_axis_title_size = 12, x_axis_text_size = 12,
                           y_axis_title_size = 12, y_axis_text_size = 12,
                           legend_pos = "bottom", legend_cols = 3,
                           legend_title_size = 12, legend_text_size = 12,
                           legend_key_size = 0.3, colour_palette = colours_vec)

barplots2 <- barplots2 + xlab("Time") + # for the x axis label
  ylab("Relative abundance") + labs(fill = "Species")

barplots2


barplots <- cowplot::plot_grid(barplots1, barplots2,
                               align = "v",
                               ncol = 1,
                               rel_heights = c(46/100, 54/100))

barplots

# ---------- Figure 4. Bacterial diversity and Metabolites PCoA  ----------
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

syncom_pallette <- c(get_palette(20))
syncom_pallette <- c("indianred1", "#6279B8", "lavenderblush3", "#DA6A00",
                     "#738564", "purple4", "#56B4E9", "indianred4",
                     "#1a3a46", "hotpink4", "honeydew1", "hotpink",
                     "cyan3", "#cd541d", "#009E73", "#EC9704",
                     "#502F4C", "#FFBA49", "ivory3", "#9C4A1A")

clusters_pallete <- c(get_palette(3))
clusters_pallete <- c("#583E26", "#F7C815", "lawngreen")

# PCoA Bacteria
res_euc <- pcoa_flex(
  metab_df      = otu_table_sctp_filt,
  metadata_df   = meta_df,
  color_var     = "ATTRIBUTE_SynCom",
  shape_var     = "ATTRIBUTE_Cluster",
  ellipse_var   = "ATTRIBUTE_Cluster",
  color_var_leg_columns = 3,
  distance      = "bray",
  preprocess    = "hellinger",
  permanova_var = "ATTRIBUTE_Cluster",
  permutations  = 999,
  points_palette = syncom_pallette,
  ellipse_palette = clusters_pallete
)

print(res_euc$plot)
res_euc$permanova


# Get untargeted metabolomics data
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
  color_var_leg_columns = 3,
  distance      = "bray",
  preprocess    = "hellinger",
  permanova_var = "ATTRIBUTE_Cluster",
  permutations  = 999,
  points_palette = syncom_pallette,
  ellipse_palette = clusters_pallete
)

plot(res_euc$plot)

res_euc$permanova


# ---------- Figure 5. Targeted Metabolites  ----------


# ---------- Supplementary Figure 3 Metabolites + Clusters Markers Heatmap  ----------
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

otu_table_sctp <- read.csv("D:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv",
                           row.names=1, sep = ";")

#ot_scree_filtered_rel_ab <- transform_feature_table(feature_table = ot_scree_filtered,
#                                                   transform_method = "rel_abundance")
#result <- cluster_barplot_panels(ot_scree_filtered_rel_ab, colour_palette = colours_vec)

# ---------- Cocultures in SNM3, SNM10 and BHI - S. aureus vs C. propinquum  ----------
# Read data
otu_table_cucultures <- read.csv("D:/SequencingData/SynCom100/4_Cocultures/emu_results/otu_table.csv",
                           row.names=1, sep = ";")


# Create dummy dataframe
bacteria <- c("S. aureus", "C. propinquum")
barcodes <- paste0("Barcode", sprintf("%02d", 1:27))

# Create random dummy data (e.g., relative abundances)
set.seed(123)  # for reproducibility
dummy_df <- data.frame(matrix(runif(2 * 27, min = 0, max = 1), 
                              nrow = 2, ncol = 27,
                              dimnames = list(bacteria, barcodes)))

# Define real sample names
real_names <- c(
  "Cpr16Sau_SNM3_R1", "Cpr16Sau_SNM3_R2", "Cpr16Sau_SNM3_R3",
  "Cpr70Sau_SNM3_R1", "Cpr70Sau_SNM3_R2", "Cpr70Sau_SNM3_R3",
  "Cpr265Sau_SNM3_R1", "Cpr265Sau_SNM3_R2", "Cpr265Sau_SNM3_R3",
  "Cpr16Sau_SNM10_R1", "Cpr16Sau_SNM10_R2", "Cpr16Sau_SNM10_R3",
  "Cpr70Sau_SNM10_R1", "Cpr70Sau_SNM10_R2", "Cpr70Sau_SNM10_R3",
  "Cpr265Sau_SNM10_R1", "Cpr265Sau_SNM10_R2", "Cpr265Sau_SNM10_R3",
  "Cpr16Sau_BHI_R1", "Cpr16Sau_BHI_R2", "Cpr16Sau_BHI_R3",
  "Cpr70Sau_BHI_R1", "Cpr70Sau_BHI_R2", "Cpr70Sau_BHI_R3",
  "Cpr265Sau_BHI_R1", "Cpr265Sau_BHI_R2", "Cpr265Sau_BHI_R3"
)

# Replace column names
colnames(dummy_df) <- real_names

# Display the resulting dataframe
dummy_df

# Step 2. Transpose for PCA
# Transpose so samples are rows and bacteria are columns
pca_input <- t(dummy_df)

# Step 3. Extract the sample group (the first part of the name)

# Extract group name (everything before _R)
sample_groups <- sub("_R[0-9]+$", "", rownames(pca_input))


# Step 4. Perform PCA
pca_res <- prcomp(pca_input, scale. = TRUE)

# Step 5. Plot PCA with ggplot2
library(ggplot2)

# Build a dataframe for plotting
pca_df <- data.frame(pca_res$x[, 1:2],
                     Group = sample_groups)

# Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of Bacterial Abundances",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"))

##### Barplots
# 0) Packages
library(tidyverse)

# 1) Build a sample metadata table from the column names
sample_meta <- tibble(Sample = colnames(dummy_df)) %>%
  tidyr::separate(Sample, into = c("Coculture", "Medium", "Replicate"),
                  sep = "_", remove = FALSE)
# Ensure consistent factor ordering in plots
sample_meta <- sample_meta %>%
  mutate(
    Coculture = factor(Coculture, levels = c("Cpr16Sau", "Cpr70Sau", "Cpr265Sau")),
    Medium    = factor(Medium,    levels = c("SNM3", "SNM10", "BHI"))
  )

# 2) Long format: Species, Sample, Abundance (+ join metadata)
df_long <- dummy_df %>%
  rownames_to_column("Species") %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(sample_meta, by = "Sample")

# 3) Convert to relative abundance per sample (so stacked bars sum to 1)
df_rel <- df_long %>%
  group_by(Sample) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

# 4) Average replicates within each Coculture × Medium × Species
df_avg <- df_rel %>%
  group_by(Coculture, Medium, Species) %>%
  summarize(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(
    # optional: control order of species in the stack
    Species = factor(Species, levels = c("S. aureus", "C. propinquum"))
  )

# 5) Plot: stacked barplot, faceted by Coculture (3 panels), x = Medium (3 bars)
p <- ggplot(df_avg, aes(x = Medium, y = MeanRelAbund, fill = Species)) +
  geom_col() +
  facet_wrap(~ Coculture, nrow = 1, drop = TRUE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Mean species composition by coculture and medium",
    x = "Medium",
    y = "Mean relative abundance",
    fill = "Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

print(p)


# ---------- Growth curves in SNM3, SNM10 and BHI - C. propinquum  ----------