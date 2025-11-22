source("C:/Users/marce/Documents/GitHub/nasalSynComs/helper_functions.R")
#source("https://raw.githubusercontent.com/marcelo-nd/nasalSynComs/helper_functions.R")
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(stringr)
#cluster

# Set WD
setwd("C:/Users/marce/OneDrive - UT Cloud/Link Lab - NasalSynCom - NasalSynCom/Paper/Data")

# ---------- Figure 2. Screening Results with strain level information ----------
otu_table_screening <- read.csv("./1_screening_otu_table.csv", row.names=1, sep = ";")


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "SC18", "SC19", "SC20",
                                   "SC21", "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30",
                                   "SC31", "SC32", "SC33", "SC34", "SC35", "SC36", "SC37", "SC38", "SC39", "SC40",
                                   "SC41", "SC42", "SC43", "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)

### Strain inoculation data processing
strain_data <- readxl::read_excel(path = "./2_nasal_syncom_strains.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)

strain_data <- tibble::column_to_rownames(strain_data, "Species")

strain_data <- remove_feature_by_prefix(strain_data, species_to_remove)

strain_data <- tibble::rownames_to_column(strain_data, "Species")

# Merge strain level data with Otu table
strain_ot <- merge_abundance_by_strain(ot_scree_filtered, strain_data)

#otu_table <- strain_ft

##### Run only for creating barplots with strain-level data for certain species.
# Merge the strain data for all except the Species we are interested in.
strain_ot <- merge_non_target_strains(strain_ot, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

strain_level_sel = TRUE

# Getting clustering results from OTU table with out the strain data.
clustering_results <- cluster_samples(ot_scree_filtered)

clusters <- clustering_results$clusters
#rel_abundance_ordered <- clustering_results$rel_abundance_ordered
sample_order <- clustering_results$sample_order
k <- clustering_results$best_k

# Create barplot
cluster_barplot_result <- cluster_barplot_panels(abundance_df = calculate_relative_abundance(strain_ot),
                                                 cluster_df = clusters,
                                                 sample_order = sample_order,
                                                 best_k = k,
                                                 strains = TRUE,
                                                 colour_palette = colours_vec)

print(cluster_barplot_result$plot)


# Calculate the mean abundance of S. aureus and C. propinquum in each cluster
cluster_mean_abundance(calculate_relative_abundance(ot_scree_filtered), species_name = "Staphylococcus aureus", k = k)
cluster_mean_abundance(calculate_relative_abundance(ot_scree_filtered), species_name = "Corynebacterium propinquum", k = k)

# ---------- Figure 3. Selected SynComs Barplots ----------
# Barplot with strain-level information for C. propinquum and D. pigrum
###### Time-series analyses
otu_table_sctp <- read.csv("./3_timepoints_otu_table.csv",
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
strain_data2 <- zero_out_species_in_samples(df = strain_data2, species_name = "Staphylococcus aureus USA300", sample_names = colnames(strain_data2))

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
metadata <- read_metadata("./4_timepoints_metadata.csv",
                          sort_table = TRUE)
metadata <- metadata[7:nrow(metadata),]
metadata_or_names <- rownames(metadata)
rownames(metadata) <- gsub("\\.mzML$", "", rownames(metadata))

# Add cluster result to metadata of Selected SynComs

meta_df <- add_cluster_column(
  meta_df = metadata,
  clusters_df = clustering_results$clusters,
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
feature_table_tic <- read_ft("./5_untargeted_quant_table.csv",
                             sort_by_names = TRUE, p_sep = ";")

# Sort feature table by sample names
feature_table_tic <- feature_table_tic[, order(colnames(feature_table_tic))]

an_table <- read.csv("./6_sirius_annotations.csv", row.names=1)

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

# ---------- Figure 5. Repetition Experiment and Targeted Metabolites  ----------
otu_table_rep_exp <- read.csv("./7_selected_syncoms_otu_table.csv",
                           row.names=1, sep = ";")
colnames(otu_table_rep_exp) <- c("SC7_1", "SC7_2", "SC7_3", "SC12_1", "SC12_2", "SC12_3",
                                 "SC20_1", "SC20_2", "SC20_3", "SC28_1", "SC28_2", "SC28_3",
                                 "SC43_1", "SC43_3") #"SC43_2"

# df is your original matrix/data.frame with species in rownames
# Ensure numeric matrix (sometimes read-in can make them character)
df_num <- as.data.frame(lapply(otu_table_rep_exp, function(x) as.numeric(as.character(x))),
                        row.names = rownames(otu_table_rep_exp))

collapsed_means <-
  df_num |>
  rownames_to_column("Species") |>
  pivot_longer(-Species, names_to = "sample", values_to = "value") |>
  mutate(SynCom = sub("_(.*)$", "", sample)) |>            # keep the SC id before the underscore
  group_by(Species, SynCom) |>
  summarize(mean = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(SynCom = factor(SynCom, levels = c("SC7","SC12","SC20","SC28","SC43"))) |>  # desired column order
  arrange(Species, SynCom) |>
  pivot_wider(names_from = SynCom, values_from = mean) |>
  column_to_rownames("Species")

# Result: one column per SynCom (SC7, SC12, SC20, SC28, SC43) and rows = species (rownames)
collapsed_means[1:3, ]  # quick peek

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

barplot_from_feature_table(feature_table = collapsed_means[1:12,], legend_cols = 1, colour_palette = colours_vec)

# Targeted metabolomics analyses
# Read Data
syncom_metabolites <- read_excel("./8_20251030_12C_Nasal_targeted_metabolomics_data_002.xlsx", sheet = "12C")

# Simple filtering
df_quant_filtered <- syncom_metabolites %>%
  dplyr::filter(Include == 1) %>%
  dplyr::select(-c("Include"))

filtered_metabolites <- df_quant_filtered$Metabolite

# Remoev unused columns
df_quant_filtered <- df_quant_filtered %>%
  dplyr::select(-c("Metabolite", "Kegg",	"BIGG",	"CHEBI",	"Polarity",	"Ion",	"Add",
                   "DP1_1", "DP1_2",	"DP1_3",	"DP2_1", "DP2_2",	"DP2_3", "DP3_1", "DP3_2",	"DP3_3", "DP1_1", "DP1_2",	"DP1_3",
                   "CTRL_1", "CTRL_2", "CTRL_3", "SAU_1", "SAU_2", "SAU_3" ))

rownames(df_quant_filtered) <- filtered_metabolites

info <- get_sample_info(df_quant_filtered)
sample_cols    <- info$sample_cols
base_names     <- info$base_names      # named vector; names are column names
unique_samples <- info$unique_samples

# quick peek
head(sample_cols)
head(base_names)
unique_samples

mats <- build_mats_from_df(df_quant_filtered, sample_cols, base_names)
mat_raw  <- mats$mat_raw
mat_mean <- mats$mat_mean
unique_samples <- mats$unique_samples

dim(mat_raw)   # metabolites x replicate columns
dim(mat_mean)  # metabolites x sample prefixes
head(colnames(mat_mean))  # sample prefixes like "CTRL", "CPR1", ...

res <- compute_lfc_and_stars(mat_raw, mat_mean, base_names, control_prefix = "CTRL150525")

lfc   <- res$lfc
stars <- res$stars

# Quick sanity check (optional)
stopifnot(identical(dim(lfc), dim(stars)),
          identical(rownames(lfc), rownames(stars)),
          identical(colnames(lfc), colnames(stars)))


rwb <- colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))
max_abs <- max(abs(lfc[is.finite(lfc)]), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 51)


pheatmap(
  lfc,
  main = "log2 Fold-Change vs CTRL (means across replicates)",
  color = rwb(50),
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = stars,
  number_color = "black",
  fontsize_number = 10,
  border_color = NA
)

# Boxplots for some metabolites

met_list <- c("Aspartic acid", "Glutamic acid", "Tyrosine", "Riboflavin", "Alanine", "Glycine")

named_cols <- c(CTRL="#4E79A7", CPR1="#F28E2B", CPR2="#E15759", CPR3="#76B7B2", SAU150525="#EDC948",
                SynCom12="#B07AA1", SynCom20="#FF9DA7", SynCom28="#9C755F", SynCom43="#59A14F", SynCom7="red")

p <- plot_metabolites_lfc_panel(
  df = df_quant_filtered,
  metabolites = met_list,
  ctrl_prefix = "CTRL150525",
  n_rows = 2, n_cols = 3,
  palette = named_cols,
  debug = TRUE  # <- prints row counts per metabolite at each step
)

print(p)




# ---------- Supplementary Figure 1. Human Microbiome Project data analyses ----------
nose_biom_path <- "./8_hmp_asv_table.biom"

asv_table_nose <- load_biom_as_table(biom_path = nose_biom_path, strain_taxonomy = TRUE, order_table = TRUE)

asv_nose_relAb <- transform_feature_table(asv_table_nose, transform_method = "rel_abundance")

asv_nose_relAb<- filter_low_abundance(asv_nose_relAb, threshold = 0.01)

# Select only the 30 more abundant species.
asv_table_nose30 <- asv_nose_relAb[1:30,]
#asv_table_nose30 <- asv_table_nose[1:30,]

#asv_nose30_relAb <- transform_feature_table(asv_table_nose30, transform_method = "rel_abundance")

#asv_nose30_relAb<- filter_low_abundance(asv_nose30_relAb, threshold = 0.01)

#species_totals <- rowMeans(asv_nose30_relAb)

species_totals <- rowMeans(asv_table_nose30)

# Barplot
barplot_from_feature_table(feature_table = asv_table_nose30, sort_type = "similarity", legend_cols = 2)

# Top 30 most abundant species Boxplot
top_species_names <- names(sort(species_totals, decreasing = TRUE))

top_species_df <- asv_nose30_relAb[top_species_names, ] %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to="Sample", values_to="RelAbundance")

ggplot(top_species_df, aes(x=reorder(Species, RelAbundance, mean), 
                           y=RelAbundance)) +
  geom_boxplot(fill="#69b3a2") +
  coord_flip() +
  labs(x="Species", y="Relative Abundance") +
  theme_minimal(base_size=14)

# Heatmap
library(vegan)
library(ComplexHeatmap)

library(cluster)

# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")
# Try silhouette method
sil_widths <- c()
for (k in 2:10) {
  pam_fit <- pam(dist_bc, diss = TRUE, k = k)
  sil_widths[k] <- pam_fit$silinfo$avg.width
}

best_k <- which.max(sil_widths)
cat("Optimal number of clusters:", best_k, "\n")

# Final clustering
pam_best <- pam(dist_bc, diss = TRUE, k = best_k)
clusters <- pam_best$clustering


# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")

# Hierarchical clustering
hc <- hclust(dist_bc, method = "ward.D2")

# Row Z-scores for comparability
z_scores <- t(scale(t(asv_table_nose30)))

# Column annotation
ha_col <- HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = structure(
    #circlize::rand_color(best_k),
    c("#B30223FF","#530E90FF","#DCFB90FF","#8091E6FF"),
    names = levels(factor(clusters))
  ))
)

# Custom color function
col_fun = circlize::colorRamp2(c(0, 1), c("white", "#FF6464"))

#asv_table30_scaled_by_sample <- scale(asv_table_nose30, scale = TRUE)

# Final heatmap with custom dendrogram # Este si es
#Heatmap(asv_table30_scaled_by_sample,
Heatmap(asv_nose30_relAb,
        name = "Relative abundance",
        top_annotation = ha_col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = as.dendrogram(hc),
        clustering_method_rows = "ward.D2",
        col = col_fun,
        column_title = paste("Samples grouped into", best_k, "clusters"),
        row_title = "Top 30 Species")

# ---------- Supplementary Figure 2. Replicates and stabilization ----------

# ---------- Supplementary Figure 3.  Growth Curves and Cocultures ----------
############ Cocultures in SNM3, SNM10 and BHI - S. aureus vs C. propinquum
otu_table_cocultures <- read.csv("./9_cocultures_otu_table.csv",
                                 row.names=1, sep = ";")

barplot_from_feature_table(otu_table_cocultures[1:2,], legend_cols = 1)

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Unassigned")
otu_table_cocultures <- remove_feature_by_prefix(otu_table_cocultures, species_to_remove)

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
colnames(otu_table_cocultures) <- real_names

# Display the resulting dataframe
otu_table_cocultures

# 1) Build a sample metadata table from the column names
sample_meta <- tibble(Sample = colnames(otu_table_cocultures)) %>%
  tidyr::separate(Sample, into = c("Coculture", "Medium", "Replicate"),
                  sep = "_", remove = FALSE)
# Ensure consistent factor ordering in plots
sample_meta <- sample_meta %>%
  mutate(
    Coculture = factor(Coculture, levels = c("Cpr16Sau", "Cpr70Sau", "Cpr265Sau")),
    Medium    = factor(Medium,    levels = c("SNM3", "SNM10", "BHI"))
  )

# 2) Long format: Species, Sample, Abundance (+ join metadata)
df_long <- otu_table_cocultures %>%
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
    Species = factor(Species, levels = c("Staphylococcus aureus", "Corynebacterium propinquum"))
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
