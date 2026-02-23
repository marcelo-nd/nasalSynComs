########################################################
## Code: Marcelo Navarro Diaz
## Contact: marcelo.n.d@ciencias.unam.mx
########################################################

#### Load helper functions and libraries
# Install and load packages
cran_packages <- c(
  "cluster", "dplyr", "tibble",
  "pheatmap", "ggplot2", "cowplot", "pak"
)

# Install missing CRAN packages
installed <- rownames(installed.packages())

missing_cran <- setdiff(cran_packages, installed)

if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

pak::pak("github::marcelo-nd/nasalSynComsPkg")

library(nasalSynComsPkg)
library(dplyr)
library(ggplot2)
library(cowplot)
library(cluster)
library(tibble)
library(pheatmap)

# Set working directory
setwd("/Users/marcelonavarrodiaz/Documents/GitHub/nasalSynComs/Data")

# ---------- Figure 2. Screening Results with strain level information ----------
# Read otu table for the screening of all SynComs
otu_table_screening <- read.csv("./Supplementary_Table_S4_Screening_OTU_table.csv", row.names=1, sep = ";")

# Read strain inoculation table
strain_data <- readxl::read_excel(path = "./Supplementary_Table_S2_Syncom_Inocula.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)
strain_data <- tibble::column_to_rownames(strain_data, "Species")

# List of species to remove (they did not grow in any of the SynComs)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes")
strain_data <- remove_feature_by_prefix(df = strain_data, patterns = species_to_remove)

strain_data <- tibble::rownames_to_column(strain_data, "Species")

# Merge strain level data with the otu table, now the otu table contains strain info instead of only species
strain_ot <- merge_abundance_by_strain(otu_table_screening, strain_data)

# Merge the strain data for all except the Species we are interested in.
strain_ot <- merge_non_target_strains(strain_ot, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

# Save color pallette
colours_vec <- c("#EF5B5B", "#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1", "olivedrab3", "#e89d56")

# Clustering the SynComs based on compositional similarity at species level.
clustering_results <- cluster_samples(otu_table_screening)
# Get clusters dataframe for all SynComs
clusters <- clustering_results$clusters
# Order of samples acording to clustering
sample_order <- clustering_results$sample_order
# Get number of clusters
k <- clustering_results$best_k

# Create barplot for Figure 2
figure2 <- cluster_barplot_panels(abundance_df = transform_feature_table(strain_ot,transform_method = "rel_abundance"),
                                  cluster_df = clusters,
                                  sample_order = sample_order,
                                  best_k = k,
                                  strains = TRUE,
                                  colour_palette = colours_vec,
                                  species_order = c("Staphylococcus aureus",
                                                    "Corynebacterium accolens",
                                                    "Corynebacterium propinquum",
                                                    "Corynebacterium pseudodiphtheriticum",
                                                    "Corynebacterium tuberculostearicum",
                                                    "Cutibacterium avidum",
                                                    "Dolosigranulum pigrum",
                                                    "Staphylococcus epidermidis",
                                                    "Staphylococcus lugdunensis"))

print(figure2$plot)


# Calculate the mean abundance of S. aureus and C. propinquum in each cluster
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Staphylococcus aureus", k = k)
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Corynebacterium propinquum", k = k)
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Dolosigranulum pigrum", k = k)

# ---------- Figure 3. Selected SynComs Barplots ----------
# Barplot with strain-level information for C. propinquum and D. pigrum
# Read otu table containing all time points and replicates for selected SynComs
otu_table_timepoints <- read.csv("./Supplementary_Table_S5_Timepoints_OTU_table.csv",
                           row.names=1, sep = ";")

# Get inoculum data, this creates a dataframe containing wich species were inoculated in each SynCom
inoculum_spp_df <- strain_data %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
  ungroup()

inoc_spps <- inoculum_spp_df$Species
inoculum_spp_df <- select(inoculum_spp_df, -1)
rownames(inoculum_spp_df) <- inoc_spps

# Add inoculation info to barplot
strain_data2 <- as.data.frame(strain_data)

strain_data2 <- strain_data2[,3:ncol(strain_data2)]

rownames(strain_data2) <- strain_data$Species

# Merge strain data with otu table.
otu_table <- merge_abundance_by_strain(otu_table_timepoints, strain_data)

##### Run only for creating barplots with strain-level data for certain species.
otu_table <- merge_non_target_strains(otu_table, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

### If inoculation included and strain-level data for certain species is going to be used.
strain_data2 <- zero_out_species_in_samples(df = strain_data2, species_name = "Staphylococcus aureus USA300", sample_names = colnames(strain_data2))

strain_data2 <- merge_non_target_strains(strain_data2, c("Dolosigranulum pigrum", "Corynebacterium propinquum"))

time_names <- c("Inoc", "T1", "T2", "T3", "T4")

# Create a table for each SynCom
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

# Create Barplots
barplots1 <- barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                                 sc12, sc13, sc14,sc19,sc22),
                           strains = TRUE, shared_samples = FALSE,
                           experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11",
                                                 "SC12", "SC13", "SC14", "SC19","SC22"),
                           x_axis_title_size = 12, x_axis_text_size = 12,
                           y_axis_title_size = 12, y_axis_text_size = 12,
                           legend_pos = "none", legend_cols = 2,
                           legend_title_size = 12, legend_text_size = 12,
                           legend_key_size = 0.3, colour_palette = colours_vec,
                           species_order = c("Staphylococcus aureus",
                                             "Corynebacterium accolens",
                                             "Corynebacterium propinquum",
                                             "Corynebacterium pseudodiphtheriticum",
                                             "Corynebacterium tuberculostearicum",
                                             "Cutibacterium avidum",
                                             "Dolosigranulum pigrum",
                                             "Staphylococcus epidermidis",
                                             "Staphylococcus lugdunensis"))

barplots1 <- barplots1 + xlab("Time") + # for the x axis label
  ylab("Relative abundance")

barplots1

barplots2 <- barplots_grid(feature_tables = list(sc23, sc24, sc25, sc27, sc31,
                                                 sc34, sc39, sc40, sc44, sc50),
                           strains = TRUE, shared_samples = FALSE,
                           experiments_names = c("SC23", "SC24", "SC25", "SC27", "SC31",
                                                 "SC34", "SC39", "SC40", "SC44","SC50"),
                           x_axis_title_size = 12, x_axis_text_size = 12,
                           y_axis_title_size = 12, y_axis_text_size = 12,
                           legend_pos = "bottom", legend_cols = 3,
                           legend_title_size = 12, legend_text_size = 12,
                           legend_key_size = 0.3, colour_palette = colours_vec,
                           species_order = c("Staphylococcus aureus",
                                             "Corynebacterium accolens",
                                             "Corynebacterium propinquum",
                                             "Corynebacterium pseudodiphtheriticum",
                                             "Corynebacterium tuberculostearicum",
                                             "Cutibacterium avidum",
                                             "Dolosigranulum pigrum",
                                             "Staphylococcus epidermidis",
                                             "Staphylococcus lugdunensis"))

barplots2 <- barplots2 + xlab("Time") + # for the x axis label
  ylab("Relative abundance") + labs(fill = "Species")

barplots2

# Create figure 3 with both barplot figures as panels
figure3 <- cowplot::plot_grid(barplots1, barplots2,
                               align = "v",
                               ncol = 1,
                               rel_heights = c(46/100, 54/100))

figure3

# ---------- Figure 4. Bacterial diversity and Metabolites PCoA  ----------
# Read metadata for selected SynComs metabolomics samples
metadata <- read_metadata("./SynCom_timepoints_metadata.csv",
                          sort_table = TRUE)
metadata <- metadata[7:nrow(metadata),]
metadata_or_names <- rownames(metadata)
rownames(metadata) <- gsub("\\.mzML$", "", rownames(metadata))

# Add clustering results to metadata of Selected SynComs
meta_df <- add_cluster_column(
  meta_df = metadata,
  clusters_df = clustering_results$clusters,
  meta_key_col      = "ATTRIBUTE_SynCom",
  cluster_key_col   = "Sample",
  cluster_value_col = "Cluster",
  new_col_name      = "ATTRIBUTE_Cluster"
)

syncom_pallette <- c("indianred1", "#6279B8", "lavenderblush3", "#DA6A00",
                     "#738564", "purple4", "#56B4E9", "indianred4",
                     "#1a3a46", "hotpink4", "honeydew1", "hotpink",
                     "cyan3", "#cd541d", "#009E73", "#EC9704",
                     "#502F4C", "#FFBA49", "ivory3", "#9C4A1A")

clusters_pallete <- c("#583E26", "#F7C815", "lawngreen")

# PCoA Bacteria
res_euc <- pcoa_flex(
  metab_df      = otu_table_timepoints,
  metadata_df   = meta_df,
  color_var     = "ATTRIBUTE_SynCom",
  shape_var     = "ATTRIBUTE_Time",
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

# Read untargeted metabolomics data
feature_table_tic <- read_ft("./Supplementary_Table_S8_Untargeted_feature_table.csv",
                             sort_by_names = TRUE, p_sep = ";")

# Sort feature table by sample names
feature_table_tic <- feature_table_tic[, order(colnames(feature_table_tic))]

# Read sirius annotation data
an_table <- read.csv("./Supplementary_Table_S9_Sirius_annotations.csv", row.names=1)

# Get limma results
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

# Summarize results
sum_ht_sirius <- summarize_markers_and_heatmap_with_classes(
  out_file      = file.path("./markers_heatmap.pdf"), # Save here if necessary, change path accordingly
  metab_df      = feature_table_tic,
  metadata_df   = meta_df,
  sample_id_col = "Sample",
  cluster_var   = "ATTRIBUTE_Cluster",
  col_annot_vars = c("ATTRIBUTE_Cluster", "ATTRIBUTE_Time"),
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
  out_width     = 15,
  out_height    = 12,
  class_na_label = "Unclassified",
  class_na_color = "#BDBDBD",
  c_legend_ncol = 2,
  r_legend_ncol = 4,
  legend_side = "bottom"
)

sum_ht_sirius

# ---------- Figure 5. Repetition Experiment and Targeted Metabolites  ----------
# Read OTU table for repetition experiment
otu_table_rep_exp <- read.csv("./Supplementary_Table_S6_Repetition_syncoms_OTU_table.csv",
                           row.names=1, sep = ";")

# Lets get the means for the 3 replicates of each SynCom
collapsed_means <-
  otu_table_rep_exp |>
  tibble::rownames_to_column("Species") |>
  pivot_longer(-Species, names_to = "sample", values_to = "value") |>
  mutate(SynCom = sub("_(.*)$", "", sample)) |>
  group_by(Species, SynCom) |>
  summarize(mean = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(SynCom = factor(SynCom, levels = c("SC7","SC12","SC20","SC28","SC43"))) |>
  arrange(Species, SynCom) |>
  pivot_wider(names_from = SynCom, values_from = mean) |>
  tibble::column_to_rownames("Species")

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

# Create barplot for Figure 5a
barplot_from_feature_table(feature_table = collapsed_means[1:12,], legend_cols = 1, colour_palette = colours_vec)

# Targeted metabolomics analyses
# Read feature table
syncom_metabolites <- read.csv("./Supplementary_Table_S10_Targeted_metabolomics_feature_table.csv",
                              row.names=1, sep = ";")

info <- get_sample_info(syncom_metabolites)
sample_cols    <- info$sample_cols
base_names     <- info$base_names
unique_samples <- info$unique_samples

# Build matrices for computing lfc and significance
mats <- build_mats_from_df(syncom_metabolites, sample_cols, base_names)
mat_raw  <- mats$mat_raw
mat_mean <- mats$mat_mean
unique_samples <- mats$unique_samples

# Compute lfc and significance
res <- compute_lfc_and_stars(mat_raw, mat_mean, base_names, control_prefix = "CTRL")

lfc   <- res$lfc
stars <- res$stars

# Quick sanity check, making sure in both matrices samples and metabolites are in the same order
stopifnot(identical(dim(lfc), dim(stars)),
          identical(rownames(lfc), rownames(stars)),
          identical(colnames(lfc), colnames(stars)))

# Define colors
rwb <- colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))
# some data for improving graphing
max_abs <- max(abs(lfc[is.finite(lfc)]), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 51)

# Figure 5b heatmap
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
# Define which metabolites we want to include in the boxplots
met_list <- c("Aspartic acid", "Glutamic acid", "Tyrosine", "Riboflavin", "Alanine", "Glycine")

named_cols <- c(CTRL="#4E79A7", CPR1="#F28E2B", CPR2="#E15759", CPR3="#76B7B2", SAU="#EDC948",
                SynCom12="#B07AA1", SynCom20="#FF9DA7", SynCom28="#9C755F", SynCom43="#59A14F", SynCom7="red")

figure5c <- plot_metabolites_lfc_panel(
  df = syncom_metabolites,
  metabolites = met_list,
  ctrl_prefix = "CTRL",
  n_rows = 2, n_cols = 3,
  palette = named_cols
)

print(figure5c)

# ---------- Supplementary Figure 1. Human Microbiome Project data analyses ----------
# Read biom file after quality control and taxonomics assignment
nose_biom_path <- "./Supplementary_Table_S1_HMP_ASV_table.biom"
asv_table_nose <- load_biom_as_table(biom_path = nose_biom_path, strain_taxonomy = TRUE, order_table = TRUE)

asv_table_nose30 <- asv_table_nose[1:30,]

# Barplot for Supplementary Figure 1a
barplot_from_feature_table(feature_table = asv_table_nose30, sort_type = "similarity", legend_cols = 2, transform_table = FALSE)

# Transform data to relative abundance
asv_nose30_relAb <- transform_feature_table(asv_table_nose30, transform_method = "rel_abundance")

species_totals <- rowMeans(asv_table_nose30)

# Top 30 most abundant species Boxplot
top_species_names <- names(sort(species_totals, decreasing = TRUE))

top_species_df <- asv_nose30_relAb[top_species_names, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to="Sample", values_to="RelAbundance")

# Make boxplot for Supplemnetary Figure 1b
ggplot(top_species_df, aes(x=reorder(Species, RelAbundance, mean), 
                           y=RelAbundance)) +
  geom_boxplot(fill="#69b3a2") +
  coord_flip() +
  labs(x="Species", y="Relative Abundance") +
  theme_minimal(base_size=14)

# Heatmap
# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")
# Try silhouette method
sil_widths <- c()
for (k in 2:10) {
  pam_fit <- cluster::pam(dist_bc, diss = TRUE, k = k)
  sil_widths[k] <- pam_fit$silinfo$avg.width
}

best_k <- which.max(sil_widths)
cat("Optimal number of clusters:", best_k, "\n")

# Final clustering
pam_best <- cluster::pam(dist_bc, diss = TRUE, k = best_k)
clusters <- pam_best$clustering

# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")

# Hierarchical clustering
hc <- hclust(dist_bc, method = "ward.D2")

# Row Z-scores for comparability
z_scores <- t(scale(t(asv_table_nose30)))

# Column annotation
ha_col <- ComplexHeatmap::HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = structure(
    #circlize::rand_color(best_k),
    c("#B30223FF","#530E90FF","#DCFB90FF","#8091E6FF"),
    names = levels(factor(clusters))
  ))
)

# Custom color function
col_fun = circlize::colorRamp2(c(0, 1), c("white", "#FF6464"))

# Hetmap for Supplementary Figure 1c
ComplexHeatmap::Heatmap(asv_nose30_relAb,
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
# Align data and create long format data frames for calculations
prepared <- prepare_data(otu_table_timepoints, metadata)
# Compute bray curtis distances between replicates
dist_tbl <- compute_within_tp_distances(prepared$meta, prepared$X, method = "bray")

# Generate plot for Supplementary Figure 2a
sup_fig_2a <- plot_replicate_similarity(dist_tbl)
print(sup_fig_2a)

# Compute bray curtis distances between time points to final state
out <- compute_distance_to_final(prepared$meta, prepared$X, method = "bray", mode = "centroid")

# Generate plot for Supplementary Figure 2b
sup_fig_2b <- plot_distance_to_final(out$per_sample, out$summary)
print(sup_fig_2b)

# ---------- Supplementary Figure 3. Growth Curves ----------
# Load growthCurveExperiment script
source("https://raw.githubusercontent.com/marcelo-nd/growthCurveExperiment/main/growthCurveExperiment.R")

# Read data tables
aerobic_df <- read.csv("./Supplementary_Table_S7_aerobic_gcs.csv")
exp_aerobic <- GrowthCurveExperiment$new(name = "Aerobic_Run")
exp_aerobic$import_table(
  data_table = aerobic_df, 
  strains_names = c("C. acc99", "C. acc157","C. acc184",
                    "C. prop16", "C. prop70", "C. prop265",
                    "C. pseDSM", "C. pse242", "C.pse244",
                    "C. tub102", "C. tub223", "C. tubDSM",
                    "D. pig21", "D. pig61", "D. pig245",
                    "S. epi28", "S. epi231","S. epi251",
                    "S. lug81", "S. lug115", "S. lug239"), 
  replicates_per_strain = 8
)

anaerobic_df <- read.csv("./Supplementary_Table_S8_anaerobic_gcs.csv")
exp_anaerobic <- GrowthCurveExperiment$new(name = "Anaerobic_Run")
exp_anaerobic$import_table(
  data_table = anaerobic_df, 
  strains_names = c("A. oct133SNM", "A. oct211SNM", "A. oct259SNM",
                    "C. acn33_SNM", "C. acn86_SNM", "C. acnes149_SNM",
                    "C. avi32SNM", "C. avi181SNM", "C. avi208SNM"), 
  replicates_per_strain = 8
)

###### Results for each species
# A. octavius

gc_A.oct <- GrowthCurveExperiment(name = "A. octavius")

gc_A.oct$add_gco(exp_anaerobic$growthCurveObjects[1])
gc_A.oct$add_gco(exp_anaerobic$growthCurveObjects[2])
gc_A.oct$add_gco(exp_anaerobic$growthCurveObjects[3])

p1 <- gc_A.oct$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. accolens

gc_C.acc <- GrowthCurveExperiment(name = "C. accolens")

gc_C.acc$add_gco(exp_aerobic$growthCurveObjects[1])
gc_C.acc$add_gco(exp_aerobic$growthCurveObjects[2])
gc_C.acc$add_gco(exp_aerobic$growthCurveObjects[3])

p2 <- gc_C.acc$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. propinquum

gc_Cpro <- GrowthCurveExperiment(name = "C. propinquum")

gc_Cpro$add_gco(exp_aerobic$growthCurveObjects[4])
gc_Cpro$add_gco(exp_aerobic$growthCurveObjects[5])
gc_Cpro$add_gco(exp_aerobic$growthCurveObjects[6])

p3 <- gcsnm_Cpro$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. pseudodiphtheriticum

gc_C.pse <- GrowthCurveExperiment(name = "C. pseudodiphtheriticum")

gc_C.pse$add_gco(exp_aerobic$growthCurveObjects[7])
gc_C.pse$add_gco(exp_aerobic$growthCurveObjects[8])
gc_C.pse$add_gco(exp_aerobic$growthCurveObjects[9])

p4 <- gc_C.pse$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. tuberculostearicum

gc_Ctub <- GrowthCurveExperiment(name = "C. tuberculostearicum")

gc_Ctub$add_gco(exp_aerobic$growthCurveObjects[10])
gc_Ctub$add_gco(exp_aerobic$growthCurveObjects[11])
gc_Ctub$add_gco(exp_aerobic$growthCurveObjects[12])

p5 <- gc_Ctub$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. acnes

gc_C.acn <- GrowthCurveExperiment(name = "C. acnes")

gc_C.acn$add_gco(exp_anaerobic$growthCurveObjects[4])
gc_C.acn$add_gco(exp_anaerobic$growthCurveObjects[5])
gc_C.acn$add_gco(exp_anaerobic$growthCurveObjects[6])

p6 <- gc_C.acn$plot_curves(yScalemin = 0, yScalemax = 0.5)

# C. avidum

gc_C.avi <- GrowthCurveExperiment(name = "C. avidum")

gc_C.avi$add_gco(exp_anaerobic$growthCurveObjects[7])
gc_C.avi$add_gco(exp_anaerobic$growthCurveObjects[8])
gc_C.avi$add_gco(exp_anaerobic$growthCurveObjects[9])

p7 <- gc_C.avi$plot_curves(yScalemin = 0, yScalemax = 0.5)

# D. pigrum

gc_D.pig <- GrowthCurveExperiment(name = "D. pigrum")

gc_D.pig$add_gco(exp_aerobic$growthCurveObjects[13])
gc_D.pig$add_gco(exp_aerobic$growthCurveObjects[14])
gc_D.pig$add_gco(exp_aerobic$growthCurveObjects[15])

p8 <- gc_D.pig$plot_curves(yScalemin = 0, yScalemax = 1)

# S. epidermidis

gc_S.epi <- GrowthCurveExperiment(name = "S. epidermidis")

gc_S.epi$add_gco(exp_aerobic$growthCurveObjects[16])
gc_S.epi$add_gco(exp_aerobic$growthCurveObjects[17])
gc_S.epi$add_gco(exp_aerobic$growthCurveObjects[18])

p9 <- gc_S.epi$plot_curves(yScalemin = 0, yScalemax = 3)

# S. lugdunensis

gc_S.lug <- GrowthCurveExperiment(name = "S. lugdunensis")

gc_S.lug$add_gco(exp_aerobic$growthCurveObjects[19])
gc_S.lug$add_gco(exp_aerobic$growthCurveObjects[20])
gc_S.lug$add_gco(exp_aerobic$growthCurveObjects[21])

p10 <- gc_S.lug$plot_curves(yScalemin = 0, yScalemax = 2)

# Create panel plot
if (!require("patchwork", quietly = TRUE))
  install.packages("patchwork")

library(patchwork)

gc_theme <- theme_bw(base_size = 12) +
  theme(
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    strip.text   = element_text(size = 12)
  )

plots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
plots <- lapply(plots, function(p) p + gc_theme)

panel_titles <- c(
  "Anaerococcus octavius",
  "Corynebacterium accolens",
  "Corynebacterium propinquum",
  "Corynebacterium pseudodiphtheriticum",
  "Corynebacterium tuberculostearicum",
  "Cutibacterium acnes",
  "Cutibacterium avidum",
  "Dolosigranulum pigrum",
  "Staphylococcus epidermidis",
  "Staphylococcus lugdunensis"
)

plots <- Map(
  function(p, ttl) p + ggtitle(ttl),
  plots,
  panel_titles
)

sup_fig_3 <- wrap_plots(plots, ncol = 2, nrow = 5, guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 14, face = "bold")  # panel labels A–J
  )

sup_fig_3

# ---------- Supplementary Figure 4. Cocultures ----------
# Cocultures barplots in SNM3, SNM10 and BHI - S. aureus vs C. propinquum
otu_table_cocultures <- read.csv("./Supplementary_Table_S9_Cocultures_OTU_table.csv",
                                 row.names=1, sep = ";")

# Build a sample metadata table from the column names
sample_meta <- tibble(Sample = colnames(otu_table_cocultures)) %>%
  tidyr::separate(Sample, into = c("Coculture", "Medium", "Replicate"),
                  sep = "_", remove = FALSE)

# Ensure consistent factor ordering in plots
sample_meta <- sample_meta %>%
  mutate(
    Coculture = factor(Coculture, levels = c("Cpr16Sau", "Cpr70Sau", "Cpr265Sau")),
    Medium    = factor(Medium,    levels = c("SNM3", "SNM10", "BHI"))
  )

# Long format: Species, Sample, Abundance (+ join metadata)
df_long <- otu_table_cocultures %>%
  tibble::rownames_to_column("Species") %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(sample_meta, by = "Sample")

# Convert to relative abundance per sample (so stacked bars sum to 1)
df_rel <- df_long %>%
  group_by(Sample) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

# Average replicates within each Coculture × Medium × Species
df_avg <- df_rel %>%
  group_by(Coculture, Medium, Species) %>%
  summarize(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(
    # Order of species in the stack
    Species = factor(Species, levels = c("Staphylococcus aureus", "Corynebacterium propinquum"))
  )

# Create barplot panel for Supplementary Figure 3
sup_fig_4 <- ggplot(df_avg, aes(x = Medium, y = MeanRelAbund, fill = Species)) +
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

print(sup_fig_4)
