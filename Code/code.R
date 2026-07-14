########################################################
## Code: Marcelo Navarro Diaz
## Contact: marcelo.n.d@ciencias.unam.mx
########################################################

#### Load helper functions and libraries
# Install and load packages
cran_packages <- c(
  "pak", "tidyverse", "patchwork", "cluster", "pheatmap", "ggplotify", "Hmisc", "corrplot"
)

# Install missing CRAN packages
installed <- rownames(installed.packages())

missing_cran <- setdiff(cran_packages, installed)

if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

pak::pak("github::marcelo-nd/nasalSynComsPkg")

# Load packages
library(nasalSynComsPkg)
library(tidyverse)
library(patchwork)

# Set working directory
setwd("/Users/marcelonavarrodiaz/Documents/GitHub/nasalSynComs/Data")

# ---------- Figure 2. Screening of SynCom composition and assembly into clusters ----------
# Read otu table for the screening of all SynComs
otu_table_screening <- read.csv("./Supplementary_Table_S4_Screening_OTU_Table.csv", row.names=1, sep = ",")

# Read strain inoculation table
df_inoc <- readxl::read_excel(path = "./Supplementary_Table_S2_Syncom_Inocula.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)
strain_data <- tibble::column_to_rownames(df_inoc, "Species")

# List of species to remove from inoculation table (they did not grow in any of the SynComs). Remove S. aureus because all SynComs have it.
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Staphylococcus  aureus")

strain_data <- remove_feature_by_prefix(df = strain_data, patterns = species_to_remove)

strain_data <- tibble::rownames_to_column(strain_data, "Species")

# Expand otu table with  strain-level data. Now the otu table contains strain info instead of only species.
strain_ot <- merge_abundance_by_strain(otu_table_screening, strain_data)

# Save color pallette
colours_vec <- c(
  "#000000", # S. aureus
  "#E69F00", # C. accolens
  "#56B4E9", # C. propinquum
  "#882255", # C. pseudodiphtericum
  "#F0E442", # C. tuberculostearicum
  "#0072B2", # C. avidum
  "#D55E00", # D. pigrum
  "#CC79A7", # S. epidermidis
  "#44AA99" # S. lugdunensis
)

clusters_vec <- c(
  "#332288", # Cluster 1
  "#117733", # Cluster 2
  "#882255"  # Cluster 3
)

# Clustering the SynComs based on compositional similarity at species level.
clustering_results <- cluster_samples(otu_table_screening)

# Get clusters dataframe for all SynComs
clusters <- clustering_results$clusters
clusters <- clusters %>%
  mutate(Cluster = factor(paste("Cluster", Cluster)))

# Order of samples according to clustering
sample_order <- clustering_results$sample_order

# Get number of clusters
k <- clustering_results$best_k

# Create barplot for Figure 2
clustered_barplot <- cluster_barplot_panels(abundance_df = transform_feature_table(strain_ot, transform_method = "rel_abundance"),
                                  cluster_df = clusters,
                                  sample_order = sample_order,
                                  best_k = k,
                                  strains = TRUE,
                                  colour_palette = colours_vec,
                                  cluster_colors = clusters_vec,
                                  species_order = c("Staphylococcus aureus",
                                                    "Corynebacterium accolens",
                                                    "Corynebacterium propinquum",
                                                    "Corynebacterium pseudodiphtheriticum",
                                                    "Corynebacterium tuberculostearicum",
                                                    "Cutibacterium avidum",
                                                    "Dolosigranulum pigrum",
                                                    "Staphylococcus epidermidis",
                                                    "Staphylococcus lugdunensis"),
                                  strip_color = "white")

#print(clustered_barplot$plot)

### Adding grid plot to show inoculation

sample_order <- clustered_barplot$sample_order

# Define the species order in the grid plot
species_order <- c(
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

inoc_summary_strains <- df_inoc %>%
  filter(stringr::word(Species, 1, 2) != "Staphylococcus aureus") %>%
  mutate(Strain = Species, 
         Species_Parent = stringr::word(Species, 1, 2)) %>%
  
  # Set Species_Parent as a factor using species_order
  mutate(Species_Parent = factor(Species_Parent, levels = species_order)) %>%
  
  group_by(Strain, Species_Parent) %>%
  summarise(across(starts_with("SC"), ~as.numeric(any(.x == 1))), .groups = "drop") %>%
  
  pivot_longer(cols = starts_with("SC"), 
               names_to = "SynCom", 
               values_to = "Presence") %>%
  
  # Sort the dataframe so strains are grouped by the species factor
  arrange(Species_Parent, Strain) %>%
  
  # Set Strain levels.
  mutate(Strain = factor(Strain, levels = rev(unique(Strain)))) %>%
  
  # Forces the X-axis to follow the clustering results exactly
  mutate(SynCom = factor(SynCom, levels = sample_order))

# Remove X axis elements from barplot.
clustered_barplot2 <- clustered_barplot$plot + 
  theme(
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance"
  )

colours_vec <- c(
  "Anaerococcus octavius"                = "#999999",
  "Corynebacterium accolens"             = "#E69F00", 
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Corynebacterium pseudodiphtheriticum" = "#882255", 
  "Corynebacterium tuberculostearicum"   = "#F0E442",
  "Cutibacterium acnes"                  = "#661100",
  "Cutibacterium avidum"                 = "#0072B2", 
  "Dolosigranulum pigrum"                = "#D55E00", 
  "Staphylococcus epidermidis"           = "#CC79A7", 
  "Staphylococcus lugdunensis"           = "#44AA99",
  "Staphylococcus aureus"                = "#000000"
)

# Make grid_plot
grid_plot_labeled <- ggplot(inoc_summary_strains, aes(x = SynCom, y = Strain)) +
  geom_tile(aes(fill = ifelse(Presence == 1, as.character(Species_Parent), NA)), 
            color = "grey92", linewidth = 0.2) +
  scale_fill_manual(values = colours_vec, na.value = "white", na.translate = FALSE) + # "na.translate = FALSE" hides the NA from the legend logic
  scale_y_discrete(position = "right") + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, hjust = 0, face = "italic"),
    axis.title.x = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "none" 
  ) +
  labs(
    x = "SynCom",
    y = NULL
  )

#grid_plot_labeled

figure2 <- (clustered_barplot2 / grid_plot_labeled) + 
  plot_layout(heights = c(4, 5))

#figure2

#ggsave("../Graphs/Figure_2.pdf", figure2, width = 18, height = 11)
#ggsave("../Graphs/Figure_2.png", figure2, width = 18, height = 11)

# Calculate the mean abundance of S. aureus and C. propinquum in each cluster
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Staphylococcus aureus", k = k)
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Corynebacterium propinquum", k = k)
cluster_mean_abundance(transform_feature_table(otu_table_screening, transform_method = "rel_abundance"), species_name = "Dolosigranulum pigrum", k = k)

# ---------- Figure 3. Compositional changes across serial passages ----------
# Read otu table containing all time points and replicates for selected SynComs
otu_table_timepoints <- read.csv("./Supplementary_Table_S5_Timepoints_OTU_Table.csv",
                           row.names=1, sep = ",")

# Get inoculum data, this creates a dataframe containing which species were inoculated in each SynCom
inoculum_spp_df <- strain_data %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>%
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

# Remove S. aureus from inocula
strain_data2 <- zero_out_species_in_samples(df = strain_data2, species_name = "Staphylococcus aureus USA300", sample_names = colnames(strain_data2))

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
colours_vec <- c(
  "Anaerococcus octavius"                = "#999999",
  "Corynebacterium accolens"             = "#E69F00", 
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Corynebacterium pseudodiphtheriticum" = "#882255", 
  "Corynebacterium tuberculostearicum"   = "#F0E442",
  "Cutibacterium acnes"                  = "#661100",
  "Cutibacterium avidum"                 = "#0072B2", 
  "Dolosigranulum pigrum"                = "#D55E00", 
  "Staphylococcus epidermidis"           = "#CC79A7", 
  "Staphylococcus lugdunensis"           = "#44AA99",
  "Staphylococcus aureus"                = "#000000"
)

cluster_colors <- c(
  "Cluster 1" = "#332288", 
  "Cluster 2" = "#117733", 
  "Cluster 3" = "#882255"
)

all_tables <- list(sc7, sc9, sc11, sc12, sc14, sc22, sc4, sc10, sc19, sc23,
                   sc24, sc25, sc31, sc39, sc44, sc50, sc13, sc27, sc34, sc40)

all_names <- c("SC7", "SC9", "SC11", "SC12", "SC14", "SC22", "SC4", "SC10", "SC19", "SC23",
               "SC24", "SC25", "SC31", "SC39", "SC44", "SC50", "SC13", "SC27", "SC34", "SC40")

# To put S. aureus at the top of the barplots
species_order <- c(
  "Staphylococcus aureus",
  "Corynebacterium accolens",
  "Corynebacterium propinquum",
  "Corynebacterium pseudodiphtheriticum",
  "Corynebacterium tuberculostearicum",
  "Cutibacterium avidum",
  "Dolosigranulum pigrum",
  "Staphylococcus epidermidis",
  "Staphylococcus lugdunensis"
)

figure3 <- barplots_grid(
  feature_tables = all_tables,
  experiments_names = all_names,
  strains = TRUE,
  n_rows = 2,
  colour_palette = colours_vec,
  species_order = species_order,
  metadata_df = clusters,
  metadata_colors = cluster_colors,
  legend_key_size = 0.5,
  legend_cols = 4
)

figure3 <- figure3 +
  labs(x = "Synthetic Community",
       y = "Relative abundance")

#figure3

#ggsave("../Graphs/Figure_3.pdf", figure3, width = 12, height = 8)
#ggsave("../Graphs/Figure_3.png", figure3, width = 12, height = 8)

# ---------- Figure 4. Metabolic profiles of the 20 selected SynComs ----------
# Read metadata for selected SynComs metabolomics samples
metadata <- read_metadata("./Supplementary_Table_S6_SynCom_Timepoints_Metadata.csv",
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

# Read untargeted metabolomics data
feature_table_tic <- read_ft("./Supplementary_Table_S7_Untargeted_Metabolomics_Feature_Table.csv",
                             sort_by_names = TRUE, p_sep = ",")

# Sort feature table by sample names
feature_table_tic <- feature_table_tic[, order(colnames(feature_table_tic))]

# Read sirius annotation data
an_table <- read.csv("./Supplementary_Table_S8_Sirius_Annotations.csv", row.names=1)

# Get limma results
res_limma <- limma_markers_by_cluster_general(
  metab_df      = feature_table_tic,
  metadata_df   = meta_df,
  #  sample_id_col = "Sample",
  cluster_var   = "ATTRIBUTE_Cluster",
  covariates    = c("ATTRIBUTE_Time"),
  block_var     = "ATTRIBUTE_SynCom",
  log_transform = TRUE, log_offset = 1,
  do_pairwise   = TRUE
)

# Summarize results
sum_ht_sirius <- summarize_markers_and_heatmap_with_classes(
  #out_file      = file.path("./markers_heatmap.pdf"), # Save here if necessary, change path accordingly
  metab_df      = feature_table_tic,
  metadata_df   = meta_df,
  sample_id_col = "Sample",
  cluster_var   = "ATTRIBUTE_Cluster",
  col_annot_vars = c("ATTRIBUTE_Cluster", "ATTRIBUTE_Time"),
  col_annot_colors = list(
    ATTRIBUTE_Cluster = c("1" = "#332288", "2" = "#117733", "3" = "#882255"),
    ATTRIBUTE_Time = c("T1" = "grey", "T2" = "orange", "T3" = "green", "T4" = "blue")
  ),
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

#sum_ht_sirius$heatmap

limma_top_table <- sum_ht_sirius$top_table

write_csv(limma_top_table, "./Supplementary_Table_Sx_limma_top_table.csv")

# Convert ComplexHeatmap object
figure4 <- wrap_elements(grid::grid.grabExpr(
  ComplexHeatmap::draw(
    sum_ht_sirius$heatmap, 
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
))

#figure4

#ggsave("../Graphs/Figure_4.pdf", figure4, width = 15, height = 17)
#ggsave("../Graphs/Figure_4.png", figure4, width = 15, height = 17)

# ---------- Figure 5. Repetition Experiment and Targeted Metabolites  ----------
# Read OTU table for repetition experiment
otu_table_rep_exp <- read.csv("./Supplementary_Table_S12_Repetition_Syncoms_OTU_Table.csv",
                           row.names=1, sep = ",")

# Get the means for the 3 replicates of each SynCom
collapsed_means <-
  otu_table_rep_exp |>
  tibble::rownames_to_column("Species") |>
  tidyr::pivot_longer(-Species, names_to = "sample", values_to = "value") |>
  mutate(SynCom = sub("_(.*)$", "", sample)) |>
  group_by(Species, SynCom) |>
  summarize(mean = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(SynCom = factor(SynCom, levels = c("SC7","SC12","SC19","SC27","SC40"))) |>
  arrange(Species, SynCom) |>
  tidyr::pivot_wider(names_from = SynCom, values_from = mean) |>
  tibble::column_to_rownames("Species")

# Colour palette for the barplots
colours_vec <- c(
  "#E69F00", # C. accolens
  "#56B4E9", # C. propinquum
  "#882255", # C. pseudodiphtericum
  "#0072B2", # C. avidum
  "#D55E00", # D. pigrum
  "#000000", # S. aureus
  "#CC79A7", # S. epidermidis,
  "#44AA99" # S. lugdunensis
)

# Create barplot for Figure 5a
fig5a_barplot <- barplot_from_feature_table(feature_table = collapsed_means, legend_cols = 1, colour_palette = colours_vec)

# Define the specific SynComs for Figure 5A inoculation grid
fig5a_samples <- c("SC7", "SC12", "SC19", "SC27", "SC40")

# Define the colour vector for the grid plot
colours_vec_full <- c(
  "Anaerococcus octavius"                = "#999999",
  "Corynebacterium accolens"             = "#E69F00", 
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Corynebacterium pseudodiphtheriticum" = "#882255", 
  "Corynebacterium tuberculostearicum"   = "#F0E442",
  "Cutibacterium acnes"                  = "#661100",
  "Cutibacterium avidum"                 = "#0072B2", 
  "Dolosigranulum pigrum"                = "#D55E00", 
  "Staphylococcus epidermidis"           = "#CC79A7", 
  "Staphylococcus lugdunensis"           = "#44AA99"
)

species_order <- names(colours_vec_full)

# Process the Inoculum Table for these 5 SynComs
df_inoc <- readxl::read_excel("Supplementary_Table_S2_Syncom_Inocula.xlsx")

### Adding grid plot to show inoculation
inoc_summary_strains <- df_inoc %>%
  filter(stringr::word(Species, 1, 2) != "Staphylococcus aureus") %>%
  mutate(Strain = Species, 
         Species_Parent = stringr::word(Species, 1, 2)) %>%
  
  # Set Species_Parent as a factor to keep Syncom Order
  mutate(Species_Parent = factor(Species_Parent, levels = species_order)) %>%
  
  group_by(Strain, Species_Parent) %>%
  summarise(across(all_of(fig5a_samples), ~as.numeric(any(.x == 1))), .groups = "drop") %>%
  
  pivot_longer(cols = starts_with("SC"), 
               names_to = "SynCom", 
               values_to = "Presence") %>%
  
  # Sort the dataframe so strains are grouped by the species factor
  arrange(Species_Parent, Strain) %>%
  
  # Set Strain levels.
  mutate(Strain = factor(Strain, levels = rev(unique(Strain)))) %>%
  
  # Standard SynCom numerical sorting
  mutate(SynCom = factor(SynCom, 
                         levels = unique(SynCom)[order(as.numeric(gsub("SC", "", unique(SynCom))))]))

# Follow the clustering result
inoc_summary_strains <- inoc_summary_strains %>%
  mutate(SynCom = factor(SynCom, levels = fig5a_samples)) %>%
  filter(!is.na(SynCom))

# Create the grid plot
grid_plot_5a <- ggplot(inoc_summary_strains, aes(x = SynCom, y = Strain)) +
  geom_tile(aes(fill = ifelse(Presence == 1, as.character(Species_Parent), NA)), 
            color = "grey92", linewidth = 0.2) +
  scale_fill_manual(values = colours_vec_full, na.value = "white", na.translate = FALSE) +
  scale_y_discrete(position = "right") + # Put the Y-axis on the right
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14, hjust = 0, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "none" # Hide the legend
  ) +
  labs(
    x = "SynCom",
    y = NULL
  )

# clean barplot
fig5a_barplot_clean <- fig5a_barplot + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

# Combine to plots into Figure 5A
fig_5a <- (fig5a_barplot_clean / grid_plot_5a) + plot_layout(heights = c(2, 4), widths = c(5, 5))

#fig_5a

##### Targeted metabolomics analyses
# Read feature table
syncom_metabolites <- read.csv("./Supplementary_Table_S14_Targeted_metabolomics_feature_table.csv",
                              row.names=1, sep = ",")

info <- get_sample_info(syncom_metabolites)
sample_cols    <- info$sample_cols
base_names     <- info$base_names
unique_samples <- info$unique_samples

# Create matrices for computing lfc and significance
mats <- build_mats_from_df(syncom_metabolites, sample_cols, base_names)
mat_raw  <- mats$mat_raw
mat_mean <- mats$mat_mean
unique_samples <- mats$unique_samples

# Compute lfc and significance
res <- compute_lfc_and_stars(mat_raw, mat_mean, base_names, control_prefix = "CTRL")

lfc   <- res$lfc
stars <- res$stars

# Verify that in both matrices samples and metabolites are in the same order
stopifnot(identical(dim(lfc), dim(stars)),
          identical(rownames(lfc), rownames(stars)),
          identical(colnames(lfc), colnames(stars)))

# Define colors
rwb <- colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))
# some data for improve plot
max_abs <- max(abs(lfc[is.finite(lfc)]), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 51)

# Figure 5b heatmap
# Labels for plotting
new_labels <- c(
  "C. propinquum 16",
  "C. propinquum 70",
  "C. propinquum 265",
  "S. aureus",
  "SynCom 12",
  "SynCom 19",
  "SynCom 27",
  "SynCom 40",
  "SynCom 7"
)

# pheatmap as a ggplot object to combine graph with 5a
fig_5b <- ggplotify::as.ggplot(pheatmap::pheatmap(
  lfc,
  labels_col = new_labels,
  color = rwb(50),
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = stars,
  number_color = "black",
  fontsize_number = 12,
  border_color = NA,
  angle_col = 45,
  silent = TRUE
))

#fig_5b

# Create combined Figure 5
figure5 <- wrap_elements(fig_5a) / fig_5b + 
  plot_layout(heights = c(5, 3),
              widths = c(6, 5)) + 
  plot_annotation(tag_levels = 'A')

#figure5

# Save the result
#ggsave("../Graphs/Figure_5.pdf", figure5, width = 12, height = 15)
#ggsave("../Graphs/Figure_5.png", figure5, width = 10, height = 11)

# ---------- Supplementary Figure 1. Human Microbiome Project data analyses ----------
# Read biom file after quality control and taxonomics assignment
nose_biom_path <- "./Supplementary_Table_S1_HMP_ASV_table.biom"
asv_table_nose <- load_biom_as_table(biom_path = nose_biom_path, strain_taxonomy = TRUE, order_table = TRUE)
# Keep only 30 most abundant species
asv_table_nose30 <- asv_table_nose[1:30,]

# SF1a Barplot
SF1a <- barplot_from_feature_table(feature_table = asv_table_nose30, sort_type = "similarity", legend_cols = 2, transform_table = FALSE, x_axis_text_size = 0)

#SF1a

# Transform data to relative abundance
asv_nose30_relAb <- transform_feature_table(asv_table_nose30, transform_method = "rel_abundance")

species_totals <- rowMeans(asv_table_nose30)

# Top 30 most abundant species Boxplot
top_species_names <- names(sort(species_totals, decreasing = TRUE))

top_species_df <- asv_nose30_relAb[top_species_names, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  tidyr::pivot_longer(-Species, names_to="Sample", values_to="RelAbundance")

# Make boxplot for Supplemnetary Figure 1b
SF1b <- ggplot(top_species_df, aes(x=reorder(Species, RelAbundance, mean), 
                           y=RelAbundance)) +
  geom_boxplot(fill="#69b3a2") +
  coord_flip() +
  labs(x="Species", y="Relative Abundance") +
  theme_minimal(base_size=14)

# Combine and force alignment
figureSF1_top <- (SF1a / SF1b) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.margin = margin(5, 5, 5, 5)) # Adds a consistent margin to all panels

#figureSF1_top

ggsave("../Graphs/Figure_SF1_a_b.pdf", figureSF1_top, width = 16, height = 15)

# SF1c Heatmap
# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")
# Silhouette method
sil_widths <- c()
for (k in 2:10) {
  pam_fit <- cluster::pam(dist_bc, diss = TRUE, k = k)
  sil_widths[k] <- pam_fit$silinfo$avg.width
}

best_k <- which.max(sil_widths)
cat("Optimal number of clusters:", best_k, "\n")

# Final clustering
pam_best <- cluster::pam(dist_bc, diss = TRUE, k = best_k)
clusters_hmp <- pam_best$clustering

# Compute Bray-Curtis distance
dist_bc <- vegan::vegdist(t(asv_nose30_relAb), method = "bray")

# Hierarchical clustering
hc <- hclust(dist_bc, method = "ward.D2")

# Row Z-scores for comparability
z_scores <- t(scale(t(asv_table_nose30)))

# Column annotation
ha_col <- ComplexHeatmap::HeatmapAnnotation(
  Cluster = factor(clusters_hmp),
  col = list(Cluster = structure(
    #circlize::rand_color(best_k),
    c("#B30223FF","#530E90FF","#DCFB90FF","#8091E6FF"),
    names = levels(factor(clusters_hmp))
  ))
)

# Custom color function
col_fun = circlize::colorRamp2(c(0, 1), c("white", "#FF6464"))

# Hetmap for Supplementary Figure 1c
SF1c <- ComplexHeatmap::Heatmap(asv_nose30_relAb,
        name = "Relative abundance",
        top_annotation = ha_col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = as.dendrogram(hc),
        clustering_method_rows = "ward.D2",
        col = col_fun)
        #column_title = paste("Samples grouped into", best_k, "clusters"))
        #row_title = "Top 30 Species")

#SF1c

# Save heatmap
pdf("../Graphs/Figure_SF1_c.pdf", width = 10, height = 5)

ComplexHeatmap::draw(
  SF1c
)

dev.off()

# ---------- Supplementary Figure 2. Replicates and stabilization ----------
# Align data and create long format data frames for calculations
prepared <- prepare_data(otu_table_timepoints, metadata)
# Compute bray curtis distances between replicates
dist_tbl <- compute_within_tp_distances(prepared$meta, prepared$X, method = "bray")

# Cluster colors
cluster_colors <- c(
  "Cluster 1" = "#332288", 
  "Cluster 2" = "#117733", 
  "Cluster 3" = "#882255"
)

# Generate plot for Supplementary Figure 2a
sup_fig_2a <- plot_replicate_similarity(
  dist_tbl = dist_tbl, 
  metadata = clusters, 
  syncom_order = all_names, 
  n_rows = 4, 
  cluster_cols = cluster_colors
)

#print(sup_fig_2a)

# Compute bray curtis distances between time points to final state
out <- compute_distance_to_final(prepared$meta, prepared$X, method = "bray", mode = "centroid")

# Generate plot for Supplementary Figure 2b
sup_fig_2b <- plot_distance_to_final(
  per_sample = out$per_sample, 
  summary_tbl = out$summary, 
  metadata = clusters, 
  syncom_order = all_names, 
  n_rows = 4, 
  cluster_cols = cluster_colors
)
#print(sup_fig_2b)

figureSF2 <- (sup_fig_2a / sup_fig_2b) + 
  plot_layout(heights = c(2, 2))

#figureSF2

#ggsave("../Graphs/Figure_SF2.pdf", figureSF2, width = 13, height = 15)
#ggsave("../Graphs/Figure_SF2.png", figureSF2, width = 13, height = 5)

# ---------- Supplementary Figure 3. Species correlations ----------
# Normalize and average the technical replicaetes from main experiment ---
# Convert T4 counts to relative abundance
t4_rel <- apply(otu_table_timepoints[, grep("_T4_", colnames(otu_table_timepoints))], 2, function(x) x/sum(x))

# Average the technical replicates
syncom_names <- unique(gsub("_T4_.*", "", colnames(t4_rel)))
t4_averaged <- matrix(0, nrow = nrow(t4_rel), ncol = length(syncom_names))
rownames(t4_averaged) <- rownames(t4_rel)
colnames(t4_averaged) <- syncom_names

for(sc in syncom_names) {
  cols <- grep(paste0("^", sc, "_T4_"), colnames(t4_rel))
  t4_averaged[, sc] <- rowMeans(t4_rel[, cols, drop = FALSE])
}

# Convert repetition experiment counts to relative abundance
rep_rel <- apply(otu_table_rep_exp, 2, function(x) x/sum(x))

# Join main experiment and repetition experiment datasets
# Align species names and combine the columns
combined_all <- merge(t4_averaged, rep_rel, by = "row.names")
rownames(combined_all) <- combined_all$Row.names
combined_all$Row.names <- NULL

# Calculate spearman Correlation (N=35)
res_global <- Hmisc::rcorr(t(combined_all), type = "spearman")

# FDR adjustment
p_matrix <- res_global$P
p_adj_matrix <- matrix(p.adjust(p_matrix, method = "fdr"), 
                       nrow = nrow(p_matrix), 
                       ncol = ncol(p_matrix))
rownames(p_adj_matrix) <- rownames(p_matrix)
colnames(p_adj_matrix) <- colnames(p_matrix)

# Correlation heatmap for Supplementary Figure
corrplot::corrplot(res_global$r, 
                   method = "color", 
                   type = "upper", 
                   order = "hclust",
                   p.mat = p_adj_matrix, 
                   sig.level = 0.05,
                   insig = "label_sig",
                   pch.cex = 1.0,
                   tl.col = "black", 
                   tl.srt = 45, 
                   tl.cex = 0.8,
                   diag = FALSE
                   #title = "Species Co-occurrence Patterns (N=35)",
                   #mar = c(0, 0, 1, 0)
)

pdf("../Graphs/Figure_SF3.pdf", width = 8, height = 8)
corrplot::corrplot(res_global$r, 
                   method = "color", 
                   type = "upper", 
                   order = "hclust",
                   p.mat = p_adj_matrix, 
                   sig.level = 0.05,
                   insig = "label_sig",
                   pch.cex = 1.0,
                   tl.col = "black", 
                   tl.srt = 45, 
                   tl.cex = 0.8,
                   diag = FALSE
)
dev.off()
# ---------- Supplementary Figure 4. Growth measurement in plates ----------
# Load data
df <- readr::read_csv("./Supplementary_Table_S10_Growth_solid_media.csv")

colours_vec <- c(
  "Anaerococcus octavius"                = "#999999",
  "Corynebacterium accolens"             = "#E69F00", 
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Corynebacterium pseudodiphtheriticum" = "#882255", 
  "Corynebacterium tuberculostearicum"   = "#F0E442",
  "Cutibacterium acnes"                  = "#661100",
  "Cutibacterium avidum"                 = "#0072B2", 
  "Dolosigranulum pigrum"                = "#D55E00", 
  "Staphylococcus epidermidis"           = "#CC79A7", 
  "Staphylococcus lugdunensis"           = "#44AA99",
  "Staphylococcus aureus"                = "#000000"
)

# Prepare the data
plot_data <- df %>%
  mutate(Base_Species = str_extract(`Species Strain`, "^[A-Za-z]+\\s[a-z]+")) %>%
  mutate(Strain = str_remove(`Species Strain`, Base_Species)) %>%
  mutate(Strain = str_trim(Strain)) %>%
  mutate(Plotmath_Label = paste0("italic('", Base_Species, "')~'", Strain, "'")) %>% #Plotmath syntax (species names in italics only)
  mutate(Plotmath_Label = fct_reorder(Plotmath_Label, OD, .fun = mean, .na_rm = TRUE))

# Build the plot
figure_SF4 <- ggplot(plot_data, aes(x = OD, y = Plotmath_Label)) +
  
  geom_boxplot(
    aes(fill = Base_Species), 
    width = 0.4, 
    alpha = 0.4, 
    color = "gray50",
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Base_Species),
    height = 0.1, 
    width = 0, 
    size = 2.5, 
    alpha = 0.8
  ) +
  scale_color_manual(values = colours_vec) +
  scale_fill_manual(values = colours_vec) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Optical Density",
    y = "Species / Strain"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"), 
    legend.position = "none" 
  )

#figure_SF4

#ggsave("../Graphs/Figure_SF4.pdf", figure_SF4, width = 8, height = 7, dpi = 300)
#ggsave("../Graphs/Figure_SF4.png", figure_SF4, width = 8, height = 7, dpi = 300)

# ---------- Supplementary Figure 5. Cocultures ----------
# Cocultures barplots in SNM3, SNM10 and BHI - S. aureus vs C. propinquum
otu_table_cocultures <- read.csv("./Supplementary_Table_S8_Cocultures_OTU_Table.csv",
                                 row.names=1, sep = ",")

# Build a sample metadata table from the column names
sample_meta <- tibble(Sample = colnames(otu_table_cocultures)) %>%
  tidyr::separate(Sample, into = c("Coculture", "Medium", "Replicate"),
                  sep = "_", remove = FALSE)

# Make sure with factor that ordering in plots is correct and create labels in italics for species names only
sample_meta <- sample_meta %>%
  mutate(
    Coculture = factor(Coculture, 
                       levels = c("Cpr16Sau", "Cpr70Sau", "Cpr265Sau"),
                       # Format the strings for plotmath parsing
                       labels = c("italic('C. propinquum')~'16'", 
                                  "italic('C. propinquum')~'70'", 
                                  "italic('C. propinquum')~'265'")),
    Medium    = factor(Medium, levels = c("SNM3", "SNM10", "BHI"))
  )

df_long <- otu_table_cocultures %>%
  tibble::rownames_to_column("Species") %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(sample_meta, by = "Sample") # Jion metadata

# Convert to relative abundance
df_rel <- df_long %>%
  group_by(Sample) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

# Average replicates
df_avg <- df_rel %>%
  group_by(Coculture, Medium, Species) %>%
  summarize(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(
    # Order of species
    Species = factor(Species, levels = c("Staphylococcus aureus", "Corynebacterium propinquum"))
  )

# Define colors pallette as in the rest of plots
colours_vec <- c(
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Staphylococcus aureus"                = "#000000"
)

# Create barplot panel for Supplementary Figure 5
figure_SF5 <- ggplot(df_avg, aes(x = Medium, y = MeanRelAbund, fill = Species)) +
  geom_col() +
  facet_wrap(~ Coculture, nrow = 1, drop = TRUE, labeller = label_parsed) + # parse labels
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = colours_vec) +
  labs(
    x = "Medium",
    y = "Mean relative abundance",
    fill = "Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

#print(figure_SF5)

#ggsave("../Graphs/Figure_SF5.pdf", figure_SF5, width = 9, height = 4)
#ggsave("../Graphs/Figure_SF5.png", figure_SF5, width = 13, height = 5)

# ---------- Supplementary Figure 6. Deferroxamine experiments ----------
#Load data
df <- read_csv("./Supplementary_Table_S15_deferoxamine_growth.csv")
# Format dataframe to plot.
df <- df %>%
  mutate(
    Treatment = factor(Treatment, levels = c("SNM3", "SNM3 + DFO")), # order of treatments
    Base_Species = case_when(
      str_detect(Species_Strain, "propinquum") ~ "Corynebacterium propinquum",
      str_detect(Species_Strain, "aureus")     ~ "Staphylococcus aureus"
    ),
    # Create labels
    Strain_Label = case_when(
      Species_Strain == "Corynebacterium propinquum 16"  ~ "italic('C. propinquum')~'16'",
      Species_Strain == "Corynebacterium propinquum 70"  ~ "italic('C. propinquum')~'70'",
      Species_Strain == "Corynebacterium propinquum 265" ~ "italic('C. propinquum')~'265'",
      Species_Strain == "Staphylococcus aureus USA300"   ~ "italic('S. aureus')~'USA300'"
    )
  ) %>%
  mutate(
    # To sotr strains correctly
    Strain_Label = factor(Strain_Label, levels = c(
      "italic('C. propinquum')~'16'",
      "italic('C. propinquum')~'70'",
      "italic('C. propinquum')~'265'",
      "italic('S. aureus')~'USA300'"
    ))
  )

# Create the plot
figure_SF6 <- ggplot(df, aes(x = Treatment, y = OD, fill = Base_Species, color = Base_Species)) +
  scale_fill_manual(values = colours_vec) +
  scale_color_manual(values = colours_vec) +
  geom_boxplot(width = 0.5, alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
  facet_wrap(~ Strain_Label, nrow = 1, labeller = label_parsed) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Media Treatment",
    y = "Optical Density (OD)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5), # Add clean boxes around panels
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(color = "black")
  )

#print(figure_SF6)

# ---------- Supplementary Figure 7. Correlation between main and repetition experiments ----------
# Clean the timepoints Table (main experiment)
# Select only SynComs in repetition experiment
target_syncoms <- c("SC7", "SC12", "SC19", "SC27", "SC40")
sc_pattern <- paste0("^(", paste(target_syncoms, collapse = "|"), ")_T4_")

# Filter the columns need for this plot from main experiment data
t4_cols <- grep(sc_pattern, colnames(otu_table_timepoints), value = TRUE)
main_t4 <- otu_table_timepoints[, t4_cols]

# Rename to identify that they come from main epxeriment
colnames(main_t4) <- gsub("_T4_", "_Main_", colnames(main_t4))

# Rename repetition experiemnt samples
colnames(otu_table_rep_exp) <- paste0("Rep_", colnames(otu_table_rep_exp))

# Merge tables (main and repetition experiments)
combined_otu <- merge(main_t4, otu_table_rep_exp, by = "row.names", all = TRUE)
rownames(combined_otu) <- combined_otu$Row.names
combined_otu$Row.names <- NULL
combined_otu[is.na(combined_otu)] <- 0

# Convert to relative abundance
combined_otu <- transform_feature_table(combined_otu, transform_method = "rel_abundance")

# Remove species that did not grow
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes")

combined_otu <- remove_feature_by_prefix(df = combined_otu, patterns = species_to_remove)

##### Correlations between repetition experiment and main experiment.
# Calculate the mean abundance per species for each experiment
main_means <- rowMeans(combined_otu[, grepl("Main_", colnames(combined_otu))])
rep_means  <- rowMeans(combined_otu[, grepl("Rep_", colnames(combined_otu))])

# Create a data frame to correlations
corr_df <- data.frame(
  Species = names(main_means),
  Main_Exp = main_means,
  Rep_Exp = rep_means
)

# Calculate the Spearman Correlation
spearman_test <- cor.test(corr_df$Main_Exp, 
                          corr_df$Rep_Exp, 
                          method = "spearman", 
                          exact = FALSE)

rho <- round(spearman_test$estimate, 3)
p_val <- spearman_test$p.value

print(paste("Spearman Rho:", rho))
print(paste("P-value:", p_val))

# Visualize the correlation
comparison_correlations <- ggplot(corr_df, aes(x = Main_Exp, y = Rep_Exp)) +
  geom_point(size = 3, alpha = 0.6, color = "#2c3e50") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + # 1:1 line
  theme_minimal(base_size = 12) +
  labs(
    subtitle = paste0("Spearman Rho = ", rho, " (p < 0.001)"),
    x = "Mean Relative Abundance (Main Experiment)",
    y = "Mean Relative Abundance (Repetition Experiment)"
  )

##### Barplots betwen main experiment and repetition experiment.
# Create temporal metadata
metadata_combined <- data.frame(SampleID = colnames(combined_otu)) %>%
  mutate(
    Experiment = ifelse(grepl("Main_", SampleID), "Main", "Repetition"),
    # Extract SynCom ID
    SynCom = str_extract(SampleID, "SC\\d+")
  )

# Reshape
df_plot <- as.data.frame(combined_otu) %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(metadata_combined, by = "SampleID")

# Create the barplot
colours_vec_full <- c(
  "Anaerococcus octavius"                = "#999999",
  "Corynebacterium accolens"             = "#E69F00", 
  "Corynebacterium propinquum"           = "#56B4E9", 
  "Corynebacterium pseudodiphtheriticum" = "#882255", 
  "Corynebacterium tuberculostearicum"   = "#F0E442",
  "Cutibacterium acnes"                  = "#661100",
  "Cutibacterium avidum"                 = "#0072B2", 
  "Dolosigranulum pigrum"                = "#D55E00", 
  "Staphylococcus epidermidis"           = "#CC79A7", 
  "Staphylococcus lugdunensis"           = "#44AA99",
  "Staphylococcus aureus"                = "#000000"
)

comparison_barplots <- ggplot(df_plot, aes(x = SampleID, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.1) + 
  facet_grid(~ SynCom + Experiment, scales = "free_x", space = "free") +
  scale_fill_manual(values = colours_vec_full) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10), 
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Relative Abundance",
    x = "Independent Replicates"
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# Combine the plots
Figure_SF7 <- (comparison_barplots / comparison_correlations) + 
  plot_layout(heights = c(1.8, 3)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.margin = margin(5, 5, 5, 5))

#Figure_SF7

# Save plot
#ggsave("../Graphs/Figure_SF7.pdf", plot = Figure_SF7, width = 9, height = 10, dpi = 300)
#ggsave("../Graphs/Figure_SF7.pdf", plot = Figure_SF7, width = 9, height = 10, dpi = 300)