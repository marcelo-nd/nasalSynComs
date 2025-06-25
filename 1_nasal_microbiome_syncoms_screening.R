# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_graphing.R")

##### Screening Results
otu_table_screening <- read.csv("E:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1, sep = ";")


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "SC18", "SC19", "SC20",
                                   "SC21", "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30",
                                   "SC31", "SC32", "SC33", "SC34", "SC35", "SC36", "SC37", "SC38", "SC39", "SC40",
                                   "SC41", "SC42", "SC43", "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)

write.csv(x = ot_scree_filtered,
          file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/1_ot_screening_filtered.csv",
          row.names = T, quote = F, sep = ",")

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

##### Screening Barplot
screening_barplot <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "none", colour_palette = colours_vec, legend_pos = "bottom", legend_cols = 3,
                                                x_axis_text_angle = 45)
screening_barplot

##### Screening Barplot sorted by S. aureus abundance
screening_barplot_sau_sorted <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "feature_value", feature_to_sort = "Staphylococcus aureus",
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)
screening_barplot_sau_sorted

##### Screening Barplot sorted by Similarity
screening_barplot_sim_sorted <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "similarity",
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1,
                                                x_axis_title_size = 11, y_axis_title_size = 11, legend_title_size = 11,
                                                x_axis_text_size = 10, legend_text_size = 10, y_axis_text_size = 10,
                                                x_axis_text_angle = 45, x_hjust = 0.5, x_vjust = 0.5)

screening_barplot_sim_sorted <- screening_barplot_sim_sorted + xlab("SynCom") + # for the x axis label
ylab("Relative abundance") +# for the y axis label
labs(fill = "Species")

screening_barplot_sim_sorted

ggsave(
  "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results//barplot_scree.pdf",
  plot = screening_barplot_sim_sorted,
  dpi = 300, device = "pdf", width = 12, height = 6
)

##### Add dendograms to barplots
dendo_plot <- dendrogram_from_feature_table(ot_scree_filtered)

barplot_dendo <- cowplot::plot_grid(dendo_plot, screening_barplot_sim_sorted, align = "v",
                         ncol = 1,
                         rel_heights = c(2/10, 8/10),
                         axis = "lr")

barplot_dendo

ggsave(
  "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/screening_barplot_dendogram.pdf",
  plot = barplot_dendo,
  dpi = 300, device = "pdf", width = 12, height = 6
)

# Screening barplot with strain information
# Convert OTU Table to strain level table.
strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data/nasal_syncom_strains.xlsx", sheet = "nasal_syncom_strains", range = "A1:BA32", col_names = TRUE)

strain_ft <- merge_abundance_by_strain(ot_scree_filtered, strain_data)

screening_barplot_strains_sim_sorted <- barplot_from_feature_table(feature_table = strain_ft, sort_type = "similarity", strains = TRUE,
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1,
                                                x_axis_title_size = 11, y_axis_title_size = 11, legend_title_size = 11,
                                                x_axis_text_size = 10, legend_text_size = 10, y_axis_text_size = 10,
                                                x_axis_text_angle = 45, x_hjust = 0.5, x_vjust = 0.5)

screening_barplot_strains_sim_sorted

dendo_plot <- dendrogram_from_feature_table(ot_scree_filtered,
                                            margin_t = 40, margin_r = 25,
                                            margin_b = -25, margin_l = 25)

dendo_plot


barplot_strains_dendo <- cowplot::plot_grid(dendo_plot, screening_barplot_strains_sim_sorted, align = "v",
                                    ncol = 1, rel_heights = c(2/10, 8/10), axis = "b")

barplot_strains_dendo

ggsave(
  "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/screening_barplot_dendogram_strains.pdf",
  plot = barplot_strains_dendo,
  dpi = 300, device = "pdf", width = 12, height = 6
)
