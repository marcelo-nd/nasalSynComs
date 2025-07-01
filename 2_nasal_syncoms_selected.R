# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_graphing.R")

###### Time-series analyses
otu_table_sctp <- read.csv("E:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv",
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

### Strain data processing
strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data/nasal_syncom_strains.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)

strain_data <- column_to_rownames(strain_data, "Species")

strain_data <- remove_feature_by_prefix(strain_data, species_to_remove)

strain_data <- rownames_to_column(strain_data, "Species")

# For inoculum with out strain-level data
inoculum_spp_df <- strain_data %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
  ungroup()

inoc_spps <- inoculum_spp_df$Species

inoculum_spp_df <- select(inoculum_spp_df, -1)

rownames(inoculum_spp_df) <- inoc_spps

######################################

# To use species-level data
# Write species-level OTU table
write.csv(x = otu_table_sctp_filt, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/2_20sc_species_ot.csv",
          row.names = T, quote = F)

otu_table <- otu_table_sctp_filt

# To use strain-level data
strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)
# Write strain-level OTU table
write.csv(x = strain_ft, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/3_20sc_strain_ot.csv",
          row.names = T, quote = F)

otu_table <- strain_ft


colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56", "black", "grey")

barplot_from_feature_table(otu_table, sort_type = "none", colour_palette = colours_vec,
                           legend_pos = "bottom", legend_cols = 3,
                           x_axis_text_angle = 90, x_axis_text_size = 7)

### Run to include all the time points and replicates
sc4 <- otu_table[1:12]
sc7 <- otu_table[13:24]
sc9 <- otu_table[25:36]
sc10 <- otu_table[37:48]
sc11 <- otu_table[49:60]
sc12 <- otu_table[61:72]
sc13 <- otu_table[73:84]
sc14 <- otu_table[85:96]
sc19 <- otu_table[97:108]
sc22 <- otu_table[109:120]
sc23 <- otu_table[121:132]
sc24 <- otu_table[133:144]
sc25 <- otu_table[145:156]
sc27 <- otu_table[157:168]
sc31 <- otu_table[169:180]
sc34 <- otu_table[181:192]
sc39 <- otu_table[193:204]
sc40 <- otu_table[205:216]
sc44 <- otu_table[217:228]
sc50 <- otu_table[229:240]

# SynCom Example - SC27
sc27_t0 <- select(inoculum_spp_df, "SC27")
sc27_t0 <- as.data.frame(sc27_t0)
rownames(sc27_t0) <- inoc_spps
sc27_t0 <- filter_features_by_col_counts(sc27_t0, min_count = 1, col_number = 1)

sc27_t1 <- sc27[1:3]
sc27_t2 <- sc27[4:6]
sc27_t3 <- sc27[7:9]
sc27_t4 <- sc27[10:12]

colours_vec_27 <- c("#ffe599", "dodgerblue4", "blueviolet","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56", "black", "grey")

barplots_grid(feature_tables = list(sc27_t0, sc27_t1, sc27_t2, sc27_t3, sc27_t4),
                            experiments_names = c("Inoc.", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec_27, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC27",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom")

# SynCom Example - SC40
sc40_t0 <- select(inoculum_spp_df, "SC40")
sc40_t0 <- as.data.frame(sc40_t0)
rownames(sc40_t0) <- inoc_spps
sc40_t0 <- filter_features_by_col_counts(sc40_t0, min_count = 1, col_number = 1)

sc40_t1 <- sc40[1:3]
sc40_t2 <- sc40[4:6]
sc40_t3 <- sc40[7:9]
sc40_t4 <- sc40[10:12]

colours_vec_40 <- c("#ffe599", "dodgerblue4", "mediumspringgreen",
                    "#EF5B5B", "olivedrab3", "#e89d56", "black", "grey")

barplots_grid(feature_tables = list(sc40_t0, sc40_t1, sc40_t2, sc40_t3, sc40_t4),
              experiments_names = c("Inoc.", "T1", "T2", "T3", "T4"),
              colour_palette = colours_vec_40, x_axis_text_size = 8, y_axis_text_size = 8,
              legend_title_size = 10, legend_text_size = 8, plot_title = "SC40",
              x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom")

# To use strain-level data
# Convert OTU Table to strain level table.
#strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)
#otu_table <- strain_ft

# Run to include only one Replicate
sc4 <- otu_table[c(2,5,8,11)]
colnames(sc4) <- c("T1", "T2", "T3", "T4")
sc7 <- otu_table[c(14,17,20,23)]
colnames(sc7) <- c("T1", "T2", "T3", "T4")
sc9 <- otu_table[c(26,29,32,35)]
colnames(sc9) <- c("T1", "T2", "T3", "T4")
sc10 <- otu_table[c(38,41,44,47)]
colnames(sc10) <- c("T1", "T2", "T3", "T4")
sc11 <- otu_table[c(50,53,56,59)]
colnames(sc11) <- c("T1", "T2", "T3", "T4")
sc12 <- otu_table[c(62,65,68,71)]
colnames(sc12) <- c("T1", "T2", "T3", "T4")
sc13 <- otu_table[c(74,77,80,83)]
colnames(sc13) <- c("T1", "T2", "T3", "T4")
sc14 <- otu_table[c(86,89,92,95)]
colnames(sc14) <- c("T1", "T2", "T3", "T4")
sc19 <- otu_table[c(98,101,104,107)]
colnames(sc19) <- c("T1", "T2", "T3", "T4")
sc22 <- otu_table[c(110,113,116,119)]
colnames(sc22) <- c("T1", "T2", "T3", "T4")
sc23 <- otu_table[c(122,125,128,131)]
colnames(sc23) <- c("T1", "T2", "T3", "T4")
sc24 <- otu_table[c(134,137,140,143)]
colnames(sc24) <- c("T1", "T2", "T3", "T4")
sc25 <- otu_table[c(146,149,152,155)] # Fix this
colnames(sc25) <- c("T1", "T2", "T3", "T4")
sc27 <- otu_table[c(158,161,164,167)]
colnames(sc27) <- c("T1", "T2", "T3", "T4")
sc31 <- otu_table[c(170,173,176,179)]
colnames(sc31) <- c("T1", "T2", "T3", "T4")
sc34 <- otu_table[c(182,185,188,191)]
colnames(sc34) <- c("T1", "T2", "T3", "T4")
sc39 <- otu_table[c(194,197,200,203)]
colnames(sc39) <- c("T1", "T2", "T3", "T4")
sc40 <- otu_table[c(206,209,212,215)]
colnames(sc40) <- c("T1", "T2", "T3", "T4")
sc44 <- otu_table[c(218,221,224,227)]
colnames(sc44) <- c("T1", "T2", "T3", "T4")
sc50 <- otu_table[c(230,233,236,239)]
colnames(sc50) <- c("T1", "T2", "T3", "T4")


### Run to include inoculation to barplot
strain_data2 <- as.data.frame(strain_data)

strain_data2 <- strain_data2[,3:ncol(strain_data2)]

rownames(strain_data2) <- strain_data$Species

# List of species to remove (only the prefix, not the full name)
#species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes")

# Apply the function
#strain_data2 <- remove_feature_by_prefix(strain_data2, species_to_remove)
#otu_table <- remove_feature_by_prefix(strain_ft, species_to_remove)

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

ggsave(
  "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/barplots_selected_grid_strains2.pdf",
  plot = barplots,
  dpi = 300, device = "pdf", width = 15, height = 8
)
