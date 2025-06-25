# Load helper functions and packages.
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiome_graphing.R")

# Load otu table for screening results
otu_table_screening <- read.csv("E:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)

colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC17", "SC18", "SC19", "SC20", "SC21",
                                   "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30", "SC31",
                                   "SC32", "SC33", "SC35", "SC36", "SC37", "SC38", "SC40", "SC41", "SC42", "SC43",
                                   "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50", "SC51", "SC52", "SC53")

ot_scree_filtered <- otu_table_screening[-10,]

# Remove species with no counts
ot_scree_filtered <- filter_features_by_counts_col_counts(ot_scree_filtered,min_count = 10, col_number = 1)

##### Screening Barplot
# Set colour palette 
colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","brown1", "olivedrab3", "darkorange3")

screening_barplot <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "none", colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

##### Screening Barplot sorted by S. aureus abundance
screening_barplot <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "feature_value", feature_to_sort = "Staphylococcus aureus",
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

##### Screening Barplot sorted by abundance similarity
screening_barplot <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "similarity",
                                                colour_palette = colours_vec, legend_pos = "bottom", legend_cols = 4)

screening_barplot

ggsave(
  "C:/Users/marce/Desktop/bar_plot_scree.pdf",
  plot = screening_barplot,
  dpi = 10, device = "pdf", width = 15, height = 6
)

##### Add dendograms to barplots
dendrogram_plot <- dendrogram_from_feature_table(ot_scree_filtered)

barplot_w_dendrogram <- cowplot::plot_grid(dendrogram_plot, screening_barplot, align = "v",
                         ncol = 1,
                         rel_heights = c(2/10, 8/10),
                         axis = "lr")

barplot_w_dendrogram

ggsave(
  "C:/Users/marce/Desktop/bar_plot_scree.pdf",
  plot = screening_barplot,
  dpi = 10, device = "pdf", width = 15, height = 6
)

# Screening barplot with strain information
# Convert OTU Table to strain level table.

strain_inoculation_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

strain_ot <- merge_abundance_by_strain(ot_scree_filtered, strain_inoculation_data)

screening_barplot_strain <- barplot_from_feature_table(feature_table = strain_ot, sort_type = "none", strains = TRUE,
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

screening_barplot_strain


########## Time-series analyses ##########
otu_table_sctp <- read.csv("E:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv", row.names=1)

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
                                                                       "SC20_T1_R1", "SC20_T1_R2", "SC20_T1_R3",
                                                                       "SC20_T2_R1", "SC20_T2_R2", "SC20_T2_R3",
                                                                       "SC20_T3_R1", "SC20_T3_R2", "SC20_T3_R3",
                                                                       "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
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
                                                                       "SC26_T1_R1", "SC26_T1_R2", "SC26_T1_R3",
                                                                       "SC26_T2_R1", "SC26_T2_R2", "SC26_T2_R3",
                                                                       "SC26_T3_R1", "SC26_T3_R2", "SC26_T3_R3",
                                                                       "SC26_TF_R1", "SC26_TF_R2", "SC26_TF_R3",
                                                                       "SC28_T1_R1", "SC28_T1_R2", "SC28_T1_R3",
                                                                       "SC28_T2_R1", "SC28_T2_R2", "SC28_T2_R3",
                                                                       "SC28_T3_R1", "SC28_T3_R2", "SC28_T3_R3",
                                                                       "SC28_TF_R1", "SC28_TF_R2", "SC28_TF_R3",
                                                                       "SC32_T1_R1", "SC32_T1_R2", "SC32_T1_R3",
                                                                       "SC32_T2_R1", "SC32_T2_R2", "SC32_T2_R3",
                                                                       "SC32_T3_R1", "SC32_T3_R2", "SC32_T3_R3",
                                                                       "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                                                                       "SC36_T1_R1", "SC36_T1_R2", "SC36_T1_R3",
                                                                       "SC36_T2_R1", "SC36_T2_R2", "SC36_T2_R3",
                                                                       "SC36_T3_R1", "SC36_T3_R2", "SC36_T3_R3",
                                                                       "SC36_TF_R1", "SC36_TF_R2", "SC36_TF_R3",
                                                                       "SC42_T1_R1", "SC42_T1_R2", "SC42_T1_R3",
                                                                       "SC42_T2_R1", "SC42_T2_R2", "SC42_T2_R3",
                                                                       "SC42_T3_R1", "SC42_T3_R2", "SC42_T3_R3",
                                                                       "SC42_TF_R1", "SC42_TF_R2", "SC42_TF_R3",
                                                                       "SC43_T1_R1", "SC43_T1_R2", "SC43_T1_R3",
                                                                       "SC43_T2_R1", "SC43_T2_R2", "SC43_T2_R3",
                                                                       "SC43_T3_R1", "SC43_T3_R2", "SC43_T3_R3",
                                                                       "SC43_TF_R1", "SC43_TF_R2", "SC43_TF_R3",
                                                                       "SC47_T1_R1", "SC47_T1_R2", "SC47_T1_R3",
                                                                       "SC47_T2_R1", "SC47_T2_R2", "SC47_T2_R3",
                                                                       "SC47_T3_R1", "SC47_T3_R2", "SC47_T3_R3",
                                                                       "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                                                                       "SC53_T1_R1", "SC53_T1_R2", "SC53_T1_R3",
                                                                       "SC53_T2_R1", "SC53_T2_R2", "SC53_T2_R3",
                                                                       "SC53_T3_R1", "SC53_T3_R2", "SC53_T3_R3",
                                                                       "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3"))


# Remove species with no counts
otu_table_sctp_filt <- filter_features_by_counts_col_counts(otu_table_sctp_sorted,
                                                        min_count = 10,
                                                        col_number = 1)

# List of species to remove (they did not grow in any of the SynComs)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

# Convert OTU Table to strain level table.

strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)

# To use species-level data
otu_table <- otu_table_sctp_filt

# To use strain-level data
otu_table <- strain_ft

# Create a df per SynCom (all replicates per time point)

sc4 <- otu_table[1:12]
sc7 <- otu_table[13:24]
sc9 <- otu_table[25:36]
sc10 <- otu_table[37:48]
sc11 <- otu_table[49:60]
sc12 <- otu_table[61:72]
sc13 <- otu_table[73:84]
sc14 <- otu_table[85:96]
sc20 <- otu_table[97:108]
sc23 <- otu_table[109:120]
sc24 <- otu_table[121:132]
sc25 <- otu_table[133:144]
sc26 <- otu_table[145:156]
sc28 <- otu_table[157:168]
sc32 <- otu_table[169:180]
sc36 <- otu_table[181:192]
sc42 <- otu_table[193:204]
sc43 <- otu_table[205:216]
sc47 <- otu_table[217:228]
sc53 <- otu_table[229:240]

# Create a df per SynCom (all time points) Single Replicate Barplots
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
sc20 <- otu_table[c(98,101,104,107)]
colnames(sc20) <- c("T1", "T2", "T3", "T4")
sc23 <- otu_table[c(110,113,116,119)]
colnames(sc23) <- c("T1", "T2", "T3", "T4")
sc24 <- otu_table[c(122,125,128,131)]
colnames(sc24) <- c("T1", "T2", "T3", "T4")
sc25 <- otu_table[c(134,137,140,143)]
colnames(sc25) <- c("T1", "T2", "T3", "T4")
sc26 <- otu_table[c(155,149,152,146)] # Fix this
colnames(sc26) <- c("T1", "T2", "T3", "T4")
sc28 <- otu_table[c(158,161,164,167)]
colnames(sc28) <- c("T1", "T2", "T3", "T4")
sc32 <- otu_table[c(170,173,176,179)]
colnames(sc32) <- c("T1", "T2", "T3", "T4")
sc36 <- otu_table[c(182,185,188,191)]
colnames(sc36) <- c("T1", "T2", "T3", "T4")
sc42 <- otu_table[c(194,197,200,203)]
colnames(sc42) <- c("T1", "T2", "T3", "T4")
sc43 <- otu_table[c(206,209,212,215)]
colnames(sc43) <- c("T1", "T2", "T3", "T4")
sc47 <- otu_table[c(218,221,224,227)]
colnames(sc47) <- c("T1", "T2", "T3", "T4")
sc53 <- otu_table[c(230,233,236,239)]
colnames(sc53) <- c("T1", "T2", "T3", "T4")


###### Including inoculum information
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
sc20 <- cbind(strain_data2["SC20"], otu_table[c(98,101,104,107)])
colnames(sc20) <- time_names
sc23 <- cbind(strain_data2["SC23"], otu_table[c(110,113,116,119)])
colnames(sc23) <- time_names
sc24 <- cbind(strain_data2["SC24"], otu_table[c(122,125,128,131)])
colnames(sc24) <- time_names
sc25 <- cbind(strain_data2["SC25"], otu_table[c(134,137,140,143)])
colnames(sc25) <- time_names
sc26 <- cbind(strain_data2["SC26"], otu_table[c(155,149,152,146)]) # fix this
colnames(sc26) <- time_names
sc28 <- cbind(strain_data2["SC28"], otu_table[c(158,161,164,167)])
colnames(sc28) <- time_names
sc32 <- cbind(strain_data2["SC32"], otu_table[c(170,173,176,179)])
colnames(sc32) <- time_names
sc36 <- cbind(strain_data2["SC36"], otu_table[c(182,185,188,191)])
colnames(sc36) <- time_names
sc42 <- cbind(strain_data2["SC42"], otu_table[c(194,197,200,203)])
colnames(sc42) <- time_names
sc43 <- cbind(strain_data2["SC43"], otu_table[c(206,209,212,215)])
colnames(sc43) <- time_names
sc47 <- cbind(strain_data2["SC47"], otu_table[c(218,221,224,227)])
colnames(sc47) <- time_names
sc53 <- cbind(strain_data2["SC53"], otu_table[c(230,233,236,239)])
colnames(sc53) <- time_names

##### Plotting the grid barplots

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56", "black", "grey")

barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                    sc12, sc13, sc14,sc20,sc23),
              strains = FALSE, shared_samples = FALSE,
              experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11",
                                    "SC12", "SC13", "SC14", "SC20","SC23"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)

barplots_grid(feature_tables = list(sc24, sc25, sc26, sc28, sc32,
                                    sc36, sc42, sc43, sc47, sc53),
              strains = TRUE, shared_samples = FALSE,
              experiments_names = c("SC24", "SC25", "SC26", "SC28", "SC32",
                                    "SC36", "sc42", "SC43", "SC47","SC53"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)


