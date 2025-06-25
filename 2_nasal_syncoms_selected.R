# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_graphing.R")

###### Time-series analyses
otu_table_sctp <- read.csv("E:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv",
                           row.names=1, sep = ";")

strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data/nasal_syncom_strains.xlsx", sheet = "nasal_syncom_strains", range = "A1:AZ32", col_names = TRUE)

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
                                                                       "SC29_T1_R1", "SC27_T1_R2", "SC27_T1_R3",
                                                                       "SC29_T2_R1", "SC27_T2_R2", "SC27_T2_R3",
                                                                       "SC29_T3_R1", "SC27_T3_R2", "SC27_T3_R3",
                                                                       "SC29_TF_R1", "SC27_TF_R2", "SC27_TF_R3",
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


# To use species-level data
otu_table <- otu_table_sctp_filt

barplot_from_feature_table(otu_table, )

# To use strain-level data
# Convert OTU Table to strain level table.
strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)
otu_table <- strain_ft

### Run to include all the timepoints and replicates
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


# Run to inlcude only one Replicate
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
colnames(sc20) <- c("T1", "T2", "T3", "T4")
sc22 <- otu_table[c(110,113,116,119)]
colnames(sc23) <- c("T1", "T2", "T3", "T4")
sc23 <- otu_table[c(122,125,128,131)]
colnames(sc24) <- c("T1", "T2", "T3", "T4")
sc24 <- otu_table[c(134,137,140,143)]
colnames(sc25) <- c("T1", "T2", "T3", "T4")
sc25 <- otu_table[c(146,149,152,155)] # Fix this
colnames(sc26) <- c("T1", "T2", "T3", "T4")
sc26 <- otu_table[c(158,161,164,167)]
colnames(sc28) <- c("T1", "T2", "T3", "T4")
sc30 <- otu_table[c(170,173,176,179)]
colnames(sc32) <- c("T1", "T2", "T3", "T4")
sc34 <- otu_table[c(182,185,188,191)]
colnames(sc36) <- c("T1", "T2", "T3", "T4")
sc39 <- otu_table[c(194,197,200,203)]
colnames(sc42) <- c("T1", "T2", "T3", "T4")
sc40 <- otu_table[c(206,209,212,215)]
colnames(sc43) <- c("T1", "T2", "T3", "T4")
sc44 <- otu_table[c(218,221,224,227)]
colnames(sc47) <- c("T1", "T2", "T3", "T4")
sc50 <- otu_table[c(230,233,236,239)]
colnames(sc53) <- c("T1", "T2", "T3", "T4")


### Run to include inoculation to barplot
strain_data2 <- as.data.frame(strain_data)

strain_data2 <- strain_data2[,3:ncol(strain_data2)]

rownames(strain_data2) <- strain_data$Species

# List of species to remove (only the prefix, not the full name)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes")

# Apply the function
strain_data2 <- remove_feature_by_prefix(strain_data2, species_to_remove)
otu_table <- remove_feature_by_prefix(strain_ft, species_to_remove)

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
colnames(sc20) <- time_names
sc22 <- cbind(strain_data2["SC22"], otu_table[c(110,113,116,119)])
colnames(sc23) <- time_names
sc23 <- cbind(strain_data2["SC23"], otu_table[c(122,125,128,131)])
colnames(sc24) <- time_names
sc24 <- cbind(strain_data2["SC24"], otu_table[c(134,137,140,143)])
colnames(sc25) <- time_names
sc25 <- cbind(strain_data2["SC25"], otu_table[c(146,149,152,155)]) # fix this
colnames(sc26) <- time_names
sc27 <- cbind(strain_data2["SC27"], otu_table[c(158,161,164,167)])
colnames(sc28) <- time_names
sc31 <- cbind(strain_data2["SC31"], otu_table[c(170,173,176,179)])
colnames(sc32) <- time_names
sc34 <- cbind(strain_data2["SC34"], otu_table[c(182,185,188,191)])
colnames(sc36) <- time_names
sc39 <- cbind(strain_data2["SC39"], otu_table[c(194,197,200,203)])
colnames(sc42) <- time_names
sc40 <- cbind(strain_data2["SC40"], otu_table[c(206,209,212,215)])
colnames(sc43) <- time_names
sc44 <- cbind(strain_data2["SC44"], otu_table[c(218,221,224,227)])
colnames(sc47) <- time_names
sc50 <- cbind(strain_data2["SC50"], otu_table[c(230,233,236,239)])
colnames(sc53) <- time_names

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56", "black", "grey")

strain_level_sel = FALSE
strain_level_sel = TRUE

barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                    sc12, sc13, sc14,sc19,sc22),
              strains = strain_level_sel, shared_samples = FALSE,
              experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11",
                                    "SC12", "SC13", "SC14", "SC19","SC22"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)

barplots_grid(feature_tables = list(sc23, sc24, sc25, sc27, sc31,
                                    sc34, sc39, sc40, sc44, sc50),
              strains = strain_level_sel, shared_samples = FALSE,
              experiments_names = c("SC23", "SC24", "SC25", "SC27", "SC31",
                                    "SC34", "sc39", "SC40", "SC44","SC50"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)

