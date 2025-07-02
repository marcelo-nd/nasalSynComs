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


### Heatmap - Bacteria and metabolites selected by random forest
colnames(otu_table) == colnames(feature_table_rf)
feature_table_heatmap(ft1 = otu_table, ft2 = feature_table_rf, corr_type = "pearson", sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr", text_angle = 45) # can use also "by" correction

### PCAs metabolites
metadata <- read_metadata("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/SC100_metadata_noqcs_nosinStrs.csv",
                          sort_table = TRUE)
metadata <- metadata[7:nrow(metadata),]
rownames(metadata) <- gsub("\\.mzML$", "", rownames(metadata))

# PCA metabolome
ft_pca_1(ft = feature_table_imputated, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = FALSE, dist_method = "bray")

# PCA with vegan
ft_pca_2(ft = feature_table_imputated, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", dist_method = "bray")



# Remove highly variable metabolites (intrareplicate variability).
filtered_ft <- filter_by_error(feature_table = feature_table, metadata_table = metadata, grouping_var = "ATTRIBUTE_Sample", error_threshold = 25)
rownames(filtered_ft) <- rownames(feature_table)
