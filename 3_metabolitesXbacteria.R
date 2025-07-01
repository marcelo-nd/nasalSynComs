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
otu_table <- remove_feature_by_prefix(otu_table, c("Corynebacterium tuberculostearicum DSM44922"))

# Load feature tables for metabolomics data
# Feature table - scaled
feature_table <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_QuantTable_scaled.csv",
                           sort_by_names = TRUE)

# Feature table - only imputated
feature_table <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_quantTable_imp.csv",
                           sort_by_names = TRUE)

# Feature table - random forest-selected metabolites
feature_table <- read_ft("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled.csv",
                           sort_by_names = TRUE)

# Remove the file extension endings for sample names.
rownames(feature_table) <- gsub("\\.mzML$", "", rownames(feature_table))

# Heatmaps
feature_table_heatmap(ft1 = otu_table, sig_stars = TRUE, corr_type = "spearman", pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr", hm_type =  "lower")

feature_table_heatmap(ft1 = otu_table, ft2 = t(feature_table), corr_type = "pearson", sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr") #by
