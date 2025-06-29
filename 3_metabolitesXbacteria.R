
# Load OTU tables
# Species-level
otu_table <- read.csv(file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/2_20sc_strain_ot.csv",
                      row.names = 1, header = TRUE)
# Strain-level
otu_table <- read.csv(file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/3_20sc_strain_ot.csv",
                      row.names = 1, header = TRUE)

# Load feature tables
# NO analog annotation
feature_table <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled.csv",
                            sort_by_names = TRUE)

feature_table <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled.csv",
                           sort_by_names = TRUE)

feature_table <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled_an.csv",
                           sort_by_names = TRUE)

rownames(feature_table) <- gsub("\\.mzML$", "", rownames(feature_table))



# Heatmaps

feature_table_heatmap2(ft1 = otu_table_sctp_filt_rel_ab, sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr")

feature_table_heatmap2(ft1 = abundance_per_strain2, corr_type = "pearson", sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr", hm_type =  "lower")

feature_table_heatmap2(ft1 = abundance_per_strain2, ft2 = t(feature_table), sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.5, stars_size = 0.8, fdr_correction = "BY")

feature_table_heatmap2(ft1 = t(feature_table), ft2 = abundance_per_strain2, sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.3, stars_size = 0.8, fdr_correction = "BY")

feature_table_heatmap2(ft1 = abundance_per_strain2, ft2 = t(feature_table), sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.3, stars_size = 0.8, fdr_correction = "BY")



# heatmaps no scaling

feature_table_heatmap2(ft1 = abundance_per_strain_not, sig_stars = TRUE, pval_adjust = TRUE, axis_text_size = 0.6, stars_size = 0.8, fdr_correction = "fdr")




############# To delete


feature_table_heatmap(otu_table_sctp_filt_rel_ab, otu_table_sctp_filt_rel_ab)

feature_table_heatmap(otu_table_sctp_filt_rel_ab, t(feature_table2))

feature_table_heatmap(abundance_per_strain2, abundance_per_strain2, sig_stars = TRUE)

feature_table_heatmap(abundance_per_strain2, t(feature_table), sig_stars = TRUE)



feature_table_heatmap_w_sig(otu_table_sctp_filt_rel_ab, otu_table_sctp_filt_rel_ab, fdr_correction = "BY")

feature_table_heatmap_w_sig(abundance_per_strain2, abundance_per_strain2, fdr_correction = "BY")

feature_table_heatmap_w_sig(otu_table_sctp_filt_rel_ab, t(feature_table), fdr_correction = "BY")

feature_table_heatmap_w_sig(abundance_per_strain2, t(feature_table), fdr_correction = "BY")
