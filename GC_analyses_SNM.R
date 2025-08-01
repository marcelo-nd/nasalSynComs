# Load growthCurveExperiment script
source("https://raw.githubusercontent.com/marcelo-nd/growthCurveExperiment/main/growthCurveExperiment.R")

source("C:/Users/marce/Documents/Github/growthCurveExperiment/growthCurveExperiment.R")

gcsnm1 <- GrowthCurveExperiment(name = "SNM 1")

gcsnm1$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_SNM(Syn)Ctub102_Ctub223_Cpro_Sepi_Dpig_240223.xlsx",
                                    plate_reader_type = "Biotek",
                                    gc_range = "B220:CT365",
                                    strains_names = c("C. tub102", "C. tub223", "C. prop16", "C. prop70", "C. prop265", "S. epi28",
                                                      "S. epi231","S. epi251", "D. pig21", "D. pig61", "D. pig245", "Blank"),
                                    strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                    strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                    blank = TRUE, blank_col = 12, pr_correction = TRUE)

gcsnm1$strains_names

gcsnm1$plot_curves(calculate_model = TRUE)



#

gcsnm2 <- GrowthCurveExperiment(name = "SNM 2")

gcsnm2$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_SNM(Syn)Slug_CtubDSM_ Cacc_Cpse_240216.xlsx",
                                    plate_reader_type = "Biotek",
                                    gc_range = "B220:CT365",
                                    strains_names = c("Blank1", "S. lug81", "S. lug115", "S. lug239", "C. tubDSM", "C. acc99",
                                                      "C. acc157","C. acc184", "C. pseDSM", "C. pse242", "C.pse244", "Blank2"),
                                    strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                    strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                    blank = TRUE, blank_col = 1, pr_correction = TRUE)

gcsnm2$strains_names

gcsnm2$plot_curves(calculate_model = TRUE)


# results for each species

# SNM C. tuberculostearicum

gcsnm_Ctub <- GrowthCurveExperiment(name = "SNM C. tuberculostearicum")

gcsnm_Ctub$add_gco(gcsnm1$growthCurveObjects[1])
gcsnm_Ctub$add_gco(gcsnm1$growthCurveObjects[2])
gcsnm_Ctub$add_gco(gcsnm2$growthCurveObjects[5])

gcsnm_Ctub$strains_names

gcsnm_Ctub$plot_curves(calculate_model = FALSE)

# SNM C. propinquum

gcsnm_Cpro <- GrowthCurveExperiment(name = "SNM C. propinquum")

gcsnm_Cpro$add_gco(gcsnm1$growthCurveObjects[3])
gcsnm_Cpro$add_gco(gcsnm1$growthCurveObjects[4])
gcsnm_Cpro$add_gco(gcsnm1$growthCurveObjects[5])

gcsnm_Cpro$strains_names

gcsnm_Cpro$plot_curves(calculate_model = FALSE)

# SNM S. epidermidis

gcsnm_S.epi <- GrowthCurveExperiment(name = "SNM S. epidermidis")

gcsnm_S.epi$add_gco(gcsnm1$growthCurveObjects[6])
gcsnm_S.epi$add_gco(gcsnm1$growthCurveObjects[7])
gcsnm_S.epi$add_gco(gcsnm1$growthCurveObjects[8])

gcsnm_S.epi$strains_names

gcsnm_S.epi$plot_curves(calculate_model = FALSE)

# SNM D. pigrum

gcsnm_D.pig <- GrowthCurveExperiment(name = "SNM D. pigrum")

gcsnm_D.pig$add_gco(gcsnm1$growthCurveObjects[9])
gcsnm_D.pig$add_gco(gcsnm1$growthCurveObjects[10])
gcsnm_D.pig$add_gco(gcsnm1$growthCurveObjects[11])

gcsnm_D.pig$strains_names

gcsnm_D.pig$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 1)

# SNM S. lugdunensis

gcsnm_S.lug <- GrowthCurveExperiment(name = "SNM S. lugdunensis")

gcsnm_S.lug$add_gco(gcsnm2$growthCurveObjects[2])
gcsnm_S.lug$add_gco(gcsnm2$growthCurveObjects[3])
gcsnm_S.lug$add_gco(gcsnm2$growthCurveObjects[4])

gcsnm_S.lug$strains_names

gcsnm_S.lug$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 2)

# SNM C. accolens

gcsnm_C.acc <- GrowthCurveExperiment(name = "SNM C. accolens")

gcsnm_C.acc$add_gco(gcsnm2$growthCurveObjects[6])
gcsnm_C.acc$add_gco(gcsnm2$growthCurveObjects[7])
gcsnm_C.acc$add_gco(gcsnm2$growthCurveObjects[8])

gcsnm_C.acc$strains_names

gcsnm_C.acc$plot_curves(calculate_model = FALSE)

# SNM C. pseudodiphtheriticum

gcsnm_C.pse <- GrowthCurveExperiment(name = "SNM C. pseudodiphtheriticum")

gcsnm_C.pse$add_gco(gcsnm2$growthCurveObjects[9])
gcsnm_C.pse$add_gco(gcsnm2$growthCurveObjects[10])
gcsnm_C.pse$add_gco(gcsnm2$growthCurveObjects[11])

gcsnm_C.pse$strains_names

gcsnm_C.pse$plot_curves(calculate_model = FALSE)

