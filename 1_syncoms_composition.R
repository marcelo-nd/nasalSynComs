source("C:/Users/marce/Documents/GitHub/nasalSynComs/helper_functions.R")

# ---------- Screening Results ----------

otu_table_screening <- read.csv("D:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1, sep = ";")


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "SC18", "SC19", "SC20",
                                   "SC21", "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30",
                                   "SC31", "SC32", "SC33", "SC34", "SC35", "SC36", "SC37", "SC38", "SC39", "SC40",
                                   "SC41", "SC42", "SC43", "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)




colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

result <- cluster_barplot_panels(ot_scree_filtered, colour_palette = colours_vec)


#ot_scree_filtered_rel_ab <- transform_feature_table(feature_table = ot_scree_filtered,
#                                                   transform_method = "rel_abundance")
#result <- cluster_barplot_panels(ot_scree_filtered_rel_ab, colour_palette = colours_vec)


# ---------- Selected SynComs ----------


