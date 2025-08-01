library(dplyr)
library(tidyr)
library(ggplot2)

# Set paths
setwd("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data")

github_path <- "C:/Users/marce/Documents/Github/"

source(paste0(github_path, "microbiome-help/otuTableReading.R"))

source(paste0(github_path, "microbiome-help/feature_table_wrangling.R"))

source(paste0(github_path, "microbiome-help/feature_table_graphing.R"))

#source(paste0(github_path, "microbiome-help/diversity_data_helper_functions.R"))

#source(paste0(github_path, "microbiome-help/table_importers.R"))

#source(paste0(github_path, "microbiome-help/graphs.R"))

# ASVs
nose_biom_path <- "./hmp_asv_table_w_tax_strain.biom"

asv_table_nose <- load_biom_as_table(biom_path = nose_biom_path, strain_taxonomy = TRUE, order_table = TRUE)

# Select only the 30 more abundant species.
asv_table_nose30 <- asv_table_nose[1:30,]

# Barplot
barplot_from_feature_table(feature_table = asv_table_nose30, sort_type = "similarity")

write.table(asv_table_nose, "./3_resultados/nose_asv_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

write.table(asv_table_nose30, "./3_resultados/nose_asv_table30_2.csv", sep = ",", col.names = FALSE, quote = FALSE)

####################################################################################
# Graph of means for each ASV/OTU

asv_table_nose2 <- asv_table_nose30

# Calculate means for each ASV/OTU
asv_table_nose2$Mean<- rowMeans(asv_table_nose2)

# Reduce table only to Means, discard all sample values
asv_table_nose2 <- select(asv_table_nose2, Mean)

# Barplot
barplot_from_feature_table(asv_table_nose2)

# Species co-occurrence analyses
####################################################################################
library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(asv_table_nose30, method = "pa")

# Infering co-ocurrences
cooccur.otus <- cooccur(otu_table_pa,
                        type = "spp_site",
                        spp_names = TRUE)

summary(cooccur.otus)
plot(cooccur.otus)


# Getting only the significant interactions
co <- print(cooccur.otus)

# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(otu_table_pa),
                    label = rownames(otu_table_pa),
                    color = "#606482",
                    shadow = TRUE) 

# Create an edges dataframe from the significant pairwise co-occurrences.
edges <- data.frame(from = co$sp1, to = co$sp2,
                    color = ifelse(co$p_lt <= 0.05,
                                   "#B0B2C1", "#3C3F51"),
                    dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

# Plotting network
library(visNetwork)
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")

# Export networka as edges list.
write.csv(edges, "C:/Users/marce/Desktop/coocur_network.csv")

# Export graph in "graphml" format
library(igraph)

g <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
print(g, e=TRUE, v=TRUE)
plot(g)
write.graph(g, "C:/Users/Marcelo/Desktop/coocur_network.graphml", format = "graphml")

####################################################################################
# Heatmaps

col_fun = circlize::colorRamp2(c(-0.7, 2, 5.5), c("#5F8D4E", "white", "#FF6464"))

# Scaled by columm (i.e. by sample).
asv_table30_scaled_by_sample <- scale(asv_table_nose30, center = TRUE, scale = TRUE)

# Graph simple heatmap
heatmap(asv_table30_scaled_by_sample, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")#heatmap(data.matrix(asv_proportions), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

# complex heatmap
htmp <- ComplexHeatmap::Heatmap(asv_table30_scaled_by_sample,
                                #name = "Scaled ASV abundance",
                                show_column_names = FALSE,
                                col = col_fun,
                                row_names_gp = grid::gpar(fontsize = 8),
                                heatmap_legend_param = list(
                                  title = "Scaled abundance",
                                  labels_rot = 0,
                                  direction = "horizontal"
                                ))

ComplexHeatmap::draw(htmp, heatmap_legend_side="bottom")

##### Scaling by ASV
asv_table30_scaled_by_asv <- t(scale(t(asv_table_nose30), center = TRUE, scale = TRUE))

col_fun = circlize::colorRamp2(c(-2, 2, 4), c("#5F8D4E", "white", "#FF6464"))

htmp_by_asv <- ComplexHeatmap::Heatmap(asv_table30_scaled_by_asv,
                                #name = "Scaled ASV abundance",
                                show_column_names = FALSE,
                                col = col_fun,
                                heatmap_legend_param = list(
                                  title = "Scaled abundance",
                                  labels_rot = 0,
                                  direction = "horizontal"
                                ))

ComplexHeatmap::draw(htmp_by_asv, heatmap_legend_side="bottom")

##### Normalization by sample

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply Min-Max normalization
asv_normalized <- as.data.frame(lapply(asv_table_nose30, min_max_norm), row.names =  row.names(asv_table_nose30))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#5F8D4E", "white", "#FF6464"))
htmp_norm <- ComplexHeatmap::Heatmap(as.matrix(asv_normalized),
                                 #name = "Scaled ASV abundance",
                                 show_column_names = FALSE,
                                 col = col_fun,
                                 heatmap_legend_param = list(
                                   title = "Normalized abundance",
                                   labels_rot = 0,
                                   direction = "horizontal"
                                 ))
ComplexHeatmap::draw(htmp_norm, heatmap_legend_side="bottom")

### Proportions heatmap

asv_proportions <- t(t(asv_table_nose30)/rowSums(t(asv_table_nose30)))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#5F8D4E", "white", "#FF6464"))
htmp_prop <- ComplexHeatmap::Heatmap(as.matrix(asv_proportions),
                                     #name = "Scaled ASV abundance",
                                     show_column_names = FALSE,
                                     col = col_fun,
                                     heatmap_legend_param = list(
                                       title = "Normalized abundance",
                                       labels_rot = 0,
                                       direction = "horizontal"
                                     ))
ComplexHeatmap::draw(htmp_prop, heatmap_legend_side="bottom")


#### Communities without S. aureus

BiocManager::install("ComplexHeatmap")

asv_no_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]==0]

#asv_no_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]<100]

# Transforming abundance data to presence/abscence
table_no_aureus_pa <- vegan::decostand(asv_no_aureus, method = "pa")

m1 = ComplexHeatmap::make_comb_mat(table_no_aureus_pa[,1:10])

m1 = ComplexHeatmap::make_comb_mat(t(table_no_aureus_pa))

ComplexHeatmap::UpSet(m1)


asv_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]>0]

table_aureus_pa <- vegan::decostand(asv_aureus, method = "pa")

m2 = ComplexHeatmap::make_comb_mat(table_no_aureus_pa[,1:10])

m2 = ComplexHeatmap::make_comb_mat(t(table_aureus_pa))

ComplexHeatmap::UpSet(m2)

######################################################################################

asv_no_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]==0]

asv_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]>0]

#col_fun = circlize::colorRamp2(c(-0.7, 2, 5.5), c("#5F8D4E", "white", "#FF6464"))

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#5F8D4E", "white", "#FF6464"))

# Scaled by column (i.e. by sample).
asv_table30_scaled_by_sample <- scale(asv_table_nose30, center = TRUE, scale = TRUE)

# Create annotation for heatmap using p/a of S. aureus.

asv_aureus <- asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]>0]

asv_table_nose30[, asv_table_nose30["Staphylococcus_aureus",]>0]

column_sa <- as.data.frame(t(asv_table_nose30["Staphylococcus_aureus",]>0))$Staphylococcus_aureus

column_sa_anotation = ComplexHeatmap::HeatmapAnnotation(S_aureus = column_sa)

# complex heatmap
htmp <- ComplexHeatmap::Heatmap(asv_table30_scaled_by_sample,
                                name = "Scaled ASV abundance",
                                top_annotation = column_sa_anotation,
                                show_column_names = FALSE,
                                col = col_fun,
                                row_names_gp = grid::gpar(fontsize = 8),
                                heatmap_legend_param = list(
                                  title = "Scaled abundance",
                                  labels_rot = 0,
                                  direction = "horizontal"
                                ))

ComplexHeatmap::draw(htmp, heatmap_legend_side="bottom")

