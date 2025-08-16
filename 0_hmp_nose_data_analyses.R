library(dplyr)
library(tidyr)
library(ggplot2)

# Set paths
setwd("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Data")

github_path <- "C:/Users/marce/Documents/Github/"

source(paste0(github_path, "microbiome-help/feature_table_wrangling.R"))

source(paste0(github_path, "microbiome-help/feature_table_graphing.R"))

# ASVs
nose_biom_path <- "./hmp_asv_table_w_tax_strain.biom"

asv_table_nose <- load_biom_as_table(biom_path = nose_biom_path, strain_taxonomy = TRUE, order_table = TRUE)

asv_nose_relAb <- transform_feature_table(asv_table_nose, transform_method = "rel_abundance")

asv_nose_relAb<- filter_low_abundance(asv_nose_relAb, threshold = 0.01)

# Select only the 30 more abundant species.
asv_table_nose30 <- asv_table_nose[1:30,]

asv_nose30_relAb <- transform_feature_table(asv_table_nose30, transform_method = "rel_abundance")

asv_nose30_relAb<- filter_low_abundance(asv_nose30_relAb, threshold = 0.01)


species_totals <- rowMeans(asv_nose30_relAb)

# Top 20 most abundant species
top_species <- sort(species_totals, decreasing = TRUE)

barplot(top_species, las=2, cex.names=0.7,
        main="Top 20 Most Abundant Species",
        ylab="Total Abundance")


top_species_names <- names(sort(species_totals, decreasing = TRUE))

top_species_df <- asv_nose30_relAb[top_species_names, ] %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to="Sample", values_to="RelAbundance")

ggplot(top_species_df, aes(x=reorder(Species, RelAbundance, mean), 
                           y=RelAbundance)) +
  geom_boxplot(fill="#69b3a2") +
  coord_flip() +
  labs(x="Species", y="Relative Abundance") +
  theme_minimal(base_size=14)


#####################################################################################################################################################################################################################################################################

species_richness <- colSums(asv_nose_relAb > 0)
richness_df <- data.frame(Sample=names(species_richness),
                          Richness=species_richness)

ggplot(richness_df, aes(x=Richness)) +
  geom_histogram(bins=30, fill="#404080", color="white", alpha=0.7) +
  labs(title="Species Richness Distribution",
       x="Number of Species", y="Count of Samples") +
  theme_minimal(base_size=14)

#####################################################################################################################################################################################################################################################################

# Barplot
barplot_from_feature_table(feature_table = asv_table_nose30, sort_type = "similarity", legend_cols = 2)

write.table(asv_table_nose, "./3_resultados/nose_asv_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

write.table(asv_table_nose30, "./3_resultados/nose_asv_table30_2.csv", sep = ",", col.names = FALSE, quote = FALSE)

#####################################################################################################################################################################################################################################################################

library(cluster)

# Try silhouette method
sil_widths <- c()
for (k in 2:10) {
  pam_fit <- pam(dist_bc, diss = TRUE, k = k)
  sil_widths[k] <- pam_fit$silinfo$avg.width
}

best_k <- which.max(sil_widths)
cat("Optimal number of clusters:", best_k, "\n")

# Final clustering
pam_best <- pam(dist_bc, diss = TRUE, k = best_k)
clusters <- pam_best$clustering


library(ComplexHeatmap)
library(circlize)

# Normalize species (row z-scores for comparability)
z_scores <- t(scale(t(asv_table_nose30)))

# Compute Bray-Curtis on top30 species
dist_bc <- vegdist(t(asv_table_nose30), method = "bray")

# Convert to a matrix for Heatmap
dist_mat <- as.matrix(dist_bc)

# Column annotation with clusters
ha_col <- HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = structure(
    circlize::rand_color(best_k), 
    names = levels(factor(clusters))
  ))
)

# Set a color palette from blue (low) to white to red (high)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(z_scores,
        name = "Z-score",
        top_annotation = ha_col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        clustering_distance_columns = "canberra",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        col = col_fun,
        column_title = paste("Samples grouped into", best_k, "clusters"),
        row_title = "Top 30 Species")


#####################################################################################################################################################################################################################################################################
library(vegan)
library(ComplexHeatmap)

# Compute Bray-Curtis distance
dist_bc <- vegdist(t(asv_nose30_relAb), method = "bray")

# Hierarchical clustering
hc <- hclust(dist_bc, method = "ward.D2")

# Row Z-scores for comparability
z_scores <- t(scale(t(asv_table_nose30)))

# Column annotation
ha_col <- HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = structure(
    circlize::rand_color(best_k),
    names = levels(factor(clusters))
  ))
)

# Custom color function
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Final heatmap with custom dendrogram # Este si es
Heatmap(asv_table30_scaled_by_sample,
        name = "Relative abundance",
        top_annotation = ha_col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = as.dendrogram(hc),
        clustering_method_rows = "ward.D2",
        col = col_fun,
        column_title = paste("Samples grouped into", best_k, "clusters"),
        row_title = "Top 30 Species")

#####################################################################################################################################################################################################################################################################
# Bray-Curtis for species (rows)
dist_rows <- vegdist(asv_nose30_relAb, method = "bray")
hc_rows <- hclust(dist_rows, method = "ward.D2")

# Heatmap with both row and column dendrograms
Heatmap(asv_table30_scaled_by_sample,
        name = "Relative abundance",
        top_annotation = ha_col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = as.dendrogram(hc),    # samples
        cluster_rows = as.dendrogram(hc_rows),  # species
        col = col_fun,
        column_title = paste("Samples grouped into", best_k, "clusters"),
        row_title = "Top 30 Species")

#####################################################################################################################################################################################################################################################################

library(reshape2)
library(ggplot2)
library(dplyr)

# Add cluster info
# Ensure tibble
rel_abund_with_cluster <- as.data.frame(t(asv_nose30_relAb))
rel_abund_with_cluster$Cluster <- factor(clusters)
rel_abund_with_cluster <- tibble::as_tibble(rel_abund_with_cluster)

# Summarize species abundances by cluster
cluster_means <- rel_abund_with_cluster %>%
  group_by(Cluster) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Melt for ggplot
df_long <- melt(cluster_means, id.vars = "Cluster",
                variable.name = "Species", value.name = "RelAbundance")

# Stacked barplot
ggplot(df_long, aes(x = Cluster, y = RelAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Cluster", y = "Mean Relative Abundance", 
       title = "Species Composition per Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(
    values = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(df_long$Species))),
      unique(df_long$Species)
    )
  )
#####################################################################################################################################################################################################################################################################
# Decide how many top species to keep
top_n <- 15   # change to 10 if you prefer

# Compute mean abundance across all clusters
species_means <- df_long %>%
  group_by(Species) %>%
  summarise(mean_abund = mean(RelAbundance), .groups = "drop")

top_species <- species_means %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

# Collapse less abundant species into "Other"
df_long$Species <- ifelse(df_long$Species %in% top_species, 
                          as.character(df_long$Species), "Other")

# Recompute after collapsing
df_long <- df_long %>%
  group_by(Cluster, Species) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

# Create randomized color palette
species_list <- unique(df_long$Species)
set.seed(2)  # change number to reshuffle colors
color_palette <- sample(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(species_list)))

# Plot
ggplot(df_long, aes(x = Cluster, y = RelAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Cluster", y = "Mean Relative Abundance",
       title = paste("Species Composition per Cluster (Top", top_n, ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = setNames(color_palette, species_list))

#####################################################################################################################################################################################################################################################################
library(dplyr)
library(ggplot2)

plot_s_aureus_vs_richness <- function(abund_df, species_name = "S.aureus") {
  
  # Transpose so samples are rows
  df_t <- as.data.frame(t(abund_df))
  colnames(df_t) <- make.names(colnames(df_t))
  
  # Check if species exists
  sp_col <- make.names(species_name)
  if (!(sp_col %in% colnames(df_t))) {
    stop(paste("Species", species_name, "not found in dataframe"))
  }
  
  # Compute richness
  df_t$Richness <- rowSums(df_t > 0)
  
  # Extract S. aureus abundance
  df_t$S_aureus <- df_t[[sp_col]]
  
  # Pearson correlation
  cor_test <- cor.test(df_t$Richness, df_t$S_aureus, method = "pearson")
  r_val <- round(cor_test$estimate, 3)
  p_val <- signif(cor_test$p.value, 3)
  
  # Scatter plot
  ggplot(df_t, aes(x = Richness, y = S_aureus)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    theme_minimal() +
    labs(x = "Species Richness",
         y = paste("Relative Abundance of", species_name),
         title = paste("Correlation between", species_name, "and Species Richness")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = max(df_t$Richness) * 0.7,
             y = max(df_t$S_aureus) * 0.9,
             label = paste0("r = ", r_val, "\n", "p = ", p_val),
             size = 5, color = "black")
}

plot_s_aureus_vs_richness(asv_nose30_relAb, species_name = "Staphylococcus_aureus")

# Heatmaps

col_fun = circlize::colorRamp2(c(-0.7, 2, 5.5), c("#5F8D4E", "white", "#FF6464"))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#5F8D4E", "white", "#FF6464"))

col_fun = circlize::colorRamp2(c(0, 1), c("white", "#FF6464"))

# Scaled by columm (i.e. by sample).
asv_table30_scaled_by_sample <- scale(asv_table_nose30, scale = TRUE)

asv_table30_scaled_by_sample <- as.matrix(asv_nose30_relAb)

colSums(asv_table30_scaled_by_sample)

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



#####################################################################################################################################################################################################################################################################
library(cluster)
library(vegan)

# transpose so samples are rows
dist_bc <- vegdist(t(asv_nose30_relAb), method="bray")

# Try different k values and compute average silhouette width
sil_widths <- c()
for (k in 2:10) {
  pam_fit <- pam(dist_bc, diss = TRUE, k = k)
  sil_widths[k] <- pam_fit$silinfo$avg.width
}

# Best k is the one with the highest silhouette width
best_k <- which.max(sil_widths)
cat("Optimal number of clusters:", best_k, "\n")

# Perform final clustering with best_k
pam_best <- pam(dist_bc, diss = TRUE, k = best_k)
clusters <- pam_best$clustering

# Add to dataframe
cluster_df <- data.frame(Sample = names(clusters),
                         Cluster = factor(clusters))

# diagnostic plot
plot(2:10, sil_widths[2:10], type="b",
     xlab="Number of clusters (k)",
     ylab="Average Silhouette Width",
     main="Silhouette Method for Optimal k")
abline(v=best_k, col="red", lty=2)


set.seed(123) # for reproducibility
library(cluster)

gap_stat <- clusGap(as.matrix(t(asv_nose30_relAb)), FUN = pam, 
                    K.max = 10, B = 50, d.power = 2)

print(gap_stat, method = "firstmax")

library(factoextra)
fviz_gap_stat(gap_stat)

ord_nmds <- metaMDS(dist_bc, k=2, trymax=100)

scores_df <- as.data.frame(scores(ord_nmds))
scores_df$Cluster <- factor(clusters)

ggplot(scores_df, aes(x=NMDS1, y=NMDS2, color=Cluster)) +
  geom_point(size=3, alpha=0.8) +
  labs(title=paste("NMDS of Samples (Optimal k =", best_k, ")")) +
  theme_minimal(base_size=14)

#####################################################################################################################################################################################################################################################################




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

