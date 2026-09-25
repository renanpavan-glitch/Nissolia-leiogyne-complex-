# -------------------- SETUP & LIBRARIES -------------------- #
library(tidyverse)
library(NbClust)
library(factoextra)
library(FactoMineR)
library(cluster) 
library(missMDA)
library(vegan)
library(pairwiseAdonis)
options(stringsAsFactors = FALSE)

# -------------------- 1. DATA READING & PROCESSING -------------------- #
famd_raw <- read_csv2("leaf_fruit2.CSV", locale = locale(encoding = "UTF-8")) |> 
  type_convert()

species_groups <- famd_raw$Sp.[1:44] 
famd_morpho <- famd_raw[1:44, 2:32] |> mutate( L5  = log(L1 / L3), L6  = log(L2 / L4),
                                               L11 = log(L7 / L9), L12 = log(L8 / L10))
# Convert character
famd_morpho[] <- lapply(famd_morpho, function(x) {
  if (is.character(x) || (is.numeric(x) && length(unique(x)) <= 3)) as.factor(x) else x
})

# -------------------- 2. MISSING DATA IMPUTATION & FAMD -------------------- #
estim_ncp <- estim_ncpFAMD(famd_morpho, ncp.max = 3)
famd_imput <- imputeFAMD(famd_morpho, ncp = estim_ncp$ncp)$completeObs

res_famd <- FAMD(famd_imput, graph = FALSE)
summary(res_famd)

# -------------------- 3. MULTIVARIATE HYPOTHESIS TESTING -------------------- #

# Gower distance matrix
dist_gower <- daisy(famd_imput, metric = "gower")

# A. Dispersion Homogeneity Test
dispersion <- betadisper(dist_gower, species_groups)
anova(dispersion)

# B. Global PERMANOVA
set.seed(123)
res_permanova <- adonis2(dist_gower ~ species_groups, permutations = 9999)
print(res_permanova)

# C. Pairwise PERMANOVA (Post-hoc test)
post_hoc_permanova <- pairwise.adonis(as.matrix(dist_gower), species_groups, p.adjust.m = "fdr")
print(post_hoc_permanova)

# -------------------- 4. FAMD VISUALIZATIONS -------------------- #
fviz_screeplot(res_famd, addlabels = TRUE, ylim = c(0, 35))
fviz_famd_ind(res_famd, repel = TRUE)
fviz_famd_var(res_famd, col.var = "orange", repel = TRUE)

# Variable contributions to principal axes
fviz_contrib(res_famd, choice = "var", axes = 1)
fviz_contrib(res_famd, choice = "var", axes = 2)
fviz_contrib(res_famd, choice = "var", axes = 3)

# Extract and format contribution table
contrib_table <- as.data.frame(res_famd$var$contrib) |> rownames_to_column(var = "Trait") |> 
  mutate(Total_Contrib = Dim.1 + Dim.2) |> arrange(desc(Total_Contrib))
head(contrib_table, 15)

# Visualizations colored by contribution gradient
fviz_famd_var(res_famd, "quanti.var", col.var = "contrib",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

fviz_famd_var(res_famd, "quali.var", col.var = "contrib",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))



# -------------------- 5. HCPC & CLUSTER PLOTS -------------------- #
res_hcpc <- HCPC(res_famd, graph = FALSE)

custom_colors_hcpc <- c("1" = "#E41A1C", "2" = "#354823", "3" = "#A2A475", "4" = "#FAD510", "5" = "#F781BF", "6" = "#377EB8", "7" = "#01665E")
custom_shapes <- c(6, 16, 18, 17, 15, 23, 1, 4)

# Simplified plotting function
plot_famd_clusters <- function(res_famd, res_hcpc, x_axis = 1, y_axis = 2) {
  # Format coordinate dataframe
  ind_df <- as.data.frame(res_famd$ind$coord) |> 
    mutate(
      cluster = as.factor(res_hcpc$data.clust$clust),
      OTU = rownames(res_famd$ind$coord)
    )
  
  colnames(ind_df)[x_axis] <- "AxisX"
  colnames(ind_df)[y_axis] <- "AxisY"
  # Calculate convex hulls grouped by cluster
  hulls <- ind_df |> 
    group_by(cluster) |> 
    slice(chull(AxisX, AxisY))
  # Axis labels with variance percentages
  var_x_label <- paste0("Dim ", x_axis, " (", round(res_famd$eig[x_axis, 2], 1), "%)")
  var_y_label <- paste0("Dim ", y_axis, " (", round(res_famd$eig[y_axis, 2], 1), "%)")
  # Generate plot
  ggplot(ind_df, aes(x = AxisX, y = AxisY, color = cluster, shape = cluster)) +
    geom_polygon(data = hulls, aes(fill = cluster, group = cluster), alpha = 0.2, color = NA) +
    geom_point(size = 5) +
    geom_text(aes(label = OTU), size = 3, vjust = -0.8, color = "black", show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = custom_colors_hcpc) +
    scale_fill_manual(values = custom_colors_hcpc) +
    scale_shape_manual(values = custom_shapes) +
    theme_classic() +
    labs(title = paste("FAMD - Dim", x_axis, "&", y_axis),
      x = var_x_label, y = var_y_label, 
      color = "Cluster", fill = "Cluster", shape = "Cluster")}

plot_famd_clusters(res_famd, res_hcpc, x_axis = 1, y_axis = 2)
plot_famd_clusters(res_famd, res_hcpc, x_axis = 1, y_axis = 3)
plot_famd_clusters(res_famd, res_hcpc, x_axis = 2, y_axis = 3)

# -------------------- 6. UPGMA DENDROGRAM -------------------- #
dist_matrix <- dist(res_famd$ind$coord, method = "euclidean")
upgma_model <- hclust(dist_matrix, method = "average")

fviz_dend(upgma_model, 
          k = 4,
          show_labels = TRUE,
          main = "UPGMA Dendrogram")
