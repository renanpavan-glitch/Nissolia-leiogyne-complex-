# ==============================================================================
# SCRIPT: FAMD ANALYSIS, CLUSTERING, AND DENDROGRAMS
# ==============================================================================

library(readr)
library(NbClust)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(cluster) 
library(missMDA)
library(tidyverse)
options(stringsAsFactors = FALSE)

# -------------------- Data Reading and Processing -------------------- #
famd <- read_csv2("leaf_fruit.CSV", locale = locale(encoding = "UTF-8")) |> type_convert()

famd <- famd[1:26, 2:17]


famd[] <- lapply(famd, function(x) {
  if (is.character(x)) {
    as.factor(x)
  } else if (is.numeric(x) && length(unique(x)) <= 3) {
    as.factor(x)
  } else {
    x
  }
})


# -------------------- FAMD imput-------------------- #
estim_ncp <- estim_ncpFAMD(famd, ncp.max = 3)
famd_imput <- imputeFAMD(famd, ncp = estim_ncp$ncp)$completeObs
res.famd <- FAMD(famd_imput, graph = FALSE)

summary(res.famd)

# Basic FAMD visualizations
fviz_famd_ind(res.famd, repel = TRUE)
fviz_famd_var(res.famd, col.var = "orange", repel = TRUE)
fviz_screeplot(res.famd, addlabels = TRUE, ylim = c(0, 35))

# Contributions
fviz_contrib(res.famd, choice = "var", axes = 1)
fviz_contrib(res.famd, choice = "var", axes = 2)
fviz_contrib(res.famd, choice = "var", axes = 3)

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

fviz_famd_var(res.famd, "quali.var", col.var = "contrib",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# -------------------- HCPC and Cluster Plots -------------------- #
custom_colors_hcpc <- c("1" = "#E41A1C", "2" = "#354823", "3" = "#A2A475",
                        "4" = "#FAD510", "5" = "#F781BF", "6" = "#377EB8", "7" = "#01665E")
custom_shapes <- c(6, 16, 18, 17, 15, 23, 1, 4)

res.hcpc <- HCPC(res.famd, graph = FALSE)

plot_famd_clusters <- function(res.famd, res.hcpc, x_axis = 1, y_axis = 2) {
  ind_df <- as.data.frame(res.famd$ind$coord)
  ind_df$cluster <- as.factor(res.hcpc$data.clust$clust)
  ind_df$OTU <- rownames(res.famd$ind$coord)
  
  colnames(ind_df)[x_axis] <- "AxisX"
  colnames(ind_df)[y_axis] <- "AxisY"
  
  hulls_list <- lapply(split(ind_df, ind_df$cluster), function(df) df[chull(df$AxisX, df$AxisY), ])
  hulls <- do.call(rbind, hulls_list)
  hulls$cluster <- rep(levels(ind_df$cluster), sapply(hulls_list, nrow))
  
  var_x <- paste0("Dim ", x_axis, " (", round(res.famd$eig[x_axis, 2], 1), "%)")
  var_y <- paste0("Dim ", y_axis, " (", round(res.famd$eig[y_axis, 2], 1), "%)")
  
  ggplot(ind_df, aes(x = AxisX, y = AxisY, color = cluster, shape = cluster)) +
    geom_polygon(data = hulls, aes(fill = cluster, group = cluster), alpha = 0.2, color = NA) +
    geom_point(size = 5) +
    geom_text(aes(label = OTU), size = 3, vjust = -0.8, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = custom_colors_hcpc) +
    scale_fill_manual(values = custom_colors_hcpc) +
    scale_shape_manual(values = custom_shapes[1:length(unique(ind_df$cluster))]) +
    theme_classic() +
    labs(title = paste("FAMD - Dim", x_axis, "&", y_axis),
         x = var_x, y = var_y, color = "Cluster", fill = "Cluster", shape = "Cluster")
}

plot_famd_clusters(res.famd, res.hcpc, x_axis = 1, y_axis = 2)
plot_famd_clusters(res.famd, res.hcpc, x_axis = 1, y_axis = 3)
plot_famd_clusters(res.famd, res.hcpc, x_axis = 2, y_axis = 3)

# -------------------- Dendrogram (UPGMA) -------------------- #
custom_colors_dend <- c("1" = "#354823", "2" = "#E41A1C", "3" = "#A2A475",
                        "4" = "#F781BF", "5" = "#FAD510", "6" = "#377EB8")

distance_matrix <- dist(res.famd$ind$coord, method = "euclidean")
upgma_model <- hclust(distance_matrix, method = "average")

fviz_dend(upgma_model, 
          k = 3,
          show_labels = TRUE,
          palette = custom_colors_dend,
          main = "UPGMA Dendrogram")
