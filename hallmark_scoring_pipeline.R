# ============================================================
# Hallmarks of Cancer Scoring in scRNA-seq (Seurat Pipeline)
# ============================================================
# Author: Renan Simões
# Description:
# This script computes hallmark gene set activity scores in
# single-cell RNA-seq data using Seurat's AddModuleScore function.
# It also generates UMAP visualizations and boxplots with
# statistical comparisons across cell types.
# ============================================================

# -----------------------------
# Load required libraries
# -----------------------------
library(Seurat)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# -----------------------------
# INPUTS (EDIT THESE PATHS)
# -----------------------------

# Path to hallmark gene file (Excel)
hallmark_file <- "data/hallmarks_genes.xlsx"

# Seurat object (must be pre-loaded or loaded here)
# Example:
# seurat_obj <- readRDS("data/seurat_object.rds")

# Metadata column for grouping (e.g., cell types)
group_column <- "celltype"

# -----------------------------
# Load hallmark gene sets
# -----------------------------
hallmarks_genes <- read_excel(hallmark_file)

# Split genes (semicolon-separated)
genes_split <- strsplit(as.character(hallmarks_genes$Genes), ";")
genes_split <- lapply(genes_split, trimws)

# Create list of gene sets per hallmark
hallmark_list <- split(
  unlist(genes_split),
  rep(hallmarks_genes$Hallmarks, lengths(genes_split))
)

# Remove duplicated genes
hallmark_list <- lapply(hallmark_list, unique)

# -----------------------------
# Filter genes present in dataset
# -----------------------------
genes_in_object <- rownames(seurat_obj)

hallmark_list_filtered <- lapply(hallmark_list, function(g) {
  intersect(g, genes_in_object)
})

# Check number of genes per hallmark
print(lapply(hallmark_list_filtered, length))

# -----------------------------
# Compute module scores
# -----------------------------
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = hallmark_list_filtered,
  name = "HallmarkScore",
  ctrl = 100
)

# Rename generated columns to actual hallmark names
real_names <- names(hallmark_list_filtered)
generated_cols <- paste0("HallmarkScore", seq_along(real_names))

for (i in seq_along(real_names)) {
  colnames(seurat_obj@meta.data)[
    colnames(seurat_obj@meta.data) == generated_cols[i]
  ] <- real_names[i]
}

# Inspect metadata
head(seurat_obj@meta.data)

# -----------------------------
# UMAP visualization
# -----------------------------
hallmark_names <- names(hallmark_list_filtered)

png("outputs/UMAP_Hallmark_Scores.png", width = 4800, height = 3600, res = 300)

p_umap <- FeaturePlot(
  object = seurat_obj,
  features = hallmark_names,
  ncol = 4,
  cols = c("lightblue", "lightgrey", "darkred"),
  combine = TRUE
) & theme(
  legend.position = "right",
  text = element_text(size = 10)
)

print(p_umap)
dev.off()

# -----------------------------
# Boxplot with statistical testing
# -----------------------------
data_plot <- seurat_obj@meta.data %>%
  select(all_of(hallmark_names), all_of(group_column))

data_long <- data_plot %>%
  pivot_longer(
    cols = all_of(hallmark_names),
    names_to = "Hallmark",
    values_to = "Score"
  )

png("outputs/Boxplots_Hallmarks.png", width = 5000, height = 4200, res = 300)

p_box <- ggplot(data_long, aes_string(x = group_column, y = "Score", fill = group_column)) +
  geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  ) +
  facet_wrap(~Hallmark, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Hallmark Activity by Cell Type",
    subtitle = "Wilcoxon rank-sum test",
    x = "Cell Type",
    y = "Module Score"
  )

print(p_box)
dev.off()

# -----------------------------
# Save updated Seurat object
# -----------------------------
saveRDS(seurat_obj, "outputs/seurat_with_hallmark_scores.rds")
