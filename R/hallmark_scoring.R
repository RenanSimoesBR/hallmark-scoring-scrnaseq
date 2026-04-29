# ===============================
# Hallmark Scoring for scRNA-seq
# ===============================

library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#-----------------------------
# 1. Read hallmark gene sets
#-----------------------------
read_hallmark_file <- function(file_path,
                               hallmark_col = "Hallmarks",
                               gene_col = "Genes",
                               sep = ";") {
  
  df <- read_excel(file_path)
  
  # Check columns
  required_cols <- c(hallmark_col, gene_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Input file must contain columns: Hallmarks and Genes")
  }
  
  # Split genes
  genes_split <- strsplit(as.character(df[[gene_col]]), sep)
  genes_split <- lapply(genes_split, trimws)
  
  # Create list
  hallmark_list <- split(
    unlist(genes_split),
    rep(df[[hallmark_col]], lengths(genes_split))
  )
  
  hallmark_list <- lapply(hallmark_list, unique)
  
  return(list(
    gene_sets = hallmark_list,
    metadata = df 
  ))
}

#-----------------------------
# 2. Filter genes present in Seurat object
#-----------------------------
filter_hallmarks <- function(hallmark_list, seurat_obj) {
  
  genes_in_object <- rownames(seurat_obj)
  
  hallmark_filtered <- lapply(hallmark_list, function(g) {
    intersect(g, genes_in_object)
  })
  
  return(hallmark_filtered)
}

#-----------------------------
# 3. Add module scores
#-----------------------------
add_hallmark_scores <- function(seurat_obj,
                                hallmark_list,
                                prefix = "HallmarkScore",
                                ctrl = 100) {
  
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = hallmark_list,
    name = prefix,
    ctrl = ctrl
  )
  
  # Rename columns
  generated_cols <- paste0(prefix, seq_along(hallmark_list))
  real_names <- names(hallmark_list)
  
  colnames(seurat_obj@meta.data)[
    match(generated_cols, colnames(seurat_obj@meta.data))
  ] <- real_names
  
  return(seurat_obj)
}

#-----------------------------
# 4. UMAP plot
#-----------------------------
plot_hallmark_umap <- function(seurat_obj,
                               hallmarks,
                               ncol = 4,
                               output = NULL) {
  
  p <- FeaturePlot(
    object = seurat_obj,
    features = hallmarks,
    ncol = ncol,
    cols = c("lightblue", "lightgrey", "darkred")
  ) & theme(
    legend.position = "right",
    text = element_text(size = 10)
  )
  
  if (!is.null(output)) {
    ggsave(output, p, width = 16, height = 12, dpi = 300)
  }
  
  return(p)
}

#-----------------------------
# 5. Boxplot with stats
#-----------------------------
plot_hallmark_boxplot <- function(seurat_obj,
                                  hallmarks,
                                  group_col,
                                  output = NULL) {
  
  data_plot <- seurat_obj@meta.data %>%
    dplyr::select(all_of(hallmarks), all_of(group_col))
  
  data_long <- data_plot %>%
    pivot_longer(
      cols = all_of(hallmarks),
      names_to = "Hallmark",
      values_to = "Score"
    )
  
  p <- ggplot(data_long, aes_string(x = group_col, y = "Score", fill = group_col)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~Hallmark, scales = "free_y", ncol = 4) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "Hallmark Activity Comparison",
      x = "Group",
      y = "Module Score"
    )
  
  if (!is.null(output)) {
    ggsave(output, p, width = 16, height = 12, dpi = 300)
  }
  
  return(p)
}
