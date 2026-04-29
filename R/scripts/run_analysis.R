library(Seurat)

source("R/hallmark_scoring.R")

#-----------------------------
# Load Seurat object
#-----------------------------
seurat_obj <- readRDS("data/seurat_object.rds")

#-----------------------------
# Load hallmark file
#-----------------------------
hallmark_data <- read_hallmark_file("data/hallmarks_genes.xlsx")

hallmarks <- hallmark_data$gene_sets

#-----------------------------
# Filter genes
#-----------------------------
hallmarks_filtered <- filter_hallmarks(hallmarks, seurat_obj)

#-----------------------------
# Add scores
#-----------------------------
seurat_obj <- add_hallmark_scores(seurat_obj, hallmarks_filtered)

#-----------------------------
# Plots
#-----------------------------
plot_hallmark_umap(
  seurat_obj,
  hallmarks = names(hallmarks_filtered),
  output = "outputs/UMAP_hallmarks.png"
)

plot_hallmark_boxplot(
  seurat_obj,
  hallmarks = names(hallmarks_filtered),
  group_col = "celltype_v2",
  output = "outputs/boxplots_hallmarks.png"
)

#-----------------------------
# Save object
#-----------------------------
saveRDS(seurat_obj, "outputs/seurat_with_scores.rds")
