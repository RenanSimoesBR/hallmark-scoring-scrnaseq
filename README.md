# Hallmarks of Cancer Scoring in Single-Cell RNA-seq

# Scoring-HallmarksOfCancer_scRNAseq

A pipeline to quantify Hallmarks of Cancer activity in single-cell RNA-seq data.

## Overview

This repository provides a reproducible and user-friendly pipeline to compute and visualize hallmarks of cancer gene set activity in single-cell RNA-seq (scRNA-seq) data using Seurat.

The workflow enables the integration of curated hallmark gene sets (derived from published literature) into single-cell analysis, allowing the exploration of biological processes at cellular resolution.

---

## Key Features

* Import hallmark gene sets from a structured Excel file
* Compute module scores using Seurat (`AddModuleScore`)
* Visualize hallmark activity across:

  * UMAP embeddings
  * Cell-type comparisons (boxplots)
* Perform statistical testing (Wilcoxon rank-sum test)
* Fully customizable and easily adaptable to different datasets

---

## Input Format

The hallmark gene set file must be provided as an Excel file with the following structure:

| Hallmarks  | Pathways | Genes               |
| ---------- | -------- | ------------------- |
| Hallmark_A | Pathway1 | gene1; gene2; gene3 |

**Column descriptions:**

* `Hallmarks`: Name of the hallmark gene set
* `Pathways`: Associated pathways (informational, not used in computation)
* `Genes`: Gene symbols separated by semicolons (`;`)

---

## Installation

Install required R packages:

```r
install.packages(c("Seurat", "tidyverse", "readxl", "ggpubr"))
```

---

## Usage

### 1. Load functions

```r
source("R/hallmark_scoring.R")
```

---

### 2. Load your Seurat object

```r
seurat_obj <- readRDS("data/seurat_object.rds")
```

---

### 3. Load hallmark gene sets

```r
hallmark_data <- read_hallmark_file("data/hallmarks_genes.xlsx")
hallmarks <- hallmark_data$gene_sets
```

---

### 4. Filter genes present in the dataset

```r
hallmarks_filtered <- filter_hallmarks(hallmarks, seurat_obj)
```

---

### 5. Compute hallmark scores

```r
seurat_obj <- add_hallmark_scores(seurat_obj, hallmarks_filtered)
```

---

### 6. Visualization

#### UMAP projection

```r
plot_hallmark_umap(
  seurat_obj,
  hallmarks = names(hallmarks_filtered),
  output = "outputs/UMAP_hallmarks.png"
)
```

#### Boxplots with statistical testing

```r
plot_hallmark_boxplot(
  seurat_obj,
  hallmarks = names(hallmarks_filtered),
  group_col = "celltype_v2",
  output = "outputs/boxplots_hallmarks.png"
)
```

## Output

* Updated Seurat object containing hallmark scores in metadata
* High-resolution figures ready for publication
* Flexible framework for downstream analyses

---

## Notes

* Gene symbols should match those in the Seurat object (case-sensitive).
* It is recommended to standardize gene names (e.g., uppercase) to avoid mismatches.
* The `Pathways` column is preserved for annotation but not used in scoring.

---

## Gene Set Source

The hallmark gene sets used in this repository were derived from:

**Sibai M, Cervilla S, Grases D ...
The spatial landscape of cancer hallmarks reveals patterns of tumor ecological dynamics and drug sensitivity
Cell Reports, 2025; 44
DOI: 10.1016/j.celrep.2024.115229**


Please cite the original publication when using these gene sets.

---

## Citation

If you use this repository in your research, please cite:

* This repository
* The original hallmark gene set publication
* Seurat (Butler et al., 2018; Stuart et al., 2019)

---

## Repository Structure

```
hallmark-scoring-scrnaseq/
│
├── R/
│   └── hallmark_scoring.R
│
├── scripts/
│   └── run_analysis.R
│
├── data/
│   └── hallmarks_genes.xlsx
│
├── outputs/
│   ├── UMAP_hallmarks.png
│   └── boxplots_hallmarks.png
│
├── README.md
└── LICENSE
```

---

## Author

Renan Simões

---

## License

This project is licensed under the MIT License.
