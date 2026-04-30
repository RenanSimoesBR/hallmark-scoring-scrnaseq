# scCancerHallmarks

## Hallmarks of Cancer Scoring in Single-Cell RNA-seq

A simple and reproducible pipeline to quantify **Hallmarks of Cancer activity** in single-cell RNA-seq (scRNA-seq) data using Seurat.

---

## Overview

This repository provides a straightforward script to compute and visualize hallmark gene set activity at single-cell resolution.

The pipeline integrates curated hallmark gene sets (derived from published literature) with scRNA-seq data, enabling the exploration of cancer-related biological processes across cell populations.

---

## Key Features

* Import hallmark gene sets from a structured Excel file
* Compute module scores using Seurat (`AddModuleScore`)
* Visualize hallmark activity using:

  * UMAP projections
  * Boxplots across cell types
* Perform statistical comparisons (Wilcoxon rank-sum test)
* Minimal, transparent, and easily adaptable pipeline

---

## Input Format

The hallmark gene set file must be an Excel file with the following structure:

| Hallmarks  | Pathways | Genes               |
| ---------- | -------- | ------------------- |
| Hallmark_A | Pathway1 | gene1; gene2; gene3 |

**Column description:**

* `Hallmarks`: Name of each hallmark gene set
* `Pathways`: Associated pathways (informational only, not used in scoring)
* `Genes`: Gene symbols separated by semicolons (`;`)

---

## Installation

Install required R packages:

```r
install.packages(c("Seurat", "tidyverse", "readxl", "ggpubr"))
```

---

## Usage

### 1. Prepare your Seurat object

The pipeline assumes that a Seurat object named:

```r
seurat_obj
```

is already loaded in your R session.

---

### 2. Edit input paths inside the script

Open:

```r
hallmark_scoring_pipeline.R
```

and adjust:

```r
hallmark_file <- "data/hallmarks_genes.xlsx"
group_column <- "celltype_v2"
```

---

### 3. Run the pipeline

```r
source("hallmark_scoring_pipeline.R")
```

---

## Outputs

The pipeline generates:

* **UMAP visualization** of hallmark activity
* **Boxplots with statistical testing** across groups
* **Updated Seurat object** with hallmark scores in metadata

### Example outputs

#### UMAP projection

![UMAP](outputs/UMAP_Hallmark_Scores.png)

---

#### Hallmark activity across cell types

![Boxplots](outputs/Boxplots_Hallmarks.png)

---

## Notes

* Gene symbols must match those in the Seurat object (case-sensitive)
* Standardizing gene names (e.g., uppercase) is recommended
* The `Pathways` column is retained for annotation but not used in computations
* The script is designed as a **minimal pipeline**, not a package

---

## Gene Set Source

The hallmark gene sets used in this repository were derived from:

**Sibai M, Cervilla S, Grases D, et al.**
*The spatial landscape of cancer hallmarks reveals patterns of tumor ecological dynamics and drug sensitivity.*
Cell Reports, 2025; 44
DOI: https://doi.org/10.1016/j.celrep.2024.115229

Please cite the original publication when using these gene sets.

---

## Citation

If you use this repository in your research, please cite:

* This repository (Renan Simões)
* The original hallmark gene set publication (Sibai et al., 2025)
* Seurat

---

## Repository Structure

```
scCancerHallmarks/
│
├── hallmark_scoring_pipeline.R
├── README.md
│
├── data/
│   └── hallmarks_genes.xlsx
│
├── outputs/
│   ├── UMAP_Hallmark_Scores.png
│   └── Boxplots_Hallmarks.png
```

---

## Author

Renan Simões

---

## License

This project is licensed under the MIT License.
