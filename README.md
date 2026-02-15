# HCC-DEDG-SVM-Analysis


# Identification of DEDGs in Hepatocellular Carcinoma (GSE214846) Using SVM

## 📌 Project Overview

This project identifies **Discriminative Differentially Expressed Genes (DEDGs)** in Hepatocellular Carcinoma (HCC). By combining statistical rigor (**DESeq2**) with the classification power of **Support Vector Machines (SVM)**, we identify genes that are both statistically significant and highly effective at distinguishing tumor from normal tissue.

### 🧬 Dataset: GSE214846

* **Disease:** Hepatocellular Carcinoma (HCC)
* **Platform:** RNA-seq (Illumina HiSeq 2500)
* **Sample Size:** 65 HCC patients (65 Tumor + 65 Paired Adjacent Normal)
* **Original Scope:** Developed a 5-Immune-Related Gene (IRG) survival signature.

---

## 🔬 Methodology

The workflow integrates high-throughput statistical testing with machine learning classification to define the final DEDG list:

### 1️⃣ Differential Expression Analysis (DEG)

* **Preprocessing:** Raw counts were normalized using Variance Stabilizing Transformation (VST) to handle heteroscedasticity.
* **Tool:** `DESeq2`
* **Statistical Filtering:** * Adjusted p-value (FDR) < 0.05 |log_2FoldChange| > 1
* **Objective:** Identify genes with biologically significant expression shifts.
* **Results:** Identified 6,758 DEGs representing the primary transcriptomic shift between HCC tumor and adjacent normal tissues.

### 2️⃣ SVM-Based Feature Selection

* **Algorithm:** Support Vector Machine (SVM) with a Radial Basis Function (RBF) Kernel.
* **Hyperparameter Tuning:** A grid search was performed using the caret package to determine optimal parameters (C = 0.25, sigma = 0.0625) via 5-fold cross-validation.
* **Validation Strategy:** 5-fold cross-validation was applied to ensure the model's generalizability and prevent overfitting.

### 3️⃣ Single-Gene Discriminative Screening

Unlike standard feature selection, this pipeline evaluates the individual classification power of each DEG:

* Each of the 6,758 DEGs was treated as an independent feature to train an SVM model.
* **Metric:** Accuracy was calculated for every gene to determine its ability to categorize samples correctly.

### 🎯 Final Identification (DEDGs)

The final list of **Discriminative Differentially Expressed Genes (DEDGs)** was curated based on a dual-threshold:

1. **Statistical Significance:** Must satisfy $DR < 0.05 and |log_2FC| > 1.
2. **Classification Performance:** Must achieve an SVM Cross-Validation Accuracy >80%.

---

## 📂 Repository Structure

```text
HCC_DEG_DEDG_SVM_GSE214846/
├── data/
│   ├── raw_counts.csv           # Raw RNA-seq count matrix
│   └── metadata.csv             # Sample info (Tumor vs Normal)
├── scripts/
│   └── GSE214846_DEDGs_Script.R # Full pipeline: DESeq2 -> SVM -> DEDGs
├── results/
│   ├── DEG_results.csv          # Full list of differentially expressed genes
│   ├── DEDGs_SVM_selected.csv   # Final list of discriminative genes
│   └── volcano_plot.png         # Visualization of significance vs fold change
└── README.md

```

---

### Prerequisites

Ensure you have **R (>= 4.0)** installed with the following packages:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("e1071", "caret", "ggplot2", "pROC", "dplyr"))

```




**Would you like me to help you draft the "Results" section summary once you have the final number of DEDGs from your script?**
