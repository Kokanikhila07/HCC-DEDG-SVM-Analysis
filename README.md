# HCC-DEDG-SVM-Analysis
**Identification of Discriminative Differentially Expressed Genes (DEDGs) in Hepatocellular Carcinoma Using SVM**
Dataset: GSE214846

**Project Overview**

This project focuses on identifying Discriminative Differentially Expressed Genes (DEDGs) in hepatocellular carcinoma (HCC) using RNA-seq data and machine learning.

The dataset used in this study:

GSE214846

Disease: Hepatocellular Carcinoma

Platform: RNA-seq

Samples: 65 HCC tumor samples with paired adjacent normal tissues

The original study developed a 5-IRG survival prediction model.
In this work, we extend the analysis by integrating differential expression analysis with Support Vector Machine (SVM)-based feature selection to identify highly discriminative genes.

Objectives

Perform Differential Gene Expression (DEG) analysis between tumor and normal samples

Train an SVM classifier to distinguish tumor vs. normal samples

**Methodology**

1️⃣ Differential Expression Analysis

Tool: DESeq2

Input: Raw RNA-seq count matrix

Design formula: ~ condition

Filtering criteria:

Adjusted p-value (FDR) < 0.05

|log2FoldChange| > 1

Output:

List of significantly upregulated and downregulated DEGs

Volcano plot visualization

Extract important discriminative genes from the trained model

Identify overlapping genes between DEGs and SVM-selected genes (DEDGs)

2️⃣ SVM-Based Discriminative Gene Identification

Algorithm: Support Vector Machine (Linear Kernel)

Package: e1071 / caret

Data: Variance Stabilized Transformed (VST) counts

Cross-validation: 5-fold

Feature importance derived from SVM weight vector

Genes that are:

Statistically significant DEGs
AND

Strong contributors to classification

→ Defined as Discriminative Differentially Expressed Genes (DEDGs)

**Repository Structure**

HCC_DEG_DEDG_SVM_GSE214846/
│
├── data/
│   ├── raw_counts.csv
│   ├── metadata.csv
│
├── scripts/
│
│   ├── GSE214846 DEDGs Script.R
│
├── results/
|
│   ├── DEG_results.csv
│   ├── DEDGs_SVM_selected.csv
│   ├── volcano_plot.png
│  
└── README.md
