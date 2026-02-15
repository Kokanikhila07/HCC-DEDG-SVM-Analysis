# =======================================================
# 1. SETUP AND DATA LOADING (GSE214846 HCC)
# =======================================================

# Load necessary libraries
library(DESeq2)
library(umap)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(stringr)

# Load count data 
df <- read.delim("GSE214846_count.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Load metadata 
metadata <- read.delim("metadata.txt", header = TRUE, sep = "\t")

# --- Data Preprocessing ---
colnames(df) <- trimws(colnames(df))
genes <- make.unique(as.character(df[,1]))
rownames(df) <- genes
count_data <- df[ , -1] # Remove the first column (gene names)
count_data <- as.matrix(count_data)

# Match metadata to count data columns
stopifnot(all(colnames(count_data) == metadata$Filename))

# Create a group factor
group <- factor(metadata$sample)
# Define the reference level (Control is often the reference)
group <- relevel(group, ref = "Control")

# =======================================================
# QUERY ANSWER: Extract Column Names into a Single Column DF
# =======================================================

# Create a data frame with all sample names in a single column
sample_names_df <- data.frame(Sample_ID = colnames(count_data))
cat("\nSample names extracted to 'sample_names_df':\n")
print(head(sample_names_df)) # Displaying the answer to your query

# =======================================================
# 2. DESEQ2 ANALYSIS AND NORMALIZATION
# =======================================================

# Create DESeqDataSet (DESeq2 requires count data to be integers)
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = DataFrame(group = group),
                              design = ~ group)

# Run DESeq2 pipeline (CRITICAL MISSING STEP ADDED)
dds <- DESeq(dds)

# Extract non-normalized counts (log-transformed raw counts)
dat_non_normalized <- log2(counts(dds, normalized = FALSE) + 1)

# Extract VST-normalized counts (for boxplots, UMAP, clustering) (CRITICAL MISSING STEP ADDED)
# VST is Variance Stabilizing Transformation, excellent for visualization
dat_normalized <- assay(vst(dds, blind = TRUE))

# Remove duplicate genes (if any)
dat_non_normalized <- dat_non_normalized[!duplicated(rownames(dat_non_normalized)), ]
dat_normalized <- dat_normalized[!duplicated(rownames(dat_normalized)), ]

# =======================================================
# 3. QUALITY CONTROL PLOTS (Boxplots & UMAP)
# =======================================================

# --- Common Plot Setup ---
ord <- order(group)
lbl_non <- "log2(raw counts + 1)"
lbl_norm <- "VST Normalized Counts"

group_colors_map <- c("Control" = "#1B9E77", "Tumor" = "#7570B3")
col_vector <- group_colors_map[group]

# --- 3a. Box-and-Whisker Plot (Non-normalized Counts) ---
par(mar=c(7,4,2,1))
boxplot(dat_non_normalized[, ord], boxwex=0.6, notch = FALSE,
        main = "Box Plot (Non-normalized Counts)",
        ylab = lbl_non,
        outline = FALSE,
        las = 2,
        col = col_vector[ord])
legend("topleft", legend = names(group_colors_map), fill = group_colors_map, bty = "n")


# --- 3b. Box-and-Whisker Plot (Normalized Counts) ---
par(mar=c(7,4,2,1))
boxplot(dat_normalized[, ord], 
        boxwex = 0.6, 
        notch = FALSE, 
        main = "Box-and-Whisker Plot Normalized Counts (VST)",
        ylab = lbl_norm, 
        outline = FALSE, 
        las = 2, 
        col = col_vector[ord])
legend("topleft", legend = names(group_colors_map), fill = group_colors_map, bty = "n")

# =======================================================
# ✅ UMAP Plots (Normalized + Non-normalized)
# =======================================================

set.seed(123)  # For reproducibility

# --- UMAP: Non-normalized Counts ---
ump_non_normalized <- umap(t(dat_non_normalized), n_neighbors = 6)

par(mar = c(3, 3, 2, 8), xpd = TRUE)
plot(ump_non_normalized$layout,
     main = "UMAP Plot (Non-normalized Counts)",
     xlab = "UMAP 1", ylab = "UMAP 2",
     col = col_vector,
     pch = 19, cex = 1.8)
text(ump_non_normalized$layout[,1],
     ump_non_normalized$layout[,2],
     labels = colnames(dat_non_normalized),
     pos = 3, cex = 0.8)

legend("topright", inset = c(-0.22, 0), legend = levels(group),
       fill = group_colors_map, pch = 19, title = "Group")


# --- UMAP: VST Normalized Counts ---
ump_normalized <- umap(t(dat_normalized), n_neighbors = 6)

par(mar = c(3, 3, 2, 8), xpd = TRUE)
plot(ump_normalized$layout,
     main = "UMAP Plot (Normalized Counts)",
     xlab = "UMAP 1", ylab = "UMAP 2",
     col = col_vector,
     pch = 19, cex = 1.8)
text(ump_normalized$layout[,1],
     ump_normalized$layout[,2],
     labels = colnames(dat_normalized),
     pos = 3, cex = 0.8)

legend("topright", inset = c(-0.22, 0), legend = levels(group),
       fill = group_colors_map, pch = 19, title = "Group")



# =======================================================
# 4. DIFFERENTIAL EXPRESSION RESULTS & VOLCANO PLOT
# =======================================================

# Extract results 
results_deseq <- results(dds)
# Remove NA rows from results (often caused by low count genes)
results_deseq <- results_deseq[!is.na(results_deseq$padj), ]


# --- 4a. Significant Genes Filtering and Saving ---
res_df <- as.data.frame(results_deseq)

# Filter: padj < 0.05
significant_genes <- res_df[res_df$padj < 0.05, ]
write.csv(significant_genes, file = "significant_genes.csv", row.names = TRUE)
nrow(significant_genes)
#[1] 15748

# Filter: padj < 0.05 and |logFC| > 1
significant_genes_1 <- res_df[
  which(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05), ]
significant_genes_1_df <- as.data.frame(significant_genes_1)
significant_genes_1_df$gene_id <- rownames(significant_genes_1_df)
write.csv(significant_genes_1_df, file = "significant_genes_logfc1_pval_0.05.csv", row.names = FALSE)
nrow(significant_genes_1_df)
#[1] 6758

# Extract genes with unadjusted p < 0.01 and |logFC| > 1 (or 2)
significant_genes$log2FoldChange <- as.numeric(as.character(significant_genes$log2FoldChange))
significant_genes$padj <- as.numeric(as.character(significant_genes$padj))

significant_genes_2 <- significant_genes[
  which(abs(significant_genes$log2FoldChange) > 1 & significant_genes$padj < 0.01), ]
nrow(significant_genes_2) 
# [1] 5980

# Save the merged annotated results
write.csv(significant_genes_2, file = "significant_genes_logfc2_pval_0.01.csv", row.names = TRUE)

# Extract genes with adjusted p < 0.05 and |logFC| > 2 (or 3)
significant_genes$log2FoldChange <- as.numeric(as.character(significant_genes$log2FoldChange))
significant_genes$padj <- as.numeric(as.character(significant_genes$padj))

significant_genes_3 <- significant_genes[
  which(abs(significant_genes$log2FoldChange) > 2 & significant_genes$padj < 0.05), ]
nrow(significant_genes_3)
#[1] 2276

# Extract genes with adjusted p < 0.01 and |logFC| > 2 (or 3)
significant_genes$log2FoldChange <- as.numeric(as.character(significant_genes$log2FoldChange))
significant_genes$padj <- as.numeric(as.character(significant_genes$padj))

significant_genes_3 <- significant_genes[
  which(abs(significant_genes$log2FoldChange) > 2 & significant_genes$padj < 0.01), ]
nrow(significant_genes_3)
# 2143

# Save the merged annotated results
write.csv(significant_genes, file = "annotated_significant_genes_3.csv", row.names = FALSE)

# Extract genes with adjusted p < 0.05 and |logFC| > 3 (or 4)
significant_genes$log2FoldChange <- as.numeric(as.character(significant_genes$log2FoldChange))
significant_genes$padj <- as.numeric(as.character(significant_genes$padj))

significant_genes_4 <- significant_genes[
  which(abs(significant_genes$log2FoldChange) > 3 & significant_genes$padj < 0.05), ]
nrow(significant_genes_4)
# 829

# Save the merged annotated results
write.csv(significant_genes, file = "annotated_significant_genes_4.csv", row.names = FALSE)


# --- 4b. Volcano Plot ---

library(EnhancedVolcano)

# Assuming 'results_normalized' contains your DESeq2 results
results <- res_df

results$log2FoldChange

head(rownames(res_df))
lab = rownames(res_df)


# Set a significance threshold (adjust as needed)
significance_threshold <- 0.05

# Create a volcano plot with gene names as labels
volcano_plot <- EnhancedVolcano(
  toptable = results,
  lab = rownames(res_df),
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 (Fold Change)",
  ylab = "-Log10 (Adjusted P-value)",
  FCcutoff = 1,
  pCutoff = significance_threshold,
  pointSize = 1,
  labSize = 3,
  xlim = c(-10, 10),
  ylim = c(0, 100),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  title = "Volcano Plot",
  subtitle = "Differential gene expression",
  legendLabels = c("Not Significant", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)),
  legendPosition = "bottom"
)

# Display the plot
print(volcano_plot)

# =======================================================
# 6. CLUSTERING
# =======================================================

# --- Hierarchical Clustering (Non-Normalized) ---
hc_non_norm <- hclust(as.dist(1-cor(dat_non_normalized, method="spearman")), method="complete")
plot(as.dendrogram(hc_non_norm),
     main = "Sample Clustering (Non-normalized, Spearman)",
     ylab = "Height")

# --- Hierarchical Clustering (Normalized) ---
hc_norm <- hclust(dist(t(dat_normalized), method="euclidean"), method="complete")
plot(as.dendrogram(hc_norm),
     main = "Sample Clustering (Normalized, Euclidean)",
     ylab = "Height")

# --- PVCLUST (Requires pvclust package) ---
require(pvclust)
pv <- pvclust(dat_normalized, method.dist="euclidean", method.hclust="ward.D2", nboot=100)
plot(pv)

matrix = dat_normalized
distance = dist(t(matrix),method="maximum")
clusters = hclust(distance)
plot(clusters)

##########################################################
#  SVM-Based Identification of DEDGs
##########################################################

library(e1071)
library(caret)
library(dplyr)

# -------------------------------------------------
# 1️⃣ Load DEGs from Uploaded File
# -------------------------------------------------

# Load your DEG file
significant_genes <- read.csv("significant_genes_logfc1_pval_0.05.csv")

# Extract gene list
deg_gene_list <- significant_genes$gene_id  # or SYMBOL if present

# Ensure genes exist in VST-normalized matrix
deg_gene_list <- intersect(deg_gene_list, rownames(dat_normalized))

# Filter normalized matrix for only the selected DEGs
svm_input_expression_final <- dat_normalized[deg_gene_list, ]

# Transpose for ML input (samples as rows)
svm_data_all <- as.data.frame(t(svm_input_expression_final))

# Add class labels (Control/Tumor)
valid_classes <- group
levels(valid_classes) <- make.names(levels(valid_classes))
svm_data_all$Class <- valid_classes

cat("Total DEGs loaded:", length(deg_gene_list), "\n")
# Total DEGs loaded: 6758 
cat("Samples available:", nrow(svm_data_all), "\n\n")
# Samples available: 130 
# -------------------------------------------------
# 2️⃣ Hyperparameter Tuning (C & Gamma)
# -------------------------------------------------
set.seed(123)

tuneGrid <- expand.grid(
  C = 2^seq(-2, 4, 2),
  sigma = 2^seq(-4, 0, 2)
)

trControl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Use first 100 genes for tuning (speed)
tuning_subset <- svm_data_all[, 1:min(100, ncol(svm_data_all) - 1)]
tuning_subset$Class <- svm_data_all$Class

cat("🔍 Performing parameter tuning...\n")

svm_tuned <- train(
  Class ~ .,
  data = tuning_subset,
  method = "svmRadial",
  preProc = c("center", "scale"),
  tuneGrid = tuneGrid,
  trControl = trControl,
  metric = "Accuracy"
)

optimal_C <- svm_tuned$bestTune$C
optimal_sigma <- svm_tuned$bestTune$sigma

cat("\n Optimal Parameters → C:", optimal_C, ", Sigma (Gamma):", optimal_sigma, "\n\n")
# Optimal Parameters → C: 0.25 , Sigma (Gamma): 0.0625
# -------------------------------------------------
# 3️⃣ Single-Gene SVM Screening (DEG by DEG)
# -------------------------------------------------

accuracy_results <- data.frame(Gene = rownames(svm_input_expression_final), Accuracy = NA)
ctrl <- trainControl(method = "cv", number = 5)

cat("🚀 Starting single-gene evaluation...\n")

num_genes <- nrow(svm_input_expression_final)

for (i in 1:num_genes) {
  gene_symbol <- rownames(svm_input_expression_final)[i]
  
  # Create single-gene classification dataset
  single_gene_data <- data.frame(
    Expression = svm_input_expression_final[gene_symbol, ],
    Class = valid_classes
  )
  
  # Train SVM on this gene
  model <- train(
    Class ~ Expression,
    data = single_gene_data,
    method = "svmRadial",
    preProc = c("center", "scale"),
    tuneGrid = expand.grid(C = optimal_C, sigma = optimal_sigma),
    trControl = ctrl,
    metric = "Accuracy"
  )
  
  accuracy_results$Accuracy[i] <- max(model$results$Accuracy)
  
  if (i %% 50 == 0) {
    cat("Processed", i, "/", num_genes, "genes...\n")
  }
}

# -------------------------------------------------
# 4️⃣ Filter DEDGs (Accuracy ≥ 80%)
# -------------------------------------------------

accuracy_results_sorted <- accuracy_results %>% arrange(desc(Accuracy))
dedgs <- accuracy_results_sorted %>% filter(Accuracy >= 0.80)

cat("\n Number of final DEDGs with ≥80% accuracy:", nrow(dedgs), "\n")
# Number of final DEDGs with ≥80% accuracy: 1250

# Merge back logFC & padj
dedgs_final <- dedgs %>%
  left_join(significant_genes, by = c("Gene" = "gene_id")) %>%
  arrange(desc(Accuracy))

# Save final results
write.csv(dedgs_final, "GSE214846_DEDGs_SVM_Accuracy80.csv", row.names = FALSE)

cat("\n💾 Final DEDGs saved to: GSE214846_DEDGs_SVM_Accuracy80.csv\n")
