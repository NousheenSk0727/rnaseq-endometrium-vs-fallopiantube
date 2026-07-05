# ENDOMETRIUM vs FALLOPIAN TUBE RNAseq ANALYSIS PIPELINE
# Differential Expression, GO Enrichment & Biomarker Identification

# LIBRARIES

library(data.table)
library(DESeq2)
library(pheatmap)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(randomForest)
library(caret)
library(ggplot2)

# STEP 1: Load Count Files
list.files("/Users/nousheenjahanshaik/Documents/")
setwd("/Users/nousheenjahanshaik/Documents/SystemsBiology /counting")
# Read all 6 HTSeq count files
count_4a <- fread("endometrium_4a.s.count.txt", data.table = FALSE, header = FALSE)
count_4b <- fread("endometrium_4b.s.count.txt", data.table = FALSE, header = FALSE)
count_5a <- fread("endometrium_5a.s.count.txt", data.table = FALSE, header = FALSE)
count_f8e <- fread("fallopiantube_8e.s.count.txt", data.table = FALSE, header = FALSE)
count_f8b <- fread("fallopiantube_8b.s.count.txt", data.table = FALSE, header = FALSE)
count_f8d <- fread("fallopiantube_8d.s.count.txt", data.table = FALSE, header = FALSE)

# Rename columns for clarity
colnames(count_4a) <- c("gene_id", "endometrium_4a")
colnames(count_4b) <- c("gene_id", "endometrium_4b")
colnames(count_5a) <- c("gene_id", "endometrium_5a")
colnames(count_f8e) <- c("gene_id", "fallopian_8e")
colnames(count_f8b) <- c("gene_id", "fallopian_8b")
colnames(count_f8d) <- c("gene_id", "fallopian_8d")


# ============================================================
# STEP 2: Define Sample Metadata
# ============================================================

# File names for all 6 samples
files <- c(
  "endometrium_4a.s.count.txt",
  "endometrium_4b.s.count.txt",
  "endometrium_5a.s.count.txt",
  "fallopiantube_8e.s.count.txt",
  "fallopiantube_8b.s.count.txt",
  "fallopiantube_8d.s.count.txt"
)

# Condition labels — 3 endometrium, 3 fallopian tube
organ <- c(
  "endometrium", "endometrium", "endometrium",
  "fallopiantube", "fallopiantube", "fallopiantube"
)

# Sample names for plot labeling
samples <- c(
  "endometrium_4a", "endometrium_4b", "endometrium_5a",
  "fallopiantube_8e", "fallopiantube_8b", "fallopiantube_8d"
)

# DESeq2 BLOCK
# Normalization → rlog Transformation → PCA → Heatmaps
# STEP 3: Create DESeq2 Dataset

# Build sample metadata table
sampleTable <- data.frame(
  sampleName = samples,
  fileName   = files,
  condition  = organ
)

# Create DESeqDataSet from HTSeq count files
# design = ~condition compares by tissue type
project_data1 <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  design      = ~condition
)

# Run DESeq2 normalization
project_data1 <- DESeq(project_data1)

# Regularized log transformation — required before PCA and heatmaps
rld_project <- rlog(project_data1, blind = FALSE)

# STEP 4: PCA Plot — Sample Overview

# Visualize sample separation by tissue type
# Confirms endometrium and fallopian tube cluster separately
plotPCA(rld_project, intgroup = "condition")

# STEP 5: Sample Distance Heatmap
# Calculate pairwise distances between all 6 samples
sampleDists <- dist(t(assay(rld_project)))
sampleDistMatrix <- as.matrix(sampleDists)

# Label rows and columns with sample names
rownames(sampleDistMatrix) <- rownames(colData(rld_project))
colnames(sampleDistMatrix) <- rownames(colData(rld_project))

# Plot sample distance heatmap
# Shows how similar/different each sample is from every other sample
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists
)

# STEP 6: Top Variable Genes Heatmap

# Annotation for heatmap columns (condition label per sample)
df <- as.data.frame(colData(rld_project)[, c("condition"), drop = FALSE])

# Select top 20 most variable genes across all samples
topVarGenes <- head(order(rowVars(assay(rld_project)), decreasing = TRUE), 30)
mat <- assay(rld_project)[topVarGenes, ]
mat <- mat - rowMeans(mat) # center data around zero

# Plot gene expression heatmap
# Shows which genes drive the difference between the two tissues
pheatmap(mat,
  annotation_col = df,
  labels_col     = rownames(colData(rld_project))
)

# edgeR
# Differential Expression → Diagnostics → GO Enrichment → Biomarker ML
# STEP 7: Create edgeR DGEList
# Read count files into edgeR DGEList object
# group = organ assigns condition; labels = sample names
project_data2 <- readDGE(files, group = organ, labels = samples)

# STEP 8: MDS Plot — edgeR Sample Clustering

# Calculate MDS coordinates
mds <- plotMDS(project_data2, plot = FALSE)

# Expand plot boundaries so labels have room
x_range <- range(mds$x) + c(-1.5, 1.5)
y_range <- range(mds$y) + c(-0.8, 0.8)

# Plot points only
plot(mds$x, mds$y,
  col  = c(rep("blue", 3), rep("red", 3)),
  pch  = 16,
  cex  = 1.5,
  xlim = x_range,
  ylim = y_range,
  xlab = "Leading logFC dim 1",
  ylab = "Leading logFC dim 2",
  main = "MDS Plot — Endometrium vs Fallopian Tube"
)

# Add labels with offsets to prevent overlap
x_offset <- c(0.2, 0.2, 0.2, 0.1, 0.1, 0.1)
y_offset <- c(0.15, 0.15, 0.15, 0.15, -0.15, -0.15)

text(mds$x + x_offset, mds$y + y_offset,
  labels = samples,
  col    = c(rep("blue", 3), rep("red", 3)),
  cex    = 0.75
)

legend("topright",
  legend = c("Endometrium", "Fallopian Tube"),
  col    = c("blue", "red"),
  pch    = 16
)

# STEP 9: Differential Expression Analysis

# Normalize library sizes
project_data2 <- calcNormFactors(project_data2)

# Set fallopian tube as reference group
organ_factor <- relevel(factor(organ), ref = "fallopiantube")
design_matrix <- model.matrix(~organ_factor)

# Estimate dispersion and fit quasi-likelihood model
project_data2 <- estimateDisp(project_data2, design_matrix)
fit <- glmQLFit(project_data2, design_matrix)
qlf <- glmQLFTest(fit)

# Extract all significant DEGs (FDR < 0.05)
diffExpGenes_project <- topTags(qlf, n = Inf, p.value = 0.05)

# Summary of up/down/not significant genes
de <- decideTests(qlf, p.value = 0.05)
summary(de)

# View top 20 DEGs sorted by FDR
head(diffExpGenes_project$table[order(diffExpGenes_project$table$FDR), ], 20)

# STEP 9C: Volcano Plot — DEG Overview

library(ggrepel)

# Prepare volcano plot data from edgeR results
volcano_data <- as.data.frame(topTags(qlf, n = Inf))
volcano_data$gene <- rownames(volcano_data)

# Classify genes as up, down or not significant
volcano_data$significance <- ifelse(
  volcano_data$FDR < 0.05 & volcano_data$logFC > 1, "Upregulated",
  ifelse(volcano_data$FDR < 0.05 & volcano_data$logFC < -1, "Downregulated", "Non-Significant")
)

# Select top 15 genes to label by FDR
genes_to_label <- head(volcano_data[order(volcano_data$FDR), "gene"], 15)

# Plot
ggplot(volcano_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1) +
  scale_color_manual(
    values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-Significant" = "gray80")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_text_repel(
    data = subset(volcano_data, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    title = "Endometrium vs Fallopian Tube — Differential Expression",
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "FDR"),
    color = "Differential Expression"
  ) +
  theme_minimal()

# STEP 9D: HOXA11 Expression Plot — Top Biomarker Candidate

# Prepare HOXA11 expression data from normalized counts
hoxa11_data <- data.frame(
  Treatment  = organ,
  Expression = as.numeric(cpm(project_data2, log = TRUE)["HOXA11", ])
)

# Plot HOXA11 expression per tissue
ggplot(hoxa11_data, aes(x = Treatment, y = Expression)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7, color = "gray40") +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(
    title = "HOXA11 Expression — Top Statistical Biomarker Candidate",
    x     = "Tissue Type",
    y     = "log2 CPM Expression"
  ) +
  theme_minimal()

# STEP 9E: FOLR1 Expression Plot — Top Fallopian Tube Biomarker Candidate

# Prepare FOLR1 expression data from normalized counts
folr1_data <- data.frame(
  Treatment  = organ,
  Expression = as.numeric(cpm(project_data2, log = TRUE)["FOLR1", ])
)

# Plot FOLR1 expression per tissue
ggplot(folr1_data, aes(x = Treatment, y = Expression)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7, color = "gray40") +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(
    title = "FOLR1 Expression — Top Fallopian Tube Biomarker Candidate",
    x     = "Tissue Type",
    y     = "log2 CPM Expression"
  ) +
  theme_minimal()

# STEP 10: BCV and Smear Plots — Expression Diagnostics

# BCV plot — shows biological variation across gene expression levels
project_data2_bcv <- estimateDisp(project_data2)
plotBCV(project_data2_bcv)

# Smear plot — highlights significant DEGs
# Blue lines mark logFC threshold of +/-1
detags <- rownames(diffExpGenes_project$table)
plotSmear(qlf, de.tags = detags)
abline(h = c(-1, 1), col = "blue")

# STEP 11: GO Enrichment Analysis

# Run GO enrichment on significant DEG gene symbols
# ont = "BP" focuses on Biological Process terms
enrich_results_project <- enrichGO(
  gene    = rownames(diffExpGenes_project$table),
  OrgDb   = org.Hs.eg.db,
  ont     = "BP",
  keyType = "SYMBOL"
)

# View top 20 enriched biological pathways
head(enrich_results_project@result[, c("Description", "pvalue", "p.adjust", "Count")], 20)

# Dotplot of top enriched GO terms
enrichplot::dotplot(enrich_results_project)

# Save DEG list for external tools (DAVID / Enrichr) if needed
write.csv(rownames(diffExpGenes_project$table),
  file      = "DE_genes_endometrium_vs_fallopiantube.csv",
  row.names = FALSE
)
# Search GO results for your 3 hypothesis pathways

# All enriched GO terms
all_go <- enrich_results_project@result

# Search for apoptosis related terms
apoptosis <- all_go[
  grep("apoptosis", all_go$Description, ignore.case = TRUE),
  c("Description", "pvalue", "p.adjust", "Count")
]

# Search for cell migration related terms
migration <- all_go[
  grep("migration", all_go$Description, ignore.case = TRUE),
  c("Description", "pvalue", "p.adjust", "Count")
]

# Search for hormonal regulation related terms
hormonal <- all_go[
  grep("hormone", all_go$Description, ignore.case = TRUE),
  c("Description", "pvalue", "p.adjust", "Count")
]

# Print results
print(apoptosis)
print(migration)
print(hormonal)

# STEP 12: Biomarker Identification — Random Forest + LOOCV

# Use top 20 DEGs by FDR as ML features
top_genes <- rownames(
  diffExpGenes_project$table[order(diffExpGenes_project$table$FDR), ]
)[1:20]

# Get log-normalized CPM counts for top genes only
counts_matrix <- as.data.frame(t(cpm(project_data2, log = TRUE)[top_genes, ]))
counts_matrix$condition <- factor(organ)
colnames(counts_matrix) <- make.names(colnames(counts_matrix))

# Train Random Forest with Leave-One-Out Cross Validation
# LOOCV is the most appropriate CV strategy for small sample sizes (n=6)
ctrl <- trainControl(method = "LOOCV")
rf_model <- train(condition ~ .,
  data       = counts_matrix,
  method     = "rf",
  trControl  = ctrl,
  importance = TRUE
)

# Print model accuracy
print(rf_model)

# Variable importance plot — top biomarker candidates
imp <- varImp(rf_model)
plot(imp, top = 20, main = "Top Biomarker Candidates")

# Importance scores table sorted by mean importance across both classes
imp_df <- imp$importance
imp_df$mean_importance <- rowMeans(imp_df)
imp_df <- imp_df[order(imp_df$mean_importance, decreasing = TRUE), ]
print(imp_df)
