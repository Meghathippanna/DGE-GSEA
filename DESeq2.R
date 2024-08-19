# Set your working directory (update this with your actual directory path)
setwd("path/to/directory/")

# Load the count files
cts <- read.csv('count_file.csv', row.names = 1, header = TRUE)

# Examine the count matrix
head(cts)

# Remove counts that are equal to zero
counts <- cts[which(rowSums(cts) > 0),]

# Load the sample sheet (sample information table)
coldata <- read.csv('samplesheet.csv', row.names = 1, header = TRUE)
coldata

# Ensure the columns in cts match the rows in coldata
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# Load DESeq2 library
library("DESeq2")

# Convert raw count data to DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Relevel condition factor
dds$condition <- factor(dds$condition, levels = c("control", "case"))
dds$condition <- droplevels(dds$condition)

# Run DESeq analysis
dds <- DESeq(dds)
res <- results(dds)

# Specify the contrast for results (case vs control)
res <- results(dds, contrast = c("condition", "case", "control"))

# Order results by adjusted p-value (padj)
resOrdered <- res[order(res$padj),]

# Display summary of results
summary(resOrdered)

# Convert results to a dataframe and add gene symbols
res.df <- as.data.frame(resOrdered)
library(org.Hs.eg.db)
res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = "ENSEMBL", column= "SYMBOL")

# View and save the ordered results
View(res.df)
write.csv(as.data.frame(res.df), file="All_candidates_annotated.csv")

# Count the number of significantly differentially expressed genes (adjusted p-value < 0.05)
sum(res.df$padj < 0.05, na.rm = TRUE)

# Select significant results (adjusted p-value < 0.05)
resSig <- subset(resOrdered, padj < 0.05)
resSig
dim(resSig)

# Plot MA for the results
plotMA(resOrdered, ylim = c(-2, 2))

# Plot gene counts for the most significant gene
d <- plotCounts(dds, gene = which.min(resOrdered$padj), intgroup = "condition", returnData = TRUE)
library("ggplot2")
ggplot(d, aes(x = condition, y = count)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) + 
  scale_y_log10(breaks = c(25, 100, 400))

# Normalized counts for selected genes
mat <- counts(dds, normalized = TRUE)[rownames(resSig),]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)

# Generate heatmap for selected genes
library(ComplexHeatmap)
Heatmap(mat.z, cluster_rows = TRUE, cluster_columns = TRUE, 
        column_labels = colnames(mat.z), name = "Z-score", 
        row_labels = resSig[rownames(mat.z), ]$symbol, 
        row_names_gp = gpar(fontsize = 6))

# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res.df, lab = res.df$symbol, x = 'log2FoldChange', y = 'padj',
                title = 'Case versus Control',
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 6.0)
