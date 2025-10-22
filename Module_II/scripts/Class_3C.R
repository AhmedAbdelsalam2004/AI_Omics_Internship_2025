if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load required libraries
library(AnnotationDbi)
library(hgu133plus2.db)  # Annotation package for GPL570 (HG-U133_Plus_2)
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(GEOquery)

# Create Results directory
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Load and process the GEO series matrix file
cat("Loading GEO series matrix file...\n")

# Read the series matrix file
gse <- getGEO(filename = "C:/Users/a3abd/OneDrive/Desktop/AI_Omics_Internship_2025/Module_II/raw_data/GSE79973_series_matrix.txt")

# Extract expression data and phenotype information
expr_matrix <- exprs(gse)
pheno_data <- pData(gse)

# Clean up sample names and create proper phenotype data
cat("Processing phenotype data...\n")

# Extract sample information from titles
sample_titles <- pheno_data$title
tumor_samples <- grepl("tumor", sample_titles, ignore.case = TRUE)
normal_samples <- grepl("normal", sample_titles, ignore.case = TRUE)

# Create clean phenotype data
pheno_data_clean <- data.frame(
  sample = colnames(expr_matrix),
  type = ifelse(tumor_samples, "tumor", 
                ifelse(normal_samples, "normal", NA)),
  patient = gsub(".*patient_([0-9]+).*", "\\1", sample_titles),
  row.names = colnames(expr_matrix)
)

# Remove any samples that couldn't be classified
valid_samples <- !is.na(pheno_data_clean$type)
expr_matrix <- expr_matrix[, valid_samples]
pheno_data_clean <- pheno_data_clean[valid_samples, ]

cat("Data loaded successfully!\n")
cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
cat("Tumor samples:", sum(pheno_data_clean$type == "tumor"), "\n")
cat("Normal samples:", sum(pheno_data_clean$type == "normal"), "\n")

# 1. Map probe IDs to gene symbols
cat("\n=== PROBE TO GENE SYMBOL MAPPING ===\n")

# Get probe IDs
probe_ids <- rownames(expr_matrix)

# Map probe IDs to gene symbols
gene_symbols <- mapIds(hgu133plus2.db, 
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Check mapping results
mapped_probes <- !is.na(gene_symbols)
cat("Probes successfully mapped to gene symbols:", sum(mapped_probes), "/", length(probe_ids), 
    "(", round(sum(mapped_probes)/length(probe_ids)*100, 2), "%)\n")

# 2. Handle multiple probes mapping to the same gene - FIXED APPROACH
cat("\n=== HANDLING DUPLICATE PROBES ===\n")

# Create mapping dataframe
mapping_df <- data.frame(ProbeID = probe_ids, GeneSymbol = gene_symbols, stringsAsFactors = FALSE)
mapping_df <- mapping_df[!is.na(mapping_df$GeneSymbol), ]

# Check for multiple probes per gene
genes_with_multiple_probes <- mapping_df %>% 
  group_by(GeneSymbol) %>% 
  filter(n() > 1) %>%
  ungroup()

cat("Genes with multiple probes:", length(unique(genes_with_multiple_probes$GeneSymbol)), "\n")
cat("Total probe-gene mappings with duplicates:", nrow(genes_with_multiple_probes), "\n")

# Strategy: Keep the probe with highest mean expression for each gene - IMPROVED
# First, get expression data for mapped probes only
mapped_probe_ids <- mapping_df$ProbeID
expr_mapped <- expr_matrix[mapped_probe_ids, ]

# Calculate mean expression for each probe
probe_means <- rowMeans(expr_mapped)

# Add means to mapping dataframe
mapping_df$MeanExpression <- probe_means

# For each gene, keep the probe with highest mean expression
mapping_df_final <- mapping_df %>%
  group_by(GeneSymbol) %>%
  arrange(desc(MeanExpression)) %>%
  slice(1) %>%  # Keep only the first (highest mean) for each gene
  ungroup()

# Create final expression matrix with unique genes
expr_final <- expr_matrix[mapping_df_final$ProbeID, ]
rownames(expr_final) <- mapping_df_final$GeneSymbol

cat("Final number of unique genes:", nrow(expr_final), "\n")

# 3. Differential Expression Analysis with Limma
cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Define design matrix
groups <- factor(pheno_data_clean$type, levels = c("normal", "tumor"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("normal", "tumor")

# Fit linear model
fit <- lmFit(expr_final, design)

# Define contrast: tumor vs normal
contrast_matrix <- makeContrasts(tumor_vs_normal = tumor - normal, levels = design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
de_results <- topTable(fit2, number = Inf, adjust.method = "BH")
de_results$Gene <- rownames(de_results)

# Add significance flags
de_results$Significant <- ifelse(de_results$adj.P.Val < 0.05 & abs(de_results$logFC) > 1, 
                                 ifelse(de_results$logFC > 1, "Up", "Down"), "Not Sig")

# Summary statistics
upregulated <- sum(de_results$Significant == "Up")
downregulated <- sum(de_results$Significant == "Down")
total_sig <- upregulated + downregulated

cat("Differential Expression Results:\n")
cat("  - Upregulated genes (adj.P.Val < 0.05, logFC > 1):", upregulated, "\n")
cat("  - Downregulated genes (adj.P.Val < 0.05, logFC < -1):", downregulated, "\n")
cat("  - Total significant genes:", total_sig, "\n")

# 4. Create Volcano Plot
cat("\n=== CREATING VOLCANO PLOT ===\n")

volcano_plot <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
  labs(title = "Volcano Plot: Tumor vs Normal Gastric Tissue",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save volcano plot
png("Results/volcano_plot.png", width = 1000, height = 800, res = 150)
print(volcano_plot)
dev.off()

# 5. Create Heatmap of Top 25 DEGs - FIXED VERSION
cat("\n=== CREATING HEATMAP OF TOP 25 DEGS ===\n")

# Get top 25 most significant genes (by adjusted p-value)
top_genes <- de_results %>%
  filter(Significant != "Not Sig") %>%
  arrange(adj.P.Val) %>%
  head(25)

if(nrow(top_genes) > 0) {
  # Check which genes actually exist in expr_final
  available_genes <- top_genes$Gene[top_genes$Gene %in% rownames(expr_final)]
  cat("Genes available for heatmap:", length(available_genes), "/", nrow(top_genes), "\n")
  
  if(length(available_genes) > 0) {
    # Extract expression data for available top genes
    heatmap_data <- expr_final[available_genes, , drop = FALSE]
    
    # Scale the data (z-score by row)
    heatmap_data_scaled <- t(scale(t(heatmap_data)))
    
    # Create annotation for samples
    sample_annotation <- HeatmapAnnotation(
      Type = pheno_data_clean$type,
      col = list(Type = c("tumor" = "red", "normal" = "blue")),
      annotation_name_side = "left"
    )
    
    # Create heatmap
    heatmap_plot <- Heatmap(heatmap_data_scaled,
                            name = "Z-score",
                            top_annotation = sample_annotation,
                            column_title = paste("Top", length(available_genes), 
                                                 "Differentially Expressed Genes\nTumor vs Normal Gastric Tissue"),
                            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                            show_row_names = TRUE,
                            show_column_names = FALSE,
                            row_names_gp = gpar(fontsize = 8),
                            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                            cluster_columns = TRUE,
                            cluster_rows = TRUE)
    
    # Save heatmap
    png("Results/heatmap_top25_deg.png", width = 1200, height = 800, res = 150)
    draw(heatmap_plot)
    dev.off()
    cat("Heatmap created successfully with", length(available_genes), "genes.\n")
  } else {
    cat("No significant genes available for heatmap after filtering.\n")
  }
} else {
  cat("No significant genes found for heatmap.\n")
}

# 6. Save DEG results to CSV files
cat("\n=== SAVING RESULTS ===\n")

# Complete results
write.csv(de_results, "Results/complete_deg_results.csv", row.names = FALSE)

# Upregulated genes
upregulated_genes <- de_results %>% 
  filter(Significant == "Up") %>%
  arrange(desc(logFC))
write.csv(upregulated_genes, "Results/upregulated_genes.csv", row.names = FALSE)

# Downregulated genes
downregulated_genes <- de_results %>% 
  filter(Significant == "Down") %>%
  arrange(logFC)
write.csv(downregulated_genes, "Results/downregulated_genes.csv", row.names = FALSE)

cat("Results saved to CSV files in Results folder.\n")

# 7. Result Summary
cat("\n=== RESULT SUMMARY ===\n")
cat("PROBE MAPPING SUMMARY:\n")
cat("- Multiple probes can map to the same gene due to alternative splicing, \n")
cat("  different transcript variants, or redundant probe design.\n")
cat("- Handled duplicates by keeping the probe with highest mean expression\n")
cat("  for each gene to ensure one unique gene per row.\n\n")

cat("DIFFERENTIAL EXPRESSION ANALYSIS:\n")
cat("- Comparison performed: Tumor vs Normal gastric tissue\n")
cat("- Contrast: tumor_vs_normal (tumor samples compared to normal samples)\n\n")

cat("DEG RESULTS SUMMARY:\n")
cat("- Upregulated genes (logFC > 1, adj.P.Val < 0.05):", upregulated, "\n")
cat("- Downregulated genes (logFC < -1, adj.P.Val < 0.05):", downregulated, "\n")
cat("- Total significant differentially expressed genes:", total_sig, "\n\n")

cat("FILES GENERATED:\n")
cat("- Results/volcano_plot.png: Volcano plot visualization\n")
cat("- Results/heatmap_top25_deg.png: Heatmap of top DEGs\n")
cat("- Results/complete_deg_results.csv: Complete DEG results\n")
cat("- Results/upregulated_genes.csv: Upregulated genes only\n")
cat("- Results/downregulated_genes.csv: Downregulated genes only\n")

# Additional diagnostic information
cat("\n=== DIAGNOSTIC INFORMATION ===\n")
cat("First few genes in de_results:", head(de_results$Gene), "\n")
cat("First few genes in expr_final:", head(rownames(expr_final)), "\n")
cat("Number of common genes:", sum(de_results$Gene %in% rownames(expr_final)), "\n")