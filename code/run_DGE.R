# load libraries 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)

# library(presto)
# library(scCustomize)
# library(SummarizedExperiment)
# library(RColorBrewer)
# library(circlize)
# library(tidyr)


seu_obj <- readRDS("rds_files/seu_for_DGE.rds")

# ---- Consolidated Pseudobulking and DESeq2 Function ----
run_pseudobulk_deg <- function(seu_obj, cell_type, min_counts = 10, alpha = 0.05,
                               plot_title = NULL, save_results = FALSE) {
  
  # Set assay to RNA for raw counts
  DefaultAssay(seu_obj) <- "RNA"
  
  # Subset to specified cell type
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) {
    stop(paste("No cells found for cell type:", cell_type))
  }
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n") # look at, maybe remove
  
  # Pseudobulk Matrix by sample (orig.ident)
  cts <- AggregateExpression(seu_subset, assays = "RNA", group.by = "orig.ident",
                             slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  # Prepare Metadata
  colData <- data.frame(samples = colnames(cts)) %>%
    mutate(
      condition = case_when(
        str_detect(samples, "Blumeria") ~ "MCT-Blumeria",
        str_detect(samples, "Water") ~ "MCT-Water",
        str_detect(samples, "Control") ~ "Control"
      )
    ) %>%
    column_to_rownames(var = "samples")
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = colData,
    design = ~ condition)
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= min_counts
  dds <- dds[keep, ]
  cat("Filtered to", nrow(dds), "genes (min total counts =", min_counts, ")\n")
  
  # Set reference level 
  dds$condition <- relevel(dds$condition, ref = "Control")
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract contrasts
  res_water_vs_ctrl <- results(dds, contrast = c("condition", "MCT-Water", "Control"),
                               alpha = alpha)
  res_blum_vs_ctrl <- results(dds, contrast = c("condition", "MCT-Blumeria", "Control"),
                              alpha = alpha)
  res_blum_vs_water <- results(dds, contrast = c("condition", "MCT-Blumeria", "MCT-Water"), 
                               alpha = alpha)
  res_water_vs_blum <- results(dds, contrast = c("condition", "MCT-Water", "MCT-Blumeria"))
  
  # Print summaries
  cat("\n=== DESeq2 Results Summary for", cell_type, "===\n")
  cat("Water vs Control:\n"); summary(res_water_vs_ctrl)
  cat("\nBlumeria vs Control:\n"); summary(res_blum_vs_ctrl)
  cat("\nBlumeria vs Water:\n"); summary(res_blum_vs_water)
  cat("\nWater vs Blumeria:\n"); summary(res_water_vs_blum)
  
  # Normalize data for PCA
  rld <- rlog(dds, blind = FALSE)
  
  # Generate plots
  title <- ifelse(is.null(plot_title), cell_type, plot_title)
  
  p_pca <- plotPCA(rld, intgroup = "condition") +
    stat_ellipse(level = 0.95, alpha = 0.3, geom = "polygon") +
    ggtitle(paste("Pseudoulk PCA:", title))
  
  p_ma_water <- plotMA(res_water_vs_ctrl) + ggtitle("MA: Water vs Control")
  p_ma_blum_water <- plotMA(res_blum_vs_water) + ggtitle("MA: Blumeria vs Water")
  p_ma_water_blum <- plotMA(res_water_vs_blum) + ggtitle("MA: Water vs Blumeria")
  
  # Display plots
  print(p_pca)
  print(p_ma_water)
  # print(p_ma_blum_water)
  print(p_ma_water_blum)
  
  # Save results if requested
  if (save_results) {
    results_list <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum,
      dds = dds,
      rld = rld
    )
    
    filename <- paste0("DEG_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = paste0("rds_files/", filename))
    cat("Results saved to:", filename, "\n")
    
    # Save plots
    # ggsave(paste0("PCA_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
    # plot = p_pca, width = 8, height = 6, dpi = 300)
    # ggsave(paste0("MA_Water_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
    # plot = p_ma_water, width = 8, height = 6, dpi = 300)
    # ggsave(paste0("MA_BlumWater_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
    # plot = p_ma_blum_water, width = 8, height = 6, dpi = 300)
  }
  
  # Return results
  return(list(
    results = list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum
    ),
    dds = dds,
    rld = rld,
    plots = list(pca = p_pca
                 # , ma_water = p_ma_water, 
                 # ma_blum_water = p_ma_blum_water
    )
  ))
}

cm_results <- run_pseudobulk_deg(seu_obj, "Cardiomyocytes", alpha = 0.1, 
                                 save_results = TRUE)

run_pseudobulk_deg(seu_obj, "Pericytes", alpha = 0.1, 
                   save_results = TRUE)
# ---- Bar Plots of up and down regulated gene counts ----
# Run DSEq2 on all cell types
deg_results <- list()
cell_types <- levels(seu_obj)  
for (cell_type in cell_types) {
  cat("\\\\n=== Processing", cell_type, "===\\\\n")
  deg_results[[cell_type]] <- run_pseudobulk_deg(seu_obj, cell_type, 
                                                 alpha = 0.1, save_results = FALSE)
}

# Extract significant DE gene counts
sig_counts <- data.frame(
  cell_type = character(),
  contrast = character(),
  upregulated = numeric(),
  downregulated = numeric()
)

for (cell_type in cell_types) {
  res <- deg_results[[cell_type]]$results
  
  # Water vs Control
  sig_up_water <- sum(res$water_vs_ctrl$padj < 0.1 & res$water_vs_ctrl$log2FoldChange > 0, na.rm = TRUE)
  sig_down_water <- sum(res$water_vs_ctrl$padj < 0.1 & res$water_vs_ctrl$log2FoldChange < 0, na.rm = TRUE)
  
  # Blumeria vs Control  
  # sig_up_blum_ctrl <- sum(res$blum_vs_ctrl$padj < 0.05 & res$blum_vs_ctrl$log2FoldChange > 0, na.rm = TRUE)
  # sig_down_blum_ctrl <- sum(res$blum_vs_ctrl$padj < 0.05 & res$blum_vs_ctrl$log2FoldChange < 0, na.rm = TRUE)
  
  # Blumeria vs Water
  # sig_up_blum_water <- sum(res$blum_vs_water$padj < 0.05 & res$blum_vs_water$log2FoldChange > 0, na.rm = TRUE)
  # sig_down_blum_water <- sum(res$blum_vs_water$padj < 0.05 & res$blum_vs_water$log2FoldChange < 0, na.rm = TRUE)
  
  # Water vs Blumeria
  sig_up_water_blum <- sum(res$water_vs_blum$padj < 0.1 & res$water_vs_blum$log2FoldChange > 0, na.rm = TRUE)
  sig_down_water_blum <- sum(res$water_vs_blum$padj < 0.1 & res$water_vs_blum$log2FoldChange < 0, na.rm = TRUE)
  
  sig_counts <- rbind(sig_counts, data.frame(
    cell_type = cell_type,
    contrast = "Water vs Ctrl",
    upregulated = sig_up_water,
    downregulated = sig_down_water
  ))
  # sig_counts <- rbind(sig_counts, data.frame(
  # cell_type = cell_type,
  # contrast = "Blumeria vs Ctrl", 
  # upregulated = sig_up_blum_ctrl,
  # downregulated = sig_down_blum_ctrl
  # ))
  # sig_counts <- rbind(sig_counts, data.frame(
  # cell_type = cell_type,
  # contrast = "Blumeria vs Water",
  # upregulated = sig_up_blum_water,
  # downregulated = sig_down_blum_water
  # ))
  sig_counts <- rbind(sig_counts, data.frame(
    cell_type = cell_type,
    contrast = "Water vs Blumeria",
    upregulated = sig_up_water_blum,
    downregulated = sig_down_water_blum
  ))
}

# Create bar plot
sig_counts_long <- sig_counts %>%
  pivot_longer(cols = c(upregulated, downregulated), 
               names_to = "direction", values_to = "count")

# Use your cell type order and colors
sig_counts_long$cell_type <- factor(sig_counts_long$cell_type, levels = rev(cell_types))
sig_counts_long$direction <- factor(sig_counts_long$direction, 
                                    levels = c("upregulated", "downregulated"))
sig_counts_long$contrast <- factor(sig_counts_long$contrast, 
                                   levels = c("Water vs Ctrl", 
                                              # "Blumeria vs Ctrl", 
                                              # "Blumeria vs Water",
                                              "Water vs Blumeria"
                                   ))

p_sig <- ggplot(sig_counts_long, aes(x = cell_type, y = count, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ contrast, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("upregulated" = "#E31A1C", "downregulated" = "#1F78B4")) +
  labs(title = "Significant DE Genes (padj < 0.1) Across Cell Types",
       x = "Cell Type", y = "Number of DE Genes",
       fill = "Direction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_sig)

ggsave("output/DE_Genes_Bar_pval_1.png",plot = p_sig, width = 12, height = 6, dpi = 300)

# ----MA Plots of up and down regulated gene counts ----

# Run DSEq2 on all cell types
deg_results <- list()
for (cell_type in cell_types) {
  cat("\n=== Processing", cell_type, "===\n")
  deg_results[[cell_type]] <- run_pseudobulk_deg(seu_integrated, cell_type, 
                                                 alpha = 0.1, save_results = FALSE)
}

# Function to create one panel
create_celltype_deg_plot <- function(deg_results, cell_types, contrast_name,
                                     padj_thresh = 0.05, lfc_thresh = 0) {
  
  # Collect data for all cell types
  plot_data <- data.frame()
  up_counts <- integer(length(cell_types))
  down_counts <- integer(length(cell_types))
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    # Extract the appropriate contrast
    #if (constrast_name == "Water vs Ctrl") {
    # res <- deg_results[[cell_type]]$results$water_vs_ctrl
    # } else if (contrast_name == "Blumeria vs Ctrl") {
    # res <- deg_results[[cell_type]]$results$blum_vs_ctrl
    # } else if (contrast_name == "Blumeria vs Water") {
    # res <- deg_results[[cell_type]]$results$blum_vs_water
    # }
    
    # Convert to data frame
    
  }
}
