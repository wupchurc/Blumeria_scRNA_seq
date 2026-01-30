# load libraries ----
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(presto)
library(patchwork)
library(scCustomize)
library(DESeq2)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Load control dataset from 10X CellRanger ----
data_dirs <- c("data/Analysis/cellranger_Control_SP_11_GEX_FL-Z0041/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_12_GEX_FL-Z0044/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_14_GEX_FL-Z0047/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_15_GEX_FL-Z0054/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_1_GEX_FL-Z0043/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_5_GEX_FL-Z0046/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_6_GEX_FL-Z0053/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_9_GEX_FL-Z0056/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_1_GEX_FL-Z0042/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_3_GEX_FL-Z0045/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_4_GEX_FL-Z0048/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_7_GEX_FL-Z0055/raw_feature_bc_matrix/")
sample_ids <- c("Control11", "Control12", "Control14", "Control15",
                "Water1", "Water5", "Water6", "Water9",
                "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7")

# Create seurat objects ----
obj_list <- list()
for (i in seq_along(data_dirs)) {
  seurat_obj <- Read10X(data.dir = data_dirs[i])
  obj_list[[i]] <- CreateSeuratObject(counts = seurat_obj, min.cells = 3, 
                                      min.features = 200, project = sample_ids[i])
}

# Add condition label
for (i in seq_along(obj_list)) {
  sid <- sample_ids[i]
  
  if (grepl("^Control",  sid)) obj_list[[i]]$condition <- "Control"
  if (grepl("^Water",    sid)) obj_list[[i]]$condition <- "MCT-Water"
  if (grepl("^Blumeria", sid)) obj_list[[i]]$condition <- "MCT-Blumeria"
}
# Doublet rate estimation----
cell.recov.counts <- lapply(obj_list, ncol)
cells_ref  <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
rate_ref   <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)   # in percent
# linear interpolation; extrapolate = TRUE if some samples outside range
dbl.rate <- approx(x = cells_ref,
                   y = rate_ref,
                   xout = cell.recov.counts,
                   method = "linear",
                   rule = 2)$y
names(dbl.rate) <- names(cell.recov.counts)
dbl.rate

# Remove Doublets ----
doubletFinder_obj_list <- list()
for (i in seq_along(obj_list)) {
  seu <- obj_list[[i]]
  
  # QC
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # Preprocessing (required for DoubletFinder)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  # print(ElbowPlot(seu))
  
  # DoubletFinder pK identification (no t-SNE needed)
  sweep.res <- paramSweep(seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # select pK value that corresponds to max bcmvn
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  # double estimation
  nExp <- round(cell.recov.counts[[i]] * (dbl.rate[[i]] / 100))
  # run doubletFinder
  seu <- doubletFinder(seu, PCs = 1:20, pN = 0.25, pK = pK, 
                          nExp = nExp, reuse.pANN = NULL, sct = FALSE)

  doubletFinder_obj_list[[i]] <- seu  
}

saveRDS(doubletFinder_obj_list, file = "doubletFinder_obj.rds")
# doubletFinder_obj_list <- readRDS("doubletFinder_obj.rds")

# Filter ----
seu_singlets_list <- vector("list", length(doubletFinder_obj_list))

for (i in seq_along(doubletFinder_obj_list)) {
  seu <- doubletFinder_obj_list[[i]]
  
  df_cols <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)
  if (length(df_cols) != 1) {
    stop("Expected exactly one DF.classifications column, found: ",
         paste(df_cols, collapse = ", "))
  }
  df_col <- df_cols[1]
  
  # logical index of singlets
  singlet_idx <- seu@meta.data[[df_col]] == "Singlet"
  seu_singlets <- seu[, singlet_idx]
  
  seu_singlets_list[[i]] <- seu_singlets
}

# assume length(sample_ids) == length(seu_singlets_list)
for (i in seq_along(seu_singlets_list)) {
  seu <- seu_singlets_list[[i]]
  prefix <- sample_ids[i]
  colnames(seu) <- paste(prefix, colnames(seu), sep = "_")
  seu_singlets_list[[i]] <- seu
}

# sanity check: no duplicates across all cells
all_cells <- unlist(lapply(seu_singlets_list, colnames))
anyDuplicated(all_cells)  # should be 0

# Skipping merge/integrate layers for v4 workflow ----
# Use the list of singlet objects instead of merged
seu_singlets_list <- lapply(
  seu_singlets_list,
  function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
    x <- ScaleData(x)
    x <- RunPCA(x)
    x
  }
)

anchors <- FindIntegrationAnchors(
  object.list = seu_singlets_list,
  reduction   = "rpca"  # or "cca"
)

seu_integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(seu_integrated) <- "integrated"
seu_integrated <- ScaleData(seu_integrated)
seu_integrated <- RunPCA(seu_integrated)
seu_integrated <- FindNeighbors(seu_integrated, dims = 1:20)
seu_integrated <- FindClusters(seu_integrated, resolution = 0.03)
seu_integrated <- RunUMAP(seu_integrated, dims = 1:20)

# Prepare Data for Export to Azimuth ----
#saveRDS(seu_integrated, file = "seu_integrated_for_azimuth.rds")
seu_integrated_for_azimuth <- readRDS("seu_integrated_for_azimuth.rds")

# TRIM FOR AZIMUTH - keep only essentials (~500MB)
DefaultAssay(seu_integrated_for_azimuth) <- "RNA"  # Azimuth uses RNA counts

# 1. Keep top variable features only
top_features <- VariableFeatures(seu_integrated_for_azimuth)[1:3000]
seu_integrated_for_azimuth <- subset(seu_integrated_for_azimuth, features = top_features)

# 2. Remove scale.data (bloats size massively)
seu_integrated_for_azimuth[["RNA"]]@scale.data <- new("dgCMatrix")

# 3. Keep only key reductions (drop intermediate ones)
Reductions(seu_integrated_for_azimuth)  # see what's there
keep_reductions <- c("pca", "umap")  # Azimuth needs these
for (red in names(Reductions(seu_integrated_for_azimuth))) {
  if (!red %in% keep_reductions) {
    seu_integrated_for_azimuth[[red]] <- NULL
  }
}

# 4. Clear command history and extra metadata
seu_integrated_for_azimuth@misc <- list()
seu_integrated_for_azimuth@commands <- list()

# 5. Save slim version
saveRDS(seu_integrated_for_azimuth, file = "seu_azimuth_ready_500MB.rds")

# Verify size
print(paste("Original size:", round(file.size("seu_integrated_for_azimuth.rds")/1e9, 2), "GB"))
print(paste("Azimuth size: ", round(file.size("seu_azimuth_ready_500MB.rds")/1e9, 2), "GB"))

# Quick structure check
print(paste("Cells:", ncol(seu_integrated_for_azimuth)))
print(paste("Features:", length(rownames(seu_integrated_for_azimuth))))
print("Assays:")
print(Assays(seu_integrated_for_azimuth))
print("Reductions:")
print(Reductions(seu_integrated_for_azimuth))

# ----

seu_integrated <- readRDS("seu_integrated_for_azimuth.rds")

# use RNA counts/normalized data for DE
DefaultAssay(seu_integrated) <- "RNA"
seu_integrated <- JoinLayers(seu_integrated, assay = "RNA")

# set identities to clusters (from integrated clustering)
Idents(seu_integrated) <- "seurat_clusters"

# --- conserved markers per cluster across conditions ---
markers_cluster0  <- FindConservedMarkers(seu_integrated, ident.1 = 0,  grouping.var = "condition")
markers_cluster1  <- FindConservedMarkers(seu_integrated, ident.1 = 1,  grouping.var = "condition")
markers_cluster2  <- FindConservedMarkers(seu_integrated, ident.1 = 2,  grouping.var = "condition")
markers_cluster3  <- FindConservedMarkers(seu_integrated, ident.1 = 3,  grouping.var = "condition")
markers_cluster4  <- FindConservedMarkers(seu_integrated, ident.1 = 4,  grouping.var = "condition")
markers_cluster5  <- FindConservedMarkers(seu_integrated, ident.1 = 5,  grouping.var = "condition")
markers_cluster6  <- FindConservedMarkers(seu_integrated, ident.1 = 6,  grouping.var = "condition")
markers_cluster7  <- FindConservedMarkers(seu_integrated, ident.1 = 7,  grouping.var = "condition")
markers_cluster8  <- FindConservedMarkers(seu_integrated, ident.1 = 8,  grouping.var = "condition")
markers_cluster9  <- FindConservedMarkers(seu_integrated, ident.1 = 9,  grouping.var = "condition")
markers_cluster10 <- FindConservedMarkers(seu_integrated, ident.1 = 10, grouping.var = "condition")

# Renaming clusters 
print(Idents(seu_integrated))
# rename cluster 0
seu_integrated <- RenameIdents(seu_integrated, '0' = 'Fibroblasts')
# rename cluster 1
seu_integrated <- RenameIdents(seu_integrated, '1' = 'Capillary EC')
# rename cluster 2 
seu_integrated <- RenameIdents(seu_integrated, '2' = 'Macrophages')
# rename cluster 3
seu_integrated <- RenameIdents(seu_integrated, '3' = 'Cardiomyocytes')
# rename cluster 4
seu_integrated <- RenameIdents(seu_integrated, '4' = 'Pericytes')
# rename cluster 5 
seu_integrated <- RenameIdents(seu_integrated, '5' = 'T cells')
# rename cluster 6
seu_integrated <- RenameIdents(seu_integrated, '6' = 'Venous EC')
# rename cluster 7
seu_integrated <- RenameIdents(seu_integrated, '7' = 'B cells')
# rename cluster 8
seu_integrated <- RenameIdents(seu_integrated, '8' = 'Lymphatic EC')
# rename cluster 9
seu_integrated <- RenameIdents(seu_integrated, '9' = 'Neutrophils')
# rename cluster 10
seu_integrated <- RenameIdents(seu_integrated, '10' = 'Neuronal')

# Store renamed Idents in metadata column
seu_integrated$cell_type <- Idents(seu_integrated)


# ---- Heatmap of markers used to validate cell type identities ----



all_markers <- c(
  "Pdgfra", "Csf1r", "Col1a1", "Cd68", "Adgre1"  # Add more as needed
)

DoHeatmap(subset(seu_integrated, downsample = 500),
          features = all_markers,
          size = 3,
          hjust = 0.5,
          vjust = - 0.2,
          angle = 0, 
          group.bar.height = 0.06) + 
  guides(color = "none")



 

# --- plotting (still using integrated UMAP) ----
DefaultAssay(seu_integrated) <- "integrated"  # for plotting/clustering context

# Reorder condition for plots
seu_integrated$condition <- factor(
  seu_integrated$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)


# Plot UMAPs with identified clusters
print(DimPlot(seu_integrated, reduction = 'umap',split.by = 'condition', 
              label = TRUE, label.size = 3, repel = TRUE) + NoLegend())

# Relative quantification ----

# shared colors (same order as Idents)
cell_types <- levels(Idents(seu_integrated))
#cell_cols  <- Seurat::DiscretePalette(length(cell_types))
cell_cols <- DiscretePalette_scCustomize(
  num_colors = length(cell_types),
  palette    = "polychrome"   # or "polychrome", "alphabet2", "varibow"
)
names(cell_cols) <- cell_types


origin_x <- -15
origin_y <- -15
len      <- 5   # length of mini-axes in UMAP units

mini_axes <- data.frame(
  x    = c(origin_x, origin_x),
  y    = c(origin_y, origin_y),
  xend = c(origin_x + len, origin_x),
  yend = c(origin_y,        origin_y + len)
)
# Use same palette in UMAP
p_umap <- DimPlot(
  seu_integrated,
  reduction  = "umap",
  split.by   = "condition",
  label      = TRUE,
  # label = FALSE,
  label.size = 3,
  cols       = cell_cols          # <— key line
) +
  NoGrid() + 
  NoLegend() +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text   = element_blank(),
    axis.title  = element_blank(),
  ) +
  geom_segment(
    data = mini_axes,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(4, "pt")),
    linewidth = 0.4
  ) +
  annotate("text",
           x = origin_x + len, y = origin_y,
           label = "UMAP1", vjust = 1.5, size = 3) +
  annotate("text",
           x = origin_x - 0.25, y = origin_y + len,
           label = "UMAP2", hjust = 0.5, angle = 90, size = 3)

# counts per cell type per condition
tab_counts <- table(
  CellType  = Idents(seu_integrated),
  Condition = seu_integrated$condition
)

# relative proportions within each condition (columns sum to 1)
tab_prop <- prop.table(tab_counts, margin = 2)
df_prop <- as.data.frame(tab_prop) %>%
  dplyr::rename(Proportion = Freq)

p_bar <- ggplot(df_prop,
                aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 1) +
  coord_flip() +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  xlab(NULL) +   # now appears above the bars
  ylab(NULL) +
  theme_minimal() +
  labs(title = "Relative Abundance (%)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),  # keep digits hidden
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),   
    axis.ticks.y     = element_blank()
  )

# plot combined graphs
p_combined <- p_umap / p_bar + plot_layout(heights = c(12, 1))

print(p_combined)

ggsave(
  filename = "umap_rel_abundance_legend.png",
  plot     = p_combined,
  width    = 8,    # adjust as needed
  height   = 6,
  dpi      = 300
)

# ---- sample-level Principal Component Analysis ----

avg <- AverageExpression(
  seu_integrated, 
  assays = "RNA",
  layer = "data",
  group.by = "orig.ident"
)$RNA

mat <- t(avg)
pca <- prcomp(mat, scale. = TRUE)

# Get variance explained percentages
pc_var <- round(100 * summary(pca)$importance[2, 1:2], 1) #PC1, PC2 %

df <- data.frame(PC1 = pca$x[,1],
                 PC2 = pca$x[,2],
                 sample = rownames(pca$x))

sample_condition <- tapply(
  seu_integrated$condition,
  seu_integrated$orig.ident,
  unique
)

df$condition <- sample_condition[df$sample]

ggplot(df, aes(PC1, PC2, color = condition, fill = condition)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, alpha = 0.3, geom = "polygon") +
  labs(x = paste0("PC1 (", pc_var[1], "%)"),
       y = paste0("PC2 (", pc_var[2], "%)")) + 
  theme_classic() +
  theme(legend.title = element_blank()) # Removes legend title

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
  
  # Print summaries
  cat("\n=== DESeq2 Results Summary for", cell_type, "===\n")
  cat("Water vs Control:\n"); summary(res_water_vs_ctrl)
  cat("\nBlumeria vs Control:\n"); summary(res_blum_vs_ctrl)
  cat("\nBlumeria vs Water:\n"); summary(res_blum_vs_water)
  
  # Normalize data for PCA
  rld <- rlog(dds, blind = FALSE)
  
  # Generate plots
  title <- ifelse(is.null(plot_title), cell_type, plot_title)
  
  p_pca <- plotPCA(rld, intgroup = "condition") +
    ggtitle(paste("Pseudoulk PCA:", title))
  
  p_ma_water <- plotMA(res_water_vs_ctrl) + ggtitle("MA: Water vs Control")
  p_ma_blum_water <- plotMA(res_blum_vs_water) + ggtitle("MA: Blumeria vs Water")
  
  # Display plots
  print(p_pca)
  print(p_ma_water)
  print(p_ma_blum_water)
  
  # Save results if requested
  if (save_results) {
    results_list <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      dds = dds,
      rld = rld
    )
    
    filename <- paste0("DEG_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = filename)
    cat("Results saved to:", filename, "\n")
    
    # Save plots
    ggsave(paste0("PCA_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
           plot = p_pca, width = 8, height = 6, dpi = 300)
    ggsave(paste0("MA_Water_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
           plot = p_ma_water, width = 8, height = 6, dpi = 300)
    ggsave(paste0("MA_BlumWater_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"),
           plot = p_ma_blum_water, width = 8, height = 6, dpi = 300)
  }
  
  # Return results
  return(list(
    results = list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water
    ),
    dds = dds,
    rld = rld,
    plots = list(pca = p_pca, ma_water = p_ma_water, ma_blum_water = p_ma_blum_water)
  ))
}
   
cm_results <- run_pseudobulk_deg(seu_integrated, "Cardiomyocytes")
  




