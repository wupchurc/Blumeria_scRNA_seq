library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)

cm_results <- readRDS("rds_files/DEG_Cardiomyocytes.rds")
fb_results <- readRDS("rds_files/DEG_Fibroblasts.rds")
mac_results <- readRDS("rds_files/DEG_Macrophages.rds")
t_results <- readRDS("DEG_T_cells.rds")
b_results <- readRDS("DEG_B_cells.rds")
ntphl_results <- readRDS("DEG_Neutrophils.rds") 
per_results <- readRDS("DEG_Pericytes.rds")
cap_ec_results <- readRDS("DEG_Capillary_EC.rds")
ven_ec_results <- readRDS("DEG_Venous_EC.rds")
lymp_ec_results <- readRDS("DEG_Lymphatic_EC.rds")
neuro_results <- readRDS("rds_files/DEG_Neuronal.rds")

# Function to process one DESeqResults → up/down ENTREZ lists
get_sig_genes <- function(res, lfc_threshold = 0, padj_threshold = 0.05) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  
  up   <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange >  lfc_threshold]
  down <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange < -lfc_threshold]
  
  list(up = up, down = down, all = df$ENTREZID)
}

# ---- Analyze Cardiomyocytes ----
# Extract sig lists for ALL contrasts
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(cm_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(cm_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(cm_results$blum_vs_water),
  water_vs_blum = get_sig_genes(cm_results$water_vs_blum)
)

shared_background <- unique(c(
  sig_lists$water_vs_ctrl$all,
  sig_lists$water_vs_blum$all
))

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Water_vs_Blum = sig_lists$water_vs_blum$up
  ),
  fun = "enrichGO",
  universe = shared_background,
  OrgDb = org.Rn.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  readable = TRUE
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Water_vs_Blum = sig_lists$water_vs_blum$down
  ),
  fun = "enrichGO",
  universe = shared_background,
  OrgDb = org.Rn.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  readable = TRUE
)

dotplot(comp_up, showCategory=10, font.size = 8, size = "Count", title = "Cardiomyocytes - Upregulated")
p <- dotplot(comp_down, showCategory=10, font.size=8, size = "Count", title = "Cardiomyocytes - Downregulated")

p + scale_size_continuous(
  name = "Gene\nCount",
  breaks = c(5,10,15,20,25),
  range = c(2,10)
)

# Individual contrast
ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.2, qvalueCutoff=0.2),
  water_vs_blum = enrichGO(gene=sig_lists$water_vs_blum$down, 
                           universe=sig_lists$water_vs_blum$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=.2, qvalueCutoff=0.5)
)

dotplot(ego_down_all$water_vs_blum, showCategory=15, title = "Water vs Blum - Downregulated")

# ---- Analyze Fibroblasts ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(fb_results$water_vs_ctrl, lfc_threshold = 0),
  blum_vs_ctrl  = get_sig_genes(fb_results$blum_vs_ctrl, lfc_threshold = 0),
  blum_vs_water = get_sig_genes(fb_results$blum_vs_water, lfc_threshold = 0), 
  water_vs_blum = get_sig_genes(fb_results$water_vs_blum, lfc_threshold = 0)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.1, qvalueCutoff=0.1),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  water_vs_blum = enrichGO(gene=sig_lists$water_vs_blum$up, 
                           universe=sig_lists$water_vs_blum$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.1, qvalueCutoff=0.1)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  water_vs_blum = enrichGO(gene=sig_lists$water_vs_blum$down, 
                           universe=sig_lists$water_vs_blum$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.1, qvalueCutoff=0.1)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    # Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    # Blum_vs_Water = sig_lists$blum_vs_water$up,
    Water_vs_Blum = sig_lists$water_vs_blum$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.1
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    # Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    # Blum_vs_Water = sig_lists$blum_vs_water$down,
    Water_vs_Blum = sig_lists$water_vs_blum$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.1
)

dotplot(comp_up, showCategory=10,font.size = 8, size = "Count", title = "Fibroblasts - Upregulated")
dotplot(comp_down, showCategory=10, font.size=6, size = "Count", title = "Fibroblasts - Downregulated")

cnetplot(comp_up, showCategory=10)
# ---- Analyze Macrophages ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(mac_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(mac_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(mac_results$blum_vs_water),
  water_vs_blum = get_sig_genes(mac_results$water_vs_blum)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

# Compare UPREGULATED across all 3 contrasts
comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    # Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    # Blum_vs_Water = sig_lists$blum_vs_water$up,
    Water_vs_Blum = sig_lists$water_vs_blum$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.1
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    # Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    # Blum_vs_Water = sig_lists$blum_vs_water$down,
    Water_vs_Blum = sig_lists$water_vs_blum$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.1
)

dotplot(comp_up, showCategory=10, font.size=8, size = "Count", title = "Macrophages - Upregulated")
dotplot(comp_down, showCategory=10, font.size=6, size = "Count", title = "Macrophages - Downregulated")

# ---- Analyze T cells ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(t_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(t_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(t_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "T cells - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "T cells - Downregulated")

# ---- Analyze B cells ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(b_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(b_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(b_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "B cells - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "B cells - Downregulated")

# ---- Analyze Neutrophils ----
# Extract sig lists for ALL contrasts
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(ntphl_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(ntphl_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(ntphl_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

# Compare UPREGULATED across all 3 contrasts
comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size=8, size = "Count", title = "Neutrophils - Upregulated")

# ---- Analyze Pericytes ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(per_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(per_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(per_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "Pericytes - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "Pericytes - Downregulated")

# ---- Analyze Capillary EC ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(cap_ec_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(cap_ec_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(cap_ec_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "Capillary EC - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "Capillary EC - Downregulated")


# ---- Analyze Venous EC ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(ven_ec_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(ven_ec_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(ven_ec_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "Venous EC - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "Venous EC - Downregulated")
# ---- Analyze Lymphatic EC ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(ven_ec_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(ven_ec_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(ven_ec_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "Venous EC - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "Venous EC - Downregulated")
# ---- Analyze Neuronal ----
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(neuro_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(neuro_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(neuro_ec_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

dotplot(comp_up, showCategory=5, font.size = 8, size = "Count", title = "Neuronal - Upregulated")
dotplot(comp_down, showCategory=5, font.size=8, size = "Count", title = "Neuronal - Downregulated")
