library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)

cm_results <- readRDS("rds_files/DEG_Cardiomyocytes.rds")
fb_results <- readRDS("rds_files/DEG_Fibroblasts.rds")
mac_results <- readRDS("rds_files/DEG_Macrophages.rds")
ntphl_results <- readRDS("rds_files/DEG_Neutrophils.rds") 

t_results <- readRDS("rds_files/DEG_T_cells.rds")
b_results <- readRDS("rds_files/DEG_B_cells.rds")

per_results <- readRDS("rds_files/DEG_Pericytes.rds")
cap_ec_results <- readRDS("rds_files/DEG_Capillary_EC.rds")
ven_ec_results <- readRDS("rds_files/rds_files/DEG_Venous_EC.rds")
lymp_ec_results <- readRDS("rds_files/DEG_Lymphatic_EC.rds")
neuro_results <- readRDS("rds_files/DEG_Neuronal.rds")

#---- Function to process one DESeqResults → up/down ENTREZ lists ----
get_sig_genes <- function(res, lfc_threshold = 0.5, padj_threshold = 0.05) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  
  up   <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange >  lfc_threshold]
  down <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange < -lfc_threshold]
  
  list(up = up, down = down, all = df$ENTREZID)
}
# ---- Print up/down regulated DEGs and background to txt file ----

water_vs_ctrl = get_sig_genes(cm_results$water_vs_ctrl)

up_genes <- water_vs_ctrl$up[!is.na(water_vs_ctrl$up)]
cat(up_genes, sep = "\n")
writeLines(up_genes, "cm_water_vs_ctrl_up_genes.txt")
down_genes <- water_vs_ctrl$down[!is.na(water_vs_ctrl$down)]
cat(down_genes, sep = "\n")
writeLines(down_genes, "cm_water_vs_ctrl_down_genes.txt")
background <- water_vs_ctrl$all[!is.na(water_vs_ctrl$all)]
cat(background, sep = "\n")
writeLines(background, "cm_water_vs_ctrl_background.txt")

water_vs_blum = get_sig_genes(cm_results$water_vs_blum)
up_genes <- water_vs_blum$up[!is.na(water_vs_blum$up)]
cat(up_genes, sep = "\n")
writeLines(up_genes, "cm_water_vs_blum_up_genes.txt")
down_genes <- water_vs_blum$down[!is.na(water_vs_blum$down)]
cat(down_genes, sep = "\n")
writeLines(down_genes, "cm_water_vs_blum_down_genes.txt")
background <- water_vs_blum$all[!is.na(water_vs_blum$all)]
cat(background, sep = "\n")
writeLines(background, "cm_water_vs_blum_background.txt")


water_vs_blum = get_sig_genes(cm_results$water_vs_blum)
down_genes <- water_vs_blum$down[!is.na(water_vs_blum$down)]
cat(down_genes, sep = "\n")
blum_vs_water = get_sig_genes(cm_results$blum_vs_water)
up_genes <- blum_vs_water$up[!is.na(blum_vs_water$up)]
cat(up_genes, sep = "\n")

# ---- Refactored Code ----

analyze_celltype <- function(results_obj, cell_name, lfc_threshold = 0.5, padj_threshold = 0.05) {
  
  sig_lists <- list(
    water_vs_ctrl = get_sig_genes(results_obj$water_vs_ctrl, lfc_threshold, padj_threshold),
    water_vs_blum = get_sig_genes(results_obj$water_vs_blum, lfc_threshold, padj_threshold)
  )
  
  #Shared background
  shared_background <- unique(c(
    sig_lists$water_vs_ctrl$all,
    sig_lists$water_vs_blum$all
  ))
  
  # Compare upregulated
  comp_up <- compareCluster(
    geneClusters = lapply(sig_lists, '[[', "up"),
    fun = "enrichGO",
    # universe = shared_background,
    universe = all_genes,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 2,
    readable = TRUE
  )
  
  # Compare downregulated
  comp_down <- compareCluster(
    geneClusters = lapply(sig_lists, '[[', "down"),
    fun = "enrichGO",
    universe = shared_background,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 2,
    readable = TRUE
  )
  
  # Plots
  p_up <- dotplot(comp_up, showCategory = 10, font.size = 8, size = "Count") + 
    ggtitle(paste(cell_name, "- Upregulated (p<0.25)")) 
  # +
    # scale_size_continuous(name = "Gene\nCount", breaks = c(5,10,15,20,25,30,35,40,45,50),
                          # limits = c(0,50), range = c(2,10))
  
  p_down <- dotplot(comp_down, showCategory = 10, font.size = 8, size = "Count") +
    ggtitle(paste(cell_name, "- Downregulated (p<0.25)")) 
  # +
    # scale_size_continuous(name = "Gene\nCount", breaks = c(5,10,15,20,25,30,40,45,50),
                          # limits = c(0,50), range = c(2,10))
  
  print(p_up)
  print(p_down)
  
  # Return results
  list (up = comp_up, down = comp_down, sig_lists = sig_lists, background = shared_background)
  
}




cm_analysis <- analyze_celltype(cm_results, "Cardiomyocytes")

fb_analysis <- analyze_celltype(fb_results, "Fibroblasts")
mac_analysis <- analyze_celltype(mac_results, "Macrophages")
ntphl_analysis <- analyze_celltype(ntphl_results, "Neutrophils")


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

dotplot(comp_up, showCategory=10, font.size = 8, size = "Count", title = "Cardiomyocytes - Upregulated")

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
  pvalueCutoff = 0.25,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  readable = TRUE
)


p <- dotplot(comp_down, showCategory=10, font.size=8, size = "Count", title = "Cardiomyocytes - Downregulated")

p + scale_size_continuous(
  name = "Gene\nCount",
  breaks = c(5,10,15,20,25),
  limits = c(0, 30),
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


