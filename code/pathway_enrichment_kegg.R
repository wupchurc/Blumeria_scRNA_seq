library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)

cm_results <- readRDS("rds_files/DEG_Cardiomyocytes.rds")

# Function to process one DESeqResults → up/down ENTREZ lists
get_sig_genes <- function(res, lfc_threshold = 0, padj_threshold = 0.1) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  
  up   <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange >  lfc_threshold]
  down <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange < -lfc_threshold]
  
  list(up = up, down = down, all = df$ENTREZID)
}

# Extract sig lists for ALL contrasts
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(cm_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(cm_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(cm_results$blum_vs_water),
  water_vs_blum = get_sig_genes(cm_results$water_vs_blum)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichKEGG(
    gene      = sig_lists$water_vs_ctrl$up,
    universe  = sig_lists$water_vs_ctrl$all,
    organism  = "rno",            # rat
    keyType   = "ncbi-geneid",           # or "ncbi-geneid"/"entrezid" depending on IDs
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  ),
  water_vs_blum = enrichKEGG(
    gene      = sig_lists$water_vs_blum$up,
    universe  = sig_lists$water_vs_blum$all,
    organism  = "rno",
    keyType   = "ncbi-geneid",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
)

ego_down_all <- list(
  water_vs_ctrl = enrichKEGG(
    gene      = sig_lists$water_vs_ctrl$down,
    universe  = sig_lists$water_vs_ctrl$all,
    organism  = "rno",            # rat
    keyType   = "ncbi-geneid",           # or "ncbi-geneid"/"entrezid" depending on IDs
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  ),
  water_vs_blum = enrichKEGG(
    gene      = sig_lists$water_vs_blum$down,
    universe  = sig_lists$water_vs_blum$all,
    organism  = "rno",
    keyType   = "ncbi-geneid",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
)


dotplot(ego_up_all$water_vs_blum, showCategory=20, font.size = 6,title = "Water vs Blum - Upregulated")
dotplot(ego_down_all$water_vs_blum, showCategory=20, font.size = 6,title = "Water vs Blum - ")


# Download rno04216 genes (rat ferroptosis)
library(KEGGREST)
rno04216_genes <- keggGet("rno04216")[[1]]$GENE  # returns "entrez:gene_name"
rno_genes <- sub("rno:([0-9]+).*", "\\1", rno04216_genes)  # extract Entrez IDs

# Test your sig genes against this pathway
ferro_test <- enricher(gene = sig_lists$water_vs_ctrl$up,
                       TERM2GENE = data.frame(ID="rno04216", Gene=rno_genes),
                       pvalueCutoff=0.1)
ferro_test
