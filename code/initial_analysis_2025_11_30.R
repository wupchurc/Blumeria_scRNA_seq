# load libraries
library(Seurat)


#10X CellRanger 
data_dir <- "data/cellranger_Control_SP_11_GEX_FL-Z0041/raw_feature_bc_matrix/"
matrix <- Read10X(data.dir = data_dir)