library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)

path2seurat <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/pass_qc_seurat/pass_qc_seurat.rds'
batch_column <- 'plate_number'
out_dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/batch_correction/' 

## Reading seurat object
seurat <- readRDS(path2seurat)

## Liger integration
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, split.by = batch_column, do.center = FALSE)
seurat <- RunOptimizeALS(seurat, k = 20, lambda = 5, split.by = batch_column)
seurat <- RunQuantileNorm(seurat, split.by = batch_column)
seurat <- FindNeighbors(seurat, reduction = "iNMF", dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.3)

# Dimensional reduction and plotting
seurat <- RunUMAP(seurat, 
                  dims = 1:ncol(seurat[["iNMF"]]), 
                  reduction = "iNMF")

out_seu_file <- paste0(out_dir, '/liger_batch_corrected.rds')

## Saving results
saveRDS(seurat, out_seu_file, compress=TRUE)

## Plotting umap
umap_path <- paste0(out_dir, '/liger_batch_correction.pdf')
pdf(umap_path)
DimPlot(seurat, 
        reduction = "umap",
        group.by=batch_column)
dev.off()
