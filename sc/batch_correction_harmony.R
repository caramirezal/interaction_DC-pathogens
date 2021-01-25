library(harmony)
library(Seurat)

path2seurat <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/pass_qc_seurat/pass_qc_seurat.rds'
batch_column <- 'plate_number'
out_dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/batch_correction/' 

## Reading seurat object
seurat <- readRDS(path2seurat)

## Harmony integration
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst")
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 20)
seurat <- RunHarmony(seurat, batch_column)


## Saving results
out_seu_file <- paste0(out_dir, '/harmony_batch_corrected.rds')
saveRDS(seurat, out_seu_file, compress=TRUE)

## Plotting umap
umap_path <- paste0(out_dir, '/harmony_batch_correction.pdf')
pdf(umap_path)
DimPlot(seurat, 
        reduction = "harmony",
        group.by=batch_column)
dev.off()
