## Batch correction in DC single cells infected and non-infected
## with Yersinia

## Dependencies
library(Seurat)
library(dplyr)
library(dittoSeq)

set.seed(333)

## Setting project directory
path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

## Settings
path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/pass_qc_seurat/pass_qc_seurat.rds'
batch_column <- 'plate_number'
k_anchor <- 100
out_dir <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/'
features <- c('nCount_SCT',
              'nFeature_SCT',
              'mapped_to_intron', 
              'percent.mt')

## Reading seurat file
seurat <- readRDS(path2seurat)



## Printing batches
batches <- unique(seurat@meta.data[, batch_column])
cat('Batches: ', batches)

seurat <- SCTransform(seurat, 
                      vars.to.regress = c(batch_column))
seurat <- RunPCA(seurat, npcs = 30) %>% RunUMAP(dims=1:30)
DimPlot(seurat)

## Splitting seurat into batches for correction
batch_cells <- lapply(batches, 
                      function(x) 
                      seurat@meta.data[, batch_column] == x)
names(batch_cells) <- batches
batch_cells <- lapply(batch_cells, 
                      function(x) colnames(seurat)[x])
batches.seu.list <- lapply(batch_cells, 
                         function(x) 
                          subset(seurat, cells = x)) 
rm(seurat)

###################################################################
## Batch correction using CCA
#batches.seu.list <- lapply(batches.seu.list,
#                           function(x)
#                           SCTransform(x)
batches.seu.list

## Batch integration
batches.features <- SelectIntegrationFeatures(
        object.list = batches.seu.list, 
        nfeatures = 3000
)
batches.seu.list <- PrepSCTIntegration(
        object.list = batches.seu.list, 
        anchor.features = batches.features
)
batches.anchors <- FindIntegrationAnchors(
        object.list = batches.seu.list, 
        normalization.method = "SCT", 
        anchor.features = batches.features, 
        k.anchor = k_anchor, 
        k.filter = 5
)
batches.integrated <- IntegrateData(
        anchorset = batches.anchors, 
        normalization.method = "SCT"
)

batches.integrated <- RunPCA(batches.integrated, 
                             npcs = 30) %>%
                        RunUMAP(dims=1:30)

## Assessing batch correction
batch_corrected_dir <- paste0(out_dir, '/batch_correction/')
dir.create(batch_corrected_dir)
pdf(paste0(batch_corrected_dir, 'batch_correction.pdf'))
DimPlot(batches.integrated)
FeaturePlot(batches.integrated, 
            features = features,
            ncol = 2)
VlnPlot(batches.integrated, 
        features = features,
        ncol = 2)
dev.off()

## Visualization of categories
pdf('figures/qc_metrics_after_batch_correction.pdf',
    width = 12)
p1 <- DimPlot(batches.integrated, pt.size = 2) +
        theme(legend.position = 'top')
p2 <- DimPlot(batches.integrated, group.by = 'disease_state') +
        theme(legend.position = 'top')
p3 <- DimPlot(batches.integrated, 
              group.by = 'bacteria_attached') +
        theme(legend.position = 'top')
p4 <- DimPlot(batches.integrated, 
              group.by = 'progenitor_type') +
        theme(legend.position = 'top')
p5 <- DimPlot(batches.integrated, 
              group.by = 'treatment') +
        theme(legend.position = 'top')
p1 + p2 + p3 + p4 + p5
dev.off()

## Saving corrected data
saveRDS(batches.integrated,
        file = paste0(out_dir, 
                      'seurat_batch_corrected.rds'),
        compress = TRUE)


####################################################################
## Data integration with LIGER ** Run in cluster
## -- Use batch_correction conda environment

#library(liger)
#library(Seurat)
#library(SeuratData)
#library(SeuratWrappers)

#path2seurat <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/pass_qc_seurat/pass_qc_seurat.rds'
#batch_column <- 'plate_number'
#out_dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/batch_correction/' 

## Reading seurat object
#seurat <- readRDS(path2seurat)

## Liger integration
#seurat <- NormalizeData(seurat)
#seurat <- FindVariableFeatures(seurat)
#seurat <- ScaleData(seurat, split.by = batch_column, do.center = FALSE)
#seurat <- RunOptimizeALS(seurat, k = 20, lambda = 5, split.by = batch_column)
#seurat <- RunQuantileNorm(seurat, split.by = batch_column)
#seurat <- FindNeighbors(seurat, reduction = "iNMF", dims = 1:20)
#seurat <- FindClusters(seurat, resolution = 0.3)

# Dimensional reduction and plotting
#seurat <- RunUMAP(seurat, 
#                  dims = 1:ncol(seurat[["iNMF"]]), 
#                  reduction = "iNMF")

#out_seu_file <- paste0(out_dir, '/liger_batch_corrected.rds')

## Saving results
#saveRDS(seurat, out_seu_file, compress=TRUE)

## Loading results
seurat <- readRDS('/Users/carlosramirez/sc/interaction_DC-pathogens/data/batch_correction_liger/liger_batch_corrected.rds')

DimPlot(seurat, group.by = 'plate_number')

######################################################################
## Batch correction with Harmony. Run this part in the cluster
## -- Use batch_correction conda environment

#library(harmony)
#library(Seurat)

#path2seurat <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/pass_qc_seurat/pass_qc_seurat.rds'
#batch_column <- 'plate_number'
#out_dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/batch_correction/' 

## Reading seurat object
#seurat <- readRDS(path2seurat)

## Harmony integration
#seurat <- NormalizeData(seurat)
#seurat <- FindVariableFeatures(seurat, selection.method = "vst")
#seurat <- ScaleData(seurat)
#seurat <- RunPCA(seurat, npcs = 20)
#seurat <- RunHarmony(seurat, batch_column)


## Saving results
#out_seu_file <- paste0(out_dir, '/harmony_batch_corrected.rds')
#saveRDS(seurat, out_seu_file, compress=TRUE)

## Plotting umap
#umap_path <- paste0(out_dir, '/harmony_batch_correction.pdf')
#pdf(umap_path)
#DimPlot(seurat, 
#        reduction = "harmony",
#        group.by=batch_column)
#dev.off()

## Harmony results
seurat <- readRDS('data/batch_correction/harmony_batch_corrected.rds')

seurat <- RunUMAP(seurat, dims = 1:20, reduction = 'harmony')
DimPlot(seurat, reduction = 'harmony')
