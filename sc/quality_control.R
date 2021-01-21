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
path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/DCs-pathogen_interactions_seu.rds'
subset_column <- 'sample_type'
subset <- 'single cell'
batch_column <- 'plate_number'
out_dir <- '/Users/carlosramirez/sc/interaction_DC-pathogens/figures/'
count_th <- 100000
perc_mit_th <- 10
pass_qc_seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/'

## Setting directory with output figures
figures_dir <- paste0(out_dir, '/seurat_qc/')
dir.create(figures_dir)

seurat <- readRDS(path2seurat)

## Setting batch column
Idents(seurat) <- seurat@meta.data[,batch_column]


## Subsetting cells
if ( 0 < length(subset)) {
        cells <- sapply(seurat@meta.data[, subset_column],
                        function(x) x %in% subset)
        cells <- colnames(seurat)[cells]
        seurat <- subset(seurat, cells = cells)
}


## QC - RNA Counts
pdf(paste0(figures_dir, 'nCount_RNA_QC.pdf'))
VlnPlot(seurat, 
        features = 'nCount_RNA', 
        log = TRUE) +
        ylab('Log(Total RNA counts)')
dev.off()

## QC number of features
pdf(paste0(figures_dir, 'nFeature_RNA_QC.pdf'))
VlnPlot(seurat, 
        features = 'nFeature_RNA') +
        ylab('Total RNA counts')
dev.off()

## Mitocondrial genes
seurat[["percent.mt"]] <- PercentageFeatureSet(
        seurat, 
        pattern = "^mt-"
)
pdf(paste0(figures_dir, 'percent_mitochondrial_genes.pdf'))
VlnPlot(seurat, 
        features = "percent.mt")
dev.off()

## Mitocondrial vs number of counts
pdf(paste0(figures_dir, 'perc_mito_vs_ncount.pdf'))
FeatureScatter(seurat, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mt", 
               pt.size = 2) +
        ggtitle('') +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        geom_hline(yintercept = perc_mit_th, 
                   linetype='dashed', 
                   color='red', 
                   size=0.8) +
        geom_vline(xintercept = count_th, 
                   linetype='dashed', 
                   color='red', 
                   size=0.8)
dev.off()

## Normalization and Scaling
seurat <- NormalizeData(seurat) %>%
                        ScaleData()

## Variable genes
seurat <- FindVariableFeatures(
        seurat, 
        selection.method = 'vst', 
        nfeatures = 3000
)
top10 <- head(VariableFeatures(seurat), 10)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste0(figures_dir, 'variable_genes.pdf'))
plot2
dev.off()

## PCA
seurat <- RunPCA(seurat, npcs = 30)
DimPlot(seurat, reduction = 'pca')
pdf(paste0(figures_dir, 'pca.pdf'))
dittoDimPlot(seurat, 
             var = batch_column, 
             reduction.use = "pca", 
             do.ellipse = TRUE, 
             do.label = TRUE, 
             legend.show = FALSE) +
        theme(panel.grid = element_blank())
dev.off()

## UMAP

seurat <- RunUMAP(seurat, dims = 1:20)

pdf(paste0(figures_dir, 'umap.pdf'))
dittoDimPlot(seurat, 
             var = batch_column, 
             reduction.use = "umap", 
             do.ellipse = TRUE, 
             do.label = TRUE, 
             legend.show = FALSE) +
        theme(panel.grid = element_blank())
dev.off()

################################################################
## Filter cells
seurat.filtered <- subset(seurat,
                          nCount_RNA < count_th &
                                  percent.mt < perc_mit_th)

## Saving pass QC cells
pass_qc_dir <- paste0(pass_qc_seurat, '/pass_qc_seurat/')
dir.create(pass_qc_dir)
pass_qc_seurat <- paste0(pass_qc_dir, '/pass_qc_seurat.rds')
saveRDS(seurat.filtered,
        file = pass_qc_seurat, 
        compress = TRUE)
rm(seurat.filtered)

## After QC check
seurat.filtered <- readRDS(pass_qc_seurat)
seurat.filtered

pdf(paste0(pass_qc_dir, 'umap_pass_qc.pdf'))
dittoDimPlot(seurat.filtered, 
             var = batch_column, 
             reduction.use = "umap", 
             do.ellipse = TRUE, 
             do.label = TRUE, 
             legend.show = FALSE) +
        theme(panel.grid = element_blank())
dev.off()

pdf(paste0(pass_qc_dir, 'pca_pass_qc.pdf'))
dittoDimPlot(seurat.filtered, 
             var = batch_column, 
             reduction.use = "pca", 
             do.ellipse = TRUE, 
             do.label = TRUE, 
             legend.show = FALSE) +
        theme(panel.grid = element_blank())
dev.off()
