## Differential expression analysis of infected or uninfected-Yersinia
## Dendritic cells 

## Dependencies
library(Seurat)
library(dplyr)
library(dittoSeq)
library(pheatmap)
library(reshape2)

set.seed(333)

## Setting project directory
path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/batch_correction/batch_corrected_seu.rds'
umap1_th <- 1.5

seurat <- readRDS(path2seurat)
seurat

DimPlot(seurat, reduction = 'pca')

## UMAP showing batch corrected data
DimPlot(seurat) +
        geom_vline(xintercept = umap1_th,
                   size=1.1, 
                   colour='red',
                   linetype='dashed')
cells <- colnames(seurat)[
        seurat@reductions$umap@cell.embeddings[,'UMAP_1'] < umap1_th
]
## After removing outlier P193 cluster
DimPlot(subset(seurat, cells=cells))

## Subsetting to cells without P193 cluster
if ( 0 < length(cells)){
        seurat.filtered <-  subset(seurat, cells=cells) 
} else {
        seurat.filtered <- seurat
}
seurat
seurat.filtered
rm(seurat)

## Visualization of categories
DimPlot(seurat.filtered, pt.size = 2) +
                theme(legend.position = 'top')
DimPlot(seurat.filtered, group.by = 'disease_state') +
        theme(legend.position = 'top')
DimPlot(seurat.filtered, 
        group.by = 'bacteria_attached') +
        theme(legend.position = 'top')
DimPlot(seurat.filtered, 
        group.by = 'progenitor_type') +
        theme(legend.position = 'top')
DimPlot(seurat.filtered, 
        group.by = 'treatment') +
        theme(legend.position = 'top')

#####################################################################
## Differential Expression Analysis

## Comparing MDP vs CMP

## Settings
n_top_deg <- 40

## Comparing MDP vs CDP
progenitor.deg <- FindMarkers(seurat.filtered, 
                              ident.1 = 'MDP', 
                              group.by = 'progenitor_type', 
                              logfc.threshold = 0, 
                              min.pct = 0)
progenitor.deg <- filter(progenitor.deg, p_val_adj < 0.05) 
progenitor.deg <- arrange(progenitor.deg, desc(abs(avg_logFC)))
dim(progenitor.deg)

progenitor.mtx <- seurat.filtered@assays$integrated@scale.data[
        rownames(progenitor.deg)[1:n_top_deg],
]
dim(progenitor.mtx)
seurat.filtered

## Heatmap visualization
col_df <- select(seurat.filtered@meta.data, progenitor_type)
head(col_df)
col_anns <- list(progenitor_type=c(CDP='steelblue',
                                   MDP='salmon'))
pheatmap(scale(progenitor.mtx), 
         show_colnames = FALSE, 
         annotation_col = col_df,
         annotation_colors = col_anns)

################################################################
## Comparing MDP (1 and 24 hours)
##  * A vs B pos
##  * A vs B neg
##  * B pos vs B neg

seurat.filtered$treatment[is.na(seurat.filtered$treatment)] <- 'NA'
dittoBarPlot(seurat.filtered, 
             "bacteria_attached", 
             group.by = "treatment", 
             scale = 'count')
hist()

## Re-label infection cell status 
seurat.filtered$'treatment_bacteria' <- paste0(
        seurat.filtered$treatment, '_',
        seurat.filtered$bacteria_attached
)
table(seurat.filtered$treatment_bacteria)
table(seurat.filtered$bacteria,
      seurat.filtered$bacteria_attached)
seurat.filtered$treatment_bacteria <-
        gsub(' ', '_', seurat.filtered$treatment_bacteria)
seurat.filtered$treatment_bacteria %>% table()
seurat.filtered$'infected_vs_bystander' <- plyr::mapvalues(
        seurat.filtered$treatment_bacteria,
        from = unique(seurat.filtered$treatment_bacteria),
        to = c('Non-treated', 'Infected', '')
)
hist(seurat.filtered$Y582_15_A_RFP, breaks = 1000)
unique(seurat.filtered$treatment_bacteria)

time <- 'Yesinia 1h post-infection'
progenitor <- 'MDP'




## Helper function to perform DEA
deg <- function(seurat, 
                reference='',
                deg_column=''){
        seurat <- FindMarkers(seurat, 
                              ident.1 = '+', 
                              group.by = 'bacteria_attached', 
                              logfc.threshold = 0, 
                              min.pct = 0, 
                              min.cells.feature = 0, 
                              min.cells.group = 1)
        seurat <- filter(seurat, p_val < 0.01) 
        seurat <- arrange(seurat, desc(abs(avg_logFC)))
        
        col_df <- select(mdp@meta.data, 
                         bacteria_attached, 
                         disease_state) 
        col_df$bacteria_attached <- ifelse(col_df$bacteria_attached == '+',
                                           'Ye', 'No_Ye')
        col_df$disease_state <- ifelse(col_df$disease_state == 'Yesinia 1h post-infection',
                                       'treatment_1h', 'treatment_24h')
        col_anns <- list(bacteria_attached=c(Ye='steelblue',
                                             No_Ye='salmon'),
                         disease_state=c(treatment_1h='green',
                                         treatment_24h='orange'))
        pheatmap(scale(mdp.inf_byst.mtx), 
                 show_colnames = FALSE, 
                 annotation_col = col_df,
                 annotation_colors = col_anns)
}

