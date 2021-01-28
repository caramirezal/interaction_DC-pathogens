## Visualization of DEA results

## Dependencies
library(ComplexHeatmap)
library(dplyr)

set.seed(333)

## Setting project directory
path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/batch_correction/seurat_batch_corrected.rds'

degs.files <- list.files('data/dea/', full.names = TRUE)
degs.list <-  lapply(
        degs.files, 
        function(x) read.table(x, header = TRUE)
)
files.names <- gsub('data/dea//degs.|.tsv.gz', '', degs.files)
names(degs.list) <- files.names
lapply(degs.list, head)


## Loading seurat file
seurat <- readRDS(path2seurat)
seurat

## Settings
n_top <- 5

top_genes_by_condition <- lapply(degs.list, 
                                 function(x) 
                                         rownames(x)[1:n_top]) %>%
                                unlist() %>%
                                unique()
top_genes_by_condition


## Extracting matrix of gene expression for the top DEG markers
mtx <- FetchData(seurat,
                 vars = c(colnames(seurat@meta.data),
                          top_genes_by_condition))
mtx.ge <- t(as.matrix(mtx[, top_genes_by_condition]))
mtx.ann <- mtx[, colnames(seurat@meta.data)]
dim(mtx.ge); dim(mtx.ann)

Heatmap(scale(mtx.ge), 
        show_column_names = FALSE)
