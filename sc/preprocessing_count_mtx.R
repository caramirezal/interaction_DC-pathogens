## Analysis of DC gene expression response to
## pathogen exposure at the single cell level

## Dependencies
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(biomaRt)
library(org.Mm.eg.db)
library(rio)

## Setting project directory
path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

dendritic <- readRDS('data/C025_Stella.scPipe.SCE.rds')

count.mtx <- assay(dendritic)

## Annotating ensemble gene stable id with gene symbol
mart <- useDataset("mmusculus_gene_ensembl", 
                   useMart("ensembl"))
genes <- rownames(dendritic)
df <- data.frame(ensembl_gene_id_version=rownames(dendritic))
G_list <- getBM(filters= "ensembl_gene_id_version", 
                attributes= c("ensembl_gene_id_version", 
                              "ensembl_gene_id",
                              "mgi_symbol"), 
                values=genes, mart= mart)
head(G_list)
dim(G_list)
unique(G_list$mgi_symbol) %>% length()

## Changing gene ensembl gene id with gene symbols
ann <- G_list[!duplicated(G_list$mgi_symbol), ] ## drop duplications
count.mtx <- count.mtx[ann$ensembl_gene_id_version, ] ## Keeping annotated
any(rownames(count.mtx) != ann$ensembl_gene_id_version) ## Checking order 
rownames(count.mtx) <- ann$mgi_symbol
count.mtx <- count.mtx[rownames(count.mtx) != '', ] ## removing missing value
rownames(count.mtx) %>% head
        
seurat <- CreateSeuratObject(
        counts = count.mtx,
        project = 'Dendritic cells-pathogens interactions',
        assay = 'RNA',
        min.cells = 1,
        min.features = 1
)
md <- as.data.frame(colData(dendritic))
## Checking cells order 
any(rownames(seurat@meta.data) != 
            rownames(md[rownames(seurat@meta.data), ]))
seurat@meta.data <- cbind(seurat@meta.data,
                          md[rownames(seurat@meta.data), ])
head(seurat@meta.data)

####################################################################
## Adding additional metadata

#seurat <- readRDS('data/DCs-pathogen_interactions_seu.rds')
sheet <- import('data/20201214 C025_sample_sheet_SEA.xlsx')

## Tiding up annotations
sheet <- mutate(sheet, id=paste(plate_number, 
                                well_position, sep = '_'))
colnames(sheet) <- gsub(':.*', '', colnames(sheet))

## Checking cell order 
any(sheet$id != colnames(seurat)) ## can be merged! By column

## Checking uniqueness in both annotations
intersect(colnames(sheet), colnames(seurat@meta.data))

## Merging annotations
seurat.md <- dplyr::select(seurat@meta.data, 
                    -disease_state, 
                    -plate_number,
                    -cell_type,
                    -index_sort_condition,
                    -well_position)
anns <- cbind(seurat.md,
              sheet)
head(anns)
colnames(anns)
dim(anns)

## Adding new annotations
seurat@meta.data <- anns

## Saving results
Idents(seurat) <- seurat$plate_number
saveRDS(seurat, 
        file = 'data/DCs-pathogen_interactions_seu.rds',
        compress = TRUE)


