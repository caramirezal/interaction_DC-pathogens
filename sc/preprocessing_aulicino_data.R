## Preprocessing Aulicino scRNA-Seq dataset on MoDCs 

## Dependencies
library(GEOquery)
library(Seurat)
library(dplyr)

path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

## Downloading metadata
series_url <- 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111543/matrix/GSE111543_series_matrix.txt.gz'
series_fname <- 'GSE111543_series_matrix.txt.gz'
download.file(series_url, series_fname)
gse <- getGEO(filename = series_fname)
gse@phenoData@data %>%
        head
file.remove(series_fname)
## formatting annotations
anns <- gse@phenoData@data
anns <- select(anns, 
               title, 
               supplementary_file_1,
               geo_accession,
               organism_ch1,
               `infection:ch1`:`well:ch1`)
colnames(anns) <- gsub(':ch1|_ch1', '', colnames(anns))
head(anns)

######################################################################
## Download count matrices data
#data_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111546&format=file'
#dir.create('tmp')
#download.file(data_url, destfile = 'tmp/GSE111546_RAW.tar')
#untar(tarfile = 'tmp/GSE111546_RAW.tar', 
#      exdir = 'tmp/GSE111546_RAW')
#count.files <- list.files('tmp/GSE111546_RAW', 
#                          full.names = TRUE)
#file.remove('tmp/GSE111546_RAW.tar')

## Loading files into R
counts.list <- lapply(count.files,
                      function(x) 
                              read.table(
                                      file = x,
                                      skip = 1,
                                      header = TRUE, 
                                      sep = '\t'
                              )
)
names(counts.list) <- gsub('tmp/GSE111546_RAW/|.txt.gz',
                           '', count.files)
counts.list[1][[1]] %>% head
## Checking that all the tables have the same order
are_ordered <- sapply(
        counts.list, 
        function(x) # The gene ids from the 1rst table is compared to others
         any(counts.list[1][[1]]$Geneid != x$Geneid)
)
sum(are_ordered)

## [1] 0     # Is OK

## Extracting column 7
gene_ids <- counts.list[1][[1]]$Geneid
counts <- sapply(counts.list, function(x) x[,7], simplify = TRUE)
rownames(counts) <- gene_ids
counts[1:5, 1:5]
rm(counts.list)
dim(counts)

dir.create('data/aulicino2018')
saveRDS(counts, 
        file = 'data/aulicino2018/counts_mtx.rds',
        compress = TRUE)


## Creation of the seurat file
counts <- readRDS('data/aulicino2018/counts_mtx.rds')
seurat <- CreateSeuratObject(
        counts = counts,
        project = 'Yersinia',
        assay = 'RNA',
        min.cells = 1,
        min.features = 1
)
seurat

## Annotation of the seurat object
anns <- mutate(anns,
               cell_name= gsub('.*suppl\\/|.txt.gz', '',
                               supplementary_file_1))
rownames(anns) <- anns$cell_name

intersect(anns$cell_name, colnames(seurat))
anns.seu <- anns[colnames(seurat), ]
anns.seu <- cbind(seurat@meta.data, anns.seu)
head(anns.seu)

seurat@meta.data <- anns.seu
saveRDS(seurat, 
        'data/aulicino2018/aulicino2018_seu.rds', 
        compress = TRUE)
