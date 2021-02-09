## Preprocessing Aulicino scRNA-Seq dataset on MoDCs 

## Dependencies
library(GEOquery)


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
data_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111546&format=file'
dir.create('tmp')
download.file(data_url, destfile = 'tmp/GSE111546_RAW.tar')
untar(tarfile = 'tmp/GSE111546_RAW.tar', 
      exdir = 'tmp/GSE111546_RAW')
count.files <- list.files('tmp/GSE111546_RAW', 
                          full.names = TRUE)
file.remove('tmp/GSE111546_RAW.tar')

## Loading 
readLines(count.files[1], n = 10)
tab <- read.table(
        file = count.files[1],
        skip = 1,
        header = TRUE, 
        sep = '\t'
)
tab <- tab[, c(1,7)]
names(tab) <- c('gene_id', 'count')
head(tab)




