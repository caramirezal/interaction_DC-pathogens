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

path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/batch_correction/seurat_batch_corrected.rds'
umap1_th <- 1.5

seurat <- readRDS(path2seurat)
seurat

DimPlot(seurat, reduction = 'pca')

## UMAP showing batch corrected data
pdf('figures/umap_filter_threshold.pdf')
DimPlot(seurat) +
        geom_vline(xintercept = umap1_th,
                   size=1.1, 
                   colour='red',
                   linetype='dashed')
dev.off()


## After removing outlier P193 cluster
cells <- colnames(seurat)[
  seurat@reductions$umap@cell.embeddings[,'UMAP_1'] < umap1_th
]
DimPlot(subset(seurat, cells=cells))

## Subsetting to cells without P193 cluster
#if ( 0 < length(cells)){
#        seurat.filtered <-  subset(seurat, cells=cells) 
#} else {
#        seurat.filtered <- seurat
#}
#seurat
#seurat.filtered
#rm(seurat)


#####################################################################
## Differential Expression Analysis

## Comparing MDP vs CMP

## Auxiliary function to perform DEA
dea <- function(seurat,
                ident,
                group_by){
  deg <- FindMarkers(seurat, 
                     ident.1 = ident, 
                     group.by = group_by, 
                     logfc.threshold = 0, 
                     min.pct = 0, 
                     min.cells.feature = 1)
  #deg <- filter(deg, p_val_adj < 0.05) 
  deg <- arrange(deg, desc(abs(avg_logFC)))
  
  return(deg)
}

## Comparing MDP vs CDP
progenitor.deg <- dea(seurat = seurat,
                      ident =  'MDP',
                      group_by = 'progenitor_type')
dim(progenitor.deg)
head(progenitor.deg)

## Settings
n_top_deg <- 50

progenitor.mtx <- seurat@assays$integrated@scale.data[
        rownames(progenitor.deg)[1:n_top_deg],
]
dim(progenitor.mtx)

## Heatmap visualization
col_df <- select(seurat@meta.data, 
                 progenitor_type,
                 plate_number,
                 disease_state, 
                 bacteria_attached,
                 treatment)
head(col_df)
col_anns <- list(progenitor_type=c(CDP='steelblue',
                                   MDP='salmon'))
pheatmap(scale(progenitor.mtx), 
         show_colnames = FALSE, 
         annotation_col = col_df,
         annotation_colors = col_anns)

################################################################
## Comparing MDP (1 and 24 hours)

## Rename cells
## Create a column named infected_vs_bystander
## MOI 0 -> control
## MOI != 0 -> bystander
## bacteria attached (+) -> infected 
seurat$treatment[is.na(seurat$treatment)] <- 'NA'
seurat$'infection_status' <- 'NA'
seurat$infection_status <- sapply(
  seurat$treatment, function(x)
    ifelse(x=='MOI 0', 'Control', 'Bystander')
)
seurat$infection_status <- sapply(
  1:ncol(seurat), function(i)
    ifelse(seurat$bacteria_attached[i]=='+', 
           'Infected', seurat$infection_status[i])
)
table(seurat$treatment,
      seurat$infection_status)

################################################################
##                                                            ##
##            Differential Expression Analysis                ##
##                                                            ##
################################################################

## Directory to store results
dir.create('data/dea/')

###############################################################
## MDPs

## 1hr
##  * A vs B pos
## Infected vs control
degs.infVsCtrl.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.1h, n = 10)
dim(degs.infVsCtrl.1h)
write.table(degs.infVsCtrl.1h,
      file = gzfile('data/dea/degs.mdp.infVsCtrl.1h.tsv.gz'),
      sep = '\t', 
      quote = FALSE, 
      row.names = TRUE
)

##  * A vs B neg
## Bystander vs control
degs.bysVsCtrl.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Control', 'Bystander')),
  ident = 'Bystander',
  group_by = 'infection_status'
)
head(degs.bysVsCtrl.1h, n = 10)
dim(degs.bysVsCtrl.1h)
write.table(degs.bysVsCtrl.1h,
            file = gzfile('data/dea/degs.mdp.bysVsCtrl.1h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * B pos vs B neg
## Infected vs Bystander
degs.infVsBys.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.1h, n = 10)
dim(degs.infVsCtrl.1h)
write.table(degs.infVsCtrl.1h,
            file = gzfile('data/dea/degs.mdp.infVsBys.1h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)


## 24hr
##  * A vs B pos
## Infected vs control
degs.infVsCtrl.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Infected', 'Control')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.24h, n = 10)
dim(degs.infVsCtrl.24h)
write.table(degs.infVsCtrl.24h,
            file = gzfile('data/dea/degs.mdp.infVsCtrl.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * A vs B neg
## Bystander vs control
degs.bysVsCtrl.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Control', 'Bystander')),
  ident = 'Bystander',
  group_by = 'infection_status'
)
head(degs.bysVsCtrl.24h, n = 10)
dim(degs.bysVsCtrl.24h)
write.table(degs.bysVsCtrl.24h,
            file = gzfile('data/dea/degs.mdp.bysVsCtrl.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * B pos vs B neg
## Infected vs Bystander
degs.infVsBys.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'MDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.24h, n = 10)
dim(degs.infVsCtrl.24h)
write.table(degs.infVsCtrl.24h,
            file = gzfile('data/dea/degs.mdp.infVsBys.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

###############################################################
###############################################################
## CDPs

## 1hr
##  * A vs B pos
## Infected vs control
degs.infVsCtrl.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.1h, n = 10)
dim(degs.infVsCtrl.1h)
write.table(degs.infVsCtrl.1h,
            file = gzfile('data/dea/degs.cdp.infVsCtrl.1h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * A vs B neg
## Bystander vs control
degs.bysVsCtrl.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Control', 'Bystander')),
  ident = 'Bystander',
  group_by = 'infection_status'
)
head(degs.bysVsCtrl.1h, n = 10)
dim(degs.bysVsCtrl.1h)
write.table(degs.bysVsCtrl.1h,
            file = gzfile('data/dea/degs.cdp.bysVsCtrl.1h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * B pos vs B neg
## Infected vs Bystander
degs.infVsBys.1h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 1h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.1h, n = 10)
dim(degs.infVsCtrl.1h)
write.table(degs.infVsCtrl.1h,
            file = gzfile('data/dea/degs.cdp.infVsBys.1h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)


## 24hr
##  * A vs B pos
## Infected vs control
degs.infVsCtrl.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Infected', 'Control')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.24h, n = 10)
dim(degs.infVsCtrl.24h)
write.table(degs.infVsCtrl.24h,
            file = gzfile('data/dea/degs.cdp.infVsCtrl.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * A vs B neg
## Bystander vs control
degs.bysVsCtrl.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Control', 'Bystander')),
  ident = 'Bystander',
  group_by = 'infection_status'
)
head(degs.bysVsCtrl.24h, n = 10)
dim(degs.bysVsCtrl.24h)
write.table(degs.bysVsCtrl.24h,
            file = gzfile('data/dea/degs.cdp.bysVsCtrl.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

##  * B pos vs B neg
## Infected vs Bystander
degs.infVsBys.24h <- dea(
  seurat = subset(seurat, 
                  disease_state == 'Yesinia 24h post-infection' &
                    progenitor_type == 'CDP' & 
                    infection_status %in% c('Infected', 'Bystander')),
  ident = 'Infected',
  group_by = 'infection_status'
)
head(degs.infVsCtrl.24h, n = 10)
dim(degs.infVsCtrl.24h)
write.table(degs.infVsCtrl.24h,
            file = gzfile('data/dea/degs.cdp.infVsBys.24h.tsv.gz'),
            sep = '\t', 
            quote = FALSE, 
            row.names = TRUE
)

ch <- read.table('data/dea/degs.cdp.bysVsCtrl.1h.tsv.gz',
                 header = TRUE)
head(ch)
