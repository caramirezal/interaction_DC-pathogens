## DotPlots of DEGs

## Dependencies
library(Seurat)
library(dplyr)
library(viridis)

set.seed(333)

## Setting project directory
path2project <- '/Users/carlosramirez/sc/interaction_DC-pathogens/'
setwd(path2project)

path2seurat <- '/Users/carlosramirez/sc/interaction_DC-pathogens/data/batch_correction/seurat_batch_corrected.rds'

degs.files <- list.files('data/dea/', full.names = TRUE)
#degs.files <- grep('cdp', degs.files, value = TRUE)
#degs.files <- grep('infVsBys', degs.files, value = TRUE)
degs.list <-  lapply(
        degs.files, 
        function(x) read.table(x, header = TRUE)
)
files.names <- gsub('data/dea//degs.|.tsv.gz', '', degs.files)
names(degs.list) <- files.names
lapply(degs.list, head)


## Loading seurat file
seurat <- readRDS(path2seurat)
seurat$disease_state <- plyr::mapvalues(
        seurat$disease_state,
        from = unique(seurat$disease_state),
        to = c('Ye 1h', 'Ye 24h')
)

##########################################################
## DotPlots

## Auxiliary function
plot_dots <- function(seurat,
                      degs,
                      n_top=15,
                      plot_name=''){
        degs.ord <- degs[!grepl('Gm[0-9]|Rik', 
                                rownames(degs)), ]
        degs.ord <- degs.ord %>%
                arrange(desc(avg_logFC)) 
        top_pos <- degs.ord %>% head(n_top) %>% rownames() 
        top_neg <- degs.ord %>% tail(n_top) %>% rownames() 
        top <- c(top_neg, top_pos)
        DotPlot(
                object = seurat,
                group.by = 'infection_status',
                features = top, 
                dot.scale = 10
        ) + 
                #coord_flip() +
                scale_color_viridis() +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      axis.text.x = element_text(angle = 45,
                                                 hjust = 1)) +
                ggtitle(plot_name) +
                xlab('') + ylab('')
}

##########################################################
## Plotting

## MDP
progenitor_type <- 'MDP'
## index
n <- 6

## Bystander vs Control 1h MDP
n <- n + 1
names(degs.list)[n]
#[1] "mdp.bysVsCtrl.1h"
mdp.bysVsCtrl.1h <- plot_dots(subset(seurat,
                 progenitor_type == progenitor_type &
                    infection_status %in% c('Bystander',
                                            'Control') &
                    disease_state == 'Ye 1h'), 
          degs = degs.list[n][[1]], 
          plot_name = names(degs.list)[n]
)
mdp.bysVsCtrl.1h


## Bystander vs Control 24h MDP
n <- n + 1
names(degs.list)[n]
#[1] "mdp.bysVsCtrl.24h"
mdp.bysVsCtrl.24h <- plot_dots(subset(seurat,
                                    progenitor_type == progenitor_type &
                                            infection_status %in% c('Bystander',
                                                                    'Control') &
                                            disease_state == 'Ye 24h'), 
                             degs = degs.list[n][[1]], 
                             plot_name = names(degs.list)[n]
)
mdp.bysVsCtrl.24h

## Infected vs Bystander 1h MDP
n <- n + 1
names(degs.list)[n]
#[1] "mdp.infVsBys.1h"
mdp.infVsBys.1h <- plot_dots(subset(seurat,
                                      progenitor_type == progenitor_type &
                                              infection_status %in% c('Infected',
                                                                      'Bystander') &
                                              disease_state == 'Ye 1h'), 
                               degs = degs.list[n][[1]], 
                               plot_name = names(degs.list)[n]
)
mdp.infVsBys.1h

n <- n + 1
## Infected vs Bystander 24h MDP
names(degs.list)[n]
#[1] "mdp.infVsBys.24h"
mdp.infVsBys.24h <- plot_dots(subset(seurat,
                                    progenitor_type == progenitor_type &
                                            infection_status %in% c('Infected',
                                                                    'Bystander') &
                                            disease_state == 'Ye 24h'), 
                             degs = degs.list[n][[1]], 
                             plot_name = names(degs.list)[n]
)
mdp.infVsBys.24h


## Infected vs Ctrl 1h MDP
n <- n + 1
names(degs.list)[n]
#[1] "mdp.infVsCtrl.1h"
mdp.infVsCtrl.1h <- plot_dots(subset(seurat,
                                     progenitor_type == progenitor_type &
                                             infection_status %in% c('Infected',
                                                                     'Control') &
                                             disease_state == 'Ye 1h'), 
                              degs = degs.list[n][[1]], 
                              plot_name = names(degs.list)[n]
)
mdp.infVsCtrl.1h


## Infected vs Ctrl 1h CDP
n <- n + 1
names(degs.list)[n]
#[1] "mdp.infVsCtrl.24h"
mdp.infVsCtrl.24h <- plot_dots(subset(seurat,
                                     progenitor_type == progenitor_type &
                                             infection_status %in% c('Infected',
                                                                     'Control') &
                                             disease_state == 'Ye 24h'), 
                              degs = degs.list[n][[1]], 
                              plot_name = names(degs.list)[n]
)
mdp.infVsCtrl.24h

dot_plots <- list(
        mdp.bysVsCtrl.1h,
        mdp.bysVsCtrl.24h,
        mdp.infVsBys.1h,
        mdp.infVsBys.24h,
        mdp.infVsCtrl.1h,
        mdp.infVsCtrl.24h
)

pdf('figures/dot_plots_degs_mdp.pdf', 
    width = 20)
gridExtra::grid.arrange(
        grobs=dot_plots, 
        ncol=2
)
dev.off()

###############
## MDPs

## CDP

## index
n <- 0

## Bystander vs Control 1h CDP
n <- n + 1
names(degs.list)[n]
#[1] "cdp.bysVsCtrl.1h"
cdp.bysVsCtrl.1h <- plot_dots(subset(seurat,
                                     progenitor_type == 'CDP' &
                                             infection_status %in% c('Bystander',
                                                                     'Control') &
                                             disease_state == 'Ye 1h'), 
                              degs = degs.list[n][[1]], 
                              plot_name = names(degs.list)[n]
)
cdp.bysVsCtrl.1h


## Bystander vs Control 24h CDP
n <- n + 1
names(degs.list)[n]
#[1] "cdp.bysVsCtrl.24h"
cdp.bysVsCtrl.24h <- plot_dots(subset(seurat,
                                      progenitor_type == 'CDP' &
                                              infection_status %in% c('Bystander',
                                                                      'Control') &
                                              disease_state == 'Ye 24h'), 
                               degs = degs.list[n][[1]], 
                               plot_name = names(degs.list)[n]
)
cdp.bysVsCtrl.24h

## Infected vs Bystander 1h CDP
n <- n + 1
names(degs.list)[n]
#[1] "cdp.infVsBys.1h"
cdp.infVsBys.1h <- plot_dots(subset(seurat,
                                    progenitor_type == 'CDP' &
                                            infection_status %in% c('Infected',
                                                                    'Bystander') &
                                            disease_state == 'Ye 1h'), 
                             degs = degs.list[n][[1]], 
                             plot_name = names(degs.list)[n]
)
cdp.infVsBys.1h

n <- n + 1
## Infected vs Bystander 24h CDP
names(degs.list)[n]
#[1] "cdp.infVsBys.24h"
cdp.infVsBys.24h <- plot_dots(subset(seurat,
                                     progenitor_type == 'CDP' &
                                             infection_status %in% c('Infected',
                                                                     'Bystander') &
                                             disease_state == 'Ye 24h'), 
                              degs = degs.list[n][[1]], 
                              plot_name = names(degs.list)[n]
)
cdp.infVsBys.24h


## Infected vs Ctrl 1h CDP
n <- n + 1
names(degs.list)[n]
#[1] "cdp.infVsCtrl.1h"
cdp.infVsCtrl.1h <- plot_dots(subset(seurat,
                                     progenitor_type == 'CDP' &
                                             infection_status %in% c('Infected',
                                                                     'Control') &
                                             disease_state == 'Ye 1h'), 
                              degs = degs.list[n][[1]], 
                              plot_name = names(degs.list)[n]
)
cdp.infVsCtrl.1h


## Infected vs Ctrl 1h CDP
n <- n + 1
names(degs.list)[n]
#[1] "cdp.infVsCtrl.24h"
cdp.infVsCtrl.24h <- plot_dots(subset(seurat,
                                      progenitor_type == 'CDP' &
                                              infection_status %in% c('Infected',
                                                                      'Control') &
                                              disease_state == 'Ye 24h'), 
                               degs = degs.list[n][[1]], 
                               plot_name = names(degs.list)[n]
)
cdp.infVsCtrl.24h

dot_plots <- list(
        cdp.bysVsCtrl.1h,
        cdp.bysVsCtrl.24h,
        cdp.infVsBys.1h,
        cdp.infVsBys.24h,
        cdp.infVsCtrl.1h,
        cdp.infVsCtrl.24h
)

pdf('figures/dot_plots_degs_cdp.pdf', 
    width = 20)
gridExtra::grid.arrange(
        grobs=dot_plots, 
        ncol=2
)
dev.off()

