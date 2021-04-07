## Alignment and counting CEL-Seq2 reads using scPipe
## Modified from a previous implementation

## Dependencies
library(scPipe)
library(Rsubread)
library(SingleCellExperiment)
library(readxl)
library(janitor)
library(dplyr)

set.seed(333)

path2project <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/'
setwd(path2project)

########################################################################

file_nn99 <- 'data/NN99_C025 Stella sample info.xlsx'

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn99,
  sheet = "Sample information",
  skip = 5,
  n_max = 1)
# NOTE: FACS data in columns >= "R"
facs_data_idx <- seq(which(LETTERS == "R"), ncol(header_row))
# NOTE: A lot of faffing about to get clean column names.
header_row <- c(
  paste0(
    colnames(header_row[, -facs_data_idx]),
    ifelse(
      grepl("e\\.g|applicable", unlist(header_row[1, -facs_data_idx])),
      "",
      unlist(header_row[1, -facs_data_idx]))),
  unlist(header_row[1, facs_data_idx], use.names = FALSE))
header_row <- gsub("\\.+[0-9]+", "", header_row)
sample_sheet_nn99 <- read_excel(
  path = file_nn99,
  sheet = "Sample information",
  skip = 7,
  col_names = header_row)
# Ensure FACS columns are stored as numeric (readxl sometimes fails, presumably
# to weird pattern of empty cells).
sample_sheet_nn99 <- sample_sheet_nn99 %>%
  mutate_at(facs_data_idx, as.numeric)
sample_sheet_nn99 <- bind_cols(
  clean_names(sample_sheet_nn99[, -facs_data_idx]),
  clean_names(sample_sheet_nn99[, facs_data_idx], case = "parsed"))
sample_sheet_nn99 <- remove_empty(sample_sheet_nn99)
# Filter out those without a cell index sequence or without an Illumina index
# read or that is an empty well.
sample_sheet_nn99 <- sample_sheet_nn99 %>%
  filter(
    !is.na(plate_number),
    sample_name != "empty")
# Some final tidying.
sample_sheet_nn99 <- sample_sheet_nn99 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    sequencing_run = "NN99") %>%
  arrange(plate_number, well_position)

sample_sheet <- sample_sheet_nn99 %>%
  mutate(rowname = paste0(plate_number, "_", well_position)) %>%
  tibble::column_to_rownames("rowname") %>%
  DataFrame(., check.names = FALSE)


#######################################################################
## Getting fastq files

fastq.dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/fastq/'
file.names <- list.files(fastq.dir) 
file.names <- grep('fastq.gz$', file.names, value = TRUE)

## file name patterns
file.pattern <- gsub('_R[1-2].fastq.gz', '', file.names) %>% unique()
file.pattern

#######################################################################
## Creating folder to store results

output_dir <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/data/scpipe_updated/'
if ( ! dir.exists(output_dir)) {
  dir.create(output_dir)
}

## Creating folders for the output of each fastq 
for (pattern in file.pattern) {
  output_folder <- paste0(output_dir, pattern)
  dir.create(output_folder)
}
## Checking directories
list.dirs(output_dir)

#######################################################################
## scPipe Settings

read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
organism <- "mmusculus_gene_ensembl"
gene_id_type <- "ensembl_gene_id"
filter_settings <- list(rmlow = TRUE, 
                        rmN = FALSE, 
                        minq = 20, 
                        numbq = 2)

# FASTQ files
#r1s <- grep('R1.fastq$', list.files(fastq.dir, full.names = TRUE), value = TRUE) 
#r2s <- grep('R2.fastq$', list.files(fastq.dir, full.names = TRUE), value = TRUE)


## Getting gene annotations from gencode version m38.p6
url_anns <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.primary_assembly.annotation.gff3.gz'
dir.create('db')
dir.create('db/GRCm38/')
dir.create('db/GRCm38/annotations')
if ( ! file.exists('db/GRCm38/annotations/gencode.vM18.primary_assembly.annotation.gff3') ) {
     download.file(
         url_anns,
         destfile = 'db/GRCm38/annotations/gencode.vM18.primary_assembly.annotation.gff3'
     )
}

## Downloading fasta sequences 
dir.create('db/GRCm38/genome_fasta/')
if ( ! file.exists('db/GRCm38/genome_fasta/GRCm38.p6.genome.fa.gz') ){
        fasta_gen_url <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.p6.genome.fa.gz'
        download.file(fasta_gen_url, 
                      destfile = 'db/GRCm38/genome_fasta/GRCm38.p6.genome.fa.gz')
}

## Building the reference
dir.create('db/GRCm38/subread_index')
if ( ! file.exists('db/GRCm38/subread_index/GRCm38_p6.reads')){
     buildindex(basename="db/GRCm38/subread_index/GRCm38_p6", 
                reference= 'db/GRCm38/genome_fasta/GRCm38.p6.genome.fa.gz')
}


############################################################################
## Filtering and formatting reads

## Creating folders to store formatted fastq
ffastq_folders <- paste0(list.dirs(output_dir, recursive = F), '/formated_fastq')
sapply(ffastq_folders, dir.create)
## checking folders were created

## fastq formatting
for (i in 1:length(file.pattern)) {
  r1 <- paste0(fastq.dir, file.pattern[i], '_R1.fastq.gz')
  r2 <- paste0(fastq.dir, file.pattern[i], '_R2.fastq.gz')
  ffastq_out <- paste0(output_dir, file.pattern[i], 
                       '/formated_fastq/', 
                       file.pattern[i], 
                       '_formatted.fastq.gz')
  if ( ! file.exists(ffastq_out)){
    sc_trim_barcode(
      outfq = ffastq_out,
      r1 = r2,
      r2 = r1,
      read_structure = read_structure,
      filter_settings = filter_settings)
  }
}

##########################################################################
## Alignment

## Creating directories to store the outputs
bam_folders <- paste0(list.dirs(output_dir,recursive = F), '/bam')
sapply(bam_folders, dir.create)
## Performing alignment
for (i in 1:length(file.pattern)) {
  formatted_fastq <- paste0(output_dir, 
                            file.pattern[i],
                            '/formated_fastq/',
                            file.pattern[i],
                            '_formatted.fastq.gz')
  output_bam <- paste0(output_dir, 
                       file.pattern[i],
                       '/bam/',
                       file.pattern[i],
                       '_aligned.bam')
  if ( ! file.exists(output_bam) ){
    print('Performing alignment')
    Rsubread::align(
      index = "db/GRCm38/subread_index/GRCm38_p6",
      readfile1 = formatted_fastq,
      output_file = output_bam,
      nthreads = 10
    )
  }
} 

##########################################################################
## Performing annotation

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
for (i in 1:length(file.pattern)) {
  in_bam <- paste0(output_dir, '/', 
                   file.pattern[i], '/bam/', 
                   file.pattern[i], '_aligned.bam')
  out_bam <- paste0(output_dir, '/', 
                    file.pattern[i], '/bam/', 
                    file.pattern[i], 'aligned_annotated.bam')
  if ( ! file.exists(out_bam) ) {
    sc_exon_mapping(
      inbam = in_bam,
      outbam = out_bam,
      annofn = 'db/GRCm38/annotations/gencode.vM18.primary_assembly.annotation.gff3',
      bam_tags = bam_tags,
      bc_len = bc_len,
      barcode_vector = barcode_vector,
      UMI_len = UMI_len,
      stnd = stnd,
      fix_chr = fix_chr)
  }
}

##############################################################################

## Mapping
dplyr::select(as.data.frame(sample_sheet), 
                         illumina_rpi_number, plate_number, 
                         illumina_index_index_sequence_separate_index_read) %>% 
           unique() %>%
            write.table(file = 'data/mapping_metadata.tsv',
                        sep = '\t',
                        row.names = FALSE)

## Splitting annotations by ilumina indexing
anns_splitted <- split(sample_sheet, 
                       sample_sheet$illumina_index_index_sequence_separate_index_read) 
anns_splitted <- lapply(anns_splitted,  
                        function(x) 
                          as.data.frame(x) %>% 
                          add_rownames(var='cell_id') %>%
                          dplyr::select(cell_id, rd1_index_cell_index_index_sequence_as_in_c_rt1_primer) %>%
                          dplyr::rename(barcode=rd1_index_cell_index_index_sequence_as_in_c_rt1_primer) %>%
                          dplyr::mutate(barcode=strtrim(barcode, 7))
)
## Checking non-duplications
lapply(anns_splitted, function(x) x$'barcode' %>% duplicated %>% sum)
lapply(anns_splitted, function(x) x$'barcode' %>% sort)
lapply(anns_splitted, head)
lapply(anns_splitted, dim)

## Mapping to fastq file names
mapping <- read_excel('data/mapping_fastq_to_anns.xlsx', 
                      sheet = 'mapping')
rownames(mapping) <- mapping$ilumina_index
mapping$fastq_file <- paste0(mapping$fastq_file, '.gz') 
mapping
## Checking ilumina indexes intersect in annotations and mapping
intersect(names(anns_splitted), mapping$ilumina_index)
## Writing split annotations
for (i in 1:nrow(mapping)){
  bc <- mapping$ilumina_index[i]
  ann_bc <- anns_splitted[bc][[1]]
  pattern <- gsub('_R1.fastq.gz', '', mapping$fastq_file[i])
  barcodes_file <- paste0(output_dir, pattern, '/', pattern, '_cell_id_barcodes.csv')
  cat(barcodes_file, ' : ', bc, ' : ',dim(ann_bc), '\n')
  write.csv(ann_bc, file = barcodes_file, row.names = FALSE, quote = FALSE)
}

#########################################################################################
## Demultiplexing

max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
## Creating folders to store results
demult_foders <- paste0(output_dir, file.pattern, '/demultiplexed/') 
sapply(demult_foders, dir.create)
for (i in 1:length(file.pattern) ){
  in_bam <- paste0(output_dir, '/', 
                   file.pattern[i], '/bam/', 
                   file.pattern[i], 'aligned_annotated.bam')
  out_dir <- paste0(output_dir, '/', 
                    file.pattern[i], '/demultiplexed/')
  barcodes <- paste0(output_dir, '/', 
                     file.pattern[i], 
                     '/', file.pattern[i], 
                     '_cell_id_barcodes.csv')
  sc_demultiplex(
      inbam = in_bam,
      outdir = out_dir,
      bc_anno = barcodes, ## change this
      max_mis = max_mis,
      bam_tags = bam_tags,
      mito = mito,
      has_UMI = has_UMI)
}

#########################################################################################
## Counting
UMI_cor <- 1
gene_fl <- FALSE
for (i in 1:length(file.pattern) ){
  demultiplex_dir <- paste0(output_dir, '/', 
                    file.pattern[i], '/demultiplexed/')
  barcodes <- paste0(output_dir, '/', 
                     file.pattern[i], 
                     '/', file.pattern[i], 
                     '_cell_id_barcodes.csv')
  sc_gene_counting(demultiplex_dir, 
                   barcodes,
                   UMI_cor = UMI_cor,
                   gene_fl = gene_fl)
}

########################################################################################
## Quality control of the alignments
counts_dir <- paste0(output_dir, '/', file.pattern, '/demultiplexed/')

sce.list <- lapply(counts_dir, create_sce_by_dir)
names(sce.list) <- file.pattern

pdf('figures/barcode_demultiplexing_UPDATE.pdf')
demultiplexing.plots <- lapply(sce.list, plot_demultiplex)
demultiplexing.plots <- lapply(1:length(demultiplexing.plots), 
                               function(i) demultiplexing.plots[i][[1]] +
                                            ggtitle(
                                              names(demultiplexing.plots)[i]
                                              ))
lapply(demultiplexing.plots, plot)
dev.off()

pdf('figures/umi_deduplications_UPDATE.pdf')
deduplication.plots <- lapply(sce.list, plot_UMI_dup, log10_x = FALSE)
deduplication.plots <- lapply(1:length(deduplication.plots), 
                               function(i) deduplication.plots[i][[1]] +
                                 ggtitle(
                                   names(deduplication.plots)[i]
                                 ))
deduplication.plots
dev.off()

pdf('figures/mapping_UPDATE.pdf')
lapply(1:length(sce.list), 
       function(i)
         plot_mapping(sce.list[i][[1]], 
                      percentage = TRUE, dataname = names(sce.list)[i]))

dev.off()

cell_ids <- c()
for (i in 1:length(sce.list)) {
  ## extract cells with counts > 5
  cells_above <- apply(assay(sce.list[i][[1]]), 2, function(x) sum(x) > 0)
  cells_above <- colnames(assay(sce.list[i][[1]]))[cells_above]
  cell_ids <- c(cell_ids, cells_above)
}
length(unique(cell_ids))

sample_sheet.sel <- sample_sheet[cell_ids,]
dim(sample_sheet.sel)
table(sample_sheet.sel$sample_name)
 dim(sample_sheet)

sample_sheet %>%
  as.data.frame() %>%
  filter(sample_name == 'single cell') %>%
  dim



