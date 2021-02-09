## Alignment and counting CEL-Seq2 reads

## Dependencies
library(scPipe)
library(Rsubread)
library(SingleCellExperiment)

set.seed(333)

path2project <- '/media/ag-cherrmann/cramirez/interaction_DC-pathogens/'
setwd(path2project)

#######################################################################
## Settings

read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
organism <- "mmusculus_gene_ensembl"
gene_id_type <- "ensembl_gene_id"
filter_settings <- list(rmlow = TRUE, 
                        rmN = FALSE, 
                        minq = 20, 
                        numbq = 2)

# FASTQ files
r1 <- 'data/tmp/merged_R1.fastq'
r2 <- 'data/tmp/merged_R2.fastq'


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
## Formatting fastq reads

## Filtering and formatting reads
dir.create('data/fastq_scpipe/')
if ( ! file.exists('data/fastq_scpipe/scpipe_formatted.fastq.gz')){
	sc_trim_barcode(
        	outfq = 'data/fastq_scpipe/scpipe_formatted.fastq.gz',
        	r1 = r2,
        	r2 = r1,
        	read_structure = read_structure,
        	filter_settings = filter_settings)
}

## Performing alignment
dir.create('data/bam_scpipe')
if ( ! file.exists('data/bam_scpipe/scpipe_alignment.bam') ){
print('Performing alignment')
Rsubread::align(
  index = "db/GRCm38/subread_index/GRCm38_p6",
  readfile1 = 'data/fastq_scpipe/scpipe_formatted.fastq.gz',
  output_file = 'data/bam_scpipe/scpipe_alignment.bam',
  nthreads = 10
)
}

## Performing annotation
bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
if ( ! file.exists('data/bam_scpipe/scpipe_annotated.bam') ) {
sc_exon_mapping(
    inbam = 'data/bam_scpipe/scpipe_alignment.bam',
    outbam = 'data/bam_scpipe/scpipe_annotated.bam',
    annofn = 'db/GRCm38/annotations/gencode.vM18.primary_assembly.annotation.gff3',
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
}


## Demultiplexing
max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
dir.create('data/scpipe_demultiplexed')
sc_demultiplex(
    inbam = 'data/bam_scpipe/scpipe_annotated.bam',
    outdir = 'data/scpipe_demultiplexed',
    bc_anno = 'data/cell_id_barcodes.csv',
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)



