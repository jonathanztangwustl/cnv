# TODO: This file has been moved to the ichorCNA repo

library('R.utils')
library('tidyverse')
library('magrittr')
library('argparser')
library('parallel')

# ARGUMENT PARSER --------------------------------------------------------------
#parser <- arg_parser('Perform CNV calling on Indexcov outputs using IchorCNA.')
#parser <- add_argument(
#    parser,
#    arg = 'data_path',
#    help = 'Path to the data folder.'
#)
#parser <- add_argument(
#    parser,
#    arg = '--offsets',
#    help = 'Offsets file to correct indexcov to WIG.',
#    default = '/scratch1/fs1/timley/fusions/jonathanztang/scripts/cnv/wig_offsets'
#    short = '-x'
#)
#parser <- add_argument(
#    parser,
#    arg = '--sample_list',
#    help = 'File including all sample names (within data_path directory).'
#    default = 'sample',
#    short = '-s'
#)
#parser <- add_argument(
#    parser,
#    arg = '--output',
#    help = 'Output directory.',
#    default = 'ichor_out',
#    short = '-o'
#)
## TODO: add Ichor inputs as necessary
#inputs <- list(
#    data_path = '.',
#    offsets = '/scratch1/fs1/timley/fusions/jonathanztang/scripts/cnv/wig_offsets',
#    sample_list = '.',
#    output = './ichor_out'
#)
#inputs <- parse_args(parser)
OFFSETS_PATH <- '/scratch1/fs1/timley/fusions/jonathanztang/scripts/cnv/wig_offsets'
OFFSETS <- read_table(OFFSETS_PATH)

# FUNCTIONS --------------------------------------------------------------------
# Load data
#   - Loads data from .bed.gz file. Interprets input file name to determine if
#     loading indexcov vs. mosdepth data. Untars, then gunzips file to a temp
#     directory, then fixes column names and returns data as a tibble.
#   - Inputs
#       - data_path: path to data file in .tar.gz format
#   - Outputs
#       - set_data: tibble of raw BED coverage data
load <- function(data_path) {
    set_name <- str_extract(data_path, "NWD[0-9]*")
    temp_dir <- paste0(tempdir(), '/', set_name, '/')
    bed_file <- 'indexcov/indexcov-indexcov.bed.gz'
    set_file <- paste0(temp_dir, '/', str_sub(bed_file, 1, -4))
    untar(
        data_path,
        files = bed_file,
        exdir = temp_dir
        )
    gunzip(paste0(temp_dir, '/', bed_file), overwrite = TRUE)
    set_data <- read_tsv(
        set_file,
        col_names = TRUE,
        col_types = 'ciid'
    )
    colnames(set_data) <- c('chromosome', 'start', 'end', 'depth')
    unlink(temp_dir)
    return(set_data)
}

# Function to take a single string and append an offset based on
# chromosome lookup from offsets. Appends n entries of depth 1.
adjust_strings <- function(chr_string, offsets) {
    chromosome <- chr_string[1] %>%
    str_extract('chr[0-9, X, Y].?\\s') %>%
        str_replace(' ', '')
    off <- offsets$value[offsets$name == chromosome]
    chr_string <- c(chr_string, rep(1, off + 1))
    return(chr_string)
}

# CONVERSION -------------------------------------------------------------------
# Single file conversion
#   Converts a single indexcov.tar.gz file to a WIG file in temporary directory.
#   Designed to work with the batch function and return a path to a temp file.
index_to_wig <- function(file_path) {
    ind <- load(file_path)
    
    # Filter only autosomal and X/Y chromosomes
    ind <- ind %>% filter(
        chromosome %in% c(
            paste0('chr', c(1:22, 'X', 'Y'))
        )
    )

    # Split and convert each chromosome to string vectors
    ind_chrs <- ind %>% group_by(chromosome) %>% group_split
    ind_strings <- vector('list', length(ind_chrs))
    ind_names <- vector('character', length(ind_chrs))
    for (i in 1:length(ind_chrs)) {
        wig_line <- ind_chrs[[i]]$depth
        ind_strings[[i]] <- wig_line
        ind_names[i] <- ind_chrs[[i]]$chromosome %>% unique %>% .[1]
    }
    names(ind_strings) <- ind_names

    # Reorder strings per chromosome
    ind_order <- ind_names %>%
        str_extract('chr[0-9, X, Y].?') %>%
        str_replace('chr', '') %>%
        str_replace('X', '23') %>%
        str_replace('Y', '24') %>%
        as.numeric
    ind_strings <- ind_strings[order(ind_order)]

    # Add headers
    ind_headers <- paste0(
        'fixedStep chrom=',
        names(ind_strings),
        ' start=1 step=16384 span=16384'
    )
    for (i in 1:length(ind_strings)) {
        ind_strings[[i]] <- c(ind_headers[i], ind_strings[[i]])
    }

    # Length adjustments
    ind_strings <- ind_strings %>% lapply(adjust_strings, OFFSETS)

    # Save wig
    wig <- unlist(ind_strings)
    temp_wig <- tempfile(fileext = '.wig')
    write(wig, temp_wig)
    return(temp_wig)
}

batch_i2w <- function(in_dir, out_dir) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    sample_list <- in_dir %>% list.dirs %>% str_subset('NWD') %>% str_extract('NWD[0-9]*')
    sample_list <- sample_list %>% paste0('./', ., '/indexcov.tar.gz')

    cluster <- makeForkCluster(detectCores() - 1)
    samples <- clusterMap(
        cluster,
        index_to_wig,
        sample_list
    )
    stopCluster(cluster)

    sample_table <- samples %>% unlist %>% enframe
    sample_table$sample <- sample_table$name %>% str_extract('NWD[0-9]*')
    
    calls <- paste0(
        'cp ',
        sample_table$value,
        ' ',
        out_dir,
        '/',
        sample_table$sample,
        '.wig'
    )

    lapply(calls, system)
}