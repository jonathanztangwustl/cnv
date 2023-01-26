seg_tools_path <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/scripts/cnv/seg_tools.R'
source(seg_tools_path)

# TODO: TESTING
#inputs <- list(
#    data_path = '.',
#    sample_file = 'chunk_16',
#    gene = 'dnmt3a',
#    type = 'indexcov',
#    out_dir = '../../analysis/run_1/TEST/',
#    mean_threshold = -0.1,
#    bin_threshold = 13,
#    wav_threshold = 0.042,
#    plot = TRUE,
#    save_dels = TRUE
#)

# ARGUMENTS ====================================================================
parser <- arg_parser('Perform CNV calling using CBS.')
parser <- add_argument(
    parser,
    arg = 'data_path',
    help = 'Path to data folder.',
)
parser <- add_argument(
    parser,
    arg = '--sample_file',
    help = 'Sample file.',
    default = 'sample.txt',
    short = '-s'
)
parser <- add_argument(
    parser,
    arg = '--gene',
    help = 'Gene to pull data from.',
    default = 'DNMT3A',
    short = '-g',
)
parser <- add_argument(
    parser,
    arg = '--type',
    help = 'Mosdepth or indexcov.',
    default = 'indexcov',
    short = '-t'
)
parser <- add_argument(
    parser,
    arg = '--out_dir',
    help = 'Output directory.',
    default = '../analysis/',
    short = '-o'
)
parser <- add_argument(
    parser,
    arg = '--mean_threshold',
    help = 'Mean threshold to compare bins to.',
    default = -0.1,
    short = '-m'
)
parser <- add_argument(
    parser,
    arg = '--bin_threshold',
    help = 'Bin threshold for minimum number of bins.',
    default = 13,
    short = '-b'
)
parser <- add_argument(
    parser,
    arg = '--wav_threshold',
    help = 'Waviness threshold',
    default = 0.042,
    short = '-w'
)
parser <- add_argument(
    parser,
    arg = '--plot',
    help = 'Whether to plot samples.',
    default = TRUE,
    short = '-p'
)
#parser <- add_argument(
#    parser,
#    arg = '--save_dels',
#    help = 'Whether to save deletion outputs. Warning: do not disable unless you want to skip the outputs!',
#    default = TRUE,
#    short = '-d'
#)
inputs <- parse_args(parser)

# SCRIPT =======================================================================
# Run segmentation -------------------------------------------------------------
if (!file.exists(inputs$out_dir)) {
    dir.create(inputs$out_dir, recursive = TRUE)
}
results <- del_gene_batch(inputs)
gc()

# Log statistics ---------------------------------------------------------------
run_log <- tibble(
    total_samples = results$data %>% length,
    filt_samples = results$data %>% lapply(`[[`, 'waviness') %>% unlist %>% `<`(inputs$wav_threshold) %>% sum,
    deletions = results$deletions %>% nrow,
    deletion_samples = results$deletions$ID %>% unique %>% length,
    filt_deletions = results$deletions %>% filter(waviness < inputs$wav_threshold) %>% nrow,
    filt_deletion_samples = results$deletions %>% filter(waviness < inputs$wav_threshold) %>% .$ID %>% unique %>% length
)
run_log$pct <- run_log$filt_deletion_samples / run_log$filt_samples
write.table(
    run_log,
    paste0(inputs$out_dir, '/', inputs$sample_file, '_stats.txt'),
    sep = '\t',
    quote = F,
    row.names = F
)

# Output all deletions ---------------------------------------------------------
all_dels <- results$deletions[results$deletions$waviness < inputs$wav_threshold, ]
all_dels$size <- all_dels$num.mark * 16384
all_dels$vaf <- 1-(2^all_dels$seg.mean)
write.table(
    all_dels,
    paste0(inputs$out_dir, '/', inputs$sample_file, '_all_deletions.txt'),
    sep = '\t',
    quote = F,
    row.names = F
)

# Output deletions in genes of interest ----------------------------------------
subset_genes <- read.table(SUBSET_REGIONS, sep = '\t', header = FALSE, col.names = c('chromosome', 'start', 'end', 'name', 'zero', 'strand')) %>% as_tibble
subset_genes$name <- toupper(subset_genes$name)
subset_genes$chromosome <- str_replace(subset_genes$chromosome, 'chr', '')
subset_genes$chromosome <- subset_genes$chromosome %>% as.numeric
check_chroms <- results$deletions$chrom %>% unique
check_genes <- subset_genes %>% filter(chromosome %in% check_chroms)
all_dels <- vector('list', nrow(check_genes))
for (gene_i in 1:nrow(check_genes)) {
    gene_info <- check_genes[gene_i, ]
    gene_dels <- results$deletions[!(results$deletions$chrom == gene_info$chromosome & (results$deletions$loc.start >= gene_info$end | results$deletions$loc.end <= gene_info$start)), ]
    if (nrow(gene_dels) > 0) {
        gene_dels$name <- gene_info$name
    } else {
        gene_dels$name <- vector('character', 0)
    }
    all_dels[[gene_i]] <- gene_dels
}
all_dels <- all_dels %>% bind_rows
all_dels <- all_dels[all_dels$waviness <= inputs$wav_threshold, ]
all_dels$size <- all_dels$num.mark * 16384
all_dels$vaf <- 1-(2^all_dels$seg.mean)

write.table(
    all_dels,
    paste0(inputs$out_dir, '/', inputs$sample_file, '_gene_deletions.txt'),
    sep = '\t',
    quote = F,
    row.names = F
)

# Plot results that meet waviness threshold ------------------------------------
if (inputs$plot) {
    results$data <- results$data[results$data %>% lapply(`[[`, 'waviness') %>% unlist %>% `<`(inputs$wav_threshold)]
    plot_raw_batch(results, inputs)
}

