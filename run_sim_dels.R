sim_del_path <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/scripts/cnv/sens_tools.R'
source(sim_del_path)

# ARGUMENTS ====================================================================
parser <- arg_parser('Perform CNV calling using CBS.')
parser <- add_argument(
    parser,
    arg = 'data_path',
    help = 'Path to data file.',
)
parser <- add_argument(
    parser,
    arg = '--chromosome',
    help = 'Chromosome',
    default = 'chr2',
    short = '-c'
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
    arg = '--out_dir',
    help = 'Output directory.',
    default = '../analysis/sim_dels/100/',
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
    help = 'Minimum number of bins to call a deletion.',
    default = 13,
    short = '-b'
)
inputs <- parse_args(parser)

# SCRIPT =======================================================================
options(scipen = 999)
results <- seg_del_batch(inputs$data_path, inputs$chromosome, inputs$gene, inputs$mean_threshold, inputs$bin_threshold)
meta_data <- get_meta(inputs$data_path)
if (!file.exists(inputs$out_dir)) {
    dir.create(inputs$out_dir)
}
#plot_batch_del_seg(results,
#    out_dir = inputs$out_dir,
#    gene = inputs$gene,
#    cov_type = meta_data$type,
#    mean_thresh = inputs$mean_threshold,
#    bin_thresh = inputs$bin_threshold
#    )
write.table(
    results %>% write_sens,
    file = paste0(
        inputs$out_dir,
        meta_data$set_name, '_',
        meta_data$type, '_',
        inputs$mean_threshold, '_',
        inputs$bin_threshold, '_',
        '.tsv'
    ),
    quote = F,
    sep = '\t',
    row.names = F
)

#results <- batch_dels(inputs$data_path, inputs$gene, inputs$mean_threshold, inputs$bin_threshold)
#meta_data <- get_meta(inputs$data_path)
#print(inputs)
#write.table(results[-1],
#    file = paste0(
#        inputs$out_dir,
#        meta_data$set_name, '_',
#        meta_data$type, '_',
#        inputs$mean_threshold, '_',
#        inputs$bin_threshold,
#        '.tsv'
#        ),
#    quote = FALSE,
#    sep = '\t',
#    row.names = FALSE
#)